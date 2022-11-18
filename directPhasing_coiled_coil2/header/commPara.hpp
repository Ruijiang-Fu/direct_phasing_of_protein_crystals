/*******************************************************************************

                        include standard library

*******************************************************************************/

//standard libraries
#include <iostream> //standard library of input and output
#include <cmath>    //standard library of math
#include <complex>  //standard library of complex numbers
#include <cstdlib>  //srand and rand, exit(EXIT_FAILURE)
#include <fstream>  //standard library for file operations
#include <iomanip>  //setw()
#include <stdlib.h> //srand, rand
#include <time.h>   //time
#include <math.h>   //round, floor, ceil, trunc
#include <limits>   //infinity
#include <string>
#include <sstream>  //string stream operations
#include <vector>
#include <unistd.h> //getpid()

//header files for Intel dfti fft
#include "mkl_dfti.h" //header file of Intel Discrete Fourier Transform Interface

//constants
#define PI 3.14159265358979323846
#define HNUM ANUM
#define KNUM BNUM
#define LNUM CNUM

using namespace std;


/*******************************************************************************

                            declare variable

*******************************************************************************/

//Constants for Calculating the Resolution of Each Reflection
const float S1 = pow( ( CRYST_B * CRYST_C * sin( CRYST_ALPHA ) ), 2.0 );
const float S2 = pow( ( CRYST_A * CRYST_C * sin( CRYST_BETA ) ), 2.0 );
const float S3 = pow( ( CRYST_A * CRYST_B * sin( CRYST_GAMMA ) ), 2.0 );
const float S4 = CRYST_A * CRYST_B * CRYST_C * CRYST_C * ( cos( CRYST_ALPHA ) * cos( CRYST_BETA ) - cos( CRYST_GAMMA ) );
const float S5 = CRYST_A * CRYST_B * CRYST_B * CRYST_C * ( cos( CRYST_GAMMA ) * cos( CRYST_ALPHA ) - cos( CRYST_BETA ) );
const float S6 = CRYST_A * CRYST_A * CRYST_B * CRYST_C * ( cos( CRYST_BETA ) * cos( CRYST_GAMMA ) - cos( CRYST_ALPHA ) );
const float NU = sqrt(1.0 - pow(cos(CRYST_ALPHA),2.0) - pow(cos(CRYST_BETA),2.0) - pow(cos(CRYST_GAMMA),2.0) \
                 + 2.0 * cos(CRYST_ALPHA) * cos(CRYST_BETA) * cos(CRYST_GAMMA));
const float CRYST_VOL = CRYST_A * CRYST_B * CRYST_C * NU;

//variables for Intel Discrete Fourier Transform Interface
const float dftiForScaleFactor = 1.0 / CRYST_VOL;
const float dftiBacScaleFactor = CRYST_VOL / float( ANUM * BNUM * CNUM );

//variable for input output files
stringstream fileNameStream;
string fileNameString;

//variable for sf and atom
int numUniqStruFact;
int numAtomAsu;

//cutoff value of density or weighted average density
float densCutoff, densCutoff1, densCutoff2;

//variables for Intel DFTI
DFTI_DESCRIPTOR_HANDLE myDescriptorHandle;
MKL_LONG status;
MKL_LONG lengthEachDimension[3] = {ANUM, BNUM, CNUM};
complex<float> realSpaceArray[ANUM][BNUM][CNUM] = {};
complex<float> reciprocalSpaceArray[HNUM][KNUM][LNUM] = {};

//variables for a whole unit cell
float densUnitCell[ANUM][BNUM][CNUM] = {};
float densUnitCellTemp[ANUM][BNUM][CNUM] = {};
float weigAvgDensUnitCell[ANUM][BNUM][CNUM] = {};
float weigAvgDensUnitCellTemp[ANUM][BNUM][CNUM] = {};
int maskUnitCell[ANUM][BNUM][CNUM] = {};
int maxPackSphereRadiusUnitCell[ANUM][BNUM][CNUM] = {};
int maxPackCylinderRadiusUnitCell[ANUM][BNUM][CNUM] = {};
float ncsCentMassCandUnitCell[ANUM][BNUM][CNUM] = {};

//convolution mask, such as sphere or cylinder
float realConvMaskUnitCell[ANUM][BNUM][CNUM] = {};
complex<float> reciprocalConvMaskUnitCell[HNUM][KNUM][LNUM] = {};

//convoluted density obtained with a covolution mask
float convDensUnitCell[ANUM][BNUM][CNUM] = {};

//use a cylinder as the convolution mask
int numGridPointInCylinder;
float cylinderDens[ANUM*BNUM*CNUM][numNcsOper][4] = {}; //numGridPointInCylinder < ANUM*BNUM*CNUM

//variables to find densMin, densMax, densAvg
const float INF = numeric_limits<float>::infinity(); //infinite
float densMin, densMax, densAvg;

//variables for Intel dfti fft
complex<float> halfGridPosi[HNUM][KNUM][LNUM] = {};
complex<float> halfGridNega[HNUM][KNUM][LNUM] = {};
float allStruFactPhase[HNUM][KNUM][LNUM] = {};
float allStruFactAmpl[HNUM][KNUM][LNUM] = {};

//ncs averaging
const float vectOnNcsAxisIni[3] = {sin(lateAnglPsi)*cos(azimAnglPhi), sin(lateAnglPsi)*sin(azimAnglPhi), cos(lateAnglPsi)};
//float vectOnNcsAxis[3] = {vectOnNcsAxisIni[0], vectOnNcsAxisIni[1], vectOnNcsAxisIni[2]};

float pointAOnNcsAxisCand[numNcsAxisCand][3] = {};
float vectOnNcsAxisCand[numNcsAxisCand][3] = {}; //direction of NCS axis candidates

int numIterNcsAvg = int(numIter*0.9);
const int numAdjaCell = 27;
const int numLargeBinNcs = 1000;
//float boundLargeBinNcs[numLargeBinNcs+1] = {}; //numLargeBin have numLargeBin+1 boundaries, calculated hist
const float lengEachDime[3] = {ANUM, BNUM, CNUM};
int ncsMask[NCSANUM][NCSBNUM][NCSCNUM] = {};
float ncsMaskDensTemp[NCSANUM][NCSBNUM][NCSCNUM] = {};
float ncsMaskDens[NCSANUM][NCSBNUM][NCSCNUM] = {};
int anum0, bnum0, cnum0; //qubic region hold ncs mask
int anum1, bnum1, cnum1;

//center of mass candidates, one of which will be selected as centMassOrth
int numCentMassCand; //<1000
float centMassCand[ANUM*BNUM*CNUM][5];

//hybrid input-output
float hioAlpha = 0.4; //CHIO; HPR 0.588
float hioBeta = 0.7; //0.5 ~ 1.0
float hioDensLimitIni = 0.5;
float hioDensLimitFina = 0.1;
const int numIterLimitDens = int(numIter*0.1);

//histogram matching
const int numSmallBin = 10000000; //~ numGridPointAsu * 20;
int numGridPointSmallBin[numSmallBin] = {}; //frequency or height
const int numLargeBin = 300; //100 ~ 500
int numGridPointLargeBin[numLargeBin] = {};
float boundLargeBin[numLargeBin+1] = {}; //numLargeBin have numLargeBin+1 boundaries, calculated hist
float stdBoundLargeBin[numLargeBin+1] = {}; //reference hist
float hma[numLargeBin] = {}; //parameter a in hist matching
float hmb[numLargeBin] = {}; //parameter b in hist matching

//solvent flattening
const int numIterSolvFlat = int(numIter*0.1);

//sigma of weighted observed data
const int numIterVarySigmWeigObs = numIter - numIterLimitDens - numIterSolvFlat;
float sigmWeigObsArray[numIter];
float sigmWeigAvgArray[numIter];

//vary solvent content
const int numIterVarySolvCont = numIter - numIterLimitDens - numIterSolvFlat;

//vary sigmWeigAvg
const int numIterVarySigmWeigAvg = numIter - numIterLimitDens - numIterSolvFlat;

//metrics of convergence
int goodProtMaskFlag = 0; //good protein mask flag
int convFlag = 0; //covergence flag
float F000[numIter] = {};
float scaleFact[numIter] = {};
float Rfree[numIter] = {};
float Rwork[numIter] = {};

//27 adjacent cells, 9 cells on top, 9 cells in the middle, and 9 cells at bottom
const float adjaCell[numAdjaCell][3] = { 
    {-1.0, -1.0, -1.0}, { 0.0, -1.0, -1.0}, { 1.0, -1.0, -1.0},
    {-1.0,  0.0, -1.0}, { 0.0,  0.0, -1.0}, { 1.0,  0.0, -1.0},
    {-1.0,  1.0, -1.0}, { 0.0,  1.0, -1.0}, { 1.0,  1.0, -1.0},
    {-1.0, -1.0,  0.0}, { 0.0, -1.0,  0.0}, { 1.0, -1.0,  0.0},
    {-1.0,  0.0,  0.0}, { 0.0,  0.0,  0.0}, { 1.0,  0.0,  0.0}, 
    {-1.0,  1.0,  0.0}, { 0.0,  1.0,  0.0}, { 1.0,  1.0,  0.0}, 
    {-1.0, -1.0,  1.0}, { 0.0, -1.0,  1.0}, { 1.0, -1.0,  1.0}, 
    {-1.0,  0.0,  1.0}, { 0.0,  0.0,  1.0}, { 1.0,  0.0,  1.0}, 
    {-1.0,  1.0,  1.0}, { 0.0,  1.0,  1.0}, { 1.0,  1.0,  1.0} };

/*******************************************************************************

                            declare function

*******************************************************************************/
//the default value is specified in function declaration

//random seed generator
unsigned long mixClockTimeGetid();

//cell parameters in CRYST1 of pdb file
string writeCryst1P1();
string writeCryst1();

//count number of atoms in a pdb file
int countNumAtomModel(string fileNameString);

//input atom string/coordinates from pdb
void inputAtomString(string fileNameString, string *atomStringArray);
void inputAtomCoor(string fileNameString, float **atomCoorArray);

//compute the centroid of ATOM in the deposited model
void computeCentOrthModel(string *atomStringArray, float *centOrth);

//compute average temperature factor of ATOM and HETATM
float computeAvgTemp(string *atomStringArray, string fileNameString);

//adjust temperature factor of ATOM and HETATM
void adjustTemp(string fileNameString, float tempAvg0, float tempDevi);

//operate fractional and orthogonal coordinates
void adjustFracCoorIntoUnitCell( float *frac);
void convertOrthToFrac( float *orth, float *frac);
void convertFracToOrth( float *frac, float *orth);

//ncs axis candidates around the initial Ncs direction got from self-rotation function
void generateNcsAxisCandLocation(float ncsAxisDeviDist);
void generateNcsAxisCandDirection(float ncsAxisDeviAngl);
void outputNcsAxisCand(string fileNameString);

//allocate ncs mask using a cylinder/sphere
void allocateNcsMaskCylinder();

//apply ncs rotation on pdb model
void applyNcsRotaOnPdbModel(string *atomStringArray, float thita, string fileNameString);

//apply ncs translation on pdb model
void applyNcsTranOnPdbModel(string *atomStringArray, float tranNcsDist, string fileNameString);

//apply ncs matrix (including rotation and translation) on pdb model
void applyNcsMatrixOnPdbModel(string *atomStringArray, string fileNameString);

//apply ncs rotation on a point
void applyNcsRotaOnPoint(float *point, float thita, float *newPoint);

//rotate a point by an angle thita about any axis through a given point with a direction 
void applyRotaAboutAnyAxis(float *pointOnAxis, float *vectOnAxis, float thita, float *point, float *newPoint);

//apply ncs translation on a point
void applyNcsTranOnPoint(float *point, float dist, float *newPoint);

//compute the minima of three floats
float computeMinima(float x1, float x2, float x3);

//compute the distance between two points
float distOrth(float *orth1, float *orth2);

//compute the angle between two vectors
float computeAngleTwoVect(float *vect1, float *vect2);

//compute the distance from a point to a straight line
float computeDistFromPointToLine(float *point, float *pointAOnLine, float *pointBOnLine);

//find a point with a distance away from a straight line
void findPointAwayFromLine(float dist, float *pointAOnLine, float *vectOnLine, float *point);

//rotate a point around a straight line
void rotaPointAroundLine(float *point, float anglPhi, float *pointAOnLine, float *vectOnLine, float *newPoint);

//input NCS axis candidates
void inputNcsAxisCand(string fileNameString);

//input reference histogram
void inputRefeHist(string fileNameString);

//input pre-computed sigmWeigObs and sigmWeigAvg
void inputSigmWeigObsAndSigmWeigAvg(string fileNameString);

//find the middle points of ncs related atoms
int findMiddlePointNcsAtom(string *atomStringArray, float **middlePoint);
int excludeOutliers(int numMiddlePoint, float **middlePoint);
void outputMiddlePointNcsAtom(int numMiddlePoint, float **middlePoint, string fileNameString);

//compute and output two points on ncs axis
void outputTwoPointOnNcsAxis(int numMiddlePoint, float **middlePoint, string fileNameString);

//input max pack sphere radius
void inputMaxPackSphereRadiusUnitCell(string fileNameString);
void inputMaxPackCylinderRadiusUnitCell(string fileNameString);

//output maximum-pack sphere radius in unit cell
void outputMaxPackSphereRadiusUnitCell(string fileNameString, float cutoff);
void outputMaxPackCylinderRadiusUnitCell(string fileNameString, float cutoff);

//make cylinder/sphere for convolution
void makeCylinderNcsMask();
void makeSphereNcsMask();

//find candidate grid points of the center of mass on the ncs axis
void findCentMassCand(float distAwayFromNcsAxis);
void outputCentMassCand(string fileNameString);

//update cent of mass
void updateCentMass(float densCutoff);
void updateCentMassForMonomer(float densCutoff);
void updateNcsCentMassUnitCell(float densCutoff);
void refineNcsCentMassLocally(float densCutoff);

void updateNcsCentMassUnitCellViaConvCylinder();
void updateNcsCentMassViaCylinder(float densCutoff);

//update density in cylinder, input weigAvgDensUnitCell, output cylinderDens
void updateCylinderDens(int ia0, int ib0, int ic0);

//inverse ncs operator for improper ncs
void inverseSecoNcsOperAsFirstNcsOper();

//apply Intel dfti fft
void setHalfGridPhase();
void applyFFTDensToStruFact();
void applyFFTStruFactToDens();

//calculate weighted average density to update protein boundary
void computeWeigAvgDensUnitCell(float sigmWeigAvg);
void updateWeigAvgDensUnitCell(string updateSpeed);
float computeDensCutoffUnitCell(float protCont);
float computeWeigAvgDensCutoffUnitCell(float protCont);

//enantiomorphDensUnitCell
void enantiomorphDensUnitCell();

//output calculated density
void outputDensUnitCell(string fileNameString, float densCutoff);
void outputResoSphere(string fileNameString);
void outputWeigAvgDensUnitCell(string fileNameString, float densCutoff);
void outputMaskUnitCell(string fileNameString, float densCutoff);
void outputConvMaskUnitCell(string fileNameString);
void outputConvDensUnitCell(string fileNameString, float densCutoff);
void outputNcsCentMassCandUnitCell(string fileNameString, float densCutoff);

//output F000
void outputF000(string fileNameString);

//output R value
void outputRvalue(string fileNameString, float *Rvalue);

//output reference histogram
void outputBoundEachLargeBin(string fileNameString);
void outputLogDensFreqEachLargeBin(string fileNameString);
void outputLogHist(string fileNameString);
void outputLogHistOld(string fileNameString);

//analyze results
bool isFileExist(string fileNameString);
int inputLogConvStat(string fileNameString);
int findOriginChoiceMaxOverlapDensAsu();

/*******************************************************************************

                            define function

*******************************************************************************/


//mix clock(), time(NULL), and getpid() to get a seed of random number
unsigned long mixClockTimeGetid()
{
    unsigned long i = clock();
    unsigned long j = time(NULL);
    unsigned long n = getpid();
    i=i-j;  i=i-n;  i=i^(n >> 13);
    j=j-n;  j=j-i;  j=j^(i << 8);
    n=n-i;  n=n-j;  n=n^(j >> 13);
    i=i-j;  i=i-n;  i=i^(n >> 12);
    j=j-n;  j=j-i;  j=j^(i << 16);
    n=n-i;  n=n-j;  n=n^(j >> 5);
    i=i-j;  i=i-n;  i=i^(n >> 3);
    j=j-n;  j=j-i;  j=j^(i << 10);
    n=n-i;  n=n-j;  n=n^(j >> 15);
    return n;
}

/******************************************************************************/
//cell parameters in CRYST1 of pdb file
string writeCryst1P1()
{
    stringstream cryst1;
    //--- cryst1 will be used in pdb file --------------------------------------
    cryst1 << "CRYST1" << setfill(' ') << setw(9) << fixed << setprecision(3) << CRYST_A
        << setfill(' ') << setw(9) << fixed << setprecision(3) << CRYST_B
        << setfill(' ') << setw(9) << fixed << setprecision(3) << CRYST_C
        << setfill(' ') << setw(7) << fixed << setprecision(2) << CRYST_ALPHA/PI*180
        << setfill(' ') << setw(7) << fixed << setprecision(2) << CRYST_BETA/PI*180
        << setfill(' ') << setw(7) << fixed << setprecision(2) << CRYST_GAMMA/PI*180
        << " P 1           1" ;
    return cryst1.str();
}

/******************************************************************************/

int countNumAtomModel(string fileNameString)
{
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    int iatom = 0;
    string lineString;
    while(getline(inputFile, lineString)) {
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0) iatom++;
    }
    inputFile.close();
    return iatom;
}

/******************************************************************************/

void inputAtomString(string fileNameString, string *atomStringArray)
{
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    int iatom = 0;
    string lineString;
    while(getline(inputFile, lineString)) {
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0) {
            atomStringArray[iatom] = lineString;
            iatom++;
        }
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void inputAtomCoor(string fileNameString, float atomCoorArray[][3])
{
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    int iatom = 0;
    string lineString;
    while(getline(inputFile, lineString)) {
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0) {
            stringstream( lineString.substr(30, 24) ) >> atomCoorArray[iatom][0] 
                >> atomCoorArray[iatom][1] >> atomCoorArray[iatom][2];
            iatom++;
        }
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void computeCentOrthModel(string *atomStringArray, float *centOrth)

{
    for(int i = 0; i <= 2; i++){
        centOrth[i] = 0;
    }
    for(int iatom = 1; iatom <= numAtomAsu - 1; ++iatom){
        string lineString = atomStringArray[iatom];
        float orth[3];
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0){
            stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
            for(int i = 0; i <= 2; i++){
                centOrth[i] = centOrth[i] + orth[i];
            }
        }
    }
    for(int i = 0; i <= 2; i++){
        centOrth[i] = centOrth[i] / float(numAtomAsu);
    }
    return;
}

/******************************************************************************/

float computeAvgTemp(string *atomStringArray, string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    float tempAvg = 0, tempMin = 1000.0, tempMax = -1000.0;
    for(int iatom = 1; iatom <= numAtomAsu - 1; ++iatom){
        string lineString = atomStringArray[iatom];
        float temp;
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0){
            stringstream( lineString.substr(60, 6) ) >> temp;
            if(temp > tempMax) tempMax = temp;
            if(temp < tempMin) tempMin = temp;
            tempAvg += temp;
        }
    }
    tempAvg /= float(numAtomAsu);
    outputFile << "For all ATOM and HETATM:" << endl;
    outputFile << "Minimum temperature factor is " << tempMin << endl;
    outputFile << "Maximum temperature factor is " << tempMax << endl;
    outputFile << "Average temperature factor is " << tempAvg << endl;
    outputFile.close();
    return tempAvg;
}

/******************************************************************************/

void adjustTemp(string fileNameString, float tempAvg0, float tempDevi)
{
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }

    float temp, tempAvg = 0, tempMin = 1000.0, tempMax = -1000.0;
    string lineString;
    int i = 0;
    while(getline(inputFile, lineString)){
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0){
            stringstream( lineString.substr(60, 6) ) >> temp;
            if(temp > tempMax) tempMax = temp;
            if(temp < tempMin) tempMin = temp;
            tempAvg += temp;
            i++;
        }
    }
    tempAvg /= float(i);

    inputFile.clear();
    inputFile.seekg(0, ios::beg);

    ofstream outputFile( "after_adjust_temp.pdb" );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }

    tempMax = tempMax + (tempAvg0 - tempAvg);
    tempMin = tempMin + (tempAvg0 - tempAvg);
    while(getline(inputFile, lineString)){
        //remove anisotropic temperature factors
        if(lineString.compare(0, 6, "ANISOU") == 0) continue;
        //adjust temperature factor
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0){
            stringstream( lineString.substr(60, 6) ) >> temp;
            temp = temp + (tempAvg0 - tempAvg);
            if(temp > tempAvg0)    
                temp = tempAvg0 + tempDevi * (temp - tempAvg0) / (tempMax - tempAvg0);
            else if(temp<tempAvg0)
                temp = tempAvg0 - tempDevi * (tempAvg0 - temp) / (tempAvg0 - tempMin);
            stringstream lineStream;
            lineStream << lineString.substr(0, 60)
                << setfill(' ') << setw(6) << setprecision(2) << fixed << temp
                << lineString.substr(66, 14);
            outputFile << lineStream.str() << endl;
        } else
            outputFile << lineString << endl;
    }

    inputFile.close();
    outputFile.close();
    return;
}

/******************************************************************************/

void adjustFracCoorIntoUnitCell( float *frac)
{
    if( frac[0] >= 0 && frac[0] < 1.0 && frac[1] >= 0 && frac[1] < 1.0 && frac[2] >= 0 && frac[2] < 1.0) 
        return;
    else {
        //make tranlations
        for( int i=0; i<=2; i++ ) 
        {
            while (frac[i] < 0) frac[i] += 1.0;
            while (frac[i] >= 1.0) frac[i] -= 1.0;
        }
        return;
    }
}

/******************************************************************************/

void convertOrthToFrac( float *orth, float *frac)
{
    float matrix[3][3];

    matrix[0][0] = 1.0/CRYST_A;
    matrix[0][1] = -cos(CRYST_GAMMA)/(CRYST_A*sin(CRYST_GAMMA));
    matrix[0][2] = (cos(CRYST_ALPHA)*cos(CRYST_GAMMA)-cos(CRYST_BETA))/(CRYST_A*sin(CRYST_GAMMA)*NU);
    matrix[1][0] = 0;
    matrix[1][1] = 1.0/(CRYST_B*sin(CRYST_GAMMA));
    matrix[1][2] = -(cos(CRYST_ALPHA)-cos(CRYST_BETA)*cos(CRYST_GAMMA))/(CRYST_B*sin(CRYST_GAMMA)*NU);
    matrix[2][0] = 0;
    matrix[2][1] = 0;
    matrix[2][2] = sin(CRYST_GAMMA)/(CRYST_C*NU);

    for( int i=0; i<=2 ; i++ ){
        frac[i] = 0;
        for( int j=0; j<=2 ; j++ ){
            frac[i] = frac[i] + matrix[i][j] * orth[j];
        }
    }
    return;
}

/******************************************************************************/

void convertFracToOrth( float *frac, float *orth)
{
    float matrix[3][3];

    matrix[0][0] = CRYST_A;
    matrix[0][1] = CRYST_B*cos(CRYST_GAMMA);
    matrix[0][2] = CRYST_C*cos(CRYST_BETA);
    matrix[1][0] = 0;
    matrix[1][1] = CRYST_B*sin(CRYST_GAMMA);
    matrix[1][2] = CRYST_C*( cos(CRYST_ALPHA)-cos(CRYST_BETA)*cos(CRYST_GAMMA) )/sin(CRYST_GAMMA);
    matrix[2][0] = 0;
    matrix[2][1] = 0;
    matrix[2][2] = CRYST_C*NU/sin(CRYST_GAMMA);

    for( int i=0; i<=2 ; i++ ){
        orth[i] = 0;
        for( int j=0; j<=2 ; j++ ){
            orth[i] = orth[i] + matrix[i][j] * frac[j];
        }
    }
    return;
}

/******************************************************************************/

void generateNcsAxisCandLocation(float ncsAxisDeviDist)
{
    int icand = 0;
    int numDist = 10;

    for(int i = 0; i <= numDist; ++i){
        float dist = float(i) * (ncsAxisDeviDist/numDist);
        float pointAOnLine[3] = {pointAOnNcsAxis[0], pointAOnNcsAxis[1], pointAOnNcsAxis[2]};
        float vectOnLine[3] = {vectOnNcsAxis[0], vectOnNcsAxis[1], vectOnNcsAxis[2]};
        float point[3];

        if(i == 0) {
            vectOnNcsAxisCand[icand][0] = vectOnNcsAxis[0];
            vectOnNcsAxisCand[icand][1] = vectOnNcsAxis[1];
            vectOnNcsAxisCand[icand][2] = vectOnNcsAxis[2];
            pointAOnNcsAxisCand[icand][0] = pointAOnNcsAxis[0];
            pointAOnNcsAxisCand[icand][1] = pointAOnNcsAxis[1];
            pointAOnNcsAxisCand[icand][2] = pointAOnNcsAxis[2];
            icand++;
            continue;
        }

        findPointAwayFromLine(dist, pointAOnLine, vectOnLine, point);

        for(int j = 0; j < (numNcsAxisCand-1)/numDist; ++j){

            float anglePhi = j * 2.0 * PI / float((numNcsAxisCand-1)/numDist);
            float newPoint[3];

            rotaPointAroundLine(point, anglePhi, pointAOnLine, vectOnLine, newPoint);

            pointAOnNcsAxisCand[icand][0] = newPoint[0];
            pointAOnNcsAxisCand[icand][1] = newPoint[1];
            pointAOnNcsAxisCand[icand][2] = newPoint[2];
            vectOnNcsAxisCand[icand][0] = vectOnNcsAxis[0];
            vectOnNcsAxisCand[icand][1] = vectOnNcsAxis[1];
            vectOnNcsAxisCand[icand][2] = vectOnNcsAxis[2];

            icand++;
        }
    }

    return;
}

/******************************************************************************/

void generateNcsAxisCandDirection(float ncsAxisDeviAngl)
{
    float psi0;
    float sinPart = sqrt(pow(vectOnNcsAxis[0],2.0) + pow(vectOnNcsAxis[1],2.0));
    float cosPart = vectOnNcsAxis[2];
    if( cosPart > 0) {
        psi0 = atan(sinPart / cosPart);
    } else if( cosPart < 0) {  
        psi0 = atan(sinPart / cosPart) + PI;
    } else if( cosPart == 0) {  
        psi0 = 0.5 * PI;
    }
    float phi0;
    if(vectOnNcsAxis[1] > 0 && vectOnNcsAxis[0] > 0) {
        phi0 = atan(vectOnNcsAxis[1]/vectOnNcsAxis[0]);
    } else if(vectOnNcsAxis[1] > 0 && vectOnNcsAxis[0] < 0) {  
        phi0 = atan(vectOnNcsAxis[1]/vectOnNcsAxis[0]) + PI;
    } else if(vectOnNcsAxis[1] < 0 && vectOnNcsAxis[0] < 0) {  
        phi0 = atan(vectOnNcsAxis[1]/vectOnNcsAxis[0]) + PI;
    } else if(vectOnNcsAxis[1] < 0 && vectOnNcsAxis[0] > 0) {
        phi0 = atan(vectOnNcsAxis[1]/vectOnNcsAxis[0]) + 2.0*PI;
    } else if(vectOnNcsAxis[1] == 0 && vectOnNcsAxis[0] >= 0) {
        phi0 = 0;
    } else if(vectOnNcsAxis[1] > 0 && vectOnNcsAxis[0] == 0) {    
        phi0 = 0.50 * PI;
    } else if(vectOnNcsAxis[1] == 0 && vectOnNcsAxis[0] < 0) {
        phi0 = PI;
    } else if(vectOnNcsAxis[1] < 0 && vectOnNcsAxis[0] == 0) {
        phi0 = 1.50 * PI;
    }
    float r0 = sqrt(pow(vectOnNcsAxis[0],2.0) + pow(vectOnNcsAxis[1],2.0) + pow(vectOnNcsAxis[2],2.0));
    float ux = vectOnNcsAxis[0] / r0;
    float uy = vectOnNcsAxis[1] / r0;
    float uz = vectOnNcsAxis[2] / r0;
    float rotaMatrix[3][3] = {};
    int icand = 0;
    int numPsi = 10;

    for(int i = 0; i <= numPsi; ++i){
        float psi = psi0 + float(i) * (ncsAxisDeviAngl/numPsi);
        while(psi > PI ) psi = 2.0 * PI - psi;
        float phi = phi0;
        float r = r0;
        float x = r * sin(psi) * cos (phi);
        float y = r * sin(psi) * sin (phi);
        float z = r * cos(psi);

        if(i == 0) {
            vectOnNcsAxisCand[icand][0] = x; //vectOnNcsAxis[0];
            vectOnNcsAxisCand[icand][1] = y; //vectOnNcsAxis[1];
            vectOnNcsAxisCand[icand][2] = z; //vectOnNcsAxis[2];
            pointAOnNcsAxisCand[icand][0] = pointAOnNcsAxis[0];
            pointAOnNcsAxisCand[icand][1] = pointAOnNcsAxis[1];
            pointAOnNcsAxisCand[icand][2] = pointAOnNcsAxis[2];
            icand++;
            continue;
        }

        for(int j = 0; j < (numNcsAxisCand-1)/numPsi; ++j){

            float thita = j * 2.0 * PI / float((numNcsAxisCand-1)/numPsi);
            rotaMatrix[0][0] = cos(thita) + ux * ux * (1.0 - cos(thita));
            rotaMatrix[0][1] = ux * uy * (1.0 - cos(thita)) - uz * sin(thita);
            rotaMatrix[0][2] = ux * uz * (1.0 - cos(thita)) + uy * sin(thita);

            rotaMatrix[1][0] = uy * ux * (1.0 - cos(thita)) + uz * sin(thita);
            rotaMatrix[1][1] = cos(thita) + uy * uy * (1.0 -cos(thita));
            rotaMatrix[1][2] = uy * uz * (1.0 - cos(thita)) - ux * sin(thita);

            rotaMatrix[2][0] = uz * ux * (1.0 - cos(thita)) - uy * sin(thita);
            rotaMatrix[2][1] = uz * uy * (1.0 - cos(thita)) + ux * sin(thita);
            rotaMatrix[2][2] = cos(thita) + uz * uz * (1.0 - cos(thita));

            vectOnNcsAxisCand[icand][0] = rotaMatrix[0][0] * x + rotaMatrix[0][1] * y + rotaMatrix[0][2] * z;
            vectOnNcsAxisCand[icand][1] = rotaMatrix[1][0] * x + rotaMatrix[1][1] * y + rotaMatrix[1][2] * z;
            vectOnNcsAxisCand[icand][2] = rotaMatrix[2][0] * x + rotaMatrix[2][1] * y + rotaMatrix[2][2] * z;

            pointAOnNcsAxisCand[icand][0] = pointAOnNcsAxis[0];
            pointAOnNcsAxisCand[icand][1] = pointAOnNcsAxis[1];
            pointAOnNcsAxisCand[icand][2] = pointAOnNcsAxis[2];

            icand++;
        }
    }

    return;
}

/******************************************************************************/

void outputNcsAxisCand(string fileNameString)
{
    ofstream outputFile(fileNameString.c_str());
    if(!outputFile.is_open()){
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE); 
    }
    outputFile << writeCryst1() << endl;

    for(int icand = 0; icand < numNcsAxisCand; ++icand){ 

        float occupancy = 1;
        float tempFact = icand;
        outputFile << "ATOM  " << setfill(' ') << setw(5) << icand 
            << ' ' << " C   ARG A   1    "
            << setfill(' ') << setw(8) << fixed << setprecision(3) << pointAOnNcsAxisCand[icand][0]
            << setfill(' ') << setw(8) << fixed << setprecision(3) << pointAOnNcsAxisCand[icand][1]
            << setfill(' ') << setw(8) << fixed << setprecision(3) << pointAOnNcsAxisCand[icand][2]
            << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
            << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
            << string(11, ' ') << 'C' << endl;
        outputFile << "ATOM  " << setfill(' ') << setw(5) << numNcsAxisCand + icand 
            << ' ' << " C   ARG A   1    "
            << setfill(' ') << setw(8) << fixed << setprecision(3) << pointAOnNcsAxisCand[icand][0] + vectOnNcsAxisCand[icand][0] * 100.0
            << setfill(' ') << setw(8) << fixed << setprecision(3) << pointAOnNcsAxisCand[icand][1] + vectOnNcsAxisCand[icand][1] * 100.0
            << setfill(' ') << setw(8) << fixed << setprecision(3) << pointAOnNcsAxisCand[icand][2] + vectOnNcsAxisCand[icand][2] * 100.0
            << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
            << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
            << string(11, ' ') << 'C' << endl;
    }

    for(int icand = 0; icand < numNcsAxisCand; ++icand){
        outputFile << "CONECT" << setfill(' ') << setw(5) << icand 
            << setfill(' ') << setw(5) << numNcsAxisCand + icand << endl;
    }

    outputFile.close();

    return;
}

/******************************************************************************/

void allocateNcsMaskCylinder()
{

    for(int ia = 0; ia <= NCSANUM - 1; ++ia){
        for(int ib = 0; ib <= NCSBNUM - 1; ++ib){
            for(int ic = 0; ic <= NCSCNUM - 1; ++ic){
                ncsMask[ia][ib][ic] = 0;
                ncsMaskDensTemp[ia][ib][ic] = 0;
                ncsMaskDens[ia][ib][ic] = 0;
            }
        }
    }

    float frac[3];
    convertOrthToFrac(pointAOnNcsAxis, frac);
    int ia0 = floor(frac[0] * ANUM);
    int ib0 = floor(frac[1] * BNUM);
    int ic0 = floor(frac[2] * CNUM);
    int cylinderRadius = maxPackCylinderRadiusUnitCell[ia0][ib0][ic0];
    int cylinderHeight = 2 * cylinderRadius; //must be an even number
    int sphereRadius = ceil( sqrt( cylinderRadius*cylinderRadius + cylinderHeight/2.0*cylinderHeight/2.0 ) );

    float cylinderRadiusReal = cylinderRadius * (CRYST_A / ANUM); //potential bug
    float cylinderHeightReal = 2.0 * cylinderRadiusReal;

    float pointBOnNcsAxis[3];
    for(int i = 0; i < 3; ++i){
        pointBOnNcsAxis[i] = pointAOnNcsAxis[i] + vectOnNcsAxis[i];
    }

    //the ncsMask is from 0~anum2(=anum1-anum0), 0~bnum2(=bnum1-bnum0), 0~cnum2(=cnum1-cnum0)
    anum0 = ia0;
    bnum0 = ib0;
    cnum0 = ic0;
    anum1 = anum0; bnum1 = bnum0; cnum1 = cnum0;

    for(int ia1 = ia0-sphereRadius; ia1 <= ia0+sphereRadius; ++ia1) {
        for(int ib1 = ib0-sphereRadius; ib1 <= ib0+sphereRadius; ++ib1) {
            for(int ic1 = ic0-sphereRadius; ic1 <= ic0+sphereRadius; ++ic1) {
                if( (ia1 - ia0)*(ia1 - ia0) + (ib1 - ib0)*(ib1 - ib0) + (ic1 - ic0)*(ic1 - ic0) <= sphereRadius*sphereRadius ) {
                    float frac1[3], orth1[3];
                    frac1[0] = (ia1 + 0.5) / float(ANUM);
                    frac1[1] = (ib1 + 0.5) / float(BNUM);
                    frac1[2] = (ic1 + 0.5) / float(CNUM);
                    convertFracToOrth(frac1, orth1);
                    float dist = computeDistFromPointToLine(orth1, pointAOnNcsAxis, pointBOnNcsAxis);
                    if(dist < cylinderRadiusReal){ //must within a dist about the axis
                        float dist1 = 0;
                        for(int i = 0; i < 3; ++i){
                            dist1 = dist1 + (orth1[i] - pointAOnNcsAxis[i]) * vectOnNcsAxis[i];
                        }
                        dist1 = abs(dist1);
                        if(dist1 <= cylinderHeightReal/2.0){ //must inside a cylinder
                            //allocate a large qubic region to hold the ncs mask
                            if(ia1 < anum0) anum0 = ia1;
                            if(ib1 < bnum0) bnum0 = ib1;
                            if(ic1 < cnum0) cnum0 = ic1;
                            if(ia1 > anum1) anum1 = ia1;
                            if(ib1 > bnum1) bnum1 = ib1;
                            if(ic1 > cnum1) cnum1 = ic1;
                        } //must inside a cylinder
                    } //must within a dist about the axis
                } //must inside a big sphere
            } //ic1
        } //ib1
    } //ia1

    //settle ncs mask into a little bit loose qubic region
    anum0 = anum0 - 3;
    anum1 = anum1 + 3;
    bnum0 = bnum0 - 3;
    bnum1 = bnum1 + 3;
    cnum0 = cnum0 - 3;
    cnum1 = cnum1 + 3;
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;

    if(anum2 > NCSANUM || bnum2 > NCSBNUM || cnum2 > NCSCNUM){
        cout << "The cubic region holding ncs mask is too small !" << endl;
        cout << anum2 << ',' << bnum2 << ',' << cnum2 << " should be less than " << NCSANUM << endl;
        exit(EXIT_FAILURE);
    }

    for(int ia1 = ia0-sphereRadius; ia1 <= ia0+sphereRadius; ++ia1) {
        for(int ib1 = ib0-sphereRadius; ib1 <= ib0+sphereRadius; ++ib1) {
            for(int ic1 = ic0-sphereRadius; ic1 <= ic0+sphereRadius; ++ic1) {
                if( (ia1 - ia0)*(ia1 - ia0) + (ib1 - ib0)*(ib1 - ib0) + (ic1 - ic0)*(ic1 - ic0) <= sphereRadius*sphereRadius ) {
                    float frac1[3], orth1[3];
                    frac1[0] = (ia1 + 0.5) / float(ANUM);
                    frac1[1] = (ib1 + 0.5) / float(BNUM);
                    frac1[2] = (ic1 + 0.5) / float(CNUM);
                    convertFracToOrth(frac1, orth1);
                    float dist = computeDistFromPointToLine(orth1, pointAOnNcsAxis, pointBOnNcsAxis);
                    if(dist < cylinderRadiusReal){ //must within a dist about the axis
                        float dist1 = 0;
                        for(int i = 0; i < 3; ++i){
                            dist1 = dist1 + (orth1[i] - pointAOnNcsAxis[i]) * vectOnNcsAxis[i];
                        }
                        dist1 = abs(dist1);
                        if(dist1 <= cylinderHeightReal/2.0){ //must inside a cylinder
                            int ia2 = ia1 - anum0;
                            int ib2 = ib1 - bnum0;
                            int ic2 = ic1 - cnum0;
                            ncsMask[ia2][ib2][ic2] = 1;
                        } //must inside a cylinder
                    } //must within a dist about the axis
                } //must inside a big sphere
            } //ic1
        } //ib1
    } //ia1

    return;
}

/******************************************************************************/

void applyNcsRotaOnPdbModel(string *atomStringArray, float thita, string fileNameString)
{
    ofstream outputFile(fileNameString.c_str());
    if(!outputFile.is_open()){
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE); 
    }

    float a0 = pointAOnNcsAxis[0];
    float b0 = pointAOnNcsAxis[1];
    float c0 = pointAOnNcsAxis[2];
    float u = vectOnNcsAxis[0];
    float v = vectOnNcsAxis[1];
    float w = vectOnNcsAxis[2];
    float L0 = u*u + v*v + w*w;
    float frac[3], orth[3], orth1[3], x, y, z;
    for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom){
        string lineString = atomStringArray[iatom];
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0){
            stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
            x = orth[0];
            y = orth[1];
            z = orth[2];
            orth1[0] = 1.0/L0*(  (a0*(v*v+w*w)-u*(b0*v+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
                +L0*x*cos(thita)+sqrt(L0)*(-c0*v+b0*w-w*y+v*z)*sin(thita)  );
            orth1[1] = 1.0/L0*(  (b0*(u*u+w*w)-v*(a0*u+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
                +L0*y*cos(thita)+sqrt(L0)*( c0*u-a0*w+w*x-u*z)*sin(thita)  );
            orth1[2] = 1.0/L0*(  (c0*(u*u+v*v)-w*(a0*u+b0*v-u*x-v*y-w*z))*(1.0-cos(thita))
                +L0*z*cos(thita)+sqrt(L0)*(-b0*u+a0*v-v*x+u*y)*sin(thita)  );
            outputFile << lineString.substr(0, 30) 
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth1[0]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth1[1]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth1[2]
                << lineString.substr(54, 26) << endl;
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void applyNcsTranOnPdbModel(string *atomStringArray, float dist, string fileNameString)
{
    ofstream outputFile(fileNameString.c_str());
    if(!outputFile.is_open()){
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE); 
    }
    float u = vectOnNcsAxis[0];
    float v = vectOnNcsAxis[1];
    float w = vectOnNcsAxis[2];
    float L0 = u*u + v*v + w*w;
    float frac[3], orth[3], orth1[3];
    for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom){
        string lineString = atomStringArray[iatom];
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0){
            stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
            orth1[0] = orth[0] + dist * u/sqrt(L0);
            orth1[1] = orth[1] + dist * v/sqrt(L0);
            orth1[2] = orth[2] + dist * w/sqrt(L0);
            outputFile << lineString.substr(0, 30) 
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth1[0]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth1[1]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth1[2]
                << lineString.substr(54, 26) << endl;
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void applyNcsMatrixOnPdbModel(string *atomStringArray, string fileNameString)
{
    ofstream outputFile(fileNameString.c_str());
    if(!outputFile.is_open()){
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE); 
    }

    float rota[3][3], tran[3], frac[3], orth[3], orth1[3];
    for (int incs = 1; incs <= numNcsOper - 1; ++incs){
        for (int i = 0; i <= 2; ++i){
            for (int j = 0; j <= 2; ++j){
                rota[i][j] = ncsOper[incs][i][j]; //load rotation matrix
            } 
            tran[i] = ncsOper[incs][i][3]; //load rotation matrix
        }
        for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom){
            string lineString = atomStringArray[iatom];
            if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0){
                stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
                for( int i = 0; i <= 2 ; ++i ){
                    orth1[i] = 0;
                    for( int j = 0; j <= 2 ; ++j ){
                        orth1[i] = orth1[i] + rota[i][j] * orth[j]; //rotation
                    }
                    orth1[i] = orth1[i] + tran[i]; //translation
                }
                outputFile << lineString.substr(0, 30) 
                    << setfill(' ') << setw(8) << setprecision(3) << fixed << orth1[0]
                    << setfill(' ') << setw(8) << setprecision(3) << fixed << orth1[1]
                    << setfill(' ') << setw(8) << setprecision(3) << fixed << orth1[2]
                    << lineString.substr(54, 26) << endl;
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void applyNcsRotaOnPoint(float point[3], float thita, float newPoint[3])
{
    float a0 = pointAOnNcsAxis[0];
    float b0 = pointAOnNcsAxis[1];
    float c0 = pointAOnNcsAxis[2];
    float u = vectOnNcsAxis[0];
    float v = vectOnNcsAxis[1];
    float w = vectOnNcsAxis[2];
    float L0 = u*u + v*v + w*w;
    float x = point[0], y = point[1], z = point[2];
    newPoint[0] = 1.0/L0*(  (a0*(v*v+w*w)-u*(b0*v+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
        +L0*x*cos(thita)+sqrt(L0)*(-c0*v+b0*w-w*y+v*z)*sin(thita)  );
    newPoint[1] = 1.0/L0*(  (b0*(u*u+w*w)-v*(a0*u+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
        +L0*y*cos(thita)+sqrt(L0)*( c0*u-a0*w+w*x-u*z)*sin(thita)  );
    newPoint[2] = 1.0/L0*(  (c0*(u*u+v*v)-w*(a0*u+b0*v-u*x-v*y-w*z))*(1.0-cos(thita))
        +L0*z*cos(thita)+sqrt(L0)*(-b0*u+a0*v-v*x+u*y)*sin(thita)  );
    return;
}

/******************************************************************************/

void applyRotaAboutAnyAxis(float pointOnAxis[3], float vectOnAxis[3], float thita, float point[3], float newPoint[3])
{
    float a0 = pointOnAxis[0];
    float b0 = pointOnAxis[1];
    float c0 = pointOnAxis[2];
    float u = vectOnAxis[0];
    float v = vectOnAxis[1];
    float w = vectOnAxis[2];
    float L0 = u*u + v*v + w*w;
    float x = point[0], y = point[1], z = point[2];
    newPoint[0] = 1.0/L0*(  (a0*(v*v+w*w)-u*(b0*v+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
        +L0*x*cos(thita)+sqrt(L0)*(-c0*v+b0*w-w*y+v*z)*sin(thita)  );
    newPoint[1] = 1.0/L0*(  (b0*(u*u+w*w)-v*(a0*u+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
        +L0*y*cos(thita)+sqrt(L0)*( c0*u-a0*w+w*x-u*z)*sin(thita)  );
    newPoint[2] = 1.0/L0*(  (c0*(u*u+v*v)-w*(a0*u+b0*v-u*x-v*y-w*z))*(1.0-cos(thita))
        +L0*z*cos(thita)+sqrt(L0)*(-b0*u+a0*v-v*x+u*y)*sin(thita)  );
    return;
}

/******************************************************************************/

void applyNcsTranOnPoint(float point[3], float dist, float newPoint[3])
{
    float u = vectOnNcsAxis[0];
    float v = vectOnNcsAxis[1];
    float w = vectOnNcsAxis[2];
    float L0 = u*u + v*v + w*w;
    newPoint[0] = point[0] + dist * u/sqrt(L0);
    newPoint[1] = point[1] + dist * v/sqrt(L0);
    newPoint[2] = point[2] + dist * w/sqrt(L0);
    return;
}

/******************************************************************************/

float computeMinima(float x1, float x2, float x3)
{
    return min(min(x1, x2), min(x2, x3));
}

/******************************************************************************/

float distOrth(float *orth1, float *orth2)
{
    float dist = 0;
    for(int i = 0; i <= 2; ++i){
        dist = dist + pow((orth2[i] - orth1[i]), 2.0);
    }
    dist = sqrt(dist);
    return dist;
}

/******************************************************************************/

float computeAngleTwoVect(float vect1[3], float vect2[3])
{
    float angle;
    angle = acos( (vect1[0] * vect2[0] + vect1[1] * vect2[1] + vect1[2] * vect2[2]) /
                  sqrt(vect1[0] * vect1[0] + vect1[1] * vect1[1] + vect1[2] * vect1[2]) /
                  sqrt(vect2[0] * vect2[0] + vect2[1] * vect2[1] + vect2[2] * vect2[2]) );
    while(angle > 2.0 * PI) angle -= 2.0 * PI;
    while(angle < 0) angle += 2.0 * PI;
    return angle;
}

/******************************************************************************/

float computeDistFromPointToLine(float *point, float *pointAOnLine, float *pointBOnLine)
{
    float x = point[0], y = point[1], z = point[2];
    float x1 = pointAOnLine[0], y1 = pointAOnLine[1], z1 = pointAOnLine[2];
    float x2 = pointBOnLine[0], y2 = pointBOnLine[1], z2 = pointBOnLine[2];
    float dist;
    dist = sqrt( pow((y-y1)*(z1-z2)-(z-z1)*(y1-y2), 2.0) + pow((z-z1)*(x1-x2)-(x-x1)*(z1-z2), 2.0) \
                + pow((x-x1)*(y1-y2)-(y-y1)*(x1-x2), 2.0) ) / sqrt( pow(x1-x2, 2.0) + pow(y1-y2, 2.0) + pow(z1-z2, 2.0) );

    return dist;
}

/******************************************************************************/

void findPointAwayFromLine(float dist, float pointAOnLine[3], float vectOnLine[3], float point[3])
{
    float xA = pointAOnLine[0], yA = pointAOnLine[1], zA = pointAOnLine[2];
    float u = vectOnLine[0];
    float v = vectOnLine[1];
    float w = vectOnLine[2];
    float L0 = u*u + v*v + w*w;
    //unit vector
    u /= sqrt(L0);
    v /= sqrt(L0);
    w /= sqrt(L0);
    //find a point P(xP,yP,zP) which is located on the plane perpendicular to the 
    //straight line AB at point A(xA,yA,zA). The distance from P to A is dist.
    //temperary coefficients c1,c2,c3,c4,c5,c6,c7
    float c1 = ( -u*v*xA+(u*u+w*w)*yA-v*w*zA ) / ( (v*v+w*w)*xA-u*v*yA-u*w*zA );
    float c2 = ( (u*xA+v*yA+w*zA)*(v*xA-u*yA) ) / ( (v*v+w*w)*xA-u*v*yA-u*w*zA );
  
    //see note zP=c3*xP+c4
    float c3 = ( c1*w*xA-w*yA+(v-c1*u)*zA ) / ( v*xA-u*yA );
    float c4 = ( c2*(w*xA-u*zA) ) / ( v*xA-u*yA );
  
    float c5 = 1.0 + c1*c1 + c3*c3;
    float c6 = -2.0*xA + 2.0*c1*(c2-yA) + 2.0*c3*(c4-zA);
    float c7 = xA*xA + pow((c2-yA), 2.0) + pow((c4-zA), 2.0) - pow(dist, 2.0);
    float xP;
    if(c6*c6-4.0*c5*c7 < 0){
        cout << "c6*c6-4.0*c5*c7 < 0" << endl;
        xP = ( -c6 + 0 ) / ( 2.0*c5 );
    } else {
        xP = ( -c6 + sqrt(c6*c6-4.0*c5*c7) ) / ( 2.0*c5 );
    }
    float yP = c1*xP + c2;
    float zP = c3*xP + c4;
    point[0] = xP;
    point[1] = yP;
    point[2] = zP;
    return;
}

/******************************************************************************/

void rotaPointAroundLine(float point[3], float anglPhi, float pointAOnLine[3], float vectOnLine[3], float newPoint[3])
{
    //this function actually is the same as the function applyNcsRotaOnPoint

    float xA = pointAOnLine[0], yA = pointAOnLine[1], zA = pointAOnLine[2];
    float u = vectOnLine[0];
    float v = vectOnLine[1];
    float w = vectOnLine[2];
    float L0 = u*u + v*v + w*w;
    //unit vector
    u /= sqrt(L0);
    v /= sqrt(L0);
    w /= sqrt(L0);
    float idenMatrix[3][3] = {  { 1.0, 0.0, 0.0},
                                { 0.0, 1.0, 0.0},
                                { 0.0, 0.0, 1.0}  };

    float rotaMatrix[3][3] = {  {  0.0, -w,  v },
                                {  w,  0.0, -u },
                                { -v,  u,  0.0 }  };

    //find the coordinates of point P' which is written as Pp(xPp,yPp,zPp)
    float xP = point[0], yP = point[1], zP = point[2];
    float vectAP[3] = {xP - xA, yP - yA, zP - zA};

    //rotate AP around straight line AB to get APp
    float matrix2[3][3] = {};
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            matrix2[i][j] = 0;
            for(int k = 0; k < 3; ++k){
                matrix2[i][j] += rotaMatrix[i][k] * rotaMatrix[k][j];
            }
        }
    }
    float matrix3[3][3] = {};
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            matrix3[i][j] = idenMatrix[i][j] + sin(anglPhi)*rotaMatrix[i][j] + (1.0-cos(anglPhi))*matrix2[i][j]; 
        }
    }
    float vectAPp[3] = {};
    for(int i = 0; i < 3; ++i){
        vectAPp[i] = 0;
        for(int j = 0; j < 3; ++j){
            vectAPp[i] += matrix3[i][j] * vectAP[j];
        }
    }
    float xPp = vectAPp[0] + xA;
    float yPp = vectAPp[1] + yA;
    float zPp = vectAPp[2] + zA;
    newPoint[0] = xPp;
    newPoint[1] = yPp;
    newPoint[2] = zPp;
    return;
}

/******************************************************************************/

void inputNcsAxisCand(string fileNameString)
{
    ifstream inputFile( fileNameString.c_str() );
    if (!inputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    int icand;
    float point[3], vect[3];
    while (inputFile >> icand >> point[0] >> point[1] >> point[2] >> vect[0] >> vect[1] >> vect[2]) {
        pointAOnNcsAxisCand[icand][0] = point[0];
        pointAOnNcsAxisCand[icand][1] = point[1];
        pointAOnNcsAxisCand[icand][2] = point[2];
        vectOnNcsAxisCand[icand][0] = vect[0];
        vectOnNcsAxisCand[icand][1] = vect[1];
        vectOnNcsAxisCand[icand][2] = vect[2];
    }
    //initilization
    float maxAngle = 10;
    int numPointForEachAngle = (numNcsAxisCand - 1) / maxAngle;
    int random = rand() % numPointForEachAngle + (numNcsAxisCand-numPointForEachAngle) ;
    for(int i = 0; i < 3; ++i){
        pointAOnNcsAxis[i] = pointAOnNcsAxisCand[random][i];
        vectOnNcsAxis[i] = vectOnNcsAxisCand[random][i];
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void inputRefeHist(string fileNameString)
{
    //read in the boundary values of reference protein hist, file name is "PDB_CODE"_hist.txt
    ifstream inputFile( fileNameString.c_str() );
    if (!inputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    int ilargeBin;
    float density;
    while (inputFile >> ilargeBin >> density) {
        stdBoundLargeBin[ilargeBin] = density;
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void inputSigmWeigObsAndSigmWeigAvg(string fileNameString)
{
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    string lineString;
    int iter;
    while(getline(inputFile, lineString)){
        stringstream(lineString.substr(0,8)) >> iter;
        stringstream(lineString.substr(9,8)) >> sigmWeigObsArray[iter];
        stringstream(lineString.substr(18,8)) >> sigmWeigAvgArray[iter];
    }
    inputFile.close();
    return;
}

/******************************************************************************/

int findMiddlePointNcsAtom(string *atomStringArray, float middlePoint[][3])
{
    int flag[numNcsOper] = {};

    for(int imiddle = 0; imiddle < 100000; ++imiddle) {
        for(int i = 0; i < 3; ++i) {
            middlePoint[imiddle][i] = 0;
        }
    }

    int numMiddlePoint, imiddle = 0;

    //outerloop over all atoms on chain A
    for(int iatom0 = 0; iatom0 < numAtomAsu; ++iatom0){

        string lineString0 = atomStringArray[iatom0];
        if(lineString0.compare(0, 4, "ATOM") != 0 || lineString0.compare(21, 1, "A") != 0) continue;

        for(int ichain = 0; ichain < numNcsOper; ++ichain){
            flag[ichain] = 0;
        }

        //innerloop over all atoms on chain B, C, D, E, ...
        for(int iatom1 = 0; iatom1 < numAtomAsu; ++iatom1){

            string lineString1 = atomStringArray[iatom1];
            if(lineString1.compare(0, 4, "ATOM") != 0 || lineString1.compare(12, 8, lineString0, 12, 8) != 0
                || lineString1.compare(23, 3, lineString0, 23, 3) != 0) continue;

            int ichain;
            if(lineString1.compare(21, 1, "A") == 0) {
                ichain = 0;
            } else if(lineString1.compare(21, 1, "B") == 0) {
                ichain = 1;
            } else if(lineString1.compare(21, 1, "C") == 0) {
                ichain = 2;
            } else if(lineString1.compare(21, 1, "D") == 0) {
                ichain = 3;
            } else if(lineString1.compare(21, 1, "E") == 0) {
                ichain = 4;
            } else {
                continue;
            }
            flag[ichain] = 1;

            float orth[3];
            stringstream(lineString1.substr(30, 24)) >> orth[0] >> orth[1] >> orth[2];
            middlePoint[imiddle][0] += orth[0]; 
            middlePoint[imiddle][1] += orth[1];
            middlePoint[imiddle][2] += orth[2];
        }

        //check that all ncs related copies have been found
        for(int ichain = 0; ichain < numNcsOper; ++ichain){
            if(flag[ichain] == 0){ //not find all ncs symmetry related copies
                middlePoint[imiddle][0] = 0;
                middlePoint[imiddle][1] = 0;
                middlePoint[imiddle][2] = 0;
                flag[0] = 0;
                break;
            }
        }
        if(flag[0] == 0) continue;
        middlePoint[imiddle][0] /= float(numNcsOper);
        middlePoint[imiddle][1] /= float(numNcsOper);
        middlePoint[imiddle][2] /= float(numNcsOper);
        imiddle ++;
        if(imiddle >= 100000){
            cout << "Found more than 100000 middles points! Increase the size of middlePoint[][]." << endl;
            exit(EXIT_FAILURE);
        }
    }
    numMiddlePoint = imiddle;

    return numMiddlePoint;
}

/******************************************************************************/

int excludeOutliers(int numMiddlePoint, float middlePoint[][3])
{
    float centroid[3] = {}, point[2][3] = {};

    //find the centroid of all middlePoint on ncs axis
    for(int i = 0; i < 3; ++i){
        centroid[i] = 0;
        for(int imiddle = 0; imiddle < numMiddlePoint; ++imiddle){
            centroid[i] += middlePoint[imiddle][i];
        }
        centroid[i] /= float(numMiddlePoint);
    }

    //find two subcentroid of middlePoint on each side of the centroid
    int ipointA = 0, ipointB = 0;
    for(int imiddle = 0; imiddle < numMiddlePoint; ++imiddle){
        float innerProd = 0;
        for(int i = 0; i < 3; ++i){
            innerProd += (middlePoint[imiddle][i] - centroid[i]) * (middlePoint[0][i] - centroid[i]);
        }
        if(innerProd >= 0){
            for(int i = 0; i < 3; ++i){
                point[0][i] += middlePoint[imiddle][i];
            }
            ipointA++;
        } else {
            for(int i = 0; i < 3; ++i){
                point[1][i] += middlePoint[imiddle][i];
            }
            ipointB++;
        }
    }
    for(int i = 0; i < 3; ++i){
        point[0][i] /= float(ipointA);
        point[1][i] /= float(ipointB);
    }

    //compute the distance from each middle point to the NCS axis to find out and exclude outliers
    float pointAOnLine[3] = {}, pointBOnLine[3] = {};
    for(int i = 0; i < 3; ++i){
        pointAOnLine[i] = point[0][i];
        pointBOnLine[i] = point[1][i];
    }
    int imiddle = 0;
    while(imiddle < numMiddlePoint){
        float pointP[3] = {};
        for(int i = 0; i < 3; ++i){
            pointP[i] = middlePoint[imiddle][i];
        }
        float dist = computeDistFromPointToLine(pointP, pointAOnLine, pointBOnLine);
        if(dist > resoCutoff / 2.0) {
            for(int i = 0; i < 3; ++i){
                middlePoint[imiddle][i]= middlePoint[numMiddlePoint-1][i];
            }
            numMiddlePoint--;
            continue;
        }
        imiddle++;
    }
    return numMiddlePoint;
}

/******************************************************************************/

void outputMiddlePointNcsAtom(int numMiddlePoint, float middlePoint[][3], string fileNameString)
{
    ofstream outputFile(fileNameString.c_str());
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    outputFile << writeCryst1() << endl;
    for(int imiddle = 0; imiddle < numMiddlePoint; ++imiddle){
        float occupancy = 1.0, tempFact = 0;
        outputFile << "ATOM  " << setfill(' ') << setw(5) << imiddle
            << ' ' << " C   ARG A   1    "
            << setfill(' ') << setw(8) << fixed << setprecision(3) << middlePoint[imiddle][0]
            << setfill(' ') << setw(8) << fixed << setprecision(3) << middlePoint[imiddle][1]
            << setfill(' ') << setw(8) << fixed << setprecision(3) << middlePoint[imiddle][2]
            << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
            << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
            << string(11, ' ') << 'C' << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputTwoPointOnNcsAxis(int numMiddlePoint, float middlePoint[][3], string fileNameString)
{
    ofstream outputFile(fileNameString.c_str());
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    outputFile << writeCryst1() << endl;

    float centroid[3] = {}, point[2][3] = {};

    //find the centroid of all middlePoint on ncs axis
    for(int i = 0; i < 3; ++i){
        centroid[i] = 0;
        for(int imiddle = 0; imiddle < numMiddlePoint; ++imiddle){
            centroid[i] += middlePoint[imiddle][i];
        }
        centroid[i] /= float(numMiddlePoint);
    }

    //find two subcentroid of middlePoint on each side of the centroid
    int ipointA = 0, ipointB = 0;
    for(int imiddle = 0; imiddle < numMiddlePoint; ++imiddle){
        float innerProd = 0;
        for(int i = 0; i < 3; ++i){
            innerProd += (middlePoint[imiddle][i] - centroid[i]) * (middlePoint[0][i] - centroid[i]);
        }
        if(innerProd >= 0){
            for(int i = 0; i < 3; ++i){
                point[0][i] += middlePoint[imiddle][i];
            }
            ipointA++;
        } else {
            for(int i = 0; i < 3; ++i){
                point[1][i] += middlePoint[imiddle][i];
            }
            ipointB++;
        }
    }
    for(int i = 0; i < 3; ++i){
        point[0][i] /= float(ipointA);
        point[1][i] /= float(ipointB);
    }

//    //make the ncs axis longer
//    for(int i = 0; i < 3; ++i){
//        float temp = point[0][i];
//        temp += 6.0 * (point[0][i] - point[1][i]);
//        point[1][i] += 3.0 * (point[1][i] - point[0][i]);
//        point[0][i] = temp;
//    }

    for(int iatom = 0; iatom < 2; ++iatom){
        float occupancy = 1.0, tempFact = 0;
        outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom
            << ' ' << " C   ARG A   1    "
            << setfill(' ') << setw(8) << fixed << setprecision(3) << point[iatom][0]
            << setfill(' ') << setw(8) << fixed << setprecision(3) << point[iatom][1]
            << setfill(' ') << setw(8) << fixed << setprecision(3) << point[iatom][2]
            << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
            << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
            << string(11, ' ') << 'C' << endl;
    }
    outputFile << "CONECT    0    1" << endl;
    outputFile.close();
    return;
}

/******************************************************************************/

void inputMaxPackSphereRadiusUnitCell(string fileNameString)
{
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                inputFile >> maxPackSphereRadiusUnitCell[ia][ib][ic];
            }
        }
    }
    inputFile.close();
    return;
}

/******************************************************************************/


void inputMaxPackCylinderRadiusUnitCell(string fileNameString)
{
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                inputFile >> maxPackCylinderRadiusUnitCell[ia][ib][ic];
            }
        }
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void outputMaxPackSphereRadiusUnitCell(string fileNameString, float cutoff)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    } else {
        outputFile << writeCryst1() << endl;
    }

    int iatom = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                if(maxPackSphereRadiusUnitCell[ia][ib][ic] > cutoff){
                    float frac[3], orth[3], occupancy, tempFact;
                    iatom++;
                    if(iatom >= 9999) iatom = 0; //too many atoms
                    frac[0] = (ia + 0.5) / float(ANUM);
                    frac[1] = (ib + 0.5) / float(BNUM);
                    frac[2] = (ic + 0.5) / float(CNUM);
                    convertFracToOrth(frac, orth);
                    occupancy = 1.0;
                    tempFact = maxPackSphereRadiusUnitCell[ia][ib][ic];
                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                        << ' ' << " C   ARG A   1    "
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                        << string(11, ' ') << 'C' << endl;
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/


void outputMaxPackCylinderRadiusUnitCell(string fileNameString, float cutoff)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    } else {
        outputFile << writeCryst1() << endl;
    }

    int iatom = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                if(maxPackCylinderRadiusUnitCell[ia][ib][ic] > cutoff){
                    float frac[3], orth[3], occupancy, tempFact;
                    iatom++;
                    if(iatom >= 9999) iatom = 0; //too many atoms
                    frac[0] = (ia + 0.5) / float(ANUM);
                    frac[1] = (ib + 0.5) / float(BNUM);
                    frac[2] = (ic + 0.5) / float(CNUM);
                    convertFracToOrth(frac, orth);
                    occupancy = 1.0;
                    tempFact = maxPackCylinderRadiusUnitCell[ia][ib][ic];
                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                        << ' ' << " C   ARG A   1    "
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                        << string(11, ' ') << 'C' << endl;
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void makeCylinderNcsMask()
{
    //use the maximum cylinder radius
    int cylinderRadius = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                realConvMaskUnitCell[ia][ib][ic] = 0;
                if(maxPackCylinderRadiusUnitCell[ia][ib][ic] > cylinderRadius){
                    cylinderRadius = maxPackCylinderRadiusUnitCell[ia][ib][ic];
                }
            }
        }
    }

    int igrid = 0;
    float cylinderRadiusReal = cylinderRadius * (CRYST_A / ANUM); //potential bug
    int cylinderHeight = 2 * cylinderRadius; //must be an even number
    float cylinderHeightReal = 2.0 * cylinderRadiusReal;
    int sphereRadius = ceil( sqrt( cylinderRadius*cylinderRadius + cylinderHeight/2.0*cylinderHeight/2.0 ) );
    for(int ia1 = -sphereRadius; ia1 <= sphereRadius; ++ia1) {
        for(int ib1 = -sphereRadius; ib1 <= +sphereRadius; ++ib1) {
            for(int ic1 = -sphereRadius; ic1 <= +sphereRadius; ++ic1) {
                //must inside a big sphere
                if( (ia1 - 0)*(ia1 - 0) + (ib1 - 0)*(ib1 - 0) + (ic1 - 0)*(ic1 - 0) <= sphereRadius*sphereRadius ) {
                    float frac1[3], orth1[3];
                    frac1[0] = (ia1 + 0.5) / float(ANUM);
                    frac1[1] = (ib1 + 0.5) / float(BNUM);
                    frac1[2] = (ic1 + 0.5) / float(CNUM);
                    convertFracToOrth(frac1, orth1);
                    float pointBOnLine[3], pointOnNcsAxis[3] = {0, 0, 0};
                    for(int i = 0; i < 3; ++i){
                        pointBOnLine[i] = pointOnNcsAxis[i] + vectOnNcsAxis[i];
                    }
                    float dist = computeDistFromPointToLine(orth1, pointOnNcsAxis, pointBOnLine);
                    if(dist < cylinderRadiusReal){ //must within a dist about the axis
                        float dist1 = 0;
                        for(int i = 0; i < 3; ++i){
                            dist1 = dist1 + (orth1[i] - pointOnNcsAxis[i]) * vectOnNcsAxis[i];
                        }
                        dist1 = abs(dist1);
                        if(dist1 <= cylinderHeightReal/2.0){ //must inside a cylinder

                            for(int incs = 0; incs < numNcsOper; ++incs){
                                float thita = 2.0 * PI * float(incs) / float(numNcsOper);
                                float newPoint[3] = {};
                                applyNcsRotaOnPoint(orth1, thita, newPoint);
                                float frac2[3] = {};
                                convertOrthToFrac(newPoint, frac2);
                                int ia2 = floor(frac2[0] * ANUM);
                                int ib2 = floor(frac2[1] * BNUM);
                                int ic2 = floor(frac2[2] * CNUM);
                                cylinderDens[igrid][incs][0] = ia2;
                                cylinderDens[igrid][incs][1] = ib2;
                                cylinderDens[igrid][incs][2] = ic2;
                                igrid++;
                            }

                            int ia2 = ia1;
                            int ib2 = ib1;
                            int ic2 = ic1;
                            while(ia2 < 0) ia2 += ANUM;
                            while(ib2 < 0) ib2 += BNUM;
                            while(ic2 < 0) ic2 += CNUM;
                            while(ia2 >= ANUM) ia2 -= ANUM;
                            while(ib2 >= BNUM) ib2 -= BNUM;
                            while(ic2 >= CNUM) ic2 -= CNUM;
                            realConvMaskUnitCell[ia2][ib2][ic2] = 1.0;
                        } //if(dist1 <= cylinderHeightReal/2.0)
                    } //if(dist < cylinderRadiusReal)
                } //if(<= sphereRadius*sphereRadius)
            } //ic1
        } //ib1
    } //ia1
    numGridPointInCylinder = igrid;

    //prepare for fft
    setHalfGridPhase(); //set half grid phase shift for fft, initialize halfGridPosi and halfGridNega

    //apply backward fft
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                realSpaceArray[ia][ib][ic] = realConvMaskUnitCell[ia][ib][ic];
            }
        }
    }

    //Intel dfti backward fft from density to structure factor
    status = DftiCreateDescriptor( &myDescriptorHandle, DFTI_SINGLE,DFTI_COMPLEX, 3, lengthEachDimension);
    status = DftiSetValue( myDescriptorHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE );
    status = DftiSetValue( myDescriptorHandle, DFTI_BACKWARD_SCALE, dftiBacScaleFactor );
    status = DftiCommitDescriptor( myDescriptorHandle);
    status = DftiComputeBackward( myDescriptorHandle, realSpaceArray, reciprocalSpaceArray);
    status = DftiFreeDescriptor(&myDescriptorHandle);

    for(int ih = 0; ih <= HNUM - 1; ih++ ) {
        for(int ik = 0; ik <= KNUM - 1; ik++ ) {
            for(int il = 0; il <= LNUM - 1; il++ ) {
                reciprocalSpaceArray[ih][ik][il] *= halfGridPosi[ih][ik][il];
                reciprocalConvMaskUnitCell[ih][ik][il] = reciprocalSpaceArray[ih][ik][il];
            }
        }
    }

    return;
}

/******************************************************************************/

void makeSphereNcsMask()
{
    //use the maximum sphere radius
    int sphereRadius = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                realConvMaskUnitCell[ia][ib][ic] = 0;
                if(maxPackSphereRadiusUnitCell[ia][ib][ic] > sphereRadius){
                    sphereRadius = maxPackSphereRadiusUnitCell[ia][ib][ic];
                }
            }
        }
    }

sphereRadius=35;

    int igrid = 0;
    float sphereRadiusReal = sphereRadius * (CRYST_A / ANUM); //potential bug
    for(int ia1 = -sphereRadius; ia1 <= sphereRadius; ++ia1) {
        for(int ib1 = -sphereRadius; ib1 <= +sphereRadius; ++ib1) {
            for(int ic1 = -sphereRadius; ic1 <= +sphereRadius; ++ic1) {
                if( (ia1 - 0)*(ia1 - 0) + (ib1 - 0)*(ib1 - 0) + (ic1 - 0)*(ic1 - 0) <= sphereRadius*sphereRadius ) {
                    int ia2 = ia1;
                    int ib2 = ib1;
                    int ic2 = ic1;
                    while(ia2 < 0) ia2 += ANUM;
                    while(ib2 < 0) ib2 += BNUM;
                    while(ic2 < 0) ic2 += CNUM;
                    while(ia2 >= ANUM) ia2 -= ANUM;
                    while(ib2 >= BNUM) ib2 -= BNUM;
                    while(ic2 >= CNUM) ic2 -= CNUM;
                    realConvMaskUnitCell[ia2][ib2][ic2] = 1.0;
                } //if(<= sphereRadius*sphereRadius)
            } //ic1
        } //ib1
    } //ia1

    //prepare for fft
    setHalfGridPhase(); //set half grid phase shift for fft, initialize halfGridPosi and halfGridNega

    //apply backward fft
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                realSpaceArray[ia][ib][ic] = realConvMaskUnitCell[ia][ib][ic];
            }
        }
    }

    //Intel dfti backward fft from density to structure factor
    status = DftiCreateDescriptor( &myDescriptorHandle, DFTI_SINGLE,DFTI_COMPLEX, 3, lengthEachDimension);
    status = DftiSetValue( myDescriptorHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE );
    status = DftiSetValue( myDescriptorHandle, DFTI_BACKWARD_SCALE, dftiBacScaleFactor );
    status = DftiCommitDescriptor( myDescriptorHandle);
    status = DftiComputeBackward( myDescriptorHandle, realSpaceArray, reciprocalSpaceArray);
    status = DftiFreeDescriptor(&myDescriptorHandle);

    for(int ih = 0; ih <= HNUM - 1; ih++ ) {
        for(int ik = 0; ik <= KNUM - 1; ik++ ) {
            for(int il = 0; il <= LNUM - 1; il++ ) {
                reciprocalSpaceArray[ih][ik][il] *= halfGridPosi[ih][ik][il];
                reciprocalConvMaskUnitCell[ih][ik][il] = reciprocalSpaceArray[ih][ik][il];
            }
        }
    }

    return;
}

/******************************************************************************/

void findCentMassCand(float distAwayFromNcsAxis)
{
    int icand = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
//                if(maxPackSphereRadiusUnitCell[ia][ib][ic] < 5) continue; //too small maxPackSphere, 5 is adjustable
                float frac[3], orth[3];
                frac[0] = (ia + 0.5) / float(ANUM);
                frac[1] = (ib + 0.5) / float(BNUM);
                frac[2] = (ic + 0.5) / float(CNUM);
                convertFracToOrth(frac, orth);
                float x = orth[0], y = orth[1], z = orth[2];
                float x1 = pointAOnNcsAxis[0], y1 = pointAOnNcsAxis[1], z1 = pointAOnNcsAxis[2];
                float x2 = pointAOnNcsAxis[0]+vectOnNcsAxis[0], y2 = pointAOnNcsAxis[1]+vectOnNcsAxis[1], z2 = pointAOnNcsAxis[2]+vectOnNcsAxis[2];
                //compute distance from point to the NCS axis
                float dist = sqrt( pow((y-y1)*(z1-z2)-(z-z1)*(y1-y2), 2.0) + pow((z-z1)*(x1-x2)-(x-x1)*(z1-z2), 2.0) \
                                  + pow((x-x1)*(y1-y2)-(y-y1)*(x1-x2), 2.0) ) / sqrt( pow(x1-x2, 2.0) + pow(y1-y2, 2.0) + pow(z1-z2, 2.0) );
                bool flag = false;
                if(dist < distAwayFromNcsAxis){
                    for(int i = 0; i < icand; ++i){ //keep a distance from current centMassCand
                        float frac1[3],orth1[3];
                        int ia1 = centMassCand[i][0];
                        int ib1 = centMassCand[i][1];
                        int ic1 = centMassCand[i][2];
                        frac1[0] = (ia1 + 0.5) / float(ANUM);
                        frac1[1] = (ib1 + 0.5) / float(BNUM);
                        frac1[2] = (ic1 + 0.5) / float(CNUM);
                        convertFracToOrth(frac1, orth1);
                        if(distOrth(orth, orth1) < resoCutoff) { //the threshold can be adjusted
                            flag = true;
                            break;
                        }
                    }
                    if(flag) continue;
                    centMassCand[icand][0] = ia;
                    centMassCand[icand][1] = ib;
                    centMassCand[icand][2] = ic;
                    centMassCand[icand][3] = maxPackSphereRadiusUnitCell[ia][ib][ic];
                    centMassCand[icand][4] = 0;
                    icand++;
                    if(icand >= ANUM*BNUM*CNUM){
                        cout << "Too many candidates for center of mass!" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }
    numCentMassCand = icand;
    //cout << numCentMassCand << endl;
    return;
}

/******************************************************************************/

void outputCentMassCand(string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    outputFile << writeCryst1() << endl;
    for(int i = 0; i <= numCentMassCand - 1; ++i){
        int ia = centMassCand[i][0];
        int ib = centMassCand[i][1];
        int ic = centMassCand[i][2];
        float frac[3], orth[3];
        frac[0] = (ia + 0.5) / float(ANUM);
        frac[1] = (ib + 0.5) / float(BNUM);
        frac[2] = (ic + 0.5) / float(CNUM);
        convertFracToOrth(frac, orth);
        float tempFact = centMassCand[i][3];
        float occupancy = 1.0;
           outputFile << "ATOM" << setfill(' ') << setw(7) << i 
            << ' ' << " C   ARG A   1    "
            << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
            << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
            << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
            << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
            << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
            << string(11, ' ') << 'C' << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void updateCentMass(float densCutoff)
{
    float frac[3], orth[3];

    for(int i = 0; i <= numCentMassCand - 1; ++i){
        int ia0 = centMassCand[i][0];
        int ib0 = centMassCand[i][1];
        int ic0 = centMassCand[i][2];
        int sphereRadius = centMassCand[i][3];
        centMassCand[i][4] = 0;
        for(int ia = ia0 - sphereRadius; ia <= ia0 + sphereRadius; ++ia){
            for(int ib = ib0 - sphereRadius; ib <= ib0 + sphereRadius; ++ib){
                for(int ic = ic0 - sphereRadius; ic <= ic0 + sphereRadius; ++ic){
                    if( (ia - ia0)*(ia - ia0) + (ib - ib0)*(ib - ib0) + (ic - ic0)*(ic - ic0) < sphereRadius*sphereRadius ){
                        int ia1 = ia;
                        int ib1 = ib;
                        int ic1 = ic;
                        while(ia1 < 0) ia1 += ANUM;
                        while(ib1 < 0) ib1 += BNUM;
                        while(ic1 < 0) ic1 += CNUM;
                        while(ia1 >= ANUM) ia1 -= ANUM;
                        while(ib1 >= BNUM) ib1 -= BNUM;
                        while(ic1 >= CNUM) ic1 -= CNUM;
                        if(weigAvgDensUnitCell[ia1][ib1][ic1] > densCutoff)
                            centMassCand[i][4] += weigAvgDensUnitCell[ia1][ib1][ic1];
                    }
                }
            }
        }
    }

    for(int i = 0; i < 1; ++i){ //only need to find the maxmum
        for(int j = i + 1; j < numCentMassCand; ++j){
            if(centMassCand[j][4] > centMassCand[i][4]){
                for(int n = 0; n <= 4; ++n){
                    float temp = centMassCand[i][n];
                    centMassCand[i][n] = centMassCand[j][n];
                    centMassCand[j][n] = temp;
                }
            }
        }
    }

    int ia = centMassCand[0][0];
    int ib = centMassCand[0][1];
    int ic = centMassCand[0][2];
    frac[0] = (ia + 0.5) / float(ANUM);
    frac[1] = (ib + 0.5) / float(BNUM);
    frac[2] = (ic + 0.5) / float(CNUM);
    convertFracToOrth(frac, orth);
    //by default, numCentMass = 1
    centMassOrth[0][0] = orth[0];
    centMassOrth[0][1] = orth[1];
    centMassOrth[0][2] = orth[2];
    return;
}

/******************************************************************************/

void updateCentMassForMonomer(float densCutoff)
{
    float frac[3], orth[3];

    for(int i = 0; i <= numCentMassCand - 1; ++i){
        int ia0 = centMassCand[i][0];
        int ib0 = centMassCand[i][1];
        int ic0 = centMassCand[i][2];
        int sphereRadius = centMassCand[i][3];
        centMassCand[i][4] = 0;

        float point[3], newPoint[3];
        frac[0] = (ia0 + 0.5) / float(ANUM);
        frac[1] = (ib0 + 0.5) / float(BNUM);
        frac[2] = (ic0 + 0.5) / float(CNUM);
        convertFracToOrth(frac, point);
        for(int incs = 0; incs <= numNcsOper - 1; ++incs){
            float thita, dist;
            if(tranNcsFlag){
                thita = incs * rotaNcsAngl;
                float tempPoint[3] = {};
                applyNcsRotaOnPoint(point, thita, tempPoint);
                dist = incs * tranNcsDist;
                applyNcsTranOnPoint(tempPoint, dist, newPoint);
            } else {
                thita = 2.0 * PI * float(incs) / float(numNcsOper);
                applyNcsRotaOnPoint(point, thita, newPoint);
            }
            convertOrthToFrac(newPoint, frac);
            int ia1 = floor(frac[0] * ANUM);
            int ib1 = floor(frac[1] * BNUM);
            int ic1 = floor(frac[2] * CNUM);
            for(int ia = ia1 - sphereRadius; ia <= ia1 + sphereRadius; ++ia){
                for(int ib = ib1 - sphereRadius; ib <= ib1 + sphereRadius; ++ib){
                    for(int ic = ic1 - sphereRadius; ic <= ic1 + sphereRadius; ++ic){
                        if( (ia - ia1)*(ia - ia1) + (ib - ib1)*(ib - ib1) + (ic - ic1)*(ic - ic1) < sphereRadius*sphereRadius ){
                            int ia2 = ia;
                            int ib2 = ib;
                            int ic2 = ic;
                            while(ia2 < 0) ia2 += ANUM;
                            while(ib2 < 0) ib2 += BNUM;
                            while(ic2 < 0) ic2 += CNUM;
                            while(ia2 >= ANUM) ia2 -= ANUM;
                            while(ib2 >= BNUM) ib2 -= BNUM;
                            while(ic2 >= CNUM) ic2 -= CNUM;
                            if(weigAvgDensUnitCell[ia2][ib2][ic2] > densCutoff)
                                centMassCand[i][4] += weigAvgDensUnitCell[ia2][ib2][ic2];
                        }
                    }
                }
            }
        } //incs
    }

    for(int i = 0; i < 1; ++i){ //only need to find the minimum
        for(int j = i + 1; j < numCentMassCand; ++j){
            if(centMassCand[j][4] > centMassCand[i][4]){
                for(int n = 0; n <= 4; ++n){
                    float temp = centMassCand[i][n];
                    centMassCand[i][n] = centMassCand[j][n];
                    centMassCand[j][n] = temp;
                }
            }
        }
    }

    int ia = centMassCand[0][0];
    int ib = centMassCand[0][1];
    int ic = centMassCand[0][2];
    frac[0] = (ia + 0.5) / float(ANUM);
    frac[1] = (ib + 0.5) / float(BNUM);
    frac[2] = (ic + 0.5) / float(CNUM);
    convertFracToOrth(frac, orth);

    //by default, numCentMass = 2
    for(int incs = 0; incs <= numNcsOper - 1; ++incs){
        float newPoint[3];
        float thita, dist;
        if(tranNcsFlag){
            thita = incs * rotaNcsAngl;
            float tempPoint[3] = {};
            applyNcsRotaOnPoint(orth, thita, tempPoint);
            dist = incs * tranNcsDist;
            applyNcsTranOnPoint(tempPoint, dist, newPoint);
        } else {
            thita = 2.0 * PI * float(incs) / float(numNcsOper);
            applyNcsRotaOnPoint(orth, thita, newPoint);
        }
        centMassOrth[incs][0] = newPoint[0];
        centMassOrth[incs][1] = newPoint[1];
        centMassOrth[incs][2] = newPoint[2];
    }
    return;
}

/******************************************************************************/

void updateNcsCentMassUnitCell(float densCutoff)
{
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {

                ncsCentMassCandUnitCell[ia][ib][ic] = 0;

                int random = rand() % 10000;
                if(random > 100) continue;
                
                float frac[3], pointOnNcsAxis[3];
                frac[0] = (ia + 0.5) / float(ANUM);
                frac[1] = (ib + 0.5) / float(BNUM);
                frac[2] = (ic + 0.5) / float(CNUM);
                convertFracToOrth(frac, pointOnNcsAxis);
                
                int cylinderRadius = maxPackCylinderRadiusUnitCell[ia][ib][ic];
                float cylinderRadiusReal = cylinderRadius * (CRYST_A / ANUM); //potential bug
                int cylinderHeight = 2 * cylinderRadius; //must be an even number
                float cylinderHeightReal = 2.0 * cylinderRadiusReal;
                int sphereRadius = ceil( sqrt( cylinderRadius*cylinderRadius + cylinderHeight/2.0*cylinderHeight/2.0 ) );
                for(int ia1 = ia-sphereRadius; ia1 <= ia+sphereRadius; ++ia1) {
                    for(int ib1 = ib-sphereRadius; ib1 <= ib+sphereRadius; ++ib1) {
                        for(int ic1 = ic-sphereRadius; ic1 <= ic+sphereRadius; ++ic1) {
                            //mush inside a big sphere
                            if( (ia1 - ia)*(ia1 - ia) + (ib1 - ib)*(ib1 - ib) + (ic1 - ic)*(ic1 - ic) <= sphereRadius*sphereRadius ) {
                                float frac1[3], orth1[3];
                                frac1[0] = (ia1 + 0.5) / float(ANUM);
                                frac1[1] = (ib1 + 0.5) / float(BNUM);
                                frac1[2] = (ic1 + 0.5) / float(CNUM);
                                convertFracToOrth(frac1, orth1);
                                float pointBOnLine[3];
                                for(int i = 0; i < 3; ++i){
                                    pointBOnLine[i] = pointOnNcsAxis[i] + vectOnNcsAxis[i];
                                }
                                float dist = computeDistFromPointToLine(orth1, pointOnNcsAxis, pointBOnLine);
                                if(dist < cylinderRadiusReal){ //must within a dist about the axis
                                    float dist1 = 0;
                                    for(int i = 0; i < 3; ++i){
                                        dist1 = dist1 + (orth1[i] - pointOnNcsAxis[i]) * vectOnNcsAxis[i];
                                    }
                                    dist1 = abs(dist1);
                                    if(dist1 <= cylinderHeightReal/2.0){ //must inside a cylinder

                                        float ncsAvgMass = 0;
                                        for(int incs = 0; incs <= numNcsOper - 1; ++incs){
                                            float thita = 2.0 * PI * float(incs) / float(numNcsOper);
                                            float frac2[3], newPoint[3];
                                            applyRotaAboutAnyAxis(pointOnNcsAxis, vectOnNcsAxis, thita, orth1, newPoint);
                                            convertOrthToFrac(newPoint, frac2);
                                            int ia2 = floor(frac2[0] * ANUM);
                                            int ib2 = floor(frac2[1] * BNUM);
                                            int ic2 = floor(frac2[2] * CNUM);
                                            while(ia2 < 0) ia2 += ANUM;
                                            while(ib2 < 0) ib2 += BNUM;
                                            while(ic2 < 0) ic2 += CNUM;
                                            while(ia2 >= ANUM) ia2 -= ANUM;
                                            while(ib2 >= BNUM) ib2 -= BNUM;
                                            while(ic2 >= CNUM) ic2 -= CNUM;
                                            if(weigAvgDensUnitCell[ia2][ib2][ic2] > densCutoff) {
                                                ncsAvgMass += weigAvgDensUnitCell[ia2][ib2][ic2];
                                            }
                                        }
                                        ncsAvgMass = ncsAvgMass / float(numNcsOper);
                                        int ia2 = ia1;
                                        int ib2 = ib1;
                                        int ic2 = ic1;
                                        while(ia2 < 0) ia2 += ANUM;
                                        while(ib2 < 0) ib2 += BNUM;
                                        while(ic2 < 0) ic2 += CNUM;
                                        while(ia2 >= ANUM) ia2 -= ANUM;
                                        while(ib2 >= BNUM) ib2 -= BNUM;
                                        while(ic2 >= CNUM) ic2 -= CNUM;
                                        if(weigAvgDensUnitCell[ia2][ib2][ic2] < ncsAvgMass){
                                            ncsCentMassCandUnitCell[ia][ib][ic] += weigAvgDensUnitCell[ia2][ib2][ic2];
                                        } else{
                                            ncsCentMassCandUnitCell[ia][ib][ic] += ncsAvgMass;
                                        }
                                    } //if(dist1 <= cylinderHeightReal/2.0)
                                } //if(dist < cylinderRadiusReal)
                            } //if(dist < sphereRadius*sphereRadius)
                        } //ic1
                    } //ib1
                } //ia1
            } //ic
        } //ib
    } //ia


    //find the grid point in unitcell with maximum ncsCentMassCandUnitCell
    if(numCentMass == 1){ //3-fold and multi-fold NCS using a cylinder
        float mass = -1000.0;
        for (int ia = 0; ia <= ANUM - 1; ++ia) {
            for (int ib = 0; ib <= BNUM - 1; ++ib) {
                for (int ic = 0; ic <= CNUM - 1; ++ic) {
                    if(ncsCentMassCandUnitCell[ia][ib][ic] > mass){
                        mass = ncsCentMassCandUnitCell[ia][ib][ic];
                        float frac[3], orth[3];
                        frac[0] = (ia + 0.5) / float(ANUM);
                        frac[1] = (ib + 0.5) / float(BNUM);
                        frac[2] = (ic + 0.5) / float(CNUM);
                        convertFracToOrth(frac, orth);
                        centMassOrth[0][0] = orth[0];
                        centMassOrth[0][1] = orth[1];
                        centMassOrth[0][2] = orth[2];
                        //this part about pointAOnNcsAxisCand is not necessary, just for minimizeNcsDensDiffNcsAvgAxis()
                        for(int icand = 0; icand < numNcsAxisCand; ++icand){
                            pointAOnNcsAxisCand[icand][0] = orth[0];
                            pointAOnNcsAxisCand[icand][1] = orth[1];
                            pointAOnNcsAxisCand[icand][2] = orth[2];
                        }
                    }
                }
            }
        }
    } else { //2-fold NCS using two spheres
        cout << "Currently, the code can not deal with 2-fold ncs." << endl;
        exit(EXIT_FAILURE);
    }

    return;
}

/******************************************************************************/

void refineNcsCentMassLocally(float densCutoff)
{
    //Currently, the code can not deal with 2-fold ncs.
    float frac[3], orth[3];
    orth[0] = centMassOrth[0][0];
    orth[1] = centMassOrth[0][1];
    orth[2] = centMassOrth[0][2];
    convertOrthToFrac(orth, frac);
    int ia00 = floor(frac[0] * ANUM);
    int ib00 = floor(frac[1] * BNUM);
    int ic00 = floor(frac[2] * CNUM);

    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                ncsCentMassCandUnitCell[ia][ib][ic] = 0;
            }
        }
    }
    for (int ia0 = ia00-10; ia0 <= ia00+10; ++ia0) {
        for (int ib0 = ib00-10; ib0 <= ib00+10; ++ib0) {
            for (int ic0 = ic00-10; ic0 <= ic00+10; ++ic0) {
                int ia = ia0;
                while(ia < 0) ia += ANUM;
                while(ia >= ANUM) ia -= ANUM;
                int ib = ib0;
                while(ib < 0) ib += BNUM;
                while(ib >= BNUM) ib -= BNUM;
                int ic = ic0;
                while(ic < 0) ic += CNUM;
                while(ic >= CNUM) ic -= CNUM;
                
                float frac[3], pointOnNcsAxis[3];
                frac[0] = (ia + 0.5) / float(ANUM);
                frac[1] = (ib + 0.5) / float(BNUM);
                frac[2] = (ic + 0.5) / float(CNUM);
                convertFracToOrth(frac, pointOnNcsAxis);
                ncsCentMassCandUnitCell[ia][ib][ic] = 0;

                int cylinderRadius = maxPackCylinderRadiusUnitCell[ia][ib][ic];
                float cylinderRadiusReal = cylinderRadius * (CRYST_A / ANUM); //potential bug
                int cylinderHeight = 2 * cylinderRadius; //must be an even number
                float cylinderHeightReal = 2.0 * cylinderRadiusReal;
                int sphereRadius = ceil( sqrt( cylinderRadius*cylinderRadius + cylinderHeight/2.0*cylinderHeight/2.0 ) );
                for(int ia1 = ia-sphereRadius; ia1 <= ia+sphereRadius; ++ia1) {
                    for(int ib1 = ib-sphereRadius; ib1 <= ib+sphereRadius; ++ib1) {
                        for(int ic1 = ic-sphereRadius; ic1 <= ic+sphereRadius; ++ic1) {
                            //mush inside a big sphere
                            if( (ia1 - ia)*(ia1 - ia) + (ib1 - ib)*(ib1 - ib) + (ic1 - ic)*(ic1 - ic) <= sphereRadius*sphereRadius ) {
                                float frac1[3], orth1[3];
                                frac1[0] = (ia1 + 0.5) / float(ANUM);
                                frac1[1] = (ib1 + 0.5) / float(BNUM);
                                frac1[2] = (ic1 + 0.5) / float(CNUM);
                                convertFracToOrth(frac1, orth1);
                                float pointBOnLine[3];
                                for(int i = 0; i < 3; ++i){
                                    pointBOnLine[i] = pointOnNcsAxis[i] + vectOnNcsAxis[i];
                                }
                                float dist = computeDistFromPointToLine(orth1, pointOnNcsAxis, pointBOnLine);
                                if(dist < cylinderRadiusReal){ //must within a dist about the axis
                                    float dist1 = 0;
                                    for(int i = 0; i < 3; ++i){
                                        dist1 = dist1 + (orth1[i] - pointOnNcsAxis[i]) * vectOnNcsAxis[i];
                                    }
                                    dist1 = abs(dist1);
                                    if(dist1 <= cylinderHeightReal/2.0){ //must inside a cylinder

                                        float ncsAvgMass = 0;
                                        for(int incs = 0; incs <= numNcsOper - 1; ++incs){
                                            float thita = 2.0 * PI * float(incs) / float(numNcsOper);
                                            float frac2[3], newPoint[3];
                                            applyRotaAboutAnyAxis(pointOnNcsAxis, vectOnNcsAxis, thita, orth1, newPoint);
                                            convertOrthToFrac(newPoint, frac2);
                                            int ia2 = floor(frac2[0] * ANUM);
                                            int ib2 = floor(frac2[1] * BNUM);
                                            int ic2 = floor(frac2[2] * CNUM);
                                            while(ia2 < 0) ia2 += ANUM;
                                            while(ib2 < 0) ib2 += BNUM;
                                            while(ic2 < 0) ic2 += CNUM;
                                            while(ia2 >= ANUM) ia2 -= ANUM;
                                            while(ib2 >= BNUM) ib2 -= BNUM;
                                            while(ic2 >= CNUM) ic2 -= CNUM;
                                            if(weigAvgDensUnitCell[ia2][ib2][ic2] > densCutoff) {
                                                ncsAvgMass += weigAvgDensUnitCell[ia2][ib2][ic2];
                                            }
                                        }
                                        ncsAvgMass = ncsAvgMass / float(numNcsOper);
                                        int ia2 = ia1;
                                        int ib2 = ib1;
                                        int ic2 = ic1;
                                        while(ia2 < 0) ia2 += ANUM;
                                        while(ib2 < 0) ib2 += BNUM;
                                        while(ic2 < 0) ic2 += CNUM;
                                        while(ia2 >= ANUM) ia2 -= ANUM;
                                        while(ib2 >= BNUM) ib2 -= BNUM;
                                        while(ic2 >= CNUM) ic2 -= CNUM;
                                        if(weigAvgDensUnitCell[ia2][ib2][ic2] < ncsAvgMass){
                                            ncsCentMassCandUnitCell[ia][ib][ic] += weigAvgDensUnitCell[ia2][ib2][ic2];
                                        } else{
                                            ncsCentMassCandUnitCell[ia][ib][ic] += ncsAvgMass;
                                        }
                                    } //if(dist1 <= cylinderHeightReal/2.0)
                                } //if(dist < cylinderRadiusReal)
                            } //if(dist < sphereRadius*sphereRadius)
                        } //ic1
                    } //ib1
                } //ia1
            } //ic0
        } //ib0
    } //ia0


    //find the grid point in unitcell with maximum ncsCentMassCandUnitCell
    if(numCentMass == 1){ //3-fold and multi-fold NCS using a cylinder
        float mass = -1000.0;
        for (int ia = 0; ia <= ANUM - 1; ++ia) {
            for (int ib = 0; ib <= BNUM - 1; ++ib) {
                for (int ic = 0; ic <= CNUM - 1; ++ic) {
                    if(ncsCentMassCandUnitCell[ia][ib][ic] > mass){
                        mass = ncsCentMassCandUnitCell[ia][ib][ic];
                        float frac[3], orth[3];
                        frac[0] = (ia + 0.5) / float(ANUM);
                        frac[1] = (ib + 0.5) / float(BNUM);
                        frac[2] = (ic + 0.5) / float(CNUM);
                        convertFracToOrth(frac, orth);
                        centMassOrth[0][0] = orth[0];
                        centMassOrth[0][1] = orth[1];
                        centMassOrth[0][2] = orth[2];
                        //this part about pointAOnNcsAxisCand is not necessary, just for minimizeNcsDensDiffNcsAvgAxis()
                        for(int icand = 0; icand < numNcsAxisCand; ++icand){
                            pointAOnNcsAxisCand[icand][0] = orth[0];
                            pointAOnNcsAxisCand[icand][1] = orth[1];
                            pointAOnNcsAxisCand[icand][2] = orth[2];
                        }
                    }
                }
            }
        }
    } else { //2-fold NCS using two spheres
        cout << "Currently, the code can not deal with 2-fold ncs." << endl;
        exit(EXIT_FAILURE);
    }

    return;
}

/******************************************************************************/

void updateNcsCentMassUnitCellViaConvCylinder()
{
    int numCand = 100;
    float centMassCand[numCand][3];

    for(int icand = 0; icand < numCand; ++icand) {
        float cylinderConvDensMax = 0;
        for (int ia = 0; ia <= ANUM - 1; ++ia) {
            for (int ib = 0; ib <= BNUM - 1; ++ib) {
                for (int ic = 0; ic <= CNUM - 1; ++ic) {
                    if(convDensUnitCell[ia][ib][ic] > cylinderConvDensMax) {
                        cylinderConvDensMax = convDensUnitCell[ia][ib][ic];
                        centMassCand[icand][0] = ia;
                        centMassCand[icand][1] = ib;
                        centMassCand[icand][2] = ic;
                    }
                }
            }
        }
        int ia0 = centMassCand[icand][0];
        int ib0 = centMassCand[icand][1];
        int ic0 = centMassCand[icand][2];
        int sphereRadius = 2; //can be modifed
        for(int ia1 = ia0-sphereRadius; ia1 <= ia0+sphereRadius; ++ia1) {
            for(int ib1 = ib0-sphereRadius; ib1 <= ib0+sphereRadius; ++ib1) {
                for(int ic1 = ic0-sphereRadius; ic1 <= ic0+sphereRadius; ++ic1) {
                    //inside a sphere, convDensUnitCell = 0
                    if( (ia1 - ia0)*(ia1 - ia0) + (ib1 - ib0)*(ib1 - ib0) + (ic1 - ic0)*(ic1 - ic0) <= sphereRadius*sphereRadius ) {
                        int ia2 = ia1;
                        int ib2 = ib1;
                        int ic2 = ic1;
                        while(ia2 < 0) ia2 += ANUM;
                        while(ib2 < 0) ib2 += BNUM;
                        while(ic2 < 0) ic2 += CNUM;
                        while(ia2 >= ANUM) ia2 -= ANUM;
                        while(ib2 >= BNUM) ib2 -= BNUM;
                        while(ic2 >= CNUM) ic2 -= CNUM;
                        convDensUnitCell[ia2][ib2][ic2] = 0;
                    }
                }
            }
        }
    }

    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                ncsCentMassCandUnitCell[ia][ib][ic] = 0;
            }
        }
    }

    //apply ncs, find the centMassCand which best satisfies NCS operations
    for(int icand = 0; icand < numCand; ++icand) {
        int ia = centMassCand[icand][0];
        int ib = centMassCand[icand][1];
        int ic = centMassCand[icand][2];
        float frac[3], pointOnNcsAxis[3];
        frac[0] = (ia + 0.5) / float(ANUM);
        frac[1] = (ib + 0.5) / float(BNUM);
        frac[2] = (ic + 0.5) / float(CNUM);
        convertFracToOrth(frac, pointOnNcsAxis);
    
        int cylinderRadius = maxPackCylinderRadiusUnitCell[ia][ib][ic]; //potential bug
        float cylinderRadiusReal = cylinderRadius * (CRYST_A / ANUM); //potential bug
        int cylinderHeight = 2 * cylinderRadius; //must be an even number
        float cylinderHeightReal = 2.0 * cylinderRadiusReal;
        int sphereRadius = ceil( sqrt( cylinderRadius*cylinderRadius + cylinderHeight/2.0*cylinderHeight/2.0 ) );
        for(int ia1 = ia-sphereRadius; ia1 <= ia+sphereRadius; ++ia1) {
            for(int ib1 = ib-sphereRadius; ib1 <= ib+sphereRadius; ++ib1) {
                for(int ic1 = ic-sphereRadius; ic1 <= ic+sphereRadius; ++ic1) {
                    //must inside a big sphere
                    if( (ia1 - ia)*(ia1 - ia) + (ib1 - ib)*(ib1 - ib) + (ic1 - ic)*(ic1 - ic) <= sphereRadius*sphereRadius ) {
                        float frac1[3], orth1[3];
                        frac1[0] = (ia1 + 0.5) / float(ANUM);
                        frac1[1] = (ib1 + 0.5) / float(BNUM);
                        frac1[2] = (ic1 + 0.5) / float(CNUM);
                        convertFracToOrth(frac1, orth1);
                        float pointBOnLine[3];
                        for(int i = 0; i < 3; ++i){
                            pointBOnLine[i] = pointOnNcsAxis[i] + vectOnNcsAxis[i];
                        }
                        float dist = computeDistFromPointToLine(orth1, pointOnNcsAxis, pointBOnLine);
                        if(dist < cylinderRadiusReal){ //must within a dist about the axis
                            float dist1 = 0;
                            for(int i = 0; i < 3; ++i){
                                dist1 = dist1 + (orth1[i] - pointOnNcsAxis[i]) * vectOnNcsAxis[i];
                            }
                            dist1 = abs(dist1);
                            if(dist1 <= cylinderHeightReal/2.0){ //must inside a cylinder
                                float ncsAvgMass = 0;
                                for(int incs = 0; incs <= numNcsOper - 1; ++incs){
                                    float thita = 2.0 * PI * float(incs) / float(numNcsOper);
                                    float frac2[3], newPoint[3];
                                    applyRotaAboutAnyAxis(pointOnNcsAxis, vectOnNcsAxis, thita, orth1, newPoint);
                                    convertOrthToFrac(newPoint, frac2);
                                    int ia2 = floor(frac2[0] * ANUM);
                                    int ib2 = floor(frac2[1] * BNUM);
                                    int ic2 = floor(frac2[2] * CNUM);
                                    while(ia2 < 0) ia2 += ANUM;
                                    while(ib2 < 0) ib2 += BNUM;
                                    while(ic2 < 0) ic2 += CNUM;
                                    while(ia2 >= ANUM) ia2 -= ANUM;
                                    while(ib2 >= BNUM) ib2 -= BNUM;
                                    while(ic2 >= CNUM) ic2 -= CNUM;
                                    if(weigAvgDensUnitCell[ia2][ib2][ic2] > densCutoff) {
                                        ncsAvgMass += weigAvgDensUnitCell[ia2][ib2][ic2];
                                    }
                                }
                                ncsAvgMass = ncsAvgMass / float(numNcsOper);
                                int ia2 = ia1;
                                int ib2 = ib1;
                                int ic2 = ic1;
                                while(ia2 < 0) ia2 += ANUM;
                                while(ib2 < 0) ib2 += BNUM;
                                while(ic2 < 0) ic2 += CNUM;
                                while(ia2 >= ANUM) ia2 -= ANUM;
                                while(ib2 >= BNUM) ib2 -= BNUM;
                                while(ic2 >= CNUM) ic2 -= CNUM;
//                                if(weigAvgDensUnitCell[ia2][ib2][ic2] < ncsAvgMass){ //using the smaller one
//                                    ncsCentMassCandUnitCell[ia][ib][ic] += weigAvgDensUnitCell[ia2][ib2][ic2];
//                                } else{
//                                    ncsCentMassCandUnitCell[ia][ib][ic] += ncsAvgMass;
//                                }
                                ncsCentMassCandUnitCell[ia][ib][ic] += ncsAvgMass * weigAvgDensUnitCell[ia2][ib2][ic2]; //using the inner product
                            } //if(dist1 <= cylinderHeightReal/2.0)
                        } //if(dist < cylinderRadiusReal)
                    } //if(dist < sphereRadius*sphereRadius)
                } //ic1
            } //ib1
        } //ia1
    } //icand


    //find the grid point in unitcell with maximum ncsCentMassCandUnitCell
    if(numCentMass == 1){ //3-fold and multi-fold NCS using a cylinder
        float mass = 0;
        for (int ia = 0; ia <= ANUM - 1; ++ia) {
            for (int ib = 0; ib <= BNUM - 1; ++ib) {
                for (int ic = 0; ic <= CNUM - 1; ++ic) {
                    if(ncsCentMassCandUnitCell[ia][ib][ic] > mass){
                        mass = ncsCentMassCandUnitCell[ia][ib][ic];
                        float frac[3], orth[3];
                        frac[0] = (ia + 0.5) / float(ANUM);
                        frac[1] = (ib + 0.5) / float(BNUM);
                        frac[2] = (ic + 0.5) / float(CNUM);
                        convertFracToOrth(frac, orth);
                        centMassOrth[0][0] = orth[0];
                        centMassOrth[0][1] = orth[1];
                        centMassOrth[0][2] = orth[2];
                        //this part about pointAOnNcsAxisCand is not necessary, just for minimizeNcsDensDiffNcsAvgAxis()
                        for(int icand = 0; icand < numNcsAxisCand; ++icand){
                            pointAOnNcsAxisCand[icand][0] = orth[0];
                            pointAOnNcsAxisCand[icand][1] = orth[1];
                            pointAOnNcsAxisCand[icand][2] = orth[2];
                        }
                    }
                }
            }
        }
    } else { //2-fold NCS using two spheres
        cout << "Currently, the code can not deal with 2-fold ncs." << endl;
        exit(EXIT_FAILURE);
    }

    return;
}

/******************************************************************************/

void updateNcsCentMassViaCylinder(float densCutoff)
{

    float frac[3];
    convertOrthToFrac(pointAOnNcsAxis, frac);
    int ia0 = floor(frac[0] * ANUM);
    int ib0 = floor(frac[1] * BNUM);
    int ic0 = floor(frac[2] * CNUM);
    int cylinderRadius = maxPackCylinderRadiusUnitCell[ia0][ib0][ic0];
    int cylinderHeight = 2 * cylinderRadius; //must be an even number
    int sphereRadius = ceil( sqrt( cylinderRadius*cylinderRadius + cylinderHeight/2.0*cylinderHeight/2.0 ) );

    float cylinderRadiusReal = cylinderRadius * (CRYST_A / ANUM); //potential bug
    float cylinderHeightReal = 2.0 * cylinderRadiusReal;

    float pointBOnNcsAxis[3];
    for(int i = 0; i < 3; ++i){
        pointBOnNcsAxis[i] = pointAOnNcsAxis[i] + vectOnNcsAxis[i];
    }
    float centNcsMask[3] = {};
    float sumDens = 0;

    for(int ia1 = ia0-sphereRadius; ia1 <= ia0+sphereRadius; ++ia1) {
        for(int ib1 = ib0-sphereRadius; ib1 <= ib0+sphereRadius; ++ib1) {
            for(int ic1 = ic0-sphereRadius; ic1 <= ic0+sphereRadius; ++ic1) {
                
                if( (ia1 - ia0)*(ia1 - ia0) + (ib1 - ib0)*(ib1 - ib0) + (ic1 - ic0)*(ic1 - ic0) <= sphereRadius*sphereRadius ) {
                    float frac1[3], orth1[3];
                    frac1[0] = (ia1 + 0.5) / float(ANUM);
                    frac1[1] = (ib1 + 0.5) / float(BNUM);
                    frac1[2] = (ic1 + 0.5) / float(CNUM);
                    convertFracToOrth(frac1, orth1);
                    float dist = computeDistFromPointToLine(orth1, pointAOnNcsAxis, pointBOnNcsAxis);
                    if(dist < cylinderRadiusReal){ //must within a dist about the axis
                        float dist1 = 0;
                        for(int i = 0; i < 3; ++i){
                            dist1 = dist1 + (orth1[i] - pointAOnNcsAxis[i]) * vectOnNcsAxis[i];
                        }
                        dist1 = abs(dist1);
                        if(dist1 <= cylinderHeightReal/2.0){ //must inside a cylinder
                            int ia2 = ia1;
                            int ib2 = ib1;
                            int ic2 = ic1;
                            while(ia2 < 0) ia2 += ANUM;
                            while(ib2 < 0) ib2 += BNUM;
                            while(ic2 < 0) ic2 += CNUM;
                            while(ia2 >= ANUM) ia2 -= ANUM;
                            while(ib2 >= BNUM) ib2 -= BNUM;
                            while(ic2 >= CNUM) ic2 -= CNUM;
                            if(weigAvgDensUnitCell[ia2][ib2][ic2] > densCutoff){
                                centNcsMask[0] += orth1[0] * weigAvgDensUnitCell[ia2][ib2][ic2];
                                centNcsMask[1] += orth1[1] * weigAvgDensUnitCell[ia2][ib2][ic2];
                                centNcsMask[2] += orth1[2] * weigAvgDensUnitCell[ia2][ib2][ic2];
                                sumDens += weigAvgDensUnitCell[ia2][ib2][ic2];
                            } //within the protein mask
                        } //must inside a cylinder
                    } //must within a dist about the axis
                } //must inside a big sphere
            } //ic1
        } //ib1
    } //ia1
    centNcsMask[0] /= sumDens;
    centNcsMask[1] /= sumDens;
    centNcsMask[2] /= sumDens;

    //find the grid point in unitcell with maximum ncsCentMassCandUnitCell
    if(numCentMass == 1){ //3-fold and multi-fold NCS using a cylinder
        centMassOrth[0][0] = centNcsMask[0];
        centMassOrth[0][1] = centNcsMask[1];
        centMassOrth[0][2] = centNcsMask[2];
        pointAOnNcsAxis[0] = centNcsMask[0];
        pointAOnNcsAxis[1] = centNcsMask[1];
        pointAOnNcsAxis[2] = centNcsMask[2];;
    } else { //2-fold NCS using two spheres
        cout << "Currently, the code can not deal with 2-fold ncs." << endl;
        exit(EXIT_FAILURE);
    }

    return;
}

/******************************************************************************/

void updateCylinderDens(int ia0, int ib0, int ic0)
{
    for(int incs = 0; incs < numNcsOper; ++incs) {
        for(int igrid = 0; igrid < numGridPointInCylinder; ++igrid) {
            int ia = round(cylinderDens[igrid][incs][0]) + ia0;
            int ib = round(cylinderDens[igrid][incs][1]) + ib0;
            int ic = round(cylinderDens[igrid][incs][2]) + ic0;
            while(ia < 0) ia += ANUM;
            while(ib < 0) ib += BNUM;
            while(ic < 0) ic += CNUM;
            while(ia >= ANUM) ia -= ANUM;
            while(ib >= BNUM) ib -= BNUM;
            while(ic >= CNUM) ic -= CNUM;
            cylinderDens[igrid][incs][3] = weigAvgDensUnitCell[ia][ib][ic];
        }
    }
    return;
}

/******************************************************************************/

void inverseSecoNcsOperAsFirstNcsOper()
{
    float rota[3][3], tran[3];
    float dete = 0;
    for(int i = 0; i <= 2; ++i){
        for(int j = 0; j <= 2; ++j){
            rota[i][j] = ncsOper[1][i][j];
        }
        tran[i] = ncsOper[1][i][3];
    }

    //find determinant of rotation matrix
    for(int i = 0; i <= 2; ++i){
        dete = dete + ( rota[0][i] * (rota[1][(i+1)%3] * rota[2][(i+2)%3] 
            - rota[1][(i+2)%3] * rota[2][(i+1)%3]) );
    }
    
    //compute inverse of rotation matrix
    for(int i = 0; i <= 2; ++i){
        for(int j = 0; j <= 2; ++j){
            ncsOper[0][i][j] = (  ( rota[(j+1)%3][(i+1)%3] * rota[(j+2)%3][(i+2)%3] ) 
                - ( rota[(j+1)%3][(i+2)%3] * rota[(j+2)%3][(i+1)%3] )  )/ dete;
        }
    }

    //compute new translation
    for(int i = 0; i <= 2; ++i){
        ncsOper[0][i][3] = 0;
        for(int j = 0; j <= 2; ++j){
            ncsOper[0][i][3] = ncsOper[0][i][3] + (-1.0) * ( ncsOper[0][i][j] * tran[j] );
        }
    }
    return;
}

/******************************************************************************/

void setHalfGridPhase()
{
    for(int ih = 0; ih <= HNUM - 1; ih++ ) {
        for(int ik = 0; ik <= KNUM - 1; ik++ ) {
            for(int il = 0; il <= LNUM - 1; il++ ) {
                int h, k, l;
                if(ih <= HNUM / 2 - 1 ) {
                    h=ih;
                } else {
                    h=ih-HNUM;
                }
                if(ik <= KNUM / 2 - 1 ) {
                    k=ik;
                } else {
                    k=ik-KNUM;
                }
                if(il <= LNUM / 2 - 1 ) {
                    l=il;
                } else {
                    l=il-LNUM;
                }
                halfGridPosi[ih][ik][il] = complex<float> ( cos(2.0*PI*(h*0.5/float(ANUM)+k*0.5/float(BNUM)+l*0.5/float(CNUM)) ), 
                                                            sin(2.0*PI*(h*0.5/float(ANUM)+k*0.5/float(BNUM)+l*0.5/float(CNUM)) ) );
                halfGridNega[ih][ik][il] = complex<float> ( cos(2.0*PI*(h*0.5/float(ANUM)+k*0.5/float(BNUM)+l*0.5/float(CNUM)) ), 
                                                           -sin(2.0*PI*(h*0.5/float(ANUM)+k*0.5/float(BNUM)+l*0.5/float(CNUM)) ) );
            }
        }
    }
    return;
}

/******************************************************************************/

void applyFFTDensToStruFact()
{
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                realSpaceArray[ia][ib][ic] = densUnitCell[ia][ib][ic];
            }
        }
    }

    //Intel dfti backward fft from density to structure factor
    status = DftiCreateDescriptor( &myDescriptorHandle, DFTI_SINGLE,DFTI_COMPLEX, 3, lengthEachDimension);
    status = DftiSetValue( myDescriptorHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE );
    status = DftiSetValue( myDescriptorHandle, DFTI_BACKWARD_SCALE, dftiBacScaleFactor );
    status = DftiCommitDescriptor( myDescriptorHandle);
    status = DftiComputeBackward( myDescriptorHandle, realSpaceArray, reciprocalSpaceArray);
    status = DftiFreeDescriptor(&myDescriptorHandle);

    for(int ih = 0; ih <= HNUM - 1; ih++ ) {
        for(int ik = 0; ik <= KNUM - 1; ik++ ) {
            for(int il = 0; il <= LNUM - 1; il++ ) {

                reciprocalSpaceArray[ih][ik][il] = reciprocalSpaceArray[ih][ik][il] * halfGridPosi[ih][ik][il];

                float realPart = real(reciprocalSpaceArray[ih][ik][il]);
                float imagPart = imag(reciprocalSpaceArray[ih][ik][il]);
                float calsf = sqrt(pow(realPart, 2.0) + pow(imagPart, 2.0));

                allStruFactAmpl[ih][ik][il]=calsf;

                //if(calsf < pow(10.0, -3.0)) continue;    //Amplitude is too small, phase won't make sense
                float phase;
                if(realPart > 0 && imagPart > 0) {
                    phase = atan(imagPart / realPart);
                } else if(realPart < 0 && imagPart > 0) {
                    phase = atan(imagPart / realPart) + PI;
                } else if(realPart < 0 && imagPart < 0) {
                    phase = atan(imagPart / realPart) + PI;
                } else if(realPart > 0 && imagPart < 0) {
                    phase = atan(imagPart / realPart) + 2.0*PI;
                } else if(realPart >= 0 && imagPart == 0) {
                    phase = 0;
                } else if(realPart == 0 && imagPart > 0) {
                    phase = 0.50 * PI;
                } else if(realPart < 0 && imagPart == 0) {
                    phase = PI;
                } else if(realPart == 0 && imagPart < 0) {
                    phase = 1.50 * PI;
                }
                while(phase > 2.0 * PI) phase = phase - 2.0 * PI;
                while(phase < 0) phase = phase + 2.0 * PI;

                allStruFactPhase[ih][ik][il] = phase;

            }
        }
    }
    return;
}

/******************************************************************************/

void applyFFTStruFactToDens()
{
    for (int ih = 0; ih <= HNUM - 1; ih++) {
        for (int ik = 0; ik <= KNUM - 1; ik++) {
            for (int il = 0; il <= LNUM - 1; il++) {
                float obssf = allStruFactAmpl[ih][ik][il];
                float phase = allStruFactPhase[ih][ik][il];
                reciprocalSpaceArray[ih][ik][il] = obssf * complex<float>(cos(phase), sin(phase)) * halfGridNega[ih][ik][il];
            }
        }
    }

    //Intel dfti forward fft from structure factor to density
    status = DftiCreateDescriptor( &myDescriptorHandle, DFTI_SINGLE,DFTI_COMPLEX, 3, lengthEachDimension);
    status = DftiSetValue( myDescriptorHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE );
    status = DftiSetValue( myDescriptorHandle, DFTI_FORWARD_SCALE, dftiForScaleFactor );
    status = DftiCommitDescriptor( myDescriptorHandle);
    status = DftiComputeForward( myDescriptorHandle, reciprocalSpaceArray, realSpaceArray);
    status = DftiFreeDescriptor(&myDescriptorHandle);

    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                densUnitCell[ia][ib][ic] = real(realSpaceArray[ia][ib][ic]);
            }
        }
    }
    return;
}

/******************************************************************************/

void computeWeigAvgDensUnitCell(float sigmWeigAvg)
{
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                if(densUnitCell[ia][ib][ic] >= 0) {
                    realSpaceArray[ia][ib][ic] = densUnitCell[ia][ib][ic];
                } else {
                    realSpaceArray[ia][ib][ic] = 0; //-densUnitCell[ia][ib][ic];
                }
            }
        }
    }

    //Intel dfti backward fft from density to structure factor
    status = DftiCreateDescriptor( &myDescriptorHandle, DFTI_SINGLE,DFTI_COMPLEX, 3, lengthEachDimension);
    status = DftiSetValue( myDescriptorHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE );
    status = DftiSetValue( myDescriptorHandle, DFTI_BACKWARD_SCALE, dftiBacScaleFactor );
    status = DftiCommitDescriptor( myDescriptorHandle);
    status = DftiComputeBackward( myDescriptorHandle, realSpaceArray, reciprocalSpaceArray);
    status = DftiFreeDescriptor(&myDescriptorHandle);

    for(int ih = 0; ih <= HNUM - 1; ih++ ) {
        for(int ik = 0; ik <= KNUM - 1; ik++ ) {
            for(int il = 0; il <= LNUM - 1; il++ ) {
                int h, k, l;
                if(ih <= HNUM / 2 - 1 ) {
                    h=ih;
                } else {
                    h=ih-HNUM;
                }
                if(ik <= KNUM / 2 - 1 ) {
                    k=ik;
                } else {
                    k=ik-KNUM;
                }
                if(il <= LNUM / 2 - 1 ) {
                    l=il;
                } else {
                    l=il-LNUM;
                }
                   float S = sqrt( h * h * S1 + k * k * S2 + l * l * S3 +
                          2.0 * h * k * S4 + 2.0 * h * l * S5 + 2.0 * k * l * S6 ) / CRYST_VOL;
                reciprocalSpaceArray[ih][ik][il] = reciprocalSpaceArray[ih][ik][il]
                * complex<float>( exp( - 2.0 * pow(PI*sigmWeigAvg*S, 2.0) ) );
            }
        }
    }

    //Intel dfti forward fft from structure factor to density
    status = DftiCreateDescriptor( &myDescriptorHandle, DFTI_SINGLE,DFTI_COMPLEX, 3, lengthEachDimension);
    status = DftiSetValue( myDescriptorHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE );
    status = DftiSetValue( myDescriptorHandle, DFTI_FORWARD_SCALE, dftiForScaleFactor );
    status = DftiCommitDescriptor( myDescriptorHandle);
    status = DftiComputeForward( myDescriptorHandle, reciprocalSpaceArray, realSpaceArray);
    status = DftiFreeDescriptor(&myDescriptorHandle);

    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                weigAvgDensUnitCellTemp[ia][ib][ic] = weigAvgDensUnitCell[ia][ib][ic]; //reserve previous value
                weigAvgDensUnitCell[ia][ib][ic] = real(realSpaceArray[ia][ib][ic]);
            }
        }
    }
    return;
}

/******************************************************************************/

void updateWeigAvgDensUnitCell(string updateSpeed)
{
    int gear;

    if(updateSpeed.find("fixed") != string::npos){
        gear = 0;
    } else if(updateSpeed.find("extremely slowly") != string::npos){
        gear = 1;
    } else if(updateSpeed.find("very slowly") != string::npos){
        gear = 2;
    } else if(updateSpeed.find("slowly") != string::npos){
        gear = 3;
    } else if(updateSpeed.find("normally") != string::npos){
        gear = 4;
    } else {
        cout << "Unrecognized updateSpeedWeigAvgDens";
        exit(EXIT_FAILURE);
    }

    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                switch(gear){
                    case 4:
                        weigAvgDensUnitCell[ia][ib][ic] = weigAvgDensUnitCell[ia][ib][ic];
                        break;
                    case 3:
                        weigAvgDensUnitCell[ia][ib][ic] = 0.5 * weigAvgDensUnitCell[ia][ib][ic] 
                                                        + 0.5 * weigAvgDensUnitCellTemp[ia][ib][ic];
                        break;
                    case 2:
                        weigAvgDensUnitCell[ia][ib][ic] = 0.05 * weigAvgDensUnitCell[ia][ib][ic] 
                                                        + 0.95 * weigAvgDensUnitCellTemp[ia][ib][ic];
                        break;
                    case 1:
                        weigAvgDensUnitCell[ia][ib][ic] = 0.01 * weigAvgDensUnitCell[ia][ib][ic] 
                                                        + 0.99 * weigAvgDensUnitCellTemp[ia][ib][ic];
                        break;
                    case 0:
                        weigAvgDensUnitCell[ia][ib][ic] = 0.00 * weigAvgDensUnitCell[ia][ib][ic] 
                                                        + 1.00 * weigAvgDensUnitCellTemp[ia][ib][ic];
                        break;
                }
            }
        }
    }

    return;
}

/******************************************************************************/

float computeDensCutoffUnitCell(float protCont)
{
    float densMin = densUnitCell[0][0][0];
    float densMax = densMin;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                float density = densUnitCell[ia][ib][ic];
                if (density < densMin) densMin = density;
                if (density > densMax) densMax = density;
            }
        }
    }

    int numGridPointHighDens = round(ANUM * BNUM * CNUM * protCont);
    float density1 = densMin;
    float density2 = densMax;
    float densCutoff;
    do    {
        densCutoff = (density2 + density1) / 2.0;
        int i=0;
        for (int ia = 0; ia <= ANUM - 1; ++ia) {
            for (int ib = 0; ib <= BNUM  - 1; ++ib) {
                for (int ic = 0; ic <= CNUM - 1; ++ic) {
                    if(densUnitCell[ia][ib][ic] > densCutoff)
                        i++;
                }
            }
        }
        if(i > numGridPointHighDens) {
            density1 = densCutoff;
        } else if(i < numGridPointHighDens) {
            density2 = densCutoff;
        } else {
            density2 = density1;
        }
    } while(  abs(density2 - density1) > pow(10.0, -4.0) );

    return densCutoff;
}

/******************************************************************************/

float computeWeigAvgDensCutoffUnitCell(float protCont)
{
    float densMin = weigAvgDensUnitCell[0][0][0];
    float densMax = densMin;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                float density = weigAvgDensUnitCell[ia][ib][ic];
                if (density < densMin) densMin = density;
                if (density > densMax) densMax = density;
            }
        }
    }

    int numGridPointHighDens = round(ANUM * BNUM * CNUM * protCont);

    float density1 = densMin;
    float density2 = densMax;
    float densCutoff;

    do    {
        densCutoff = (density2 + density1) / 2.0;
        int i = 0;
        for (int ia = 0; ia <= ANUM - 1; ++ia) {
            for (int ib = 0; ib <= BNUM  - 1; ++ib) {
                for (int ic = 0; ic <= CNUM - 1; ++ic) {
                    if(weigAvgDensUnitCell[ia][ib][ic] > densCutoff)
                        i++;
                }
            }
        }
        if(i > numGridPointHighDens) {
            density1 = densCutoff;
        } else if(i < numGridPointHighDens) {
            density2 = densCutoff;
        } else {
            density2 = density1;
        }
    } while(  abs(density2 - density1) > pow(10.0, -4.0) );

    return densCutoff;
}

/******************************************************************************/

void enantiomorphDensUnitCell()
{
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                densUnitCellTemp[ia][ib][ic] = densUnitCell[ia][ib][ic];
            }
        }
    }
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                float frac[3];
                frac[0] = - (ia + 0.5) / float(ANUM);
                frac[1] = - (ib + 0.5) / float(BNUM);
                frac[2] = - (ic + 0.5) / float(CNUM);
                adjustFracCoorIntoUnitCell(frac);
                int ia0 = floor(frac[0] * ANUM); //round down value (interger)
                 int ib0 = floor(frac[1] * BNUM);
                 int ic0 = floor(frac[2] * CNUM);
                densUnitCell[ia0][ib0][ic0] = densUnitCellTemp[ia][ib][ic];
            }
        }
    }
    //reset densUnitCell to zero
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                densUnitCellTemp[ia][ib][ic] = 0;
            }
        }
    }
    return;
}

/******************************************************************************/

void outputDensUnitCell(string fileNameString, float densCutoff)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    } else {
        outputFile << writeCryst1() << endl;
    }

    int iatom = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                float orth[3], frac[3];
                float tempFact;
                float occupancy;
                if(densUnitCell[ia][ib][ic] > densCutoff) { //protein region
                    iatom++;
                    if(iatom >= 9999) iatom = 0; //too many atoms

                    frac[0] = (ia + 0.5) / float(ANUM);
                    frac[1] = (ib + 0.5) / float(BNUM);
                    frac[2] = (ic + 0.5) / float(CNUM);
                    convertFracToOrth(frac, orth);

                    tempFact = densUnitCell[ia][ib][ic];
                    occupancy = 1.0;

                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                    << ' ' << " C   ARG A   1    " 
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                    << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                    << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                    << string(11, ' ') << 'C' << endl;
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputResoSphere(string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    
    int iatom = 0;
    for(int h = -HNUM/2; h <= HNUM/2 - 1; h++ ) {
        for(int k = -KNUM/2; k <= KNUM/2 - 1; k++ ) {
            for(int l = -LNUM/2; l <= LNUM/2 - 1; l++ ) {
                float random = rand()/(double)RAND_MAX;
                if(random > 0.1) continue; //keep 10% grid points
                float S = sqrt( h * h * S1 + k * k * S2 + l * l * S3 +
                          2.0 * h * k * S4 + 2.0 * h * l * S5 + 2.0 * k * l * S6 ) / CRYST_VOL;
                float tempFact = 1./S;
                if(tempFact > resoCutoff) tempFact = resoCutoff;
                else tempFact = 0;
                if(iatom > 9999) iatom = 0;
                outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom
                    << ' ' << " H  " << string(14, ' ')
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << (h+0.5)
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << (k+0.5)
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << (l+0.5)
                    << setfill(' ') << setw(6) << fixed << setprecision(2) << 1.0
                    << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact << endl;
                iatom++;
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputWeigAvgDensUnitCell(string fileNameString, float densCutoff)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    } else {
        outputFile << writeCryst1() << endl;
    }

    int iatom = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                float orth[3], frac[3];
                float tempFact;
                float occupancy;
                if(weigAvgDensUnitCell[ia][ib][ic] > densCutoff) {  //protein region
                    iatom++;
                    if(iatom >= 9999) iatom = 0; //too many atoms

                    frac[0] = (ia + 0.5) / float(ANUM);
                    frac[1] = (ib + 0.5) / float(BNUM);
                    frac[2] = (ic + 0.5) / float(CNUM);
                    convertFracToOrth(frac, orth);

                    tempFact = weigAvgDensUnitCell[ia][ib][ic];
                    occupancy = 1.0;

                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                    << ' ' << " C   ARG A   1    "
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                    << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                    << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                    << string(11, ' ') << 'C' << endl;
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputMaskUnitCell(string fileNameString, float cutoff)

{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    } else {
        outputFile << writeCryst1() << endl;
    }

    int iatom = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                float orth[3], frac[3];
                float tempFact;
                float occupancy;
                if(maskUnitCell[ia][ib][ic] > cutoff) {  //mask region
                    iatom++;
                    if(iatom >= 9999) iatom = 0; //too many atoms

                    frac[0] = (ia + 0.5) / float(ANUM);
                    frac[1] = (ib + 0.5) / float(BNUM);
                    frac[2] = (ic + 0.5) / float(CNUM);
                    convertFracToOrth(frac, orth);

                    tempFact = maskUnitCell[ia][ib][ic];
                    occupancy = 1.0;

                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                    << ' ' << " C   ARG A   1    "
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                    << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                    << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                    << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                    << string(11, ' ') << 'C' << endl;
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputConvMaskUnitCell(string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    } else {
        outputFile << writeCryst1() << endl;
    }

    int iatom = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                if(realConvMaskUnitCell[ia][ib][ic] > 0){
                    float frac[3], orth[3], occupancy, tempFact;
                    iatom++;
                    if(iatom >= 9999) iatom = 0; //too many atoms
                    frac[0] = (ia + 0.5) / float(ANUM);
                    frac[1] = (ib + 0.5) / float(BNUM);
                    frac[2] = (ic + 0.5) / float(CNUM);
                    convertFracToOrth(frac, orth);
                    occupancy = 1.0;
                    tempFact = realConvMaskUnitCell[ia][ib][ic];
                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                        << ' ' << " C   ARG A   1    "
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                        << string(11, ' ') << 'C' << endl;
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputConvDensUnitCell(string fileNameString, float densCutoff)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    } else {
        outputFile << writeCryst1() << endl;
    }

    float densAvg = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                densAvg += convDensUnitCell[ia][ib][ic];
            }
        }
    }
    densAvg /= ANUM*BNUM*CNUM;


    int iatom = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                if(convDensUnitCell[ia][ib][ic] > densCutoff){
                    float frac[3], orth[3], occupancy, tempFact;
                    iatom++;
                    if(iatom >= 9999) iatom = 0; //too many atoms
                    frac[0] = (ia + 0.5) / float(ANUM);
                    frac[1] = (ib + 0.5) / float(BNUM);
                    frac[2] = (ic + 0.5) / float(CNUM);
                    convertFracToOrth(frac, orth);
                    occupancy = 1.0;
                    tempFact = convDensUnitCell[ia][ib][ic]/densAvg*100.0;
                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                        << ' ' << " C   ARG A   1    "
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                        << string(11, ' ') << 'C' << endl;
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputNcsCentMassCandUnitCell(string fileNameString, float densCutoff)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    } else {
        outputFile << writeCryst1() << endl;
    }

    float densAvg = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                densAvg += ncsCentMassCandUnitCell[ia][ib][ic];
            }
        }
    }
    densAvg /= ANUM*BNUM*CNUM;


    int iatom = 0;
    for (int ia = 0; ia <= ANUM - 1; ++ia) {
        for (int ib = 0; ib <= BNUM  - 1; ++ib) {
            for (int ic = 0; ic <= CNUM - 1; ++ic) {
                if(ncsCentMassCandUnitCell[ia][ib][ic] > densCutoff){
                    float frac[3], orth[3], occupancy, tempFact;
                    iatom++;
                    if(iatom >= 9999) iatom = 0; //too many atoms
                    frac[0] = (ia + 0.5) / float(ANUM);
                    frac[1] = (ib + 0.5) / float(BNUM);
                    frac[2] = (ic + 0.5) / float(CNUM);
                    convertFracToOrth(frac, orth);
                    occupancy = 1.0;
                    tempFact = ncsCentMassCandUnitCell[ia][ib][ic]/densAvg*100.0;
                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                        << ' ' << " C   ARG A   1    "
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                        << string(11, ' ') << 'C' << endl;
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputF000(string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int iter = 0; iter <= numIter - 1; iter++) {
        outputFile << setfill(' ') << setw(8) << iter << ','
                << setfill(' ') << setw(11) << setprecision(1) << fixed << F000[iter] << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputScaleFact(string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int iter = 0; iter <= numIter - 1; iter++) {
        outputFile << setfill(' ') << setw(8) << iter << ','
                << setfill(' ') << setw(8) << setprecision(3) << fixed << scaleFact[iter] << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputRvalue(string fileNameString, float *Rvalue)

{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int iter = 0; iter <= numIter - 1; iter++) {
        outputFile << setfill(' ') << setw(8) << iter << ','
                << setfill(' ') << setw(8) << setprecision(3) << fixed << Rvalue[iter] << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputBoundEachLargeBin(string fileNameString)

{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int ilargeBin = 0; ilargeBin <= numLargeBin; ++ilargeBin){
        outputFile << setfill(' ') << setw(12) << fixed << ilargeBin 
                   << setfill(' ') << setw(15) << setprecision(7) << fixed << boundLargeBin[ilargeBin] << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputLogDensFreqEachLargeBin(string fileNameString)

{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int ilargeBin = 0; ilargeBin <= numLargeBin-1; ++ilargeBin){
        outputFile << setfill(' ') << setw(12) << fixed << ilargeBin << ','
                   << setfill(' ') << setw(12) << fixed << numGridPointLargeBin[ilargeBin] << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputLogHist(string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }

    int i = 0;
    float area = 0;
    float density; 
    float deltaDens; // e/A^3
    deltaDens = float( int(floor(densMin*100)) - int(floor(densMin*100)) % 2 ) / 100.0;
    for(int ismallBin = 0; ismallBin <= numSmallBin-1; ++ismallBin){
        i = i + numGridPointSmallBin[ismallBin];
        density = densMin + float(ismallBin) / float(numSmallBin) * (densMax-densMin);
        if(density > deltaDens || ismallBin == numSmallBin-1) {
            area = area + float(i)*0.02;  //draw N cylinders with width 0.02e/A^3 on output histgram
            i = 0;
            deltaDens = deltaDens + 0.02;
        }
    }
    i = 0;
    deltaDens = float( int(floor(densMin*100)) - int(floor(densMin*100)) % 2 ) / 100.0;
    for(int ismallBin = 0; ismallBin <= numSmallBin-1; ++ismallBin){
        i = i + numGridPointSmallBin[ismallBin];
        density = densMin + float(ismallBin) / float(numSmallBin) * (densMax-densMin);
        if(density > deltaDens || ismallBin == numSmallBin-1) {
            outputFile << setfill(' ') << setw(15) << fixed << setprecision(7)
                << deltaDens - 0.01 << ','
                << setfill(' ') << setw(15) << fixed << setprecision(7)
                << float(i)/area << endl;
                //note: 1/A integral(x*dx) = 1, A = integral(xdx), 
                //x is float(i), dx is (densMax-densMin)/float(100)
                //integral( (x/A)*dx ) = 1, x is float(i), A is area
            i = 0;
            deltaDens = deltaDens + 0.02;
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputLogHistOld(string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }

    int i = 0;
    float area = 0;
    for(int ismallBin = 0; ismallBin <= numSmallBin-1; ++ismallBin){
        i = i + numGridPointSmallBin[ismallBin];
        if((ismallBin % 100000 == 0 && ismallBin > 0) || ismallBin == numSmallBin-1){
            area = area + float(i)*(densMax-densMin)/float(100);  //draw 100 cylinders on output histgram
            i = 0;
        }
    }
    i = 0;
    for(int ismallBin = 0; ismallBin <= numSmallBin-1; ++ismallBin){
        i = i + numGridPointSmallBin[ismallBin];
        if((ismallBin % 100000 == 0 && ismallBin > 0) || ismallBin == numSmallBin-1){
            outputFile << setfill(' ') << setw(15) << fixed << setprecision(7)
                << densMin+float(ismallBin-100000/2)/float(numSmallBin)*(densMax-densMin) << ','
                << setfill(' ') << setw(15) << fixed << setprecision(7)
                << float(i)/area << endl;
                //note: 1/A integral(x*dx) = 1, A = integral(xdx), 
                //x is float(i), dx is (densMax-densMin)/float(100)
                //integral( (x/A)*dx ) = 1, x is float(i), A is area
            i = 0;
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

bool isFileExist(string fileNameString)
{
    ifstream infile(fileNameString.c_str());
    return infile.good();
    //ifstream destructor will be called upon exiting isFileExist and it will close the stream.
}

/******************************************************************************/

int inputLogConvStat(string fileNameString)
{
    ifstream inputFile( fileNameString.c_str() );
    if (!inputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    string lineString;
    int iorigin;
    while(getline(inputFile, lineString)){
        stringstream( lineString.substr(9, 5) ) >> iorigin;
    }
    inputFile.close();
    return iorigin;
}

/******************************************************************************/
