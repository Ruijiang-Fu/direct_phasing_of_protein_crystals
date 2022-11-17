#define IALOOPASU ia=0;ia<ANUM;++ia
#define IBLOOPASU ib=0;ib<BNUM;++ib
#define ICLOOPASU ic=0;ic<CNUM/8;++ic

#define ASUCOND (x>=0)&&(x<1.0)&&(y>=0)&&(y<1.0)&&(z>=0)&&(z<0.125)

#define NOTUNIQREFLCOND (h<k||(h==0&&k==0&&l%4!=0))

#define IHLOOPUNIQREFL ih=0;ih<HNUM/2;++ih
#define IKLOOPUNIQREFL ik=0;ik<KNUM/2;++ik
#define ILLOOPUNIQREFL il=0;il<LNUM/2;++il


/*******************************************************************************

                            declare variable

*******************************************************************************/

const int numSymmOper = 8; //number of symmetry transformations
const int numOrigTran = 4; //number of allowed origin translations, enantiomorph P4122
const int numGridPointAsu = ANUM*BNUM*CNUM/numSymmOper;

//variables of asymmetric unit
float densAsu[ANUM][BNUM][CNUM/8] = {}; //dens for previous iteration
float densAsuResv[ANUM][BNUM][CNUM/8] = {}; //dens for previous iteration
float densAsuTemp[ANUM][BNUM][CNUM/8] = {}; //dens for current iteration
float densAsuSum[ANUM][BNUM][CNUM/8] = {}; //sum of dens for several successful runs
float densAsuAvg[ANUM][BNUM][CNUM/8] = {}; //avg of dens for several successful runs
float variDensAsuTemp[ANUM][BNUM][CNUM/8] = {}; //variance of density
float weigAvgDensAsu[ANUM][BNUM][CNUM/8] = {};
int protMaskAsu[ANUM][BNUM][CNUM/8] = {};
int protMaskAsuTemp[ANUM][BNUM][CNUM/8] = {};
int protMaskAsuAvg[ANUM][BNUM][CNUM/8] = {};
int maskAsu[ANUM][BNUM][CNUM/8] = {};

float listGridPointAsu[numGridPointAsu][4] = {};
float listGridPointNcsMask[numGridPointAsu][5] = {};
float listGridPointNcsMaskSurf[numGridPointAsu][5] = {};

//variables of unique structure factor
int uniqStruFactStat[HNUM/2][KNUM/2][LNUM/2] = {};
float uniqStruFactReso[HNUM/2][KNUM/2][LNUM/2] = {};
float uniqStruFactAmplObs[HNUM/2][KNUM/2][LNUM/2] = {};
float uniqStruFactAmplObsOrig[HNUM/2][KNUM/2][LNUM/2] = {};
float uniqStruFactAmplObsSigma[HNUM/2][KNUM/2][LNUM/2] = {};
float uniqStruFactAmplCal[HNUM/2][KNUM/2][LNUM/2] = {};
float uniqStruFactPhase[HNUM/2][KNUM/2][LNUM/2] = {};
float missStruFactAmpl[HNUM/2][KNUM/2][LNUM/2] = {};

//metrics of convergence
const int numResoShell = 9;
float resoShell[numResoShell] = {resoCutoff, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0};
float truePhaseOrigTran[HNUM/2][KNUM/2][LNUM/2][numOrigTran] = {};
int trueProtMaskOrigTran[ANUM][BNUM][CNUM/8][numOrigTran] = {};
float RworkResoShell[numIter][numResoShell] = {};
float RfreeResoShell[numIter][numResoShell] = {};
float meanPhaseErrorOrigTranResoShell[numIter][numResoShell][numOrigTran] = {};
float meanPhaseErrorOrigTranResoSphere[numIter][numResoShell][numOrigTran] = {};
float corrCoefOrigTran[numIter][numResoShell][numOrigTran] = {};
float protMaskMatchOrigTran[numIter][numOrigTran] = {};

//crystallographic symmetry transformations
const float symmOper[numSymmOper][3][4] =
{  { { 1.0, 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0, 0.0 } },
   { {-1.0, 0.0, 0.0, 0.0 }, { 0.0,-1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0, 0.5 } },
   { { 0.0,-1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0, 0.75} },
   { { 0.0, 1.0, 0.0, 0.0 }, {-1.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0, 0.25} },
   { {-1.0, 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0, 0.0 }, { 0.0, 0.0,-1.0, 0.0 } },
   { { 1.0, 0.0, 0.0, 0.0 }, { 0.0,-1.0, 0.0, 0.0 }, { 0.0, 0.0,-1.0, 0.5 } },
   { { 0.0, 1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0,-1.0, 0.25} },
   { { 0.0,-1.0, 0.0, 0.0 }, {-1.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0,-1.0, 0.75} }  };

//allowed origin translations
const float origTran[numOrigTran][3][4]=
{  { { 1.0, 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0, 0.0 } },
   { { 1.0, 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0, 0.5 } },
   { { 1.0, 0.0, 0.0, 0.5 }, { 0.0, 1.0, 0.0, 0.5 }, { 0.0, 0.0, 1.0, 0.0 } },
   { { 1.0, 0.0, 0.0, 0.5 }, { 0.0, 1.0, 0.0, 0.5 }, { 0.0, 0.0, 1.0, 0.5 } }   };


/*******************************************************************************

                            declare function

*******************************************************************************/

//cell parameters in CRYST1 of pdb file
string writeCryst1();

//compute resolution of each unique structure factor
void computeUniqStruFactReso();

//sort and output unique structure factor amplitudes observed in experiment
void sortAndOutputUniqStruFactAmplObs(string fileNameString);

//generate crystallographic equivalent molecules in unit cell
void generateCrysEquiMoleInUnitCell(string fileNameString, string *atomStringArray);

//generate unit cell for all allowed origin translations
void generateUnitCellForAllOrigTran(string fileNameString, string *atomStringArray);
void generateAllOrigTranForCoot(string fileNameString, string *atomStringArray);

//input files
void inputNcsMaskFromTxt(string fileNameString);
void inputProtMaskAsu(string fileNameString);
void inputprotMaskAsuTempFromPdb(string fileNameString);
void inputDensAsuFromPdb(string fileNameString, float percent);
void inputUniqStruFactAmplObs(string fileNameString);
void inputFmodel(string fileNameString);
void inputUniqStruFactAmplCalAndPhaseFromCif(string fileNameString);
void inputUniqStruFactAmplCalAndPhaseFromTxt(string fileNameString);

void generateRandDensAsu(unsigned long seed);
void generateRandDensAsuInProtMask(unsigned long seed);
void addRandDensAsuInProtMask(unsigned long seed, float randDensMax);
void generateProtMaskAsuFromNcsMask();

//vary solvent content
float varySolvContLinearly(int iter);
float shakeSolvCont(int iter);

//vary sigmWeigAvg
float varySigmWeigAvgLinearly(int iter);
float varySigmWeigAvgNonLinearly(int iter);

//vary sigmWeigObs
float varySigmWeigObsLinearly(int iter);
float varySigmWeigObsNonLinearly(int iter);

//vary weight of observed data
void varyWeightUniqStruFactAmplObs(float sigmWeigObs);

//operate densAsu
void reserveDensAsu();
void restoreDensAsu();

//apply symmetric operations to get densUnitCell from densAsu
void applySymmOperDensAsu();
void applySymmOperWeigAvgDensAsu();
void applySymmOperVariDensAsuTemp();

//adjust fractional coordinate into the default asu of a specific spacegroup
void adjustFracCoorIntoAsu(float *frac);

//get unique calsf and phase from all calsf and phase
void catchUniqStruFactAmplCalAndPhase();

//calculate R value, meanPhaseError, corrCoef, protMaskMatch
void scaleUniqStruFactAmplCal();
float computeScaleFact();
float computeRvalue(int obsStatus);
void computeRvalueResoShell(float *RvalueResoShell, int obsStatus0);
void computeMeanPhaseErrorResoSphere(float **meanPhaseError);
void computeMeanPhaseErrorResoShell(float **meanPhaseError);
void computeCorrCoefInSphere(float **corrCoef);
void computeCorrCoefInShell(float **corrCoef);
void computeProtMaskMatch(float *protMaskMatch);

//fill in missing reflections
void fillMissRefl(int iter);
//clear missing reflections
void clearMissRefl(int iter);
//flipping phases of missing reflections
void flipPhaseMissRefl(int iter);
//clear high-reso reflections
void flipPhaseHighResoRefl(int iter);

//generate whole reciprocal space sf from unique sf
void computeAllStruFactAmplFromUniqStruFactAmplObs();
void computeAllStruFactAmplFromUniqStruFactAmplCal();
void computeAllStruFactPhaseFromUniqStruFactPhase();

//get densAsu and make protMaskAsu
void catchDensAsuFromDensUnitCell();
void catchDensAsuTempFromDensUnitCell();
void computeVariDensAsuTemp();
void catchWeigAvgDensAsuFromWeigAvgDensUnitCell();
float computeWeigAvgDensCutoffAsu(float protCont);
void makeProtMaskAsu(float densCutoff);
void makeProtMaskAsuAvg(float densCutoff);

//compute reference histogram
float findActuSolvCont();
void findDensAsuMinMaxAvg(string region);
void findWeigAvgDensAsuMinMaxAvg(string region);
float excludeExorDensAsu(int numGridPoint);
void adjustExorDensAsu(float densMin0, float densMax0);
void assignDensAsuIntoSmallBin(string region);
void assignWeigAvgDensAsuIntoSmallBin(string region);
void combineSmallBinIntoLargeBin(string region);

//compute reference histogram
void computeRefeHist(int iter);

//compute maximum-pack sphere radius in unit cell
void computeMaxPackSphereRadiusUnitCell(int ia0, float maxOverlapAllow);
void computeMaxPackSphereRadiusUnitCellForMonomer(int ia0, float maxOverlapAllow);
void computeMaxPackCylinderRadiusUnitCell(int ia0, float maxOverlapAllow);

//ncs averaging
void allocateNcsMask(int iter);
void outputNcsCore(string fileNameString);
void outputNcsMask(string fileNameString);
void seperateNcsMask(int iter);
void outputSepeNcsMask(string fileNameString);
void outputProtMask(string fileNameString);

//optimize NCS axis by minimizing NCS-density difference
void minimizeNcsDensDiffNcsAvgAxis();

//density modification
//proper(closed) ncs, e.g., rotation
void applyNcsAvgMatrix();
void applyNcsAvgAxis();
//improper(open) ncs, e.g., rotation+translation
void applyImprNcsAvgMatrix();
void applyImprNcsAvgAxis();
void applyHIO(int iter);
void applyCHIO(int iter);
void applyHistMatch();
void applySolvFlat();

//true phases if available
void generateTruePhaseForAllOrigTran(string outputStatus, string fileNameString);

//true prot mask if available
void generateTrueProtMaskForAllOrigTran();

//output Rvalue in resolution shell
void outputRvalueResoShell(string fileNameString, float **RvalueResoShell);

//output mean phase error if available
void outputMeanPhaseErrorResoShell(string fileNameString, int iorigin0);
void outputMeanPhaseErrorResoSphere(string fileNameString, int iorigin0);

//output correlation coefficient if available
void outputCorrCoef(string fileNameString, int iorig0);

//output protein mask match if available
void outputProtMaskMatch(string fileNameString, int iorigin0);

void outputUniqStruFactPhaseAndMissAmpl(string fileNameString);

//display calculated result in Coot
void outputCifFileForCoot(string fileNameString);

//output final density
void outputDensAsu(string fileNameString, float densCutoff);
void outputDensAsuTxt(string fileNameString);
void outputWeigAvgDensAsu(string fileNameString, float densCutoff);
void outputProtMaskAsuTxt(string fileNameString);
void outputRvalueInResoShell(string fileNameString);

//check convergence
int checkOrigChoiceByProtMaskMatch(float *protMaskMatch);
int checkOrigChoiceByMeanPhaseError(float **meanPhaseError);
int checkConvByMeanPhaseError(float **meanPhaseError);

//analyze results
void translateOrigDensAsu(int iorigin0);
void translateOrigProtMaskAsu(int iorigin0);
void accumulateDensAsu();
void computeAvgDensAsu(int num);
void computeDensAsuAvg(int num);
void translateOrigDensAsuMinDiff();
void translateOrigProtMaskAsuMinDiff();

//cut and output single-molecule mask from calculated protein mask in unitcell
void cutAndOutputSingMoleMask(float **atomCoorArray, string fileNameString);


/*******************************************************************************

                            define function

*******************************************************************************/


//cell parameters in CRYST1 of pdb file
string writeCryst1()
{
    stringstream cryst1;
    //--- cryst1 will be used in pdb file --------------------------------------
    cryst1 << "CRYST1" << setfill(' ') << setw(9) << fixed << setprecision(3) << CRYST_A
        << setfill(' ') << setw(9) << fixed << setprecision(3) << CRYST_B
        << setfill(' ') << setw(9) << fixed << setprecision(3) << CRYST_C
        << setfill(' ') << setw(7) << fixed << setprecision(2) << CRYST_ALPHA/PI*180
        << setfill(' ') << setw(7) << fixed << setprecision(2) << CRYST_BETA/PI*180
        << setfill(' ') << setw(7) << fixed << setprecision(2) << CRYST_GAMMA/PI*180
        << " P 43 2 2    " << numSymmOper;
    return cryst1.str();
}

/******************************************************************************/

//compute resolution of each unique structure factor
void computeUniqStruFactReso()

{
    int irefl = 0;
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                if(h == 0 && k == 0 && l == 0) {
                    uniqStruFactReso[ih][ik][il] = 999.0; //F000
                } else if( NOTUNIQREFLCOND ) {
                    uniqStruFactReso[ih][ik][il] = 0; //missing orders of diffraction due to interference
                } else {
                    float S=sqrt( h * h * S1 + k * k * S2 + l * l * S3 +
                        2.0 * h * k * S4 + 2.0 * h * l * S5 + 2.0 * k * l * S6 ) / CRYST_VOL;
                    float resolution = 1.0 / S ;
                    uniqStruFactReso[ih][ik][il] = resolution;
                    if(resolution >= resoCutoff) {
                        irefl++;
                    }
                }
            }
        }
    }
    numUniqStruFact = irefl;
    return;
    //note: zero the resolutions of F000, missing orders of diffraction due to interference, high resolution reflections
    //As a result, all zero resolution will be excluded from phasing process. 
}

/******************************************************************************/

void sortAndOutputUniqStruFactAmplObs(string fileNameString)

{
    int obsStatus;
    float obssf, obssfSigma;

    //count the number of unique structure factors
    int iuniq = 0;
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                if(h == 0 && k ==0 && l ==0)continue;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution >= resoCutoff) {
                    iuniq++;
                }
            }
        }
    }
    int numUniqStruFact = iuniq;

    //sort unique structure factor in descending order of its resolution
    float uniqStruFact[numUniqStruFact][7];
    iuniq = 0;
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                if(h == 0 && k ==0 && l ==0)continue;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution >= resoCutoff) {
                    int obsStatus = uniqStruFactStat[ih][ik][il];
                    float obssf = uniqStruFactAmplObs[ih][ik][il];
                    float obssfSigma = uniqStruFactAmplObsSigma[ih][ik][il];
                    uniqStruFact[iuniq][0] = h;
                    uniqStruFact[iuniq][1] = k;
                    uniqStruFact[iuniq][2] = l;
                    uniqStruFact[iuniq][3] = obsStatus;
                    uniqStruFact[iuniq][4] = obssf;
                    uniqStruFact[iuniq][5] = obssfSigma;
                    uniqStruFact[iuniq][6] = resolution;
                    iuniq++;
                }
            }
        }
    }
    for(int i = 0; i < numUniqStruFact - 1; ++i){
        for(int j = i + 1; j <= numUniqStruFact - 1; ++j){
            if(uniqStruFact[j][6] > uniqStruFact[i][6]){
                for(int k = 0; k <= 6; k++){
                    float temp = uniqStruFact[i][k];
                    uniqStruFact[i][k] = uniqStruFact[j][k];
                    uniqStruFact[j][k] = temp;
                }

            }
        }
    }

    //output unique structure factor in descending order of its resolution
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()){
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int iuniq = 0; iuniq <= numUniqStruFact - 1; iuniq++){
        int h = int(uniqStruFact[iuniq][0]);
        int k = int(uniqStruFact[iuniq][1]);
        int l = int(uniqStruFact[iuniq][2]);
        int obsStatus = int(uniqStruFact[iuniq][3]);
        float obssf = uniqStruFact[iuniq][4];
        float obssfSigma = uniqStruFact[iuniq][5];
        float resolution = uniqStruFact[iuniq][6];
        outputFile  << setfill(' ') << setw(5) << h
                    << setfill(' ') << setw(5) << k
                    << setfill(' ') << setw(5) << l
                    << setfill(' ') << setw(2) << obsStatus
                    << setfill(' ') << setw(9) << setprecision(1) << fixed << obssf
                    << setfill(' ') << setw(7) << setprecision(1) << fixed << obssfSigma
                    << setfill(' ') << setw(8) << setprecision(3) << fixed << resolution << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void generateCrysEquiMoleInUnitCell(string fileNameString, string *atomStringArray)
{
    float rota[3][3], tran[3];
    float orth[3], frac[3], temp[3];

    for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
        //catch symmetry operation matrix
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = symmOper[isymm][i][j]; //load rotation matrix
            }
            tran[i] = symmOper[isymm][i][3]; //load translation vector
        }
        //apply symmetry operation
        string *atomStringArrayTemp = new string[numAtomAsu]; //temparary
        for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom) {
            string lineString = atomStringArray[iatom];
            stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
            //converse orth to frac coordinate
            convertOrthToFrac(orth, temp);
            for( int i=0; i<=2 ; ++i )
            {
                frac[i] = 0;
                for( int j=0; j<=2 ; ++j )
                {
                    frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                }
                frac[i] = frac[i] + tran[i]; //translation
            }
            //converse frac to orth coordinate
            convertFracToOrth(frac, orth);
            stringstream lineStream;
            lineStream << lineString.substr(0, 6)
                << setfill(' ') << setw(5) << iatom
                << lineString.substr(11, 19)
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[0]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[1]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[2]
                << lineString.substr(54, 26);
            atomStringArrayTemp[iatom] = lineStream.str();
        }

        //translate the center of the new molecule into unit cell
        float centOrth[3];
        computeCentOrthModel(atomStringArrayTemp, centOrth);
        float centFrac[3];
        convertOrthToFrac(centOrth, centFrac);
        for(int i = 0; i <= 2; ++i){
            while(centFrac[i] < 0 || centFrac[i] >= 1.0){
                float translation = 0;
                if(centFrac[i] < 0) translation = 1.0;
                if(centFrac[i] >= 1.0) translation = -1.0;
                centFrac[i] = centFrac[i] + translation;
                for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom){
                    string lineString = atomStringArrayTemp[iatom];
                    stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
                    convertOrthToFrac(orth, frac);
                    frac[i] = frac[i] + translation;
                    convertFracToOrth(frac, orth);
                    stringstream lineStream;
                    lineStream << lineString.substr(0, 6)
                        << setfill(' ') << setw(5) << iatom
                        << lineString.substr(11, 19)
                        << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[0]
                        << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[1]
                        << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[2]
                        << lineString.substr(54, 26);
                    atomStringArrayTemp[iatom] = lineStream.str();
                }
            }
        }
        //output pdb file
        fileNameStream.str( string() );
        fileNameStream.clear();
        fileNameStream << fileNameString << isymm << ".pdb";
        ofstream outputFile( fileNameStream.str().c_str() );
        if(!outputFile.is_open()){
            cout << "Can't open file " << fileNameStream.str() << endl;
            exit(EXIT_FAILURE);
        }
        //write in cell paremeters in pdb header
        outputFile << writeCryst1() << endl;
        for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom) {
            outputFile << atomStringArrayTemp[iatom] << endl;
        }
        outputFile.close();    
    }

    return;
}

/******************************************************************************/

void generateUnitCellForAllOrigTran(string fileNameString, string *atomStringArray)

{

    string *atomStringArrayUnitCell = new string[ numSymmOper * numAtomAsu ];

    float rota[3][3], tran[3];
    float orth[3], frac[3], temp[3];
    //apply symmtry operation to get the first unit cell
    for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = symmOper[isymm][i][j]; //load rotation matrix
            }
            tran[i] = symmOper[isymm][i][3]; //load translation vector
        }
        for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom) {
            string lineString = atomStringArray[iatom];
            stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
            //converse orth to frac coordinate
            convertOrthToFrac(orth, temp);
            for( int i=0; i<=2 ; ++i )
            {
                frac[i] = 0;
                for( int j=0; j<=2 ; ++j )
                {
                    frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                }
                frac[i] = frac[i] + tran[i]; //translation
            }
            adjustFracCoorIntoUnitCell(frac);
            //converse frac to orth coordinate
            convertFracToOrth(frac, orth);
            stringstream lineStream;
            lineStream << lineString.substr(0, 6)
                << setfill(' ') << setw(5) << (isymm * numAtomAsu + iatom)
                << lineString.substr(11, 10)
                << char('A'+isymm)
                << lineString.substr(22, 8)
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[0]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[1]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[2]
                << lineString.substr(54, 26);
            atomStringArrayUnitCell[isymm * numAtomAsu + iatom] = lineStream.str();
        }
    }

    //apply allowed origin translation to unit cell for each origin translation
    for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
        fileNameStream.str( string() );
        fileNameStream.clear();
        fileNameStream << fileNameString << iorig << ".pdb";
        ofstream outputFile( fileNameStream.str().c_str() );
        if(!outputFile.is_open()){
            cout << "Can't open file " << fileNameStream.str() << endl;
            exit(EXIT_FAILURE);
        }
        //write in cell paremeters in pdb header
        outputFile << writeCryst1P1() << endl;
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = origTran[iorig][i][j]; //load rotation matrix
            }
            tran[i] = origTran[iorig][i][3]; //load rotation matrix
        }
        for(int iatom = 0; iatom <= numSymmOper * numAtomAsu - 1; ++iatom) {
            string lineString = atomStringArrayUnitCell[iatom];
            stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
            //converse orth to frac coordinate
            convertOrthToFrac(orth, temp);
            for( int i=0; i<=2 ; ++i )
            {
                frac[i] = 0;
                for( int j=0; j<=2 ; ++j )
                {
                    frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                }
                frac[i] = frac[i] + tran[i]; //translation
            }
            adjustFracCoorIntoUnitCell(frac);
            //converse frac to orth coordinate
            convertFracToOrth(frac, orth);
            stringstream lineStream;
            lineStream << lineString.substr(0, 6)
                << setfill(' ') << setw(5) << iatom
                << lineString.substr(11, 19)
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[0]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[1]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[2]
                << lineString.substr(54, 26);
            outputFile << lineStream.str() << endl;
        }
        outputFile.close();    
    }
    return;
}

/******************************************************************************/

void generateAllOrigTranForCoot(string fileNameString, string *atomStringArray)

{
    float rota[3][3], tran[3];
    float orth[3], frac[3], temp[3];

    for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
        //catch allowed origin translation matrix
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = origTran[iorig][i][j]; //load rotation matrix
            }
            tran[i] = origTran[iorig][i][3]; //load translation vector
        }
        //apply allowed origin translation
        string *atomStringArrayTemp = new string[numAtomAsu]; //temparary
        for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom) {
            string lineString = atomStringArray[iatom];
            stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
            //converse orth to frac coordinate
            convertOrthToFrac(orth, temp);
            for( int i=0; i<=2 ; ++i )
            {
                frac[i] = 0;
                for( int j=0; j<=2 ; ++j )
                {
                    frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                }
                frac[i] = frac[i] + tran[i]; //translation
            }
            //converse frac to orth coordinate
            convertFracToOrth(frac, orth);
            stringstream lineStream;
            lineStream << lineString.substr(0, 6)
                << setfill(' ') << setw(5) << iatom
                << lineString.substr(11, 19)
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[0]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[1]
                << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[2]
                << lineString.substr(54, 26);
            atomStringArrayTemp[iatom] = lineStream.str();
        }

        //translate the center of the new molecule into unit cell
        float centOrth[3];
        computeCentOrthModel(atomStringArrayTemp, centOrth);
        float centFrac[3];
        convertOrthToFrac(centOrth, centFrac);
        for(int i = 0; i <= 2; ++i){
            while(centFrac[i] < 0 || centFrac[i] >= 1.0){
                float translation = 0;
                if(centFrac[i] < 0) translation = 1.0;
                if(centFrac[i] >= 1.0) translation = -1.0;
                centFrac[i] = centFrac[i] + translation;
                for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom){
                    string lineString = atomStringArrayTemp[iatom];
                    stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
                    convertOrthToFrac(orth, frac);
                    frac[i] = frac[i] + translation;
                    convertFracToOrth(frac, orth);
                    stringstream lineStream;
                    lineStream << lineString.substr(0, 6)
                        << setfill(' ') << setw(5) << iatom
                        << lineString.substr(11, 19)
                        << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[0]
                        << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[1]
                        << setfill(' ') << setw(8) << setprecision(3) << fixed << orth[2]
                        << lineString.substr(54, 26);
                    atomStringArrayTemp[iatom] = lineStream.str();
                }
            }
        }
        //output pdb file
        fileNameStream.str( string() );
        fileNameStream.clear();
        fileNameStream << fileNameString << iorig << ".pdb";
        ofstream outputFile( fileNameStream.str().c_str() );
        if(!outputFile.is_open()){
            cout << "Can't open file " << fileNameStream.str() << endl;
            exit(EXIT_FAILURE);
        }
        //write in cell paremeters in pdb header
        outputFile << writeCryst1() << endl;
        for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom) {
            outputFile << atomStringArrayTemp[iatom] << endl;
        }
        outputFile.close();    
    }

    return;
}

/******************************************************************************/

void inputNcsMaskFromTxt(string fileNameString)
{
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    anum0 = 1000;
    bnum0 = 1000;
    cnum0 = 1000;
    anum1 = -1000;
    bnum1 = -1000;
    cnum1 = -1000;
    int ia, ib, ic;
    while(inputFile >> ia >> ib >> ic){
        if(ia < anum0) anum0 = ia;
        if(ib < bnum0) bnum0 = ib;
        if(ic < cnum0) cnum0 = ic;
        if(ia > anum1) anum1 = ia;
        if(ib > bnum1) bnum1 = ib;
        if(ic > cnum1) cnum1 = ic;
    }
    //settle ncs mask into a little bit loose qubic region
    anum0 = anum0 - 3;
    anum1 = anum1 + 3;
    bnum0 = bnum0 - 3;
    bnum1 = bnum1 + 3;
    cnum0 = cnum0 - 3;
    cnum1 = cnum1 + 3;

    inputFile.clear();
    inputFile.seekg(0, ios::beg);

    while(inputFile >> ia >> ib >> ic){
        ia = ia - anum0;
        ib = ib - bnum0;
        ic = ic - cnum0;
        ncsMask[ia][ib][ic] = 1;
    }

    inputFile.close();
    return;
}

/******************************************************************************/

void inputProtMaskAsu(string fileNameString)

{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    protMaskAsu[ia][ib][ic] = 0;
                }
            }
        }
    }
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    inputFile >> protMaskAsu[ia][ib][ic];
                }
            }
        }
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void inputprotMaskAsuTempFromPdb(string fileNameString)
{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    protMaskAsuTemp[ia][ib][ic] = 0;
                }
            }
        }
    }
    //read in density in asu from deposited pdb file, file name is "PDB_CODE".pdb
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    string lineString;
    float orth[3], frac[3];
    int ia, ib, ic;
    while(getline(inputFile, lineString)) {
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0 ){
            stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
               //converse orth to frac coordinate
            convertOrthToFrac(orth, frac);
            //when frac is not inside default asu, move it into default asu
            adjustFracCoorIntoAsu(frac);
            //each atom is represented by density 1.0
            ia = floor(frac[0] * ANUM);
            ib = floor(frac[1] * BNUM);
            ic = floor(frac[2] * CNUM);
            protMaskAsuTemp[ia][ib][ic] = 1;
        }
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void inputDensAsuFromPdb(string fileNameString, float percent)

{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsu[ia][ib][ic] = 0;
                }
            }
        }
    }
    //read in density in asu from deposited pdb file, file name is "PDB_CODE".pdb
    ifstream inputFile( fileNameString.c_str() );
    if(!inputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    string lineString;
    float orth[3], frac[3];
    int ia, ib, ic;
    float random;
    while(getline(inputFile, lineString)) {
        if(lineString.compare(0, 4, "ATOM") == 0 || lineString.compare(0, 6, "HETATM") == 0 ){

            random = ( rand() % 100 ) * 0.01 ; //0 to 0.99
            if(random < percent){

                stringstream( lineString.substr(30, 24) ) >> orth[0] >> orth[1] >> orth[2];
                //converse orth to frac coordinate
                convertOrthToFrac(orth, frac);
                //when frac is not inside default asu, move it into default asu
                adjustFracCoorIntoAsu(frac);
                //each atom is represented by density 1.0
                ia = floor(frac[0] * ANUM);
                ib = floor(frac[1] * BNUM);
                ic = floor(frac[2] * CNUM);
                densAsu[ia][ib][ic] = 1.0;
            }
        }
        
    }
    inputFile.close();

    return;
}


/******************************************************************************/

void inputUniqStruFactAmplObs(string fileNameString)

{
    //read in observerd diffraction data in experiment, file name is "PDB_CODE"_uniq_sf.txt
    ifstream inputFile( fileNameString.c_str() );
    if (!inputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }

    int h, k, l, obsStatus;
    float obssf, obssfSigma, resolution;
    //data format (3%5d, %2d, %9.1f, %7.1, %8.3f) for h, k, l, obsStatus, obssf, obssfSigma, resolution
    while (inputFile >> h >> k >> l >> obsStatus >> obssf >> obssfSigma >> resolution) {
        if(resolution >= resoCutoff) { 
            if(resolution >= resoCutoffLow) obsStatus = 0;
            uniqStruFactStat[h][k][l] = obsStatus;  
            uniqStruFactAmplObs[h][k][l] = obssf;
            uniqStruFactAmplObsOrig[h][k][l] = obssf;
            uniqStruFactAmplObsSigma[h][k][l] = obssfSigma;
        }
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void inputFmodel(string fileNameString)
{
    ifstream inputFile( fileNameString.c_str() );
    if (!inputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    string lineString;
    while (getline(inputFile, lineString)) {
        if(lineString.compare(0, 6, " 1 1 1") == 0){
            int h, k, l, obsStatus;
            float obssf, phase;
            stringstream(lineString.substr(9, 15)) >> h >> k >> l;
            stringstream(lineString.substr(27, 23)) >> obssf >> phase;
            float S, resolution;
            S=sqrt( h * h * S1 + k * k * S2 + l * l * S3 +
                2.0 * h * k * S4 + 2.0 * h * l * S5 + 2.0 * k * l * S6 ) / CRYST_VOL;
            resolution = 1.0 / S ;
            phase = phase / 180.0 * PI;
            if(resolution >= resoCutoff) {  
                uniqStruFactAmplCal[h][k][l] = obssf;
                uniqStruFactPhase[h][k][l] = phase;
            }
        }
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void inputUniqStruFactAmplCalAndPhaseFromCif(string fileNameString)
{
    ifstream inputFile( fileNameString.c_str() );
    if (!inputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    string lineString;
    while (getline(inputFile, lineString)) {
        if(lineString.compare(0, 5, "1 1 1") == 0){
            int h, k, l, obsStatus;
            float obssf, calsf, phase;
            const char *status = lineString.substr(6, 1).c_str();
            stringstream(lineString.substr(7, 15)) >> h >> k >> l;
            stringstream(lineString.substr(22, 9)) >> obssf;
            stringstream(lineString.substr(38, 11)) >> calsf;
            stringstream(lineString.substr(55, 8)) >> phase;
            float S, resolution;
            S=sqrt( h * h * S1 + k * k * S2 + l * l * S3 +
                2.0 * h * k * S4 + 2.0 * h * l * S5 + 2.0 * k * l * S6 ) / CRYST_VOL;
            resolution = 1.0 / S ;
            phase = phase / 180.0 * PI;
            if(resolution >= resoCutoff) {  
                uniqStruFactAmplCal[h][k][l] = calsf;
                uniqStruFactPhase[h][k][l] = phase;
                if(status[0] == 'x' || status[0] == 'f') missStruFactAmpl[h][k][l] = obssf;
            }
        }
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void inputUniqStruFactAmplCalAndPhaseFromTxt(string fileNameString)

{
    //read in calculated structure factor phases and missing amplitudes
    ifstream inputFile( fileNameString.c_str() );
    if (!inputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    //data format (3%5d, %2d, %11.1f, %8.1f) for h, k, l, obsStatus, calsf, phase
    int h, k, l, obsStatus;
    float calsf, phase;
    while (inputFile >> h >> k >> l >> obsStatus >> calsf >> phase) {
        phase = phase / 180.0 * PI;
        uniqStruFactAmplCal[h][k][l] = calsf;
        uniqStruFactPhase[h][k][l] = phase;
        if(obsStatus == 0 || obsStatus == 2) missStruFactAmpl[h][k][l] = calsf;
    }
    inputFile.close();
    return;
}

/******************************************************************************/

void generateRandDensAsu(unsigned long seed)

{
    srand(seed); //initialize random seed

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsu[ia][ib][ic] = ( rand() % 10000 ) * 0.0001;
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void addRandDensAsuInProtMask(unsigned long seed, float randDensMax)

{
    srand(seed); //initialize random seed

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(protMaskAsu[ia][ib][ic]==1) {
                        densAsu[ia][ib][ic] += randDensMax * ( rand() % 10000 ) * 0.0001;
                    }
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void generateRandDensAsuInProtMask(unsigned long seed)

{
    srand(seed); //initialize random seed

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(protMaskAsu[ia][ib][ic]==1) {
                        densAsu[ia][ib][ic] = ( rand() % 10000 ) * 0.0001;
                    }
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void generateProtMaskAsuFromNcsMask()
{
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;

    int ia1, ib1, ic1;
    float frac[3];
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    protMaskAsu[ia1][ib1][ic1] = 1;
                }
            }
        }
    }

    return;
}

/******************************************************************************/

float varySolvContLinearly(int iter)
{
    if(iter <= numIterVarySolvCont) {
        solvCont = solvContIni + (solvContFina - solvContIni) * float(iter) / float(numIterVarySolvCont);
    } else{
        solvCont = solvContFina;
    }

    return solvCont;
}

/******************************************************************************/

float shakeSolvCont(int iter)
{
    if(iter <= numIterVarySolvCont) {
        if(solvContFina > solvContIni){
            solvCont += ((iter % 3) - 1)* 0.005 * float(numIterVarySolvCont - iter) / float(numIterVarySolvCont); 
        } else { // -1%, 0, +1%
            solvCont = solvContIni + ((iter % 3) - 1)* 0.005 * float(numIterVarySolvCont - iter) / float(numIterVarySolvCont); 
        }
    } else{
        solvCont = solvContFina;
    }

    return solvCont;
}

/******************************************************************************/

float varySigmWeigAvgLinearly(int iter)
{
    float sigmWeigAvg;
    if(iter <= numIterVarySigmWeigAvg) {
        sigmWeigAvg = sigmWeigAvgIni + (sigmWeigAvgFina - sigmWeigAvgIni) * float(iter) / float(numIterVarySigmWeigAvg);
    } else{
        sigmWeigAvg = sigmWeigAvgFina;
    }
    return sigmWeigAvg;
}

/******************************************************************************/

float varySigmWeigAvgNonLinearly(int iter)
{
    sigmWeigAvg = sigmWeigAvgArray[iter]; //pre-computed value
    return sigmWeigAvg;
}

/******************************************************************************/

float varySigmWeigObsLinearly(int iter)
{
    if(iter <= numIterVarySigmWeigObs) {
        sigmWeigObs = sigmWeigObsIni * (1.0 - float(iter) / float(numIterVarySigmWeigObs) );
    } else{
        sigmWeigObs = 0;
    }
    return sigmWeigObs;
}

/******************************************************************************/

float varySigmWeigObsNonLinearly(int iter)
{
    sigmWeigObs = sigmWeigObsArray[iter]; //pre-computed value
    return sigmWeigObs;
}

/******************************************************************************/

void varyWeightUniqStruFactAmplObs(float sigmWeigObs)
{
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                float S = sqrt( h * h * S1 + k * k * S2 + l * l * S3 +
                    2.0 * h * k * S4 + 2.0 * h * l * S5 + 2.0 * k * l * S6 ) / CRYST_VOL;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution >= resoCutoff) { 
                    uniqStruFactAmplObs[ih][ik][il] = uniqStruFactAmplObsOrig[ih][ik][il] * exp( - 2.0 * pow(PI*sigmWeigObs*S, 2.0) );
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void reserveDensAsu()
{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsuResv[ia][ib][ic] = densAsu[ia][ib][ic];
                    densAsu[ia][ib][ic] = densAsuTemp[ia][ib][ic];
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void restoreDensAsu()
{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsu[ia][ib][ic] = densAsuResv[ia][ib][ic];
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void applySymmOperDensAsu()

{
    float rota[3][3], tran[3];
    float frac[3], temp[3];

    for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = symmOper[isymm][i][j]; //load rotation matrix
            }
            tran[i] = symmOper[isymm][i][3];    //load rotation matrix
        }
                
        for( int IALOOPASU ) {
            float x = (ia + 0.5) / float(ANUM);
            for( int IBLOOPASU ) {
                float y = (ib + 0.5) / float(BNUM);
                for( int ICLOOPASU ) {
                    float z = (ic + 0.5) / float(CNUM);
                    if( ASUCOND ) {
                        temp[0] = (ia + 0.5) / float(ANUM);
                        temp[1] = (ib + 0.5) / float(BNUM);
                        temp[2] = (ic + 0.5) / float(CNUM);
                        for( int i=0; i<=2 ; ++i )
                        {
                            frac[i] = 0;
                            for( int j=0; j<=2 ; ++j )
                            {
                                frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                            }
                            frac[i] = frac[i] + tran[i]; //translation
                        }
                        adjustFracCoorIntoUnitCell(frac);
                        int ia0 = floor(frac[0] * ANUM); //round down value (interger)
                         int ib0 = floor(frac[1] * BNUM);
                         int ic0 = floor(frac[2] * CNUM);
                        densUnitCell[ia0][ib0][ic0] = densAsu[ia][ib][ic];
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void applySymmOperWeigAvgDensAsu()

{
    float rota[3][3], tran[3];
    float frac[3], temp[3];

    for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = symmOper[isymm][i][j]; //load rotation matrix
            }
            tran[i] = symmOper[isymm][i][3];    //load rotation matrix
        }
        for( int IALOOPASU ) {
            float x = (ia + 0.5) / float(ANUM);
            for( int IBLOOPASU ) {
                float y = (ib + 0.5) / float(BNUM);
                for( int ICLOOPASU ) {
                    float z = (ic + 0.5) / float(CNUM);
                    if( ASUCOND ) {
                        temp[0] = (ia + 0.5) / float(ANUM);
                        temp[1] = (ib + 0.5) / float(BNUM);
                        temp[2] = (ic + 0.5) / float(CNUM);
                        for( int i=0; i<=2 ; ++i )
                        {
                            frac[i] = 0;
                            for( int j=0; j<=2 ; ++j )
                            {
                                frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                            }
                            frac[i] = frac[i] + tran[i]; //translation
                        }
                        adjustFracCoorIntoUnitCell(frac);
                        int ia0 = floor(frac[0] * ANUM); //round down value (interger)
                        int ib0 = floor(frac[1] * BNUM);
                        int ic0 = floor(frac[2] * CNUM);
                        weigAvgDensUnitCell[ia0][ib0][ic0] = weigAvgDensAsu[ia][ib][ic];
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void applySymmOperVariDensAsuTemp()
{
    float rota[3][3], tran[3];
    float frac[3], temp[3];

    for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = symmOper[isymm][i][j]; //load rotation matrix
            }
            tran[i] = symmOper[isymm][i][3]; //load rotation matrix
        }
        for( int IALOOPASU ) {
            float x = (ia + 0.5) / float(ANUM);
            for( int IBLOOPASU ) {
                float y = (ib + 0.5) / float(BNUM);
                for( int ICLOOPASU ) {
                    float z = (ic + 0.5) / float(CNUM);
                    if( ASUCOND ) {
                        temp[0] = (ia + 0.5) / float(ANUM);
                        temp[1] = (ib + 0.5) / float(BNUM);
                        temp[2] = (ic + 0.5) / float(CNUM);
                        for( int i=0; i<=2 ; ++i )
                        {
                            frac[i] = 0;
                            for( int j=0; j<=2 ; ++j )
                            {
                                frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                            }
                            frac[i] = frac[i] + tran[i]; //translation
                        }
                        adjustFracCoorIntoUnitCell(frac);
                        int ia0 = floor(frac[0] * ANUM); //round down value (interger)
                        int ib0 = floor(frac[1] * BNUM);
                        int ic0 = floor(frac[2] * CNUM);
                        densUnitCell[ia0][ib0][ic0] = variDensAsuTemp[ia][ib][ic];
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void adjustFracCoorIntoAsu(float *frac) 

{
    //when frac is outside the unit cell, translate it into unit cell
    adjustFracCoorIntoUnitCell(frac);

    //if frac is already inside the default asu, return frac
    float x = frac[0]; 
    float y = frac[1];
    float z = frac[2];
    if( ASUCOND )
        return;

    float symmFrac[numSymmOper][3];

    //when frac is not inside the default asu, find all of its symmetric copies and store them in symmFrac
    for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
        float rota[3][3], tran[3], temp[3];
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = symmOper[isymm][i][j]; //pick up rotation matrix
            }
            tran[i] = symmOper[isymm][i][3]; //pick up transaltion vector
        }
        for( int i=0; i<=2 ; ++i )
        {
            temp[i] = 0;
            for( int j=0; j<=2 ; ++j )
            {
                temp[i] = temp[i] + rota[i][j] * frac[j]; //rotate
            }
            temp[i] = temp[i] + tran[i]; //translate
        }

        adjustFracCoorIntoUnitCell(temp);

        for(int i = 0; i <= 2 ; ++i)
        {
            symmFrac[isymm][i] = temp[i];
        }
    }

    //select the one and the only one which is inside the default asymmetric unit
    int flag = 0;
    for(int isymm = 0; isymm <= numSymmOper - 1 ; ++isymm) {
        float x = symmFrac[isymm][0];
        float y = symmFrac[isymm][1];
        float z = symmFrac[isymm][2];
        if( ASUCOND) {
            for(int i = 0; i <=2 ; ++i) {
                frac[i] = symmFrac[isymm][i];
            }
            flag++;
        }
    }
    
    if(flag != 1){
        cout << "Error, not only one in adjustFracCoorIntoAsu" << endl;
        cout << flag << endl;
        exit (EXIT_FAILURE);
    }

    return;
}

/******************************************************************************/

void catchUniqStruFactAmplCalAndPhase()

{
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                uniqStruFactAmplCal[ih][ik][il] = allStruFactAmpl[ih][ik][il];
                uniqStruFactPhase[ih][ik][il] = allStruFactPhase[ih][ik][il];
            }
        }
    }

    return;

}

/******************************************************************************/

void scaleUniqStruFactAmplCal()
{
    float sumFobs = 0;
    float sumFcal = 0;
    float sumFcalTimeFobs=0;
    float sumFcalSqu=0;
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution < resoCutoff) continue;
                int obsStatus = uniqStruFactStat[ih][ik][il];
                float obssf = uniqStruFactAmplObs[ih][ik][il];
                float calsf = uniqStruFactAmplCal[ih][ik][il];
                if(obsStatus == 1) {
                    sumFobs = sumFobs + obssf;
                    sumFcal = sumFcal + calsf;
                    sumFcalTimeFobs = sumFcalTimeFobs + calsf * obssf;
                    sumFcalSqu = sumFcalSqu + pow(calsf, 2.0);
                }
            }
        }
    }
    float scaleFact = sumFcalTimeFobs / sumFcalSqu;

    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution < resoCutoff) continue;
                if(h == 0 && k == 0 && l == 0) continue; //keep calculated F000
                float calsf = uniqStruFactAmplCal[ih][ik][il];
                uniqStruFactAmplCal[ih][ik][il] = scaleFact * calsf;
            }
        }
    }

    return;
}

/******************************************************************************/

float computeScaleFact()
{
    float sumFcalTimeFobs=0;
    float sumFcalSqu=0;
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution < resoCutoff) continue;
                int obsStatus = uniqStruFactStat[ih][ik][il];
                float obssf = uniqStruFactAmplObsOrig[ih][ik][il];
                float calsf = uniqStruFactAmplCal[ih][ik][il];
                if(obsStatus == 1) {
                    sumFcalTimeFobs = sumFcalTimeFobs + calsf * obssf;
                    sumFcalSqu = sumFcalSqu + pow(calsf, 2.0);
                }
            }
        }
    }
    float scaleFact = sumFcalTimeFobs / sumFcalSqu;
    return scaleFact;
}

/******************************************************************************/

float computeRvalue(int obsStatus0)

{
    float sumFobs = 0;
    float sumFcal = 0;
    float sumFcalTimeFobs=0;
    float sumFcalSqu=0;
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution < resoCutoff) continue;
                int obsStatus = uniqStruFactStat[ih][ik][il];
                float obssf = uniqStruFactAmplObsOrig[ih][ik][il];
                float calsf = uniqStruFactAmplCal[ih][ik][il];
                if(obsStatus == obsStatus0) {
                    sumFobs = sumFobs + obssf;
                    sumFcal = sumFcal + calsf;
                    sumFcalTimeFobs = sumFcalTimeFobs + calsf * obssf;
                    sumFcalSqu = sumFcalSqu + pow(calsf, 2.0);
                }
            }
        }
    }
    float scaleFact = sumFcalTimeFobs / sumFcalSqu;

    float Rvalue = 0;
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution < resoCutoff) continue;
                int obsStatus = uniqStruFactStat[ih][ik][il];
                float obssf = uniqStruFactAmplObsOrig[ih][ik][il];
                float calsf = uniqStruFactAmplCal[ih][ik][il];
                if(obsStatus == obsStatus0) {
                    Rvalue = Rvalue + abs( obssf - scaleFact * calsf );
                }
            }
        }
    }
    Rvalue = Rvalue / sumFobs;
    return Rvalue;
}

/******************************************************************************/

void computeRvalueResoShell(float RvalueResoShell[numResoShell], int obsStatus0)
{
    float sumFobs = 0;
    float sumFcal = 0;
    float sumFcalTimeFobs=0;
    float sumFcalSqu=0;
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution < resoCutoff) continue;
                int obsStatus = uniqStruFactStat[ih][ik][il];
                float obssf = uniqStruFactAmplObsOrig[ih][ik][il];
                float calsf = uniqStruFactAmplCal[ih][ik][il];
                if(obsStatus == obsStatus0) {
                    sumFobs = sumFobs + obssf;
                    sumFcal = sumFcal + calsf;
                    sumFcalTimeFobs = sumFcalTimeFobs + calsf * obssf;
                    sumFcalSqu = sumFcalSqu + pow(calsf, 2.0);
                }
            }
        }
    }
    float scaleFact = sumFcalTimeFobs / sumFcalSqu;

    for(int ishell = 0; ishell <= numResoShell - 1; ++ishell){
        RvalueResoShell[ishell] = 0;
        for( int IHLOOPUNIQREFL ) {
            for( int IKLOOPUNIQREFL ) {
                for( int ILLOOPUNIQREFL ) {
                    int h = ih;
                    int k = ik;
                    int l = il;
                    float resolution = uniqStruFactReso[ih][ik][il];
                    if(resolution < resoCutoff) continue;
                    if( (ishell < numResoShell - 1 && resolution >= resoShell[ishell] && resolution < resoShell[ishell+1])
                        || (ishell == numResoShell - 1 && resolution >= resoShell[ishell]) ){
                        int obsStatus = uniqStruFactStat[ih][ik][il];
                        float obssf = uniqStruFactAmplObsOrig[ih][ik][il];
                        float calsf = uniqStruFactAmplCal[ih][ik][il];
                        if(obsStatus == obsStatus0) {
                            RvalueResoShell[ishell] = RvalueResoShell[ishell] + abs( obssf - scaleFact * calsf );
                        }
                    }
                }
            }
        }
        RvalueResoShell[ishell] = RvalueResoShell[ishell] / sumFobs;
    }
    return;
}

/******************************************************************************/

void computeMeanPhaseErrorResoSphere(float meanPhaseError[numResoShell][numOrigTran])

{
    for(int ishell = 0; ishell <= numResoShell - 1; ++ishell){
        for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
            meanPhaseError[ishell][iorig] = 0;
        }
        int i = 0;
        for( int IHLOOPUNIQREFL ) {
            for( int IKLOOPUNIQREFL ) {
                for( int ILLOOPUNIQREFL ) {
                    int h = ih;
                    int k = ik;
                    int l = il;
                    float resolution = uniqStruFactReso[ih][ik][il];
                    if(resolution < resoCutoff) continue;
                    if(resolution >= resoShell[ishell]){
                        int obsStatus = uniqStruFactStat[ih][ik][il];
                        float phase = uniqStruFactPhase[ih][ik][il];
                        if(obsStatus == 1) {
                            i++;
                            for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
                                meanPhaseError[ishell][iorig] += acos(cos(truePhaseOrigTran[ih][ik][il][iorig] - phase));
                            }
                        }
                    }
                }
            }
        }        
        for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
            meanPhaseError[ishell][iorig] = meanPhaseError[ishell][iorig] / float(i) / PI * 180.0;
        }
    }
    return;
}

/******************************************************************************/

void computeMeanPhaseErrorResoShell(float meanPhaseError[numResoShell][numOrigTran])

{
    for(int ishell = 0; ishell <= numResoShell - 1; ++ishell){
        for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
            meanPhaseError[ishell][iorig] = 0;
        }
        int i = 0;
        for( int IHLOOPUNIQREFL ) {
            for( int IKLOOPUNIQREFL ) {
                for( int ILLOOPUNIQREFL ) {
                    int h = ih;
                    int k = ik;
                    int l = il;
                    float resolution = uniqStruFactReso[ih][ik][il];
                    if(resolution < resoCutoff) continue;
                    if( (ishell < numResoShell - 1 && resolution >= resoShell[ishell] && resolution < resoShell[ishell+1])
                        || (ishell == numResoShell - 1 && resolution >= resoShell[ishell]) ){
                        int obsStatus = uniqStruFactStat[ih][ik][il];
                        float phase = uniqStruFactPhase[ih][ik][il];
                        if(obsStatus == 1) {
                            i++;
                            for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
                                meanPhaseError[ishell][iorig] += acos(cos(truePhaseOrigTran[ih][ik][il][iorig] - phase));
                            }
                        }
                    }
                }
            }
        }
        for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
            meanPhaseError[ishell][iorig] = meanPhaseError[ishell][iorig] / float(i) / PI * 180.0;
        }
    }
    return;
}

/******************************************************************************/

void computeCorrCoefInSphere(float corrCoef[numResoShell][numOrigTran])

{
    for(int ishell = 0; ishell <= numResoShell - 1; ++ishell){
        for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
            corrCoef[ishell][iorig] = 0;
        }
        float sumFobsSqu = 0;
        float sumFcalSqu = 0;
        float sumFobsTimeFcalTimeCos[numOrigTran] = {};
        for( int IHLOOPUNIQREFL ) {
            for( int IKLOOPUNIQREFL ) {
                for( int ILLOOPUNIQREFL ) {
                    int h = ih;
                    int k = ik;
                    int l = il;
                    float resolution = uniqStruFactReso[ih][ik][il];
                    if(resolution < resoCutoff) continue;
                    if(resolution >= resoShell[ishell]){
                        int obsStatus = uniqStruFactStat[ih][ik][il];
                        float obssf = uniqStruFactAmplObsOrig[ih][ik][il];
                        float calsf = uniqStruFactAmplCal[ih][ik][il];
                        float phase = uniqStruFactPhase[ih][ik][il];
                        if(obsStatus == 1) {
                            sumFobsSqu += pow(obssf, 2.0);
                            sumFcalSqu += pow(calsf, 2.0);
                            for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
                                sumFobsTimeFcalTimeCos[iorig] += obssf * calsf * cos(truePhaseOrigTran[ih][ik][il][iorig] - phase);
                            }
                        }
                    }
                }
            }
        }
        for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
            corrCoef[ishell][iorig] = sumFobsTimeFcalTimeCos[iorig] / sqrt(sumFobsSqu * sumFcalSqu);
        }
    }
    return;
}

/******************************************************************************/

void computeCorrCoefInShell(float corrCoef[numResoShell][numOrigTran])

{
    for(int ishell = 0; ishell <= numResoShell - 1; ++ishell){
        for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
            corrCoef[ishell][iorig] = 0;
        }
        float sumFobsSqu = 0;
        float sumFcalSqu = 0;
        float sumFobsTimeFcalTimeCos[numOrigTran] = {};
        for( int IHLOOPUNIQREFL ) {
            for( int IKLOOPUNIQREFL ) {
                for( int ILLOOPUNIQREFL ) {
                    int h = ih;
                    int k = ik;
                    int l = il;
                    float resolution = uniqStruFactReso[ih][ik][il];
                    if(resolution < resoCutoff) continue;
                    if( (ishell < numResoShell - 1 && resolution >= resoShell[ishell] && resolution < resoShell[ishell+1])
                        || (ishell == numResoShell - 1 && resolution >= resoShell[ishell]) ){
                        int obsStatus = uniqStruFactStat[ih][ik][il];
                        float obssf = uniqStruFactAmplObsOrig[ih][ik][il];
                        float calsf = uniqStruFactAmplCal[ih][ik][il];
                        float phase = uniqStruFactPhase[ih][ik][il];
                        if(obsStatus == 1) {
                            sumFobsSqu += pow(obssf, 2.0);
                            sumFcalSqu += pow(calsf, 2.0);
                            for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
                                sumFobsTimeFcalTimeCos[iorig] += obssf * calsf * cos(truePhaseOrigTran[ih][ik][il][iorig] - phase);
                            }
                        }
                    }
                }
            }
        }
        for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
            corrCoef[ishell][iorig] = sumFobsTimeFcalTimeCos[iorig] / sqrt(sumFobsSqu * sumFcalSqu);
        }
    }
    return;
}

/******************************************************************************/

void computeProtMaskMatch(float protMaskMatch[numOrigTran])

{
    for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
        protMaskMatch[iorig] = 0;
    }
    int i = 0;
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(protMaskAsu[ia][ib][ic] == 0){
                        continue;
                    } else{
                        i++;
                        for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
                            if( trueProtMaskOrigTran[ia][ib][ic][iorig] == 1) protMaskMatch[iorig]++;
                        }
                    }
                }
            }
        }
    }

    for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {
        protMaskMatch[iorig] = protMaskMatch[iorig] / float(i);
    }
    return;
}

/******************************************************************************/

void fillMissRefl(int iter)

{
    int h, k, l, obsStatus;
    float obssf, calsf;
    float S, resolution;

    float sumFobs = 0;
    float sumFcal = 0;
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                h = ih;
                k = ik;
                l = il;
                resolution = uniqStruFactReso[ih][ik][il];
                if(resolution < resoCutoff) continue;
                obsStatus = uniqStruFactStat[ih][ik][il];
                obssf = uniqStruFactAmplObs[ih][ik][il];
                calsf = uniqStruFactAmplCal[ih][ik][il];
                if(obsStatus == 1) { //work data set
                    sumFobs = sumFobs + obssf;
                    sumFcal = sumFcal + calsf;
                }
            }
        }
    }

    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                h = ih;
                k = ik;
                l = il;
                if(iter > numIter - numIterSolvFlat){
                    resolution = uniqStruFactReso[ih][ik][il];
                    if(resolution < resoCutoff){
                        missStruFactAmpl[ih][ik][il] = 0;
                        continue;
                    }
                }
                if( NOTUNIQREFLCOND ) continue;
                obsStatus = uniqStruFactStat[ih][ik][il];
                calsf = uniqStruFactAmplCal[ih][ik][il];
                if(obsStatus != 1) { //F000, missing data set, and test data set
                    missStruFactAmpl[ih][ik][il] = calsf * (sumFobs / sumFcal);
                }
            }
        }
    }
    F000[iter] = uniqStruFactAmplCal[0][0][0];
    missStruFactAmpl[0][0][0] = uniqStruFactAmplCal[0][0][0]; //special F000
    return;

}

/******************************************************************************/

void clearMissRefl(int iter)

{
    int h, k, l, obsStatus;
    float obssf, calsf;
    float S, resolution;

    if(iter > numIter - numIterSolvFlat) return;

    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                h = ih;
                k = ik;
                l = il;
                if( NOTUNIQREFLCOND  ) continue;
                if(h == 0 && k == 0 && l == 0) continue; //F000
                obsStatus = uniqStruFactStat[ih][ik][il];
                if(obsStatus != 1) { //missing data set, test data set, and extra high-resolution refls involved in FFT
                    missStruFactAmpl[ih][ik][il] = missStruFactAmpl[ih][ik][il] * (float(iter) / float(numIter - numIterSolvFlat));
                }
            }
        }
    }

    return;

}

/******************************************************************************/

void flipPhaseMissRefl(int iter)

{
    int h, k, l, obsStatus;
    float obssf, calsf;
    float S, resolution;

    if(iter > numIter - numIterSolvFlat) return;

    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                h = ih;
                k = ik;
                l = il;
                if( NOTUNIQREFLCOND  ) continue;
                if(h == 0 && k == 0 && l == 0) continue; //F000
                obsStatus = uniqStruFactStat[ih][ik][il];
                if(obsStatus != 1) { //missing data set, test data set, and extra high-resolution refls involved in FFT
                    uniqStruFactPhase[ih][ik][il] = -uniqStruFactPhase[ih][ik][il];
                }
            }
        }
    }

    return;

}

/******************************************************************************/

void flipPhaseHighResoRefl(int iter)

{
    int h, k, l, obsStatus;
    float obssf, calsf;
    float S, resolution;

    if(iter > numIter - numIterSolvFlat) return;

    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                h = ih;
                k = ik;
                l = il;
                if( NOTUNIQREFLCOND  ) continue;
                if(h == 0 && k == 0 && l == 0) continue; //F000
                resolution = sigmWeigObs * PI / sqrt(log(2.0)/2.0);
                if(resolution < uniqStruFactReso[ih][ik][il]) {
                    uniqStruFactPhase[ih][ik][il] = -uniqStruFactPhase[ih][ik][il];
                }
            }
        }
    }

    return;

}

/******************************************************************************/

void computeAllStruFactAmplFromUniqStruFactAmplObs()

{
    int h, k, l, obsStatus;
    float obssf;
    float rota[3][3], tran[3];
    float primitive_hkl[3];
    float equivalent_hkl[3];
    
    for( int ih=0; ih<HNUM; ++ih ) {
        for( int ik=0; ik<KNUM; ++ik ) {
            for( int il=0; il<LNUM; ++il ) {
                allStruFactAmpl[ih][ik][il] = 0; //clear previous value
            }
        }
    }

    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                h = ih;
                k = ik;
                l = il;
                if( NOTUNIQREFLCOND ||
                    uniqStruFactReso[ih][ik][il] < resoCutoff  ) {
                    continue;  //missing diffraction due to interference
                }
                obsStatus = uniqStruFactStat[ih][ik][il];
                if(obsStatus == 1) {
                    obssf = uniqStruFactAmplObs[ih][ik][il]; //observed data
                } else {
                    obssf = missStruFactAmpl[ih][ik][il]; //missing and test data
                }

                primitive_hkl[0] = h;
                primitive_hkl[1] = k;
                primitive_hkl[2] = l;

                for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
                    for(int i = 0; i <= 2; ++i)
                    {
                        for(int j = 0; j <= 2; ++j) 
                        {
                            rota[i][j] = symmOper[isymm][i][j]; //load rotation matrix
                        }
                        tran[i] = symmOper[isymm][i][3]; //load rotation matrix
                    }

                    for( int i = 0; i <= 2 ; ++i )
                    {
                        equivalent_hkl[i] = 0;
                        for( int j = 0; j <= 2 ; ++j )
                        {                        //equivalent_hkl=matmul(transpose(rotation),primitive_hkl);
                            equivalent_hkl[i] = equivalent_hkl[i] + rota[j][i] * primitive_hkl[j]; //rotation
                        }
                    }

                    for(int i = 1; i>= -1; i = i-2) { //native F(h,k,l) i=1, Fresnel pair F(-h,-k,-l) i=-1;

                        int h0 = i * round(equivalent_hkl[0]);
                        int k0 = i * round(equivalent_hkl[1]);
                        int l0 = i * round(equivalent_hkl[2]);

                        int ih0, ik0, il0;

                        if(h0 >= 0 && h0 <= HNUM/2 - 1 ) {
                            ih0 = h0;
                        } else if(h0 < 0 && h0 >= -HNUM/2 + 1) {
                            ih0 = h0 + HNUM;
                        } else {
                            cout << "Wrong! uniqStruFact to allStruFact, h0 = " << h0 
                                << ">=" << -HNUM/2 + 1 << "<=" << HNUM/2 -1 << endl;
                            exit(EXIT_FAILURE);
                        }
                        if(k0 >= 0 && k0 <= KNUM/2 - 1) {
                            ik0 = k0;
                        } else if(k0 < 0 && k0 >= -KNUM/2 + 1) {
                            ik0 = k0 + KNUM;
                        } else {
                            cout << "Wrong! uniqStruFact to allStruFact, k0 = " << k0 
                                << ">=" << -KNUM/2 + 1 << "<=" << KNUM/2 -1 << endl;
                            exit(EXIT_FAILURE);
                        }
                        if(l0 >= 0 && l0 <= LNUM/2 - 1) {
                            il0 = l0;
                        } else if(l0 < 0 && l0 >= -LNUM/2 + 1) {
                            il0 = l0 + LNUM;
                        } else {
                            cout << "Wrong! uniqStruFact to allStruFact, l0 = " << l0 
                                << ">=" << -LNUM/2 + 1 << "<=" << LNUM/2 -1 << endl;
                            exit(EXIT_FAILURE);
                        }
                        allStruFactAmpl[ih0][ik0][il0] = obssf;
                    }
                }
            }
        }
    }

    return;
}

/******************************************************************************/


void computeAllStruFactAmplFromUniqStruFactAmplCal()

{
    int h, k, l, obsStatus;
    float calsf;
    float rota[3][3], tran[3];
    float primitive_hkl[3];
    float equivalent_hkl[3];

    for( int ih=0; ih<HNUM; ++ih ) {
        for( int ik=0; ik<KNUM; ++ik ) {
            for( int il=0; il<LNUM; ++il ) {
                allStruFactAmpl[ih][ik][il] = 0; //clear previous value
            }
        }
    }

    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                h = ih;
                k = ik;
                l = il;
                if( NOTUNIQREFLCOND ||
                    uniqStruFactReso[ih][ik][il] < resoCutoff  ) {
                    continue;  //missing diffraction due to interference
                }
                calsf = uniqStruFactAmplCal[ih][ik][il]; //calculated synthetic data

                primitive_hkl[0] = h;
                primitive_hkl[1] = k;
                primitive_hkl[2] = l;

                for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
                    for(int i = 0; i <= 2; ++i)
                    {
                        for(int j = 0; j <= 2; ++j) 
                        {
                            rota[i][j] = symmOper[isymm][i][j]; //load rotation matrix
                        }
                        tran[i] = symmOper[isymm][i][3]; //load rotation matrix
                    }

                    for( int i = 0; i <= 2 ; ++i )
                    {
                        equivalent_hkl[i] = 0;
                        for( int j = 0; j <= 2 ; ++j )
                        {                        //equivalent_hkl=matmul(transpose(rotation),primitive_hkl);
                            equivalent_hkl[i] = equivalent_hkl[i] + rota[j][i] * primitive_hkl[j]; //rotation
                        }
                    }

                    for(int i = 1; i>= -1; i = i-2) { //native F(h,k,l) i=1, Fresnel pair F(-h,-k,-l) i=-1;

                        int h0 = i * round(equivalent_hkl[0]);
                        int k0 = i * round(equivalent_hkl[1]);
                        int l0 = i * round(equivalent_hkl[2]);

                        int ih0, ik0, il0;

                        if(h0 >= 0 && h0 <= HNUM/2 - 1 ) {
                            ih0 = h0;
                        } else if(h0 < 0 && h0 >= -HNUM/2 + 1) {
                            ih0 = h0 + HNUM;
                        } else {
                            cout << "Wrong! uniqStruFact to allStruFact, h0 = " << h0 
                                << ">=" << -HNUM/2 + 1 << "<=" << HNUM/2 -1 << endl;
                            exit(EXIT_FAILURE);
                        }
                        if(k0 >= 0 && k0 <= KNUM/2 - 1) {
                            ik0 = k0;
                        } else if(k0 < 0 && k0 >= -KNUM/2 + 1) {
                            ik0 = k0 + KNUM;
                        } else {
                            cout << "Wrong! uniqStruFact to allStruFact, k0 = " << k0 
                                << ">=" << -KNUM/2 + 1 << "<=" << KNUM/2 -1 << endl;
                            exit(EXIT_FAILURE);
                        }
                        if(l0 >= 0 && l0 <= LNUM/2 - 1) {
                            il0 = l0;
                        } else if(l0 < 0 && l0 >= -LNUM/2 + 1) {
                            il0 = l0 + LNUM;
                        } else {
                            cout << "Wrong! uniqStruFact to allStruFact, l0 = " << l0 
                                << ">=" << -LNUM/2 + 1 << "<=" << LNUM/2 -1 << endl;
                            exit(EXIT_FAILURE);
                        }
                        allStruFactAmpl[ih0][ik0][il0] = calsf;
                    }
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void computeAllStruFactPhaseFromUniqStruFactPhase()

{
    int h, k, l, obsStatus;
    float phase;
    float rota[3][3], tran[3];
    float primitive_hkl[3];
    float equivalent_hkl[3];

    for( int ih=0; ih<HNUM; ++ih ) {
        for( int ik=0; ik<KNUM; ++ik ) {
            for( int il=0; il<LNUM; ++il ) {
                allStruFactPhase[ih][ik][il] = 0; //clear previous value
            }
        }
    }

    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                h = ih;
                k = ik;
                l = il;
                if(    NOTUNIQREFLCOND  ||
                    uniqStruFactReso[ih][ik][il] < resoCutoff ) {
                    continue;  //missing diffraction due to interference
                }
                phase = uniqStruFactPhase[ih][ik][il];

                primitive_hkl[0] = h;
                primitive_hkl[1] = k;
                primitive_hkl[2] = l;

                for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
                    for(int i = 0; i <= 2; ++i)
                    {
                        for(int j = 0; j <= 2; ++j) 
                        {
                            rota[i][j] = symmOper[isymm][i][j]; //load rotation matrix
                        }
                        tran[i] = symmOper[isymm][i][3]; //load rotation matrix
                    }

                    for( int i = 0; i <= 2 ; ++i )
                    {
                        equivalent_hkl[i] = 0;
                        for( int j = 0; j <= 2 ; ++j )
                        {                        //equivalent_hkl=matmul(transpose(rotation),primitive_hkl);
                            equivalent_hkl[i] = equivalent_hkl[i] + rota[j][i] * primitive_hkl[j]; //rotation
                        }
                    }

                    float phase0;

                    for(int i = 1; i >= -1; i = i-2) { //native F(h,k,l) i=1, Fresnel pair F(-h,-k,-l) i=-1;

                        int h0 = i * round(equivalent_hkl[0]);
                        int k0 = i * round(equivalent_hkl[1]);
                        int l0 = i * round(equivalent_hkl[2]);

                        int ih0, ik0, il0;

                        if(h0 >= 0 && h0 <= HNUM/2 - 1 ) {
                            ih0 = h0;
                        } else if(h0 < 0 && h0 >= -HNUM/2 + 1) {
                            ih0 = h0 + HNUM;
                        } else {
                            cout << "Wrong! uniqStruFact to allStruFact, h0 = " << h0 
                                << ">=" << -HNUM/2 + 1 << "<=" << HNUM/2 -1 << endl;
                            exit(EXIT_FAILURE);
                        }
                        if(k0 >= 0 && k0 <= KNUM/2 - 1) {
                            ik0 = k0;
                        } else if(k0 < 0 && k0 >= -KNUM/2 + 1) {
                            ik0 = k0 + KNUM;
                        } else {
                            cout << "Wrong! uniqStruFact to allStruFact, k0 = " << k0 
                                << ">=" << -KNUM/2 + 1 << "<=" << KNUM/2 -1 << endl;
                            exit(EXIT_FAILURE);
                        }
                        if(l0 >= 0 && l0 <= LNUM/2 - 1) {
                            il0 = l0;
                        } else if(l0 < 0 && l0 >= -LNUM/2 + 1) {
                            il0 = l0 + LNUM;
                        } else {
                            cout << "Wrong! uniqStruFact to allStruFact, l0 = " << l0 
                                << ">=" << -LNUM/2 + 1 << "<=" << LNUM/2 -1 << endl;
                            exit(EXIT_FAILURE);
                        }

                        if(i == 1){
                            phase0 = 0;
                            for( int i=0; i<=2 ; ++i )
                            {                    //phase0=phase-2.0*PI*dot_product(tran,primitive_hkl)
                                phase0 = phase0 + tran[i] * primitive_hkl[i];
                            }
                            phase0 = phase - 2.0 * PI * phase0;
                            while (phase0 > 2.0 * PI) phase0 = phase0 - 2.0 * PI;
                            while (phase0 < 0) phase0 = phase0 + 2.0 * PI;
                            allStruFactPhase[ih0][ik0][il0] = phase0;
                        } else if(i == -1){
                            phase0 = 2.0*PI - phase0;
                            while (phase0 > 2.0 * PI) phase0 = phase0 - 2.0 * PI;
                            while (phase0 < 0) phase0 = phase0 + 2.0 * PI;
                            allStruFactPhase[ih0][ik0][il0] = phase0;
                        }
                    }
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void catchDensAsuFromDensUnitCell()

{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsu[ia][ib][ic] = densUnitCell[ia][ib][ic];
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void catchDensAsuTempFromDensUnitCell()

{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsuTemp[ia][ib][ic] = densUnitCell[ia][ib][ic];
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void computeVariDensAsuTemp()
{
    densAvg = 0;
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAvg = densAvg + densAsuTemp[ia][ib][ic];
                }
            }
        }
    }
    densAvg = densAvg / float(numGridPointAsu);
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    variDensAsuTemp[ia][ib][ic] = pow((densAsuTemp[ia][ib][ic] - densAvg), 2.0);
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void catchWeigAvgDensAsuFromWeigAvgDensUnitCell()

{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    weigAvgDensAsu[ia][ib][ic] = weigAvgDensUnitCell[ia][ib][ic];
                }
            }
        }
    }

    return;
}

/******************************************************************************/

float computeWeigAvgDensCutoffAsu(float protCont)
{
    float densCutoff;

    densMin = weigAvgDensAsu[0][0][0];
    densMax = densMin;
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    float density = weigAvgDensAsu[ia][ib][ic];
                    if(density < densMin) densMin = density;
                    if(density > densMax) densMax = density;
                }
            }
        }
    }

    for(int ismallBin = 0; ismallBin <= numSmallBin-1; ++ismallBin)
    {
        numGridPointSmallBin[ismallBin] = 0; //initialize density frequency
    }
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    float density = weigAvgDensAsu[ia][ib][ic];
                    int ismallBin = round( (density-densMin) / (densMax-densMin) * numSmallBin );
                    if(ismallBin < 0) ismallBin = 0;
                    if(ismallBin >= numSmallBin) ismallBin = numSmallBin - 1;
                    numGridPointSmallBin[ismallBin]++;
                }
            }
        }
    }

    int numGridPointHighDens = round(numGridPointAsu * protCont);

    int igrid = 0;
    for(int ismallBin = numSmallBin-1; ismallBin >= 0; --ismallBin){
        igrid = igrid + numGridPointSmallBin[ismallBin];
        if(igrid >= numGridPointHighDens){
            densCutoff = densMin + float(ismallBin) / float(numSmallBin) * (densMax-densMin);
            break;
        }
    }    

    for(int ismallBin = 0; ismallBin <= numSmallBin-1; ++ismallBin)
    {
        numGridPointSmallBin[ismallBin] = 0; //reinitialize
    }

    return densCutoff;

}

/******************************************************************************/

void makeProtMaskAsu(float densCutoff)

{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(weigAvgDensAsu[ia][ib][ic] >= densCutoff){
                        protMaskAsu[ia][ib][ic] = 1;
                    } else {
                        protMaskAsu[ia][ib][ic] = 0;
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void makeProtMaskAsuAvg(float densCutoff)

{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(weigAvgDensAsu[ia][ib][ic] >= densCutoff){
                        protMaskAsuAvg[ia][ib][ic] = 1;
                    } else {
                        protMaskAsuAvg[ia][ib][ic] = 0;
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

float findActuSolvCont()

{
    float solvCont;
    int i = 0;
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(protMaskAsu[ia][ib][ic] == 0) {
                        i++;
                    }
                }
            }
        }
    }
    solvCont = float(i)/float(numGridPointAsu);
    return solvCont;

}

/******************************************************************************/

void findDensAsuMinMaxAvg(string region)
{
    densMin = densAsu[0][0][0];
    densMax = densMin;
    densAvg = 0;
    int i = 0;
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(    ( region.compare("prot") == 0 && protMaskAsu[ia][ib][ic] == 1 ) 
                        || ( region.compare("solv") == 0 && protMaskAsu[ia][ib][ic] == 0 )
                        ||     region.compare("cell") == 0 ) 
                    { 
                        i++;
                        float density = densAsu[ia][ib][ic];
                        if(density < densMin) densMin = density;
                        if(density > densMax) densMax = density;
                        densAvg = densAvg + density;
                    }
                }
            }
        }
    }
    densAvg = densAvg / float(i);
    return;
}

/******************************************************************************/

void findWeigAvgDensAsuMinMaxAvg(string region)
{
    densMin = weigAvgDensAsu[0][0][0];
    densMax = densMin;
    densAvg = 0;
    int i = 0;
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(    ( region.compare("prot") == 0 && protMaskAsu[ia][ib][ic] == 1 ) 
                        || ( region.compare("solv") == 0 && protMaskAsu[ia][ib][ic] == 0 )
                        ||     region.compare("cell") == 0 ) 
                    { 
                        i++;
                        float density = weigAvgDensAsu[ia][ib][ic];
                        if(density < densMin) densMin = density;
                        if(density > densMax) densMax = density;
                        densAvg = densAvg + density;
                    }
                }
            }
        }
    }
    densAvg = densAvg / float(i);
    return;
}

/******************************************************************************/

float excludeExorDensAsu(int numGridPoint)

{
    findDensAsuMinMaxAvg("cell");

    float densCutoff;
    float density1 = densMin;
    float density2 = densMax;
    do {
        densCutoff = (density2 + density1) / 2.0;
        int i = 0;
        for( int IALOOPASU ) {
            float x = (ia + 0.5) / float(ANUM);
            for( int IBLOOPASU ) {
                float y = (ib + 0.5) / float(BNUM);
                for( int ICLOOPASU ) {
                    float z = (ic + 0.5) / float(CNUM);
                    if( ASUCOND ) {
                        if(densAsu[ia][ib][ic] > densCutoff) i++;
                    }
                }
            }
        }
        if(i > numGridPoint) {
            density1 = densCutoff;
        } else if(i < numGridPoint) {
            density2 = densCutoff;
        } else {
            density2 = density1;
        }
    } while ( abs(density2 - density1) > pow(10.0, -5.0));

    return densCutoff;

}

/******************************************************************************/

void adjustExorDensAsu(float densMin0, float densMax0)

{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(densAsu[ia][ib][ic] < densMin0) densAsu[ia][ib][ic] = densMin0;
                    if(densAsu[ia][ib][ic] > densMax0) densAsu[ia][ib][ic] = densMax0;
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void assignDensAsuIntoSmallBin(string region)

{
    findDensAsuMinMaxAvg("prot"); //always in protein region

    for(int ismallBin = 0; ismallBin <= numSmallBin-1; ++ismallBin)
    {
        numGridPointSmallBin[ismallBin] = 0; //initialize density frequency
    }
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(    ( region.compare("prot") == 0 && protMaskAsu[ia][ib][ic] == 1 ) 
                        || ( region.compare("solv") == 0 && protMaskAsu[ia][ib][ic] == 0 )
                        ||   region.compare("cell") == 0 ) 
                    {
                        float density = densAsu[ia][ib][ic];
                        int ismallBin = round( (density-densMin) / (densMax-densMin) * numSmallBin );
                        if(ismallBin < 0) ismallBin = 0;
                        if(ismallBin >= numSmallBin) ismallBin = numSmallBin - 1;
                        numGridPointSmallBin[ismallBin]++;
                    }
                }
            }
        }
    }

    return;

}

/******************************************************************************/


void assignWeigAvgDensAsuIntoSmallBin(string region)

{
    findWeigAvgDensAsuMinMaxAvg("cell"); //always in whole cell

    for(int ismallBin = 0; ismallBin <= numSmallBin-1; ++ismallBin)
    {
        numGridPointSmallBin[ismallBin] = 0; //initialize density frequency
    }
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(    ( region.compare("prot") == 0 && protMaskAsu[ia][ib][ic] == 1 ) 
                        || ( region.compare("solv") == 0 && protMaskAsu[ia][ib][ic] == 0 )
                        ||   region.compare("cell") == 0 ) 
                    {
                        float density = weigAvgDensAsu[ia][ib][ic];
                        int ismallBin = round( (density-densMin) / (densMax-densMin) * numSmallBin );
                        if(ismallBin < 0) ismallBin = 0;
                        if(ismallBin >= numSmallBin) ismallBin = numSmallBin - 1;
                        numGridPointSmallBin[ismallBin]++;
                    }
                }
            }
        }
    }

    return;

}

/******************************************************************************/

void combineSmallBinIntoLargeBin(string region)

{
    findDensAsuMinMaxAvg("prot"); //always in protein region

    int i = 0;
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(    ( region.compare("prot") == 0 && protMaskAsu[ia][ib][ic] == 1 ) 
                        || ( region.compare("solv") == 0 && protMaskAsu[ia][ib][ic] == 0 )
                        ||   region.compare("cell") == 0 ) 
                    {
                        i++;
                    }
                }
            }
        }
    }

    int densFreqEachLargeBin = round( float(i) / float(numLargeBin) ); 
    boundLargeBin[0] = densMin;
    numGridPointLargeBin[numLargeBin-1] = 0; //prepare for the last bin

    int ilargeBin = 0;
    i = 0;
    for(int ismallBin = 0; ismallBin <= numSmallBin-1; ++ismallBin){
        i = i + numGridPointSmallBin[ismallBin];
        if(i >= densFreqEachLargeBin || ismallBin == numSmallBin-1) {
            if(ilargeBin <= numLargeBin-1) {
                numGridPointLargeBin[ilargeBin] = i; 
            } else if(ilargeBin > numLargeBin-1){ 
                numGridPointLargeBin[numLargeBin-1] = numGridPointLargeBin[numLargeBin-1] + i; 
            }
            i = 0;
            ilargeBin++;
            if(ilargeBin <= numLargeBin-1) {
                boundLargeBin[ilargeBin] = boundLargeBin[0] + float(ismallBin) / float(numSmallBin) * (densMax-densMin);
            }
        }
    }
    while(ilargeBin <= numLargeBin-1) {
        boundLargeBin[ilargeBin] = densMax;
        ilargeBin++;
    }
    boundLargeBin[numLargeBin] = densMax;
    return;

}

/******************************************************************************/

void computeRefeHist(int iter)
{
    //reserve current value of densAsu
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsuTemp[ia][ib][ic] = densAsu[ia][ib][ic];
                    protMaskAsuTemp[ia][ib][ic] = protMaskAsu[ia][ib][ic];
                }
            }
        }
    }

    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                uniqStruFactAmplCal[ih][ik][il] = 0;
                uniqStruFactPhase[ih][ik][il] = 0;
            }
        }
    }
    //read in Fatom to make protein mask
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << "_fatom.hkl,.HKL,.sca";
    inputFmodel(fileNameStream.str());

    //compute allStruFactAmpl from uniqStruFactAmplCal, input uniqStruFactAmplCal, output allStruFactAmpl
    computeAllStruFactAmplFromUniqStruFactAmplCal();

    //compute allStruFactPhase from uniqStruFactPhase, input uniqStruFactPhase, output allStruFactPhase
    computeAllStruFactPhaseFromUniqStruFactPhase();

    //prepare for fft
    setHalfGridPhase(); //set half grid phase shift for fft, initialize halfGridPosi and halfGridNega

    //apply forward fft, input allStruFactAmpl and allStruFactPhase, output densUnitCell
    applyFFTStruFactToDens();

    //compute weighted average density in unit cell, input densUnitCell, output weigAvgDensUnitCell
    float sigmWeig = 0.6;
    computeWeigAvgDensUnitCell(sigmWeig);

    //get weigAvgDensAsu from weigAvgDensUnitCell, input weigAvgDensUnitCell, output weigAvgDensAsu
    catchWeigAvgDensAsuFromWeigAvgDensUnitCell();

    //compute a densCutoff value on weigAvgDensAsu, input weigAvgDensAsu, return densCutoff
    densCutoff = computeWeigAvgDensCutoffAsu(1.0-solvCont);

    //make a protein mask in asu, input weigAvgDensAsu, output protMaskAsu
    makeProtMaskAsu(densCutoff);

    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                uniqStruFactAmplCal[ih][ik][il] = 0;
                uniqStruFactPhase[ih][ik][il] = 0;
            }
        }
    }
    //read in Fmodel to compute density
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << "_fmodel.hkl,.HKL,.sca";
    inputFmodel(fileNameStream.str());

    //vary the weight of Fmodel
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                float S = sqrt( h * h * S1 + k * k * S2 + l * l * S3 +
                    2.0 * h * k * S4 + 2.0 * h * l * S5 + 2.0 * k * l * S6 ) / CRYST_VOL;
                float resolution = 1.0 / S;
                if(resolution >= resoCutoff){
                    uniqStruFactAmplCal[ih][ik][il] = uniqStruFactAmplCal[ih][ik][il] * exp( - 2.0 * pow(PI*sigmWeigObs*S, 2.0) );
                }
            }
        }
    }

    //compute allStruFactAmpl from uniqStruFactAmplCal, input uniqStruFactAmplCal, output allStruFactAmpl
    computeAllStruFactAmplFromUniqStruFactAmplCal();

    //set F000
    float F000 = 30000;
    allStruFactAmpl[0][0][0] = F000;

    //compute allStruFactPhase from uniqStruFactPhase, input uniqStruFactPhase, output allStruFactPhase
    computeAllStruFactPhaseFromUniqStruFactPhase();

    //apply forward fft, input allStruFactAmpl and allStruFactPhase, output densUnitCell
    applyFFTStruFactToDens();

    //get densAsu from densUnitCell, input densUnitCell, output densAsu
    catchDensAsuFromDensUnitCell();

    //exclude top 10 outliers with exordinary high density values 
    int numGridPoint = 10;  
    float densMax0 = excludeExorDensAsu(numGridPoint); //densMax0 excludes exordinary high density values

    //exclude bottom 10 outliers with exordinary low density values 
    numGridPoint = numGridPointAsu - 10;  
    float densMin0 = excludeExorDensAsu(numGridPoint); //densMin0 excludes exordinary low density values

    //adjust exordinary density values to densMin0 or densMax0
    adjustExorDensAsu(densMin0, densMax0);

    //assign each grid point into a specific small bin, ideally, each small bin contains one or none grid point
    assignDensAsuIntoSmallBin("prot");

    //combine small bins into large bins, ideally, each large bin contains the same number of grid points
    combineSmallBinIntoLargeBin("prot");

    for(int ilargeBin = 0; ilargeBin <= numLargeBin; ++ilargeBin){
        stdBoundLargeBin[ilargeBin] = boundLargeBin[ilargeBin];
    }

    //retrieve the value of densAsu
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsu[ia][ib][ic] = densAsuTemp[ia][ib][ic];
                    protMaskAsu[ia][ib][ic] = protMaskAsuTemp[ia][ib][ic];
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void computeMaxPackSphereRadiusUnitCell(int ia0, float maxOverlapAllow)
{
    for( int ia=ia0-1; ia<ia0; ++ia ) {
        for( int ib=0; ib<BNUM; ++ib ) {
            for( int ic=0; ic<CNUM; ++ic ) {
                int occupiedVolume = 0, occupiedVolume1;
                int ioverlap = 0;
                int sphereRadius = 0;
                do{
                    sphereRadius++;
                    //use maskUnitCell to check overlap status
                    for( int ia1=0; ia1<ANUM; ++ia1 ) {
                        for( int ib1=0; ib1<BNUM; ++ib1 ) {
                            for( int ic1=0; ic1<CNUM; ++ic1 ) {
                                maskUnitCell[ia1][ib1][ic1] = 0;
                            }
                        }
                    }
                    occupiedVolume1 = occupiedVolume;
                    occupiedVolume = 0;
                    for(int ia1 = ia-sphereRadius; ia1 <= ia+sphereRadius; ++ia1) {
                        if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                        for(int ib1 = ib-sphereRadius; ib1 <= ib+sphereRadius; ++ib1) {
                            if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                            for(int ic1 = ic-sphereRadius; ic1 <= ic+sphereRadius; ++ic1) {
                                if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                                if( (ia1 - ia)*(ia1 - ia) + (ib1 - ib)*(ib1 - ib) + (ic1 - ic)*(ic1 - ic) \
                                   <= sphereRadius*sphereRadius ) {
                                    float temp[3], frac[3], rota[3][3], tran[3];
                                    temp[0] = (ia1 + 0.5) / float(ANUM);
                                    temp[1] = (ib1 + 0.5) / float(BNUM);
                                    temp[2] = (ic1 + 0.5) / float(CNUM);
                                    for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
                                        if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                                        for(int i = 0; i <= 2; ++i){
                                            for(int j = 0; j <= 2; ++j){
                                                rota[i][j] = symmOper[isymm][i][j]; // load rotation matrix
                                            } 
                                            tran[i] = symmOper[isymm][i][3]; // load rotation matrix
                                        }
                                        for( int i=0; i<=2 ; ++i ){
                                            frac[i] = 0;
                                            for( int j=0; j<=2 ; ++j ){
                                                frac[i] += rota[i][j] * temp[j]; // rotation
                                            }
                                            frac[i] += tran[i]; // translation
                                        }
                                        adjustFracCoorIntoUnitCell(frac);
                                        int ia2 = floor(frac[0] * ANUM);
                                        int ib2 = floor(frac[1] * BNUM);
                                        int ic2 = floor(frac[2] * CNUM);
                                        if(maskUnitCell[ia2][ib2][ic2] == 0){
                                            maskUnitCell[ia2][ib2][ic2] = 1;
                                            occupiedVolume++;
                                        }else{
                                            ioverlap++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }while(ioverlap <= round(maxOverlapAllow*occupiedVolume1));
                maxPackSphereRadiusUnitCell[ia][ib][ic] = sphereRadius - 1;
            }
        }
    }
    return;
}

/******************************************************************************/

void computeMaxPackSphereRadiusUnitCellForMonomer(int ia0, float maxOverlapAllow)
{
    for( int ia=ia0-1; ia<ia0; ++ia ) {
        for( int ib=0; ib<BNUM; ++ib ) {
            for( int ic=0; ic<CNUM; ++ic ) {
                int occupiedVolume = 0, occupiedVolume1;
                int ioverlap = 0;
                int sphereRadius = 0;
                do{
                    sphereRadius++;
                    //use maskUnitCell to check overlap status
                    for( int ia1=0; ia1<ANUM; ++ia1 ) {
                        for( int ib1=0; ib1<BNUM; ++ib1 ) {
                            for( int ic1=0; ic1<CNUM; ++ic1 ) {
                                maskUnitCell[ia1][ib1][ic1] = 0;
                            }
                        }
                    }
                    occupiedVolume1 = occupiedVolume;
                    occupiedVolume = 0;
                    float frac[3], point[3], newPoint[3];
                    frac[0] = (ia + 0.5) / float(ANUM);
                    frac[1] = (ib + 0.5) / float(BNUM);
                    frac[2] = (ic + 0.5) / float(CNUM);
                    convertFracToOrth(frac, point);
                    for(int incs = 0; incs <= numNcsOper - 1; ++incs){
                        if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
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
                        int ia2 = floor(frac[0] * ANUM);
                        int ib2 = floor(frac[1] * BNUM);
                        int ic2 = floor(frac[2] * CNUM);
                        for(int ia1 = ia2-sphereRadius; ia1 <= ia2+sphereRadius; ++ia1) {
                            if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                            for(int ib1 = ib2-sphereRadius; ib1 <= ib2+sphereRadius; ++ib1) {
                                if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                                for(int ic1 = ic2-sphereRadius; ic1 <= ic2+sphereRadius; ++ic1) {
                                    if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                                    if( (ia1 - ia2)*(ia1 - ia2) + (ib1 - ib2)*(ib1 - ib2) + (ic1 - ic2)*(ic1 - ic2)
                                        <= sphereRadius*sphereRadius ) {
                                        float temp[3], frac[3], orth[3], rota[3][3], tran[3];
                                        temp[0] = (ia1 + 0.5) / float(ANUM);
                                        temp[1] = (ib1 + 0.5) / float(BNUM);
                                        temp[2] = (ic1 + 0.5) / float(CNUM);
                                        for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
                                            if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                                            for(int i = 0; i <= 2; ++i){
                                                for(int j = 0; j <= 2; ++j){
                                                    rota[i][j] = symmOper[isymm][i][j]; // load rotation matrix
                                                } 
                                                tran[i] = symmOper[isymm][i][3]; // load rotation matrix
                                            }
                                            for( int i=0; i<=2 ; ++i ){
                                                frac[i] = 0;
                                                for( int j=0; j<=2 ; ++j ){
                                                    frac[i] += rota[i][j] * temp[j]; // rotation
                                                }
                                                frac[i] += tran[i]; // translation
                                            }
                                            adjustFracCoorIntoUnitCell(frac);
                                            int ia3 = floor(frac[0] * ANUM);
                                            int ib3 = floor(frac[1] * BNUM);
                                            int ic3 = floor(frac[2] * CNUM);
                                            if(maskUnitCell[ia3][ib3][ic3] == 0){
                                                maskUnitCell[ia3][ib3][ic3] = 1;
                                                occupiedVolume++;
                                            }else{
                                                ioverlap++;
                                            }
                                        } //isymm
                                    } //if inside sphere
                                }
                            }
                        }
                    } //incs
                }while(ioverlap <= round(maxOverlapAllow*occupiedVolume1));
                maxPackSphereRadiusUnitCell[ia][ib][ic] = sphereRadius - 1;
            }
        }
    }
    return;
}

/******************************************************************************/

void computeMaxPackCylinderRadiusUnitCell(int ia0, float maxOverlapAllow)
{
    for( int ia=ia0-1; ia<ia0; ++ia ) {
        for( int ib=0; ib<BNUM; ++ib ) {
            for( int ic=0; ic<CNUM; ++ic ) {
                float frac[3], orth[3];
                frac[0] = (ia + 0.5) / float(ANUM);
                frac[1] = (ib + 0.5) / float(BNUM);
                frac[2] = (ic + 0.5) / float(CNUM);
                convertFracToOrth(frac, orth);

                int occupiedVolume = 0, occupiedVolume1;
                int ioverlap = 0;
                int cylinderRadius = 1;
                float cylinderRadiusReal = 0;
                do{
                    cylinderRadius++;
                    cylinderRadiusReal = cylinderRadius * (CRYST_A / ANUM); //
                    //use maskUnitCell to check overlap status
                    for( int ia1=0; ia1<ANUM; ++ia1 ) {
                        for( int ib1=0; ib1<BNUM; ++ib1 ) {
                            for( int ic1=0; ic1<CNUM; ++ic1 ) {
                                maskUnitCell[ia1][ib1][ic1] = 0;
                            }
                        }
                    }
                    occupiedVolume1 = occupiedVolume;
                    occupiedVolume = 0;
                    int cylinderHeight = 2 * cylinderRadius; //must be an even number
                    float cylinderHeightReal = 2.0 * cylinderRadiusReal;
                    int sphereRadius = ceil( sqrt( cylinderRadius*cylinderRadius + cylinderHeight/2.0*cylinderHeight/2.0 ) );
                    for(int ia1 = ia-sphereRadius; ia1 <= ia+sphereRadius; ++ia1) {
                        if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                        for(int ib1 = ib-sphereRadius; ib1 <= ib+sphereRadius; ++ib1) {
                            if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                            for(int ic1 = ic-sphereRadius; ic1 <= ic+sphereRadius; ++ic1) {
                                if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                                //mush inside a big sphere
                                if( (ia1 - ia)*(ia1 - ia) + (ib1 - ib)*(ib1 - ib) + (ic1 - ic)*(ic1 - ic) \
                                   <= sphereRadius*sphereRadius ) {
                                    float frac1[3], orth1[3];
                                    frac1[0] = (ia1 + 0.5) / float(ANUM);
                                    frac1[1] = (ib1 + 0.5) / float(BNUM);
                                    frac1[2] = (ic1 + 0.5) / float(CNUM);
                                    convertFracToOrth(frac1, orth1);
                                    float pointBOnLine[3];
                                    for(int i = 0; i < 3; ++i){
                                        pointBOnLine[i] = orth[i] + vectOnNcsAxisIni[i];
                                    }
                                    float dist = computeDistFromPointToLine(orth1, orth, pointBOnLine);
                                    if(dist < cylinderRadiusReal){
                                        float dist1 = 0;
                                        for(int i = 0; i < 3; ++i){
                                            dist1 = dist1 + (orth1[i] - orth[i]) * vectOnNcsAxisIni[i];
                                        }
                                        if(dist1 <= cylinderHeightReal){
                                            float frac2[3], rota[3][3], tran[3];
                                            for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
                                                if(ioverlap > round(maxOverlapAllow*occupiedVolume1)) break;
                                                for(int i = 0; i <= 2; ++i){
                                                    for(int j = 0; j <= 2; ++j){
                                                        rota[i][j] = symmOper[isymm][i][j]; // load rotation matrix
                                                    } 
                                                    tran[i] = symmOper[isymm][i][3]; // load rotation matrix
                                                }
                                                for( int i=0; i<=2 ; ++i ){
                                                    frac2[i] = 0;
                                                    for( int j=0; j<=2 ; ++j ){
                                                        frac2[i] += rota[i][j] * frac1[j]; // rotation
                                                    }
                                                    frac2[i] += tran[i]; // translation
                                                }
                                                adjustFracCoorIntoUnitCell(frac2);
                                                int ia2 = floor(frac2[0] * ANUM);
                                                int ib2 = floor(frac2[1] * BNUM);
                                                int ic2 = floor(frac2[2] * CNUM);
                                                if(maskUnitCell[ia2][ib2][ic2] == 0){
                                                    maskUnitCell[ia2][ib2][ic2] = 1;
                                                    occupiedVolume++;
                                                }else{
                                                    ioverlap++;

                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }while(ioverlap <= round(maxOverlapAllow*occupiedVolume1));
                maxPackCylinderRadiusUnitCell[ia][ib][ic] = cylinderRadius - 1;
            }
        }
    }
    return;
}

/******************************************************************************/

void allocateNcsMask(int iter)
{
    //obtain weigAvg density frequence in whole asu, input weigAvgDensAsu, output numGridPointSmallBin
    assignWeigAvgDensAsuIntoSmallBin("cell");

    //group weigAvg density in ascending order, input weigAvgDensAsu, numGridPointSmallBin, output boundLargeBinNcs
    float boundLargeBinNcs[numLargeBinNcs + 1] = {};
    int densFreqEachLargeBin = round( float(numGridPointAsu) / float(numLargeBinNcs) ); 
    findWeigAvgDensAsuMinMaxAvg("cell");
    boundLargeBinNcs[0] = densMin;
    int ilargeBin = 0;
    int i = 0;
    for(int ismallBin = 0; ismallBin <= numSmallBin-1; ++ismallBin){
        i = i + numGridPointSmallBin[ismallBin];
        if(i >= densFreqEachLargeBin || ismallBin == numSmallBin-1) {
            i = 0;
            ilargeBin++;
            if(ilargeBin <= numLargeBinNcs-1) {
                boundLargeBinNcs[ilargeBin] = boundLargeBinNcs[0] + float(ismallBin) / float(numSmallBin) * (densMax-densMin);
            }
        }
    }
    while(ilargeBin <= numLargeBinNcs-1) {
        boundLargeBinNcs[ilargeBin] = densMax;
        ilargeBin++;
    }
    boundLargeBinNcs[numLargeBinNcs] = densMax;

    //rewrite weigAvgDens in a list in descending order of weigAvgDens value
    for(int igrid = 0; igrid <= numGridPointAsu - 1; ++igrid){
        listGridPointAsu[igrid][0] = 0;
        listGridPointAsu[igrid][1] = 0;
        listGridPointAsu[igrid][2] = 0;
        listGridPointAsu[igrid][3] = 0;
        listGridPointNcsMask[igrid][0] = 0;
        listGridPointNcsMask[igrid][1] = 0;
        listGridPointNcsMask[igrid][2] = 0;
        listGridPointNcsMask[igrid][3] = 0;
        listGridPointNcsMask[igrid][4] = 0;
        listGridPointNcsMaskSurf[igrid][0] = 0;
        listGridPointNcsMaskSurf[igrid][1] = 0;
        listGridPointNcsMaskSurf[igrid][2] = 0;
        listGridPointNcsMaskSurf[igrid][3] = 0;
        listGridPointNcsMaskSurf[igrid][4] = 0;
    }
    int igrid = 0;
    float density;
    for(ilargeBin = numLargeBinNcs - 1; ilargeBin >= 0; --ilargeBin) {
        for( int IALOOPASU ) {
            float x = (ia + 0.5) / float(ANUM);
            for( int IBLOOPASU ) {
                float y = (ib + 0.5) / float(BNUM);
                for( int ICLOOPASU ) {
                    float z = (ic + 0.5) / float(CNUM);
                    if( ASUCOND ) {
                        density = weigAvgDensAsu[ia][ib][ic];
                        if(density < boundLargeBinNcs[ilargeBin+1] && density >= boundLargeBinNcs[ilargeBin]) {
                            listGridPointAsu[igrid][0] = ia;
                            listGridPointAsu[igrid][1] = ib;
                            listGridPointAsu[igrid][2] = ic;
                            listGridPointAsu[igrid][3] = density;
                            igrid++;
                        } else if((ilargeBin == numLargeBinNcs-1) && density == densMax){
                            listGridPointAsu[igrid][0] = ia;
                            listGridPointAsu[igrid][1] = ib;
                            listGridPointAsu[igrid][2] = ic;
                            listGridPointAsu[igrid][3] = density;
                            igrid++;
                        }
                    }
                }
            }
        }
    }
    //grow the region of ncs averaging
    float rota[3][3], tran[3];
    float orth[3], orth1[3], frac[3], frac1[3], frac2[3];

    int numAdjaCand = 5;
    float dist;
    int touchFlag = 0;
    int igrid2 = 0;

    for(int igrid = 0; igrid <= numGridPointAsu - 1; ++igrid){

        //find the grid points on the surface of the growing ncs mask at present
        if(touchFlag == 0 && igrid < densFreqEachLargeBin){
            igrid2 = 0;
            for(int igrid1 = 0; igrid1 <= igrid - 1; ++igrid1){
                //at the beginning, listGridPointNcsMaskSurf store grid around center of mass
                listGridPointNcsMaskSurf[igrid2][0] = listGridPointNcsMask[igrid1][0];
                listGridPointNcsMaskSurf[igrid2][1] = listGridPointNcsMask[igrid1][1];
                listGridPointNcsMaskSurf[igrid2][2] = listGridPointNcsMask[igrid1][2];
                listGridPointNcsMaskSurf[igrid2][3] = listGridPointNcsMask[igrid1][3];
                igrid2++;
            }
        } else if(touchFlag == 0 && igrid % densFreqEachLargeBin == 0){
            for( int IALOOPASU ) {
                float x = (ia + 0.5) / float(ANUM);
                for( int IBLOOPASU ) {
                    float y = (ib + 0.5) / float(BNUM);
                    for( int ICLOOPASU ) {
                        float z = (ic + 0.5) / float(CNUM);
                        if( ASUCOND ) {
                            maskAsu[ia][ib][ic] = 0;
                        }
                    }
                }
            }
            for(int igrid1 = 0; igrid1 <= igrid - 1; ++igrid1){
                int ia = listGridPointAsu[igrid1][0];
                int ib = listGridPointAsu[igrid1][1];
                int ic = listGridPointAsu[igrid1][2];
                maskAsu[ia][ib][ic] = 1;
            }
            igrid2 = 0;
            for(int igrid1 = 0; igrid1 <= igrid - 1; ++igrid1){
                int ia = listGridPointAsu[igrid1][0];
                int ib = listGridPointAsu[igrid1][1];
                int ic = listGridPointAsu[igrid1][2];
                int flag = 0;
                for(int ia1 = ia - 1; ia1 <= ia + 1; ++ia1){
                    if(flag == 1)break;
                    for(int ib1 = ib - 1; ib1 <= ib + 1; ++ib1){
                        if(flag == 1)break;
                        for(int ic1 = ic - 1; ic1 <= ic + 1; ++ic1){
                            if(flag == 1)break;
                            if(ia1 == ia && ib1 == ib && ic1 == ic)continue;
                            frac[0] = (ia1 + 0.5) / float(ANUM);
                            frac[1] = (ib1 + 0.5) / float(BNUM);
                            frac[2] = (ic1 + 0.5) / float(CNUM);
                            adjustFracCoorIntoAsu(frac);
                            int ia0 = floor(frac[0] * ANUM);
                            int ib0 = floor(frac[1] * BNUM);
                            int ic0 = floor(frac[2] * CNUM);
                            if(maskAsu[ia0][ib0][ic0] == 0){
                                flag = 1;
                                //listGridPointNcsMaskSurf stores surface grid of the current ncs core
                                listGridPointNcsMaskSurf[igrid2][0] = listGridPointNcsMask[igrid1][0];
                                listGridPointNcsMaskSurf[igrid2][1] = listGridPointNcsMask[igrid1][1];
                                listGridPointNcsMaskSurf[igrid2][2] = listGridPointNcsMask[igrid1][2];
                                listGridPointNcsMaskSurf[igrid2][3] = listGridPointNcsMask[igrid1][3];
                                igrid2++;
                            }
                        }
                    }
                }
            }
        } 

        for(int i = 0; i <= 2; ++i){
            frac[i] = ( listGridPointAsu[igrid][i] + 0.5 ) / lengEachDime[i];
        }

        //apply symmtry operation to get symmetry copies in unit cell
        float distAdjaCand[numSymmOper*numAdjaCell][4] = {};
        for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
            for(int i = 0; i <= 2; ++i){
                for(int j = 0; j <= 2; ++j){
                    rota[i][j] = symmOper[isymm][i][j]; //load rotation matrix
                } 
                tran[i] = symmOper[isymm][i][3]; //load rotation matrix
            }
            for( int i=0; i<=2 ; ++i ){
                frac1[i] = 0;
                for( int j=0; j<=2 ; ++j ){
                    frac1[i] = frac1[i] + rota[i][j] * frac[j]; //rotation
                }
                frac1[i] = frac1[i] + tran[i]; //translation
            }
            adjustFracCoorIntoUnitCell(frac1);  
            //27 adjacent cells, 9 cells on top, 9 cells in the middle, and 9 cells at bottom
            for(int icell = 0; icell <= numAdjaCell - 1; ++icell){
                for(int i = 0; i <= 2; ++i){
                    frac2[i] = frac1[i] + adjaCell[icell][i];
                }
                convertFracToOrth(frac2, orth);
                dist = 10000.0;
                for(int icent = 0; icent <= numCentMass - 1; ++icent){
                    float distTemp = 0;
                    for(int i = 0; i <= 2; ++i){
                        distTemp = distTemp + pow((orth[i] - centMassOrth[icent][i]), 2.0);
                    }
                    distTemp = sqrt(distTemp);
                    if(distTemp < dist) dist = distTemp;
                }
                distAdjaCand[isymm*numAdjaCell+icell][0] = orth[0];
                distAdjaCand[isymm*numAdjaCell+icell][1] = orth[1];
                distAdjaCand[isymm*numAdjaCell+icell][2] = orth[2];
                distAdjaCand[isymm*numAdjaCell+icell][3] = dist;
            }
        }

        //several grids among numSymmOper*numAdjaCell close to the default mass center
        for(int i = 0; i <= numAdjaCand - 1; ++i){
            for(int j = i+1; j <= numSymmOper*numAdjaCell - 1; ++j){
                if(distAdjaCand[j][3] < distAdjaCand[i][3]){
                    float temp;
                    for(int k = 0; k <= 3; ++k){
                        temp = distAdjaCand[i][k];
                        distAdjaCand[i][k] = distAdjaCand[j][k];
                        distAdjaCand[j][k] = temp;
                    }
                }
            }
        }

        for(int i = 0; i <= numAdjaCand - 1; ++i){ //grids among numSymmOper*numAdjaCell close to the default mass center
            for(int j = 0; j <= 2; ++j){
                orth[j] = distAdjaCand[i][j];
            }
            float distMin = 10000.0;
            if(igrid == 0 || touchFlag == 0){
                dist = 10000.0;
                for(int icent = 0; icent <= numCentMass - 1; ++icent){
                    float distTemp = 0;
                    for(int i = 0; i <= 2; ++i){
                        distTemp = distTemp + pow((orth[i] - centMassOrth[icent][i]), 2.0);
                    }
                    distTemp = sqrt(distTemp);
                    if(distTemp < dist) dist = distTemp;
                }
                if(dist < distMin) distMin = dist;
            } 
            if(igrid > 0){
                for(int igrid1 = 0; igrid1 <= igrid2 - 1; ++igrid1){
                    //select 300 grid points on growing-ncs-mask surface, 
                    //if there are too many grid on surface, select part of them 1/(igrid2/300), eg 1/2, 1/3
                    if(igrid2 > 300 && igrid1 % (igrid2/300) != 0) continue;
                    for(int j = 0; j <= 2; ++j){
                        frac1[j] = ( round(listGridPointNcsMaskSurf[igrid1][j]) + 0.5 ) / lengEachDime[j];
                    }
                    convertFracToOrth(frac1, orth1);
                    if(touchFlag == 0){
                        dist = distOrth(orth, orth1) + listGridPointNcsMaskSurf[igrid1][3];
                    } else{
                        dist = distOrth(orth, orth1);
                    }
                    if(dist < distMin) distMin = dist;
                }
            } 
            distAdjaCand[i][3] = distMin;
        }

        for(int i = 0; i <= numAdjaCand - 2; ++i){
            for(int j = i+1; j <= numAdjaCand - 1; ++j){
                if(distAdjaCand[j][3] < distAdjaCand[i][3]){
                    float temp;
                    for(int k = 0; k <= 3; ++k){
                        temp = distAdjaCand[i][k];
                        distAdjaCand[i][k] = distAdjaCand[j][k];
                        distAdjaCand[j][k] = temp;
                    }
                }
            }
        }

        if(abs(distAdjaCand[1][3]-distAdjaCand[0][3]) <= resoCutoff/2.0){
            touchFlag = 1; //ncs mask begin to touch its symmetric ones
        }

        for(int j = 0; j <= 2; ++j){
            orth[j] = distAdjaCand[0][j];
        }
        convertOrthToFrac(orth, frac);

        listGridPointNcsMask[igrid][0] = floor(frac[0] * ANUM);
        listGridPointNcsMask[igrid][1] = floor(frac[1] * BNUM);
        listGridPointNcsMask[igrid][2] = floor(frac[2] * CNUM);
        listGridPointNcsMask[igrid][3] = distAdjaCand[0][3];
        if(touchFlag == 0) {
            listGridPointNcsMask[igrid][4] = 1;
        } else{
            listGridPointNcsMask[igrid][4] = 2;
        }
    }

    //allocate a large qubic region to hold the ncs mask
    anum0 = listGridPointNcsMask[0][0];
    bnum0 = listGridPointNcsMask[0][1];
    cnum0 = listGridPointNcsMask[0][2];
    anum1 = anum0; bnum1 = bnum0; cnum1 = cnum0;
    int ia, ib, ic;
    for(int igrid = 0; igrid <= numGridPointAsu - 1; ++igrid){
        ia = listGridPointNcsMask[igrid][0];
        ib = listGridPointNcsMask[igrid][1];
        ic = listGridPointNcsMask[igrid][2];
        if(ia < anum0) anum0 = ia;
        if(ib < bnum0) bnum0 = ib;
        if(ic < cnum0) cnum0 = ic;
        if(ia > anum1) anum1 = ia;
        if(ib > bnum1) bnum1 = ib;
        if(ic > cnum1) cnum1 = ic;
    }
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

    for(int ia = 0; ia <= NCSANUM - 1; ++ia){
        for(int ib = 0; ib <= NCSBNUM - 1; ++ib){
            for(int ic = 0; ic <= NCSCNUM - 1; ++ic){
                ncsMask[ia][ib][ic] = 0;
                ncsMaskDensTemp[ia][ib][ic] = 0;
                ncsMaskDens[ia][ib][ic] = 0;
            }
        }
    }

    for(int igrid = 0; igrid <= numGridPointAsu - 1; ++igrid){
        int ia = listGridPointNcsMask[igrid][0] - anum0;
        int ib = listGridPointNcsMask[igrid][1] - bnum0;
        int ic = listGridPointNcsMask[igrid][2] - cnum0;
        //ncsMask[ia][ib][ic] = 1;
        ncsMask[ia][ib][ic] = round(listGridPointNcsMask[igrid][4]); //==1 ~ ncs core, >=1 ~ ncs mask 
    }

    return;
}

/******************************************************************************/

void outputNcsCore(string fileNameString)

{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    outputFile << writeCryst1() << endl;
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;
    int ia1, ib1, ic1;
    float orth[3], frac[3], frac1[3];
    int iatom = 0;
    float tempFact, occupancy = 1.0;
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                if(ncsMask[ia][ib][ic] == 1){ // ==1 ncs core, >=1 ncs mask
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    for(int i = 0; i <= 2; ++i){
                        frac1[i] = frac[i];
                    }
                    adjustFracCoorIntoAsu(frac1);
                    ia1 = floor(frac1[0] * ANUM);
                    ib1 = floor(frac1[1] * BNUM);
                    ic1 = floor(frac1[2] * CNUM);
                    convertFracToOrth(frac, orth);
                    tempFact = weigAvgDensAsu[ia1][ib1][ic1];
                    occupancy = 1.0;
                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                        << ' ' << " C   ARG A   1    " 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                        << string(11, ' ') << 'C' << endl;
                    iatom++;
                    if(iatom > 9999) iatom = 0; //too many atoms
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputNcsMask(string fileNameString)

{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    outputFile << writeCryst1() << endl;
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;
    int ia1, ib1, ic1;
    float orth[3], frac[3], frac1[3];
    int iatom = 0;
    float tempFact, occupancy = 1.0;
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    for(int i = 0; i <= 2; ++i){
                        frac1[i] = frac[i];
                    }
                    adjustFracCoorIntoAsu(frac1);
                    ia1 = floor(frac1[0] * ANUM);
                    ib1 = floor(frac1[1] * BNUM);
                    ic1 = floor(frac1[2] * CNUM);
                    convertFracToOrth(frac, orth);
                    tempFact = weigAvgDensAsu[ia1][ib1][ic1];
                    occupancy = 1.0;
                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                        << ' ' << " C   ARG A   1    " 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                        << string(11, ' ') << 'C' << endl;
                    iatom++;
                    if(iatom > 9999) iatom = 0; //too many atoms
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void seperateNcsMask(int iter)
{
    int densFreqEachLargeBin = round( float(numGridPointAsu) / float(numLargeBinNcs) ); 
    //listGridPointNcsMask previously store the ncs mask ia, ib, ic, in decending order of their weigAvgDensAsu
    for(int igrid = 0; igrid <= numGridPointAsu - 1; ++igrid){
        //keep listGridPointNcsMask[igrid][0], listGridPointNcsMask[igrid][1], listGridPointNcsMask[igrid][2]
        //previously store 1 for ncs inner core, 2 for ncs outer region
        listGridPointNcsMask[igrid][3] = 0; //now ==1 for mole A, ==2 for mole B, ...
        listGridPointNcsMask[igrid][4] = 0;
        listGridPointNcsMaskSurf[igrid][0] = 0;
        listGridPointNcsMaskSurf[igrid][1] = 0;
        listGridPointNcsMaskSurf[igrid][2] = 0;
        listGridPointNcsMaskSurf[igrid][3] = 0;
        listGridPointNcsMaskSurf[igrid][4] = 0;
    }

    float orth[3], orth1[3], frac[3], frac1[3];
    float dist;
    int touchFlag = 0;
    int igrid2 = 0;

    for(int igrid = 0; igrid <= numGridPointAsu - 1; ++igrid){

        //find the grid points on the surface of the growing ncs mask at present
        if(touchFlag == 0 && igrid < densFreqEachLargeBin){
            igrid2 = 0;
            for(int igrid1 = 0; igrid1 <= igrid - 1; ++igrid1){
                //at the beginning, listGridPointNcsMaskSurf store grid around center of mass
                listGridPointNcsMaskSurf[igrid2][0] = listGridPointNcsMask[igrid1][0];
                listGridPointNcsMaskSurf[igrid2][1] = listGridPointNcsMask[igrid1][1];
                listGridPointNcsMaskSurf[igrid2][2] = listGridPointNcsMask[igrid1][2];
                listGridPointNcsMaskSurf[igrid2][3] = listGridPointNcsMask[igrid1][3];
                listGridPointNcsMaskSurf[igrid2][4] = listGridPointNcsMask[igrid1][4];
                igrid2++;
            }
        } else if(touchFlag == 0 && igrid % densFreqEachLargeBin == 0){
            for( int IALOOPASU ) {
                float x = (ia + 0.5) / float(ANUM);
                for( int IBLOOPASU ) {
                    float y = (ib + 0.5) / float(BNUM);
                    for( int ICLOOPASU ) {
                        float z = (ic + 0.5) / float(CNUM);
                        if( ASUCOND ) {
                            maskAsu[ia][ib][ic] = 0;
                        }
                    }
                }
            }
            for(int igrid1 = 0; igrid1 <= igrid - 1; ++igrid1){
                int ia = listGridPointAsu[igrid1][0];
                int ib = listGridPointAsu[igrid1][1];
                int ic = listGridPointAsu[igrid1][2];
                maskAsu[ia][ib][ic] = 1;
            }
            igrid2 = 0;
            for(int igrid1 = 0; igrid1 <= igrid - 1; ++igrid1){
                int ia = listGridPointAsu[igrid1][0];
                int ib = listGridPointAsu[igrid1][1];
                int ic = listGridPointAsu[igrid1][2];
                int flag = 0;
                for(int ia1 = ia - 1; ia1 <= ia + 1; ++ia1){
                    if(flag == 1)break;
                    for(int ib1 = ib - 1; ib1 <= ib + 1; ++ib1){
                        if(flag == 1)break;
                        for(int ic1 = ic - 1; ic1 <= ic + 1; ++ic1){
                            if(flag == 1)break;
                            if(ia1 == ia && ib1 == ib && ic1 == ic)continue;
                            frac[0] = (ia1 + 0.5) / float(ANUM);
                            frac[1] = (ib1 + 0.5) / float(BNUM);
                            frac[2] = (ic1 + 0.5) / float(CNUM);
                            adjustFracCoorIntoAsu(frac);
                            int ia0 = floor(frac[0] * ANUM);
                            int ib0 = floor(frac[1] * BNUM);
                            int ic0 = floor(frac[2] * CNUM);
                            if(maskAsu[ia0][ib0][ic0] == 0){
                                flag = 1;
                                //listGridPointNcsMaskSurf stores surface grid of the current ncs core
                                listGridPointNcsMaskSurf[igrid2][0] = listGridPointNcsMask[igrid1][0];
                                listGridPointNcsMaskSurf[igrid2][1] = listGridPointNcsMask[igrid1][1];
                                listGridPointNcsMaskSurf[igrid2][2] = listGridPointNcsMask[igrid1][2];
                                listGridPointNcsMaskSurf[igrid2][3] = listGridPointNcsMask[igrid1][3];
                                listGridPointNcsMaskSurf[igrid2][4] = listGridPointNcsMask[igrid1][4];
                                igrid2++;
                            }
                        }
                    }
                }
            }        
        } 

        for(int i = 0; i <= 2; ++i){
            frac[i] = ( listGridPointNcsMask[igrid][i] + 0.5 ) / lengEachDime[i];
        }
        convertFracToOrth(frac, orth);

        float distMin = 10000.0;
        float distMinPlus = 10000.0;
        int centChoice = 0;
        if(igrid == 0 || touchFlag == 0){
            for(int icent = 0; icent <= numCentMass - 1; ++icent){
                dist = 0;
                for(int i = 0; i <= 2; ++i){
                    dist = dist + pow((orth[i] - centMassOrth[icent][i]), 2.0);
                }
                dist = sqrt(dist);
                if(dist < distMin) {
                    if(icent != centChoice){
                        centChoice = icent;
                        distMinPlus = distMin;
                    }
                    distMin = dist;
                }
            }
        } 
        if(igrid > 0){
            for(int igrid1 = 0; igrid1 <= igrid2 - 1; ++igrid1){
                //select 300 grid points on growing-ncs-mask surface, 
                //if there are too many grid on surface, select part of them 1/(igrid2/300), eg 1/2, 1/3
                if(igrid2 > 300 && igrid1 % (igrid2/300) != 0) continue; 
                for(int j = 0; j <= 2; ++j){
                    frac1[j] = ( round(listGridPointNcsMaskSurf[igrid1][j]) + 0.5 ) / lengEachDime[j];
                }
                convertFracToOrth(frac1, orth1);
                if(touchFlag == 0){
                    dist = distOrth(orth, orth1) + listGridPointNcsMaskSurf[igrid1][3];
                } else{
                    dist = distOrth(orth, orth1);
                }
                if(dist < distMin) {
                    if(round(listGridPointNcsMaskSurf[igrid1][4]) != centChoice){
                        centChoice = listGridPointNcsMaskSurf[igrid1][4];
                        distMinPlus = distMin;
                    }
                    distMin = dist;
                }
            }
        } 
        listGridPointNcsMask[igrid][3] = distMin;
        listGridPointNcsMask[igrid][4] = centChoice;

        if(abs(distMinPlus-distMin) <= resoCutoff/2.0){
            touchFlag = 1; //ncs mask belonging to different copies begin to touch
        }
    }

    for(int igrid = 0; igrid <= numGridPointAsu - 1; ++igrid){
        int ia = listGridPointNcsMask[igrid][0] - anum0;
        int ib = listGridPointNcsMask[igrid][1] - bnum0;
        int ic = listGridPointNcsMask[igrid][2] - cnum0;
        ncsMask[ia][ib][ic] = round(listGridPointNcsMask[igrid][4]) + 1; //0+1=1 mole A; 1+1=2 moleB; 
    }

    return;
}

/******************************************************************************/

void outputSepeNcsMask(string fileNameString)

{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    outputFile << writeCryst1() << endl;
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;
    int ia1, ib1, ic1;
    float orth[3], frac[3], frac1[3];
    int iatom = 0;
    char chain;
    float tempFact, occupancy = 1.0;
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                if(ncsMask[ia][ib][ic] >= 1){
                    if(ncsMask[ia][ib][ic] == 1) chain = 'A';
                    if(ncsMask[ia][ib][ic] == 2) chain = 'B';
                    if(ncsMask[ia][ib][ic] == 3) chain = 'C';
                    if(ncsMask[ia][ib][ic] == 4) chain = 'D';
                    if(ncsMask[ia][ib][ic] == 5) chain = 'E';
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    for(int i = 0; i <= 2; ++i){
                        frac1[i] = frac[i];
                    }
                    adjustFracCoorIntoAsu(frac1);
                    ia1 = floor(frac1[0] * ANUM);
                    ib1 = floor(frac1[1] * BNUM);
                    ic1 = floor(frac1[2] * CNUM);
                    convertFracToOrth(frac, orth);
                    tempFact = weigAvgDensAsu[ia1][ib1][ic1];
                    occupancy = 1.0;
                       outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                        << ' ' << " C   ARG " << chain << "   1    " 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                        << string(11, ' ') << 'C' << endl;
                    iatom++;
                    if(iatom > 9999) iatom = 0; //too many atoms
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputProtMask(string fileNameString)

{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    outputFile << writeCryst1() << endl;
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;
    int ia1, ib1, ic1;
    float orth[3], frac[3], frac1[3];
    int iatom = 0;
    float tempFact, occupancy = 1.0;
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    for(int i = 0; i <= 2; ++i){
                        frac1[i] = frac[i];
                    }
                    adjustFracCoorIntoAsu(frac1);
                    ia1 = floor(frac1[0] * ANUM);
                    ib1 = floor(frac1[1] * BNUM);
                    ic1 = floor(frac1[2] * CNUM);
                    if(protMaskAsu[ia1][ib1][ic1] == 1){
                        convertFracToOrth(frac, orth);
                        tempFact = weigAvgDensAsu[ia1][ib1][ic1];
                        occupancy = 1.0;
                           outputFile << "ATOM  " << setfill(' ') << setw(5) << iatom 
                            << ' ' << " C   ARG A   1    " 
                            << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[0]
                            << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[1] 
                            << setfill(' ') << setw(8) << fixed << setprecision(3) << orth[2] 
                            << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                            << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                            << string(11, ' ') << 'C' << endl;
                        iatom++;
                        if(iatom > 9999) iatom = 0; //too many atoms
                    }
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void minimizeNcsDensDiffNcsAvgAxis()
{
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;

    //apply ncs operator
    float rota[3][3], tran[3];
    float orth[3], orth1[3], frac[3], frac1[3];
    int ia1, ib1, ic1;
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                ncsMaskDensTemp[ia][ib][ic] = 0;
                ncsMaskDens[ia][ib][ic] = 0;
                //if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    ncsMaskDensTemp[ia][ib][ic] = densAsuTemp[ia1][ib1][ic1];
                    ncsMaskDens[ia][ib][ic] = densAsuTemp[ia1][ib1][ic1]; //corresponds to incs = 0
                //}
            }
        }
    }

    float ncsDensDiffMin;
    int icand0;

    for(int icand = 0; icand < numNcsAxisCand; ++icand){

        float ncsDensDiff = 0;
        
        float a0 = pointAOnNcsAxisCand[icand][0];
        float b0 = pointAOnNcsAxisCand[icand][1];
        float c0 = pointAOnNcsAxisCand[icand][2];
        float u = vectOnNcsAxisCand[icand][0];
        float v = vectOnNcsAxisCand[icand][1];
        float w = vectOnNcsAxisCand[icand][2];
    
        float L0 = u*u + v*v + w*w;
        for(int ia = 0; ia <= anum2; ++ia){
            for(int ib = 0; ib <= bnum2; ++ib){
                for(int ic = 0; ic <= cnum2; ++ic){
    
                    if(ncsMask[ia][ib][ic] == 0) continue;
                    frac[0] = (ia + anum0 + 0.5) / float(ANUM);
                    frac[1] = (ib + bnum0 + 0.5) / float(BNUM);
                    frac[2] = (ic + cnum0 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    if(protMaskAsu[ia1][ib1][ic1] == 0)continue;
    
                    frac[0] = (ia + anum0 + 0.5) / float(ANUM);
                    frac[1] = (ib + bnum0 + 0.5) / float(BNUM);
                    frac[2] = (ic + cnum0 + 0.5) / float(CNUM);
                    convertFracToOrth(frac, orth);
                    float x = orth[0];
                    float y = orth[1];
                    float z = orth[2];
                    int incsAvg = 1; //have considered incs = 0
                    for(int incs = 1; incs <= numNcsOper - 1; ++incs) { //start from incs = 1
                        float thita = 2.0*PI/float(numNcsOper)*float(incs);
                        orth1[0] = 1.0/L0*(  (a0*(v*v+w*w)-u*(b0*v+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
                            +L0*x*cos(thita)+sqrt(L0)*(-c0*v+b0*w-w*y+v*z)*sin(thita)  );
                        orth1[1] = 1.0/L0*(  (b0*(u*u+w*w)-v*(a0*u+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
                            +L0*y*cos(thita)+sqrt(L0)*( c0*u-a0*w+w*x-u*z)*sin(thita)  );
                        orth1[2] = 1.0/L0*(  (c0*(u*u+v*v)-w*(a0*u+b0*v-u*x-v*y-w*z))*(1.0-cos(thita))
                            +L0*z*cos(thita)+sqrt(L0)*(-b0*u+a0*v-v*x+u*y)*sin(thita)  );
                        convertOrthToFrac(orth1, frac1);
    
                        ia1 = floor( frac1[0] * ANUM) - anum0;
                        ib1 = floor( frac1[1] * BNUM) - bnum0;
                        ic1 = floor( frac1[2] * CNUM) - cnum0;
    
                        if(ia1 < 0 || ia1 > anum2 || ib1 < 0 || ib1 > bnum2 || ic1 < 0 || ic1 > cnum2){
                            frac[0] = (ia1 + anum0 + 0.5) / float(ANUM);
                            frac[1] = (ib1 + bnum0 + 0.5) / float(BNUM);
                            frac[2] = (ic1 + cnum0 + 0.5) / float(CNUM);
                            adjustFracCoorIntoAsu(frac);
                            int ia2 = floor(frac[0] * ANUM);
                            int ib2 = floor(frac[1] * BNUM);
                            int ic2 = floor(frac[2] * CNUM);
                            ncsMaskDens[ia][ib][ic] = ncsMaskDens[ia][ib][ic] + densAsuTemp[ia2][ib2][ic2];
                            incsAvg++;
                            continue;
                        }
    
                        //if(ncsMask[ia1][ib1][ic1] >= 1){
                            ncsMaskDens[ia][ib][ic] = ncsMaskDens[ia][ib][ic] + ncsMaskDensTemp[ia1][ib1][ic1];
                            incsAvg++;
                        //}
    
                    }
                    ncsMaskDens[ia][ib][ic] = ncsMaskDens[ia][ib][ic] / float(incsAvg);
                }
            }
        }
        //compute density correlation coefficient between densAsuTemp and ncsMaskDens
        for(int ia = 0; ia <= anum2; ++ia){
            for(int ib = 0; ib <= bnum2; ++ib){
                for(int ic = 0; ic <= cnum2; ++ic){
                    if(ncsMask[ia][ib][ic] >= 1){
                        ia1 = ia + anum0;
                        ib1 = ib + bnum0;
                        ic1 = ic + cnum0;
                        frac[0] = (ia1 + 0.5) / float(ANUM);
                        frac[1] = (ib1 + 0.5) / float(BNUM);
                        frac[2] = (ic1 + 0.5) / float(CNUM);
                        adjustFracCoorIntoAsu(frac);
                        ia1 = floor(frac[0] * ANUM);
                        ib1 = floor(frac[1] * BNUM);
                        ic1 = floor(frac[2] * CNUM);
                        if(protMaskAsu[ia1][ib1][ic1] == 1){
                            ncsDensDiff += abs(densAsuTemp[ia1][ib1][ic1] - ncsMaskDens[ia][ib][ic]);
                        }
                    }
                }
            }
        }
        if(icand == 0 || ncsDensDiff < ncsDensDiffMin){
            ncsDensDiffMin = ncsDensDiff;
            icand0 = icand;
        }
        //cout << icand << '\t' << ncsDensDiff << ',' <<'\t'<< icand0 << '\t'<< ncsDensDiffMin << endl;
    } //icand

    for(int i = 0; i < 3; ++i){
        pointAOnNcsAxis[i] = pointAOnNcsAxisCand[icand0][i];
        vectOnNcsAxis[i] = vectOnNcsAxisCand[icand0][i];
        pointBOnNcsAxis[i] = pointAOnNcsAxis[i] + vectOnNcsAxis[i];
    }

    return;
}

/******************************************************************************/

void applyNcsAvgMatrix()
{
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;

    //apply ncs operator
    float rota[3][3], tran[3];
    float orth[3], orth1[3], frac[3], frac1[3];
    int ia1, ib1, ic1;
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                ncsMaskDensTemp[ia][ib][ic] = 0;
                ncsMaskDens[ia][ib][ic] = 0;
                //if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    ncsMaskDensTemp[ia][ib][ic] = densAsuTemp[ia1][ib1][ic1];
                    ncsMaskDens[ia][ib][ic] = densAsuTemp[ia1][ib1][ic1]; //corresponds to incs = 0
                //}
            }
        }
    }

    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){

                if(ncsMask[ia][ib][ic] == 0) continue;
                frac[0] = (ia + anum0 + 0.5) / float(ANUM);
                frac[1] = (ib + bnum0 + 0.5) / float(BNUM);
                frac[2] = (ic + cnum0 + 0.5) / float(CNUM);
                adjustFracCoorIntoAsu(frac);
                ia1 = floor(frac[0] * ANUM);
                ib1 = floor(frac[1] * BNUM);
                ic1 = floor(frac[2] * CNUM);
                if(protMaskAsu[ia1][ib1][ic1] == 0)continue;

                frac[0] = (ia + anum0 + 0.5) / float(ANUM);
                frac[1] = (ib + bnum0 + 0.5) / float(BNUM);
                frac[2] = (ic + cnum0 + 0.5) / float(CNUM);
                convertFracToOrth(frac, orth);
                int incsAvg = 1; //have considered incs = 0
                for(int incs = 1; incs <= numNcsOper - 1; ++incs) { //start from incs = 1
                    for(int i = 0; i <= 2; ++i){
                        for(int j = 0; j <= 2; ++j){
                            rota[i][j] = ncsOper[incs][i][j]; //load rotation matrix
                        } 
                        tran[i] = ncsOper[incs][i][3]; //load rotation matrix
                    }
                    for( int i = 0; i <= 2 ; ++i ){
                        orth1[i] = 0;
                        for( int j = 0; j <= 2 ; ++j ){
                            orth1[i] = orth1[i] + rota[i][j] * orth[j]; //rotation
                        }
                        orth1[i] = orth1[i] + tran[i]; //translation
                    }
                    convertOrthToFrac(orth1, frac1);
                    ia1 = floor( frac1[0] * ANUM) - anum0;
                    ib1 = floor( frac1[1] * BNUM) - bnum0;
                    ic1 = floor( frac1[2] * CNUM) - cnum0;

                    if(ia1 < 0 || ia1 > anum2 || ib1 < 0 || ib1 > bnum2 || ic1 < 0 || ic1 > cnum2){
                        frac[0] = (ia1 + anum0 + 0.5) / float(ANUM);
                        frac[1] = (ib1 + bnum0 + 0.5) / float(BNUM);
                        frac[2] = (ic1 + cnum0 + 0.5) / float(CNUM);
                        adjustFracCoorIntoAsu(frac);
                        int ia2 = floor(frac[0] * ANUM);
                        int ib2 = floor(frac[1] * BNUM);
                        int ic2 = floor(frac[2] * CNUM);
                        ncsMaskDens[ia][ib][ic] = ncsMaskDens[ia][ib][ic] + densAsuTemp[ia2][ib2][ic2];
                        incsAvg++;
                        continue;
                    }

                    //if(ncsMask[ia1][ib1][ic1] >= 1){
                        ncsMaskDens[ia][ib][ic] = ncsMaskDens[ia][ib][ic] + ncsMaskDensTemp[ia1][ib1][ic1];
                        incsAvg++;
                    //}
                }
                ncsMaskDens[ia][ib][ic] = ncsMaskDens[ia][ib][ic] / float(incsAvg);
            }
        }
    }
    //update densAsuTemp
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    if(protMaskAsu[ia1][ib1][ic1] == 1){
                        densAsuTemp[ia1][ib1][ic1] = ncsMaskDens[ia][ib][ic];
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void applyNcsAvgAxis()
{
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;

    //apply ncs operator
    float rota[3][3], tran[3];
    float orth[3], orth1[3], frac[3], frac1[3];
    int ia1, ib1, ic1;
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                ncsMaskDensTemp[ia][ib][ic] = 0;
                ncsMaskDens[ia][ib][ic] = 0;
                //if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    ncsMaskDensTemp[ia][ib][ic] = densAsuTemp[ia1][ib1][ic1];
                    ncsMaskDens[ia][ib][ic] = densAsuTemp[ia1][ib1][ic1]; //corresponds to incs = 0
                //}
            }
        }
    }
    float a0 = pointAOnNcsAxis[0];
    float b0 = pointAOnNcsAxis[1];
    float c0 = pointAOnNcsAxis[2];
    float u = vectOnNcsAxis[0];
    float v = vectOnNcsAxis[1];
    float w = vectOnNcsAxis[2];
    float L0 = u*u + v*v + w*w;
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){

                if(ncsMask[ia][ib][ic] == 0) continue;
                frac[0] = (ia + anum0 + 0.5) / float(ANUM);
                frac[1] = (ib + bnum0 + 0.5) / float(BNUM);
                frac[2] = (ic + cnum0 + 0.5) / float(CNUM);
                adjustFracCoorIntoAsu(frac);
                ia1 = floor(frac[0] * ANUM);
                ib1 = floor(frac[1] * BNUM);
                ic1 = floor(frac[2] * CNUM);
                if(protMaskAsu[ia1][ib1][ic1] == 0)continue;

                frac[0] = (ia + anum0 + 0.5) / float(ANUM);
                frac[1] = (ib + bnum0 + 0.5) / float(BNUM);
                frac[2] = (ic + cnum0 + 0.5) / float(CNUM);
                convertFracToOrth(frac, orth);
                float x = orth[0];
                float y = orth[1];
                float z = orth[2];
                int incsAvg = 1; //have considered incs = 0
                for(int incs = 1; incs <= numNcsOper - 1; ++incs) { //start from incs = 1
                    float thita = 2.0 * PI / float(numNcsOper) * float(incs);
                    orth1[0] = 1.0/L0*(  (a0*(v*v+w*w)-u*(b0*v+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
                        +L0*x*cos(thita)+sqrt(L0)*(-c0*v+b0*w-w*y+v*z)*sin(thita)  );
                    orth1[1] = 1.0/L0*(  (b0*(u*u+w*w)-v*(a0*u+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
                        +L0*y*cos(thita)+sqrt(L0)*( c0*u-a0*w+w*x-u*z)*sin(thita)  );
                    orth1[2] = 1.0/L0*(  (c0*(u*u+v*v)-w*(a0*u+b0*v-u*x-v*y-w*z))*(1.0-cos(thita))
                        +L0*z*cos(thita)+sqrt(L0)*(-b0*u+a0*v-v*x+u*y)*sin(thita)  );
                    convertOrthToFrac(orth1, frac1);

                    /*
                    //8 point trilinear interpolation
                    float x1 = frac1[0];
                    float y1 = frac1[1];
                    float z1 = frac1[2];

                    int ia_0 = floor(x1 * ANUM - 0.5);
                    int ib_0 = floor(y1 * BNUM - 0.5);
                    int ic_0 = floor(z1 * CNUM - 0.5);
                    if(ia_0 < 0 || ia_0 > anum2 || ib_0 < 0 || ib_0 > bnum2 || ic_0 < 0 || ic_0 > cnum2)continue;

                    int ia_1 = ia_0 + 1;
                    int ib_1 = ib_0 + 1;
                    int ic_1 = ic_0 + 1;
                    if(ia_1 < 0 || ia_1 > anum2 || ib_1 < 0 || ib_1 > bnum2 || ic_1 < 0 || ic_1 > cnum2)continue;

                    float x_0 = (ia_0 + 0.5) / float(ANUM);
                    float y_0 = (ib_0 + 0.5) / float(BNUM);
                    float z_0 = (ic_0 + 0.5) / float(CNUM);
                    float x_1 = (ia_1 + 0.5) / float(ANUM);
                    float y_1 = (ib_1 + 0.5) / float(BNUM);
                    float z_1 = (ic_1 + 0.5) / float(CNUM);

                    float xd = (x1 - x_0) / (x_1 - x_0);
                    float yd = (y1 - y_0) / (y_1 - y_0);
                    float zd = (z1 - z_0) / (z_1 - z_0);

                    float density_000 = ncsMaskDensTemp[ia_0][ib_0][ic_0];
                    float density_100 = ncsMaskDensTemp[ia_1][ib_0][ic_0];
                    float density_010 = ncsMaskDensTemp[ia_0][ib_1][ic_0];
                    float density_110 = ncsMaskDensTemp[ia_1][ib_1][ic_0];
                    float density_001 = ncsMaskDensTemp[ia_0][ib_0][ic_1];
                    float density_101 = ncsMaskDensTemp[ia_1][ib_0][ic_1];
                    float density_011 = ncsMaskDensTemp[ia_0][ib_1][ic_1];
                    float density_111 = ncsMaskDensTemp[ia_1][ib_1][ic_1];
                    
                    float density_00 = density_000 * (1.0 - xd) + density_100 * xd;
                    float density_10 = density_010 * (1.0 - xd) + density_110 * xd;
                    float density_01 = density_001 * (1.0 - xd) + density_101 * xd;
                    float density_11 = density_011 * (1.0 - xd) + density_111 * xd;
           
                    float density_0 = density_00 * (1.0 - yd) + density_10 * yd;
                    float density_1 = density_01 * (1.0 - yd) + density_11 * yd;
                        
                    float density = density_0 * (1.0 - zd) + density_1 * zd;

                    ncsMaskDens[ia][ib][ic] = ncsMaskDens[ia][ib][ic] + density;
                    incsAvg++;
                    */

                    ia1 = floor( frac1[0] * ANUM) - anum0;
                    ib1 = floor( frac1[1] * BNUM) - bnum0;
                    ic1 = floor( frac1[2] * CNUM) - cnum0;

                    if(ia1 < 0 || ia1 > anum2 || ib1 < 0 || ib1 > bnum2 || ic1 < 0 || ic1 > cnum2){
                        frac[0] = (ia1 + anum0 + 0.5) / float(ANUM);
                        frac[1] = (ib1 + bnum0 + 0.5) / float(BNUM);
                        frac[2] = (ic1 + cnum0 + 0.5) / float(CNUM);
                        adjustFracCoorIntoAsu(frac);
                        int ia2 = floor(frac[0] * ANUM);
                        int ib2 = floor(frac[1] * BNUM);
                        int ic2 = floor(frac[2] * CNUM);
                        ncsMaskDens[ia][ib][ic] = ncsMaskDens[ia][ib][ic] + densAsuTemp[ia2][ib2][ic2];
                        incsAvg++;
                        continue;
                    }

                    //if(ncsMask[ia1][ib1][ic1] >= 1){
                        ncsMaskDens[ia][ib][ic] = ncsMaskDens[ia][ib][ic] + ncsMaskDensTemp[ia1][ib1][ic1];
                        incsAvg++;
                    //}

                }
                ncsMaskDens[ia][ib][ic] = ncsMaskDens[ia][ib][ic] / float(incsAvg);
            }
        }
    }
    //update densAsuTemp
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    if(protMaskAsu[ia1][ib1][ic1] == 1){
                        densAsuTemp[ia1][ib1][ic1] = ncsMaskDens[ia][ib][ic];
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void applyImprNcsAvgMatrix()
{
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;

    //apply ncs operator
    float rota[3][3], tran[3];
    float orth[3], orth1[3], frac[3], frac1[3];
    int ia1, ib1, ic1;
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                ncsMaskDensTemp[ia][ib][ic] = 0;
                ncsMaskDens[ia][ib][ic] = 0;
                //if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    ncsMaskDensTemp[ia][ib][ic] = densAsuTemp[ia1][ib1][ic1];
                    ncsMaskDens[ia][ib][ic] = densAsuTemp[ia1][ib1][ic1]; //before applying improper NCS avg
                //}
            }
        }
    }

    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){

                if(ncsMask[ia][ib][ic] == 0) continue;
                frac[0] = (ia + anum0 + 0.5) / float(ANUM);
                frac[1] = (ib + bnum0 + 0.5) / float(BNUM);
                frac[2] = (ic + cnum0 + 0.5) / float(CNUM);
                adjustFracCoorIntoAsu(frac);
                ia1 = floor(frac[0] * ANUM);
                ib1 = floor(frac[1] * BNUM);
                ic1 = floor(frac[2] * CNUM);
                if(protMaskAsu[ia1][ib1][ic1] == 0)continue;

                //ncsMask[ia][ib][ic] = 1 for mole A, = 2 for mole B
                int incs = ncsMask[ia][ib][ic] - 1;
                frac[0] = (ia + anum0 + 0.5) / float(ANUM);
                frac[1] = (ib + bnum0 + 0.5) / float(BNUM);
                frac[2] = (ic + cnum0 + 0.5) / float(CNUM);
                convertFracToOrth(frac, orth);
                for(int i = 0; i <= 2; ++i){
                    for(int j = 0; j <= 2; ++j){
                        rota[i][j] = ncsOper[incs][i][j]; //load rotation matrix
                    } 
                    tran[i] = ncsOper[incs][i][3]; //load rotation matrix
                }
                for( int i = 0; i <= 2 ; ++i ){
                    orth1[i] = 0;
                    for( int j = 0; j <= 2 ; ++j ){
                        orth1[i] = orth1[i] + rota[i][j] * orth[j]; //rotation
                    }
                    orth1[i] = orth1[i] + tran[i]; //translation
                }
                convertOrthToFrac(orth1, frac1);
                ia1 = floor( frac1[0] * ANUM) - anum0;
                ib1 = floor( frac1[1] * BNUM) - bnum0;
                ic1 = floor( frac1[2] * CNUM) - cnum0;

                if(ia1 < 0 || ia1 > anum2 || ib1 < 0 || ib1 > bnum2 || ic1 < 0 || ic1 > cnum2){
                    frac[0] = (ia1 + anum0 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + bnum0 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + cnum0 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    int ia2 = floor(frac[0] * ANUM);
                    int ib2 = floor(frac[1] * BNUM);
                    int ic2 = floor(frac[2] * CNUM);
                    ncsMaskDens[ia][ib][ic] = (ncsMaskDensTemp[ia][ib][ic] + densAsuTemp[ia2][ib2][ic2]) / 2.0;
                    continue;
                }

                //if(ncsMask[ia1][ib1][ic1] == incs + 1 + 1){
                    ncsMaskDens[ia][ib][ic] = (ncsMaskDensTemp[ia][ib][ic] + ncsMaskDensTemp[ia1][ib1][ic1]) / 2.0;
                //}
            }
        }
    }

    //update densAsuTemp
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    if(protMaskAsu[ia1][ib1][ic1] == 1){
                        densAsuTemp[ia1][ib1][ic1] = ncsMaskDens[ia][ib][ic];
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void applyImprNcsAvgAxis()
{
    int anum2 = anum1 - anum0;
    int bnum2 = bnum1 - bnum0;
    int cnum2 = cnum1 - cnum0;

    //apply ncs operator
    float rota[3][3], tran[3];
    float orth[3], orth1[3], frac[3], frac1[3];
    int ia1, ib1, ic1;
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                ncsMaskDensTemp[ia][ib][ic] = 0;
                ncsMaskDens[ia][ib][ic] = 0;
                //if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    ncsMaskDensTemp[ia][ib][ic] = densAsuTemp[ia1][ib1][ic1];
                    ncsMaskDens[ia][ib][ic] = densAsuTemp[ia1][ib1][ic1]; //before applying improper NCS avg
                //}
            }
        }
    }
    float a0 = pointAOnNcsAxis[0];
    float b0 = pointAOnNcsAxis[1];
    float c0 = pointAOnNcsAxis[2];
    float u = vectOnNcsAxis[0];
    float v = vectOnNcsAxis[1];
    float w = vectOnNcsAxis[2];
    float L0 = u*u + v*v + w*w;

    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){

                if(ncsMask[ia][ib][ic] == 0) continue;
                frac[0] = (ia + anum0 + 0.5) / float(ANUM);
                frac[1] = (ib + bnum0 + 0.5) / float(BNUM);
                frac[2] = (ic + cnum0 + 0.5) / float(CNUM);
                adjustFracCoorIntoAsu(frac);
                ia1 = floor(frac[0] * ANUM);
                ib1 = floor(frac[1] * BNUM);
                ic1 = floor(frac[2] * CNUM);
                if(protMaskAsu[ia1][ib1][ic1] == 0)continue;

                //ncsMask[ia][ib][ic] = 1 for mole A, = 2 for mole B
                int incs = ncsMask[ia][ib][ic] - 1;
                frac[0] = (ia + anum0 + 0.5) / float(ANUM);
                frac[1] = (ib + bnum0 + 0.5) / float(BNUM);
                frac[2] = (ic + cnum0 + 0.5) / float(CNUM);
                convertFracToOrth(frac, orth);
                float x = orth[0];
                float y = orth[1];
                float z = orth[2];

                float thita = pow(-1.0, float(incs)) * rotaNcsAngl;

                orth1[0] = 1.0/L0*(  (a0*(v*v+w*w)-u*(b0*v+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
                    +L0*x*cos(thita)+sqrt(L0)*(-c0*v+b0*w-w*y+v*z)*sin(thita)  );
                orth1[1] = 1.0/L0*(  (b0*(u*u+w*w)-v*(a0*u+c0*w-u*x-v*y-w*z))*(1.0-cos(thita))
                    +L0*y*cos(thita)+sqrt(L0)*( c0*u-a0*w+w*x-u*z)*sin(thita)  );
                orth1[2] = 1.0/L0*(  (c0*(u*u+v*v)-w*(a0*u+b0*v-u*x-v*y-w*z))*(1.0-cos(thita))
                    +L0*z*cos(thita)+sqrt(L0)*(-b0*u+a0*v-v*x+u*y)*sin(thita)  );

                orth1[0] += pow(-1.0, float(incs)) * tranNcsDist * u / sqrt(L0);
                orth1[1] += pow(-1.0, float(incs)) * tranNcsDist * v / sqrt(L0);
                orth1[2] += pow(-1.0, float(incs)) * tranNcsDist * w / sqrt(L0);

                convertOrthToFrac(orth1, frac1);
                ia1 = floor( frac1[0] * ANUM) - anum0;
                ib1 = floor( frac1[1] * BNUM) - bnum0;
                ic1 = floor( frac1[2] * CNUM) - cnum0;

                if(ia1 < 0 || ia1 > anum2 || ib1 < 0 || ib1 > bnum2 || ic1 < 0 || ic1 > cnum2){
                    frac[0] = (ia1 + anum0 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + bnum0 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + cnum0 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    int ia2 = floor(frac[0] * ANUM);
                    int ib2 = floor(frac[1] * BNUM);
                    int ic2 = floor(frac[2] * CNUM);
                    ncsMaskDens[ia][ib][ic] = (ncsMaskDensTemp[ia][ib][ic] + densAsuTemp[ia2][ib2][ic2]) / 2.0;
                    continue;
                }

                //if(ncsMask[ia1][ib1][ic1] == incs + 1 + 1){
                    ncsMaskDens[ia][ib][ic] = (ncsMaskDensTemp[ia][ib][ic] + ncsMaskDensTemp[ia1][ib1][ic1]) / 2.0;
                //}
            }
        }
    }

    //update densAsuTemp
    for(int ia = 0; ia <= anum2; ++ia){
        for(int ib = 0; ib <= bnum2; ++ib){
            for(int ic = 0; ic <= cnum2; ++ic){
                if(ncsMask[ia][ib][ic] >= 1){
                    ia1 = ia + anum0;
                    ib1 = ib + bnum0;
                    ic1 = ic + cnum0;
                    frac[0] = (ia1 + 0.5) / float(ANUM);
                    frac[1] = (ib1 + 0.5) / float(BNUM);
                    frac[2] = (ic1 + 0.5) / float(CNUM);
                    adjustFracCoorIntoAsu(frac);
                    ia1 = floor(frac[0] * ANUM);
                    ib1 = floor(frac[1] * BNUM);
                    ic1 = floor(frac[2] * CNUM);
                    if(protMaskAsu[ia1][ib1][ic1] == 1){
                        densAsuTemp[ia1][ib1][ic1] = ncsMaskDens[ia][ib][ic];
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void applyHIO(int iter)
{
    float hioDensLimit;
    if(iter < (numIter-numIterLimitDens-numIterSolvFlat)) {
        hioDensLimit = hioDensLimitIni;
    } else if(iter < (numIter-numIterSolvFlat)) {
        hioDensLimit = hioDensLimitIni - (hioDensLimitIni-hioDensLimitFina)
           *float( iter-(numIter-numIterLimitDens-numIterSolvFlat) )/float(numIterLimitDens);
    } else {
        hioDensLimit = hioDensLimitFina;
    }

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(protMaskAsu[ia][ib][ic] == 1) {
                        densAsu[ia][ib][ic] = densAsuTemp[ia][ib][ic];
                    } else {
                        densAsu[ia][ib][ic] = densAsu[ia][ib][ic] - hioBeta * densAsuTemp[ia][ib][ic];
                        //turn off HIO at the end of iterations and before solvent flattening
                        if(densAsu[ia][ib][ic] > hioDensLimit) {
                            densAsu[ia][ib][ic] = hioDensLimit;
                        } else if(densAsu[ia][ib][ic] < -hioDensLimit) {
                            densAsu[ia][ib][ic] = -hioDensLimit;
                        }
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void applyCHIO(int iter)
{
    float hioDensLimit;
    if(iter < (numIter-numIterLimitDens-numIterSolvFlat)) {
        hioDensLimit = hioDensLimitIni;
    } else if(iter < (numIter-numIterSolvFlat)) {
        hioDensLimit = hioDensLimitIni - (hioDensLimitIni-hioDensLimitFina)
           *float( iter-(numIter-numIterLimitDens-numIterSolvFlat) )/float(numIterLimitDens);
    } else {
        hioDensLimit = hioDensLimitFina;
    }

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(protMaskAsu[ia][ib][ic] == 1 
                        && weigAvgDensAsu[ia][ib][ic] >= 1.1 * densCutoff) {
                        densAsu[ia][ib][ic] = densAsuTemp[ia][ib][ic];
                    } else if(protMaskAsu[ia][ib][ic] == 1
                        && weigAvgDensAsu[ia][ib][ic] < 1.1 * densCutoff){
                        densAsu[ia][ib][ic] = densAsu[ia][ib][ic] - (1.0 - hioAlpha) / hioAlpha * densAsuTemp[ia][ib][ic];
                    } else {
                        densAsu[ia][ib][ic] = densAsu[ia][ib][ic] - hioBeta * densAsuTemp[ia][ib][ic];
                        //turn off HIO at the end of iterations and before solvent flattening
                        if(densAsu[ia][ib][ic] > hioDensLimit) {
                            densAsu[ia][ib][ic] = hioDensLimit;
                        } else if(densAsu[ia][ib][ic] < -hioDensLimit) {
                            densAsu[ia][ib][ic] = -hioDensLimit;
                        }
                    }
                }
            }
        }
    }
    return;
}


/******************************************************************************/

void applyHistMatch()

{
    assignDensAsuIntoSmallBin("prot");

    combineSmallBinIntoLargeBin("prot");

    for(int ilargeBin = 0; ilargeBin <= numLargeBin - 1; ++ilargeBin) {
        hma[ilargeBin] = ( stdBoundLargeBin[ilargeBin+1] - stdBoundLargeBin[ilargeBin] )
                       /( boundLargeBin[ilargeBin+1] - boundLargeBin[ilargeBin] );

        hmb[ilargeBin] = ( stdBoundLargeBin[ilargeBin] * boundLargeBin[ilargeBin+1]
                         -stdBoundLargeBin[ilargeBin+1] * boundLargeBin[ilargeBin] )
                       /( boundLargeBin[ilargeBin+1] - boundLargeBin[ilargeBin] );
    }

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if( protMaskAsu[ia][ib][ic] == 1 ) {
                        float density = densAsu[ia][ib][ic];
                        for(int ilargeBin=0; ilargeBin <= numLargeBin - 1; ++ilargeBin) {
                            if(density >= boundLargeBin[ilargeBin] && density <= boundLargeBin[ilargeBin+1]) {
                                densAsu[ia][ib][ic] = hma[ilargeBin] * density+hmb[ilargeBin];
                                break; //go to next grid point on protMaskAsu
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void applySolvFlat()
{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(protMaskAsu[ia][ib][ic] == 0) { //solvent
                        densAsu[ia][ib][ic] = 0;
                    }
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void generateTruePhaseForAllOrigTran(string outputStatus, string fileNameString)

{
    //compute allStruFactAmpl from uniqStruFactAmplCal
    computeAllStruFactAmplFromUniqStruFactAmplCal();

    //compute allStruFactPhase from uniqStruFactPhase
    computeAllStruFactPhaseFromUniqStruFactPhase();

    //prepare for fft
    setHalfGridPhase(); //set half grid phase shift for fft, initialize halfGridPosi and halfGridNega

    //apply forward fft, input allStruFactAmpl and allStruFactPhase, output densUnitCell
    applyFFTStruFactToDens();

    //reserve initial densUnitCell
    for( int ia=0; ia<ANUM; ++ia ) {
        for( int ib=0; ib<BNUM; ++ib ) {
            for( int ic=0; ic<CNUM; ++ic ) {
                densUnitCellTemp[ia][ib][ic] = densUnitCell[ia][ib][ic];
            }
        }
    }

    for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {

        float rota[3][3], tran[3];
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = origTran[iorig][i][j]; //load rotation matrix
            }
            tran[i] = origTran[iorig][i][3]; //load rotation matrix
        }

        for( int ia=0; ia<ANUM; ++ia ) {
            for( int ib=0; ib<BNUM; ++ib ) {
                for( int ic=0; ic<CNUM; ++ic ) {
                    float frac[3], temp[3];
                    temp[0] = (ia + 0.5) / float(ANUM);
                    temp[1] = (ib + 0.5) / float(BNUM);
                    temp[2] = (ic + 0.5) / float(CNUM);
                    for( int i=0; i<=2 ; ++i )
                    {
                        frac[i] = 0;
                        for( int j=0; j<=2 ; ++j )
                        {                        //frac = matmul(rotation,temp) + translation;
                            frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                        }
                        frac[i] = frac[i] + tran[i];
                    }
                    adjustFracCoorIntoUnitCell(frac);
                    int ia0 = floor(frac[0] * ANUM);
                    int ib0 = floor(frac[1] * BNUM);
                    int ic0 = floor(frac[2] * CNUM);
                    densUnitCell[ia0][ib0][ic0] = densUnitCellTemp[ia][ib][ic];
                }
            }
        }

        //apply backward fft, input densUnitCell, output allStruFactAmpl and allStruFactPhase
        applyFFTDensToStruFact();

        //catch uniqStruFactAmplCal and uniqStruFactPhase from allStruFactAmpl and allStruFactPhase
        catchUniqStruFactAmplCalAndPhase();

        ofstream outputFile;
        size_t found = outputStatus.find("not");
        if(found == string::npos) {    //output unique sf phases for current allowd origin translation
            stringstream fileNameStream;
            fileNameStream << fileNameString << iorig << ".txt"; 
            outputFile.open( fileNameStream.str().c_str() );
            if(!outputFile.is_open()){
                cout << "Can't open file " << fileNameStream.str() << endl;
                exit(EXIT_FAILURE);
            }
        } else { //find "not", don't output true phases
        }

        for( int IHLOOPUNIQREFL ) {
            for( int IKLOOPUNIQREFL ) {
                for( int ILLOOPUNIQREFL ) {
                    int h = ih;
                    int k = ik;
                    int l = il;
                    float resolution = uniqStruFactReso[ih][ik][il];
                    if(resolution >= resoCutoff) { 
                        truePhaseOrigTran[ih][ik][il][iorig] = uniqStruFactPhase[ih][ik][il];
                        if(found == string::npos){
                            outputFile << setfill(' ') << setw(5) << h
                            << setfill(' ') << setw(5) << k    << setfill(' ') << setw(5) << l
                            << setfill(' ') << setw(8) << setprecision(1)  << fixed 
                            << uniqStruFactPhase[ih][ik][il] / PI * 180.0 << endl;
                        }
                    }
                }
            }
        }
        outputFile.close();
    }

    //reinitialize densUnitCell and densUnitCellTemp
    for( int ia=0; ia<ANUM; ++ia ) {
        for( int ib=0; ib<BNUM; ++ib ) {
            for( int ic=0; ic<CNUM; ++ic ) {
                densUnitCell[ia][ib][ic] = 0;
                densUnitCellTemp[ia][ib][ic] = 0;  
            }
        }
    }

    //reinitialize allStruFactAmpl and allStruFactPhase to zero
    for( int ih=0; ih<HNUM; ++ih ) {
        for( int ik=0; ik<KNUM; ++ik ) {
            for( int il=0; il<LNUM; ++il ) {
                allStruFactAmpl[ih][ik][il] = 0;
                allStruFactPhase[ih][ik][il] = 0;
            }
        }
    }

    //reinitialize variables of unique structure factor
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                uniqStruFactAmplCal[ih][ik][il] = 0;
                uniqStruFactPhase[ih][ik][il] = 0;
            }
        }
    }
    return;
}

/******************************************************************************/

void generateTrueProtMaskForAllOrigTran()

{
    //compute allStruFactAmpl from uniqStruFactAmplCal
    computeAllStruFactAmplFromUniqStruFactAmplCal();

    //compute allStruFactPhase from uniqStruFactPhase
    computeAllStruFactPhaseFromUniqStruFactPhase();

    //prepare for fft
    setHalfGridPhase(); //set half grid phase shift for fft, initialize halfGridPosi and halfGridNega

    //apply forward fft, input allStruFactAmpl and allStruFactPhase, output densUnitCell
    applyFFTStruFactToDens();

    //compute weighted average density in unit cell, input densUnitCell, output weigAvgDensUnitCell
    computeWeigAvgDensUnitCell(0.6);

    //get weigAvgDensAsu from weigAvgDensUnitCell, input weigAvgDensUnitCell, output weigAvgDensAsu
    catchWeigAvgDensAsuFromWeigAvgDensUnitCell();

    //compute a densCutoff value on weigAvgDensAsu by bins, input weigAvgDensAsu, return densCutoff
    densCutoff = computeWeigAvgDensCutoffAsu(1.0-solvCont);

    //make a protein mask in asu, input weigAvgDensAsu, output protMaskAsu
    makeProtMaskAsu(densCutoff);

    for(int iorig = 0; iorig <= numOrigTran - 1; ++iorig) {

        float rota[3][3], tran[3];
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = origTran[iorig][i][j]; //load rotation matrix
            }
            tran[i] = origTran[iorig][i][3]; //load rotation matrix
        }

        for( int IALOOPASU ) {
            float x = (ia + 0.5) / float(ANUM);
            for( int IBLOOPASU ) {
                float y = (ib + 0.5) / float(BNUM);
                for( int ICLOOPASU ) {
                    float z = (ic + 0.5) / float(CNUM);
                    if( ASUCOND ) {
                        float frac[3], temp[3];
                        temp[0] = (ia + 0.5) / float(ANUM);
                        temp[1] = (ib + 0.5) / float(BNUM);
                        temp[2] = (ic + 0.5) / float(CNUM);
                        for( int i=0; i<=2 ; ++i )
                        {
                            frac[i] = 0;
                            for( int j=0; j<=2 ; ++j )
                            {                        //frac = matmul(rotation,temp) + translation;
                                frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                            }
                            frac[i] = frac[i] + tran[i];
                        }
                        adjustFracCoorIntoAsu(frac);
                        int ia0 = floor(frac[0] * ANUM);
                        int ib0 = floor(frac[1] * BNUM);
                        int ic0 = floor(frac[2] * CNUM);
                        trueProtMaskOrigTran[ia0][ib0][ic0][iorig] = protMaskAsu[ia][ib][ic];
                    }
                }
            }
        }
    }

    //reinitialize densUnitCell and weigAvgDensUnitCell
    for( int ia=0; ia<ANUM; ++ia ) {
        for( int ib=0; ib<BNUM; ++ib ) {
            for( int ic=0; ic<CNUM; ++ic ) {
                densUnitCell[ia][ib][ic] = 0;
                weigAvgDensUnitCell[ia][ib][ic] = 0;
            }
        }
    }

    //reinitialize weigAvgDensAsu and protMaskAsu
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    weigAvgDensAsu[ia][ib][ic] = 0;
                    protMaskAsu[ia][ib][ic] = 0;
                }
            }
        }
    }

    //reinitialize allStruFactAmpl and allStruFactPhase to zero
    for( int ih=0; ih<HNUM; ++ih ) {
        for( int ik=0; ik<KNUM; ++ik ) {
            for( int il=0; il<LNUM; ++il ) {
                allStruFactAmpl[ih][ik][il] = 0;
                allStruFactPhase[ih][ik][il] = 0;
            }
        }
    }

    //reinitialize variables of unique structure factor
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                uniqStruFactAmplCal[ih][ik][il] = 0;
                uniqStruFactPhase[ih][ik][il] = 0;
            }
        }
    }

    return;

}

/******************************************************************************/

void outputRvalueResoShell(string fileNameString, float RvalueResoShell[numIter][numResoShell])
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int iter = 0; iter <= numIter - 1; ++iter) {
        outputFile << setfill(' ') << setw(8) << iter;
        for(int ishell = 0; ishell <= numResoShell-1; ++ishell) {
            outputFile << ',' << setfill(' ') << setw(8) << setprecision(3) << fixed << RvalueResoShell[iter][ishell];
        }
        outputFile << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputMeanPhaseErrorResoShell(string fileNameString, int iorig0)
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int iter = 0; iter <= numIter - 1; ++iter) {
        outputFile << setfill(' ') << setw(8) << iter;
        for(int ishell = 0; ishell <= numResoShell-1; ++ishell) {
            outputFile << ',' << setfill(' ') << setw(8) << setprecision(1) << fixed 
                << meanPhaseErrorOrigTranResoShell[iter][ishell][iorig0];
        }
        outputFile << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputMeanPhaseErrorResoSphere(string fileNameString, int iorig0)
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int iter = 0; iter <= numIter - 1; ++iter) {
        outputFile << setfill(' ') << setw(8) << iter;
        for(int ishell = 0; ishell <= numResoShell-1; ++ishell) {
            outputFile << ',' << setfill(' ') << setw(8) << setprecision(1) << fixed 
                << meanPhaseErrorOrigTranResoSphere[iter][ishell][iorig0];
        }
        outputFile << endl;
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputCorrCoef(string fileNameString, int iorig0)
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int iter = 0; iter <= numIter - 1; ++iter) {
        outputFile << setfill(' ') << setw(8) << iter;
        for(int ishell = 0; ishell <= numResoShell-1; ++ishell) {
            outputFile << ',' << setfill(' ') << setw(8) << setprecision(3) << fixed << corrCoefOrigTran[iter][ishell][iorig0];
        }
        outputFile << endl;
    }
    outputFile.close();

    return;
}

/******************************************************************************/

void outputProtMaskMatch(string fileNameString, int iorig0)
{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for(int iter = 0; iter <= numIter - 1; ++iter) {
        outputFile << setfill(' ') << setw(8) << iter
        << ',' << setfill(' ') << setw(8) << setprecision(3) << fixed << protMaskMatchOrigTran[iter][iorig0] << endl;
    }
    outputFile.close();

    return;
}

/******************************************************************************/

void outputUniqStruFactPhaseAndAmplCal(string fileNameString)

{
    ofstream outputFile( fileNameString.c_str() );
    if(!outputFile.is_open()) {
        cout << "Can't open file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution >= resoCutoff) {
                    int obsStatus = uniqStruFactStat[ih][ik][il];
                    float calsf = uniqStruFactAmplCal[ih][ik][il];
                    float phase = uniqStruFactPhase[ih][ik][il];
                    outputFile << setfill(' ') << setw(5) << h
                        << setfill(' ') << setw(5) << k
                        << setfill(' ') << setw(5) << l
                        << setfill(' ') << setw(2) << obsStatus
                        << setfill(' ') << setw(11) << setprecision(1) << fixed << calsf
                        << setfill(' ') << setw(8) << setprecision(1) << fixed << phase/PI*180.0
                        << endl;
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputCifFileForCoot(string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    //header of cif file
    outputFile << "data_r" << PDB_CODE << "sf" << endl;
    outputFile << '#' << endl;
    outputFile << "_cell.length_a     " << setfill(' ') << setw(9) << setprecision(3) << fixed << CRYST_A << endl;
    outputFile << "_cell.length_b     " << setfill(' ') << setw(9) << setprecision(3) << fixed << CRYST_B << endl;
    outputFile << "_cell.length_c     " << setfill(' ') << setw(9) << setprecision(3) << fixed << CRYST_C << endl;
    outputFile << "_cell.angle_alpha  " << setfill(' ') << setw(9) << setprecision(2) << fixed << CRYST_ALPHA/PI*180 << endl;
    outputFile << "_cell.angle_beta   " << setfill(' ') << setw(9) << setprecision(2) << fixed << CRYST_BETA/PI*180 << endl;
    outputFile << "_cell.angle_gamma  " << setfill(' ') << setw(9) << setprecision(2) << fixed << CRYST_GAMMA/PI*180 << endl;
    outputFile << '#' << endl;
    outputFile << "loop_" << endl;
    outputFile << "_symmetry_equiv.id" << endl;
    outputFile << "_symmetry_equiv.pos_as_xyz" << endl;
    outputFile << "1  ' X, Y,  Z'" << endl;
    outputFile << "2  '-X,-Y,  Z+1/2'" << endl;
    outputFile << "3  '-Y, X,  Z+3/4'" << endl;
    outputFile << "4  ' Y,-X,  Z+1/4''" << endl;
    outputFile << "5  '-X, Y, -Z'" << endl;
    outputFile << "6  ' X,-Y, -Z+1/2'" << endl;
    outputFile << "7  ' Y, X, -Z+1/4'" << endl;
    outputFile << "8  '-Y,-X, -Z+3/4''" << endl;
    outputFile << '#' << endl;
    outputFile << "loop_" << endl;
    outputFile << "_refln.crystal_id" << endl;
    outputFile << "_refln.wavelength_id" << endl;
    outputFile << "_refln.scale_group_code" << endl;
    outputFile << "_refln.status" << endl;
    outputFile << "_refln.index_h" << endl;
    outputFile << "_refln.index_k" << endl;
    outputFile << "_refln.index_l" << endl;
    outputFile << "_refln.F_meas_au" << endl;
    outputFile << "_refln.F_meas_sigma_au" << endl;
    outputFile << "_refln.F_calc" << endl;
    outputFile << "_refln.fom" << endl;
    outputFile << "_refln.phase_calc" << endl;
    //output data
    for( int IHLOOPUNIQREFL ) {
        for( int IKLOOPUNIQREFL ) {
            for( int ILLOOPUNIQREFL ) {
                int h = ih;
                int k = ik;
                int l = il;
                if(h == 0 && k == 0 && l == 0)continue;
                float resolution = uniqStruFactReso[ih][ik][il];
                if(resolution >= resoCutoff) {
                    int obsStatus = uniqStruFactStat[ih][ik][il];
                    char status;
                    float F_meas_au, F_meas_sigma_au, F_calc, phase_calc;
                    F_calc = uniqStruFactAmplCal[ih][ik][il];
                    phase_calc = uniqStruFactPhase[ih][ik][il];
                    if(obsStatus == 1) {
                        status = 'o';
                        F_meas_au = uniqStruFactAmplObs[ih][ik][il];
                        F_meas_sigma_au = uniqStruFactAmplObsSigma[ih][ik][il];
                    } else if(obsStatus == 2){
                        status = 'f';
                        F_meas_au = uniqStruFactAmplObs[ih][ik][il];
                        F_meas_sigma_au = uniqStruFactAmplObsSigma[ih][ik][il];
                    } else {
                        status = 'x';
                        F_meas_au = missStruFactAmpl[ih][ik][il];
                        F_meas_sigma_au = 0;
                    }
                    int crystal_id = 1, wavelength_id = 1, scale_group_code = 1;
                    float fom = 1.0;
                    outputFile << crystal_id << ' ' << wavelength_id << ' ' << scale_group_code << ' ' << status 
                        << setfill(' ') << setw(5) << h 
                        << setfill(' ') << setw(5) << k
                        << setfill(' ') << setw(5) << l 
                        << setfill(' ') << setw(9) << setprecision(1) << fixed << F_meas_au
                        << setfill(' ') << setw(7) << setprecision(1) << fixed << F_meas_sigma_au
                        << setfill(' ') << setw(11) << setprecision(1) << fixed << F_calc
                        << setfill(' ') << setw(6) << setprecision(2) << fixed << fom
                        << setfill(' ') << setw(8) << setprecision(1) << fixed << phase_calc/PI*180.0
                        << endl;
                }
            }
        }
    }

    outputFile.close();
    return;
}

/******************************************************************************/

void outputDensAsu(string fileNameString, float densCutoff)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    } else {
        outputFile << writeCryst1() << endl;
    }
    
    int iatom = 0;
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    float orth[3], frac[3];
                    float tempFact = densAsu[ia][ib][ic];
                    float occupancy = 1.0;
                    if(densAsu[ia][ib][ic] > densCutoff) {
                        iatom++;
                        if(iatom >= 9999) iatom = 0; //too many atoms

                        frac[0] = (ia + 0.5) / float(ANUM);
                        frac[1] = (ib + 0.5) / float(BNUM);
                        frac[2] = (ic + 0.5) / float(CNUM);
                        convertFracToOrth(frac, orth);

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
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputDensAsuTxt(string fileNameString)

{
    
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    outputFile << setfill(' ') << setw(11) << setprecision(5) << fixed << densAsu[ia][ib][ic] << endl;
                }
            }
        }
    }

    outputFile.close();
    return;
}

/******************************************************************************/

void outputWeigAvgDensAsu(string fileNameString, float densCutoff)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    } else {
        outputFile << writeCryst1() << endl;
    }
    
    int iatom = 0;
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    float orth[3], frac[3];
                    float tempFact = weigAvgDensAsu[ia][ib][ic];
                    float occupancy = 1.0;
                    if(weigAvgDensAsu[ia][ib][ic] > densCutoff) {
                        iatom++;
                        if(iatom >= 9999) iatom = 0; //too many atoms

                        frac[0] = (ia + 0.5) / float(ANUM);
                        frac[1] = (ib + 0.5) / float(BNUM);
                        frac[2] = (ic + 0.5) / float(CNUM);
                        convertFracToOrth(frac, orth);

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
    }
    outputFile.close();
    return;
}

/******************************************************************************/

void outputProtMaskAsuTxt(string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    outputFile << setfill(' ') << setw(1) << fixed << protMaskAsu[ia][ib][ic] << endl;
                }
            }
        }
    }

    outputFile.close();
    return;
}

/******************************************************************************/

void outputRvalueInResoShell(string fileNameString)
{
    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }

    for(int ishell = 29; ishell >= int(resoCutoff); --ishell ){
        float Rvalue = 0;
        float sumFobs = 0;
        for( int IHLOOPUNIQREFL ) {
            for( int IKLOOPUNIQREFL ) {
                for( int ILLOOPUNIQREFL ) {
                    int h = ih;
                    int k = ik;
                    int l = il;
                    float resolution = uniqStruFactReso[ih][ik][il];
                    if(resolution < resoCutoff) continue;
                    if(resolution >= ishell && resolution < ishell+1){
                        int obsStatus = uniqStruFactStat[ih][ik][il];
                        float obssf = uniqStruFactAmplObs[ih][ik][il];
                        float calsf = uniqStruFactAmplCal[ih][ik][il];
                        if(obsStatus == 1) {
                            Rvalue = Rvalue + abs( obssf - calsf );
                            sumFobs = sumFobs + obssf;
                        }
                    }
                }
            }
        }
        if(sumFobs < 0.01) continue; //no reflections in this shell
        Rvalue = Rvalue / sumFobs;
        outputFile << setfill(' ') << setw(5) << setprecision(1) << fixed << ishell+0.5
            << ',' << setfill(' ') << setw(8) << setprecision(2) << fixed << Rvalue << endl;
    }

    outputFile.close();
    return;
}

/******************************************************************************/

int checkOrigChoiceByProtMaskMatch(float protMaskMatch[numOrigTran])
{
    int iorig0 = 0;
    float protMaskMatchMax = 0;
    for(int iorig = 0; iorig <= numOrigTran-1; ++iorig) {
        if(protMaskMatch[iorig] > protMaskMatchMax){
            iorig0 = iorig;
            protMaskMatchMax = protMaskMatch[iorig];
        }
    }
    if(protMaskMatchMax > goodProtMaskThreshold && goodProtMaskFlag == 0) {
        fileNameString = "_WE_GOT_A_GOOD_PROT_MASK!!!";
        ofstream outputFile(fileNameString.c_str());
        if (!outputFile.is_open()) {
            cout << "Can't open the file " << fileNameString << endl;
            exit(EXIT_FAILURE);
        }
        outputFile.close();
        goodProtMaskFlag = 1;
    }
    return iorig0;
}

/******************************************************************************/

int checkOrigChoiceByMeanPhaseError(float meanPhaseError[numResoShell][numOrigTran])
{
    int iorig0;
    float meanPhaseErrorMin;
    //meanPhaseError[0][numOrigTran] contains the mean phase error of all observed data
    for(int iorig = 0; iorig <= numOrigTran-1; ++iorig) {
        if(iorig == 0){
            iorig0 = iorig;
            meanPhaseErrorMin = meanPhaseError[0][iorig];
        } else if(meanPhaseError[0][iorig] < meanPhaseErrorMin) {
            iorig0 = iorig;
            meanPhaseErrorMin = meanPhaseError[0][iorig];
        }
    }
    return iorig0;
}

/******************************************************************************/

void checkConvByMeanPhaseError(float meanPhaseError[numResoShell][numOrigTran])
{
    float meanPhaseErrorMin = 100.0;
    //meanPhaseError[0][numOrigTran] contains the mean phase error of all observed data
    for(int iorig = 0; iorig <= numOrigTran-1; ++iorig) {
        if(meanPhaseError[0][iorig] < meanPhaseErrorMin)
            meanPhaseErrorMin = meanPhaseError[0][iorig];
    }
    if(meanPhaseErrorMin < convThreshold && convFlag == 0) {
        fileNameString = "_WE_GOT_IT!!!";
        ofstream outputFile(fileNameString.c_str());
        if (!outputFile.is_open()) {
            cout << "Can't open the file " << fileNameString << endl;
            exit(EXIT_FAILURE);
        }
        outputFile.close();
        convFlag = 1;
    }
    return;
}

/******************************************************************************/

void translateOrigDensAsu(int iorig)
{
    float rota[3][3], tran[3];

    for(int i = 0; i <= 2; ++i)
    {
        for(int j = 0; j <= 2; ++j) 
        {
            rota[i][j] = origTran[iorig][i][j]; //load rotation matrix
        }
        tran[i] = origTran[iorig][i][3]; //load translation matrix
    }

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsuTemp[ia][ib][ic] = densAsu[ia][ib][ic];
                }
            }
        }
    }

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    float frac[3], temp[3];
                    temp[0] = (ia + 0.5) / float(ANUM);
                    temp[1] = (ib + 0.5) / float(BNUM);
                    temp[2] = (ic + 0.5) / float(CNUM);
                    for( int i=0; i<=2 ; ++i )
                    {
                        frac[i] = 0;
                        for( int j=0; j<=2 ; ++j )
                        {                        //frac = matmul(rotation,temp) + translation;
                            frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                        }
                        frac[i] = frac[i] + tran[i];
                    }
                    adjustFracCoorIntoAsu(frac);
                    int ia0 = floor(frac[0] * ANUM);
                    int ib0 = floor(frac[1] * BNUM);
                    int ic0 = floor(frac[2] * CNUM);
                    densAsu[ia0][ib0][ic0] = densAsuTemp[ia][ib][ic];
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void translateOrigProtMaskAsu(int iorig)
{
    float rota[3][3], tran[3];

    for(int i = 0; i <= 2; ++i)
    {
        for(int j = 0; j <= 2; ++j) 
        {
            rota[i][j] = origTran[iorig][i][j]; //load rotation matrix
        }
        tran[i] = origTran[iorig][i][3]; //load translation matrix
    }

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    protMaskAsuTemp[ia][ib][ic] = protMaskAsu[ia][ib][ic];
                }
            }
        }
    }

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    float frac[3], temp[3];
                    temp[0] = (ia + 0.5) / float(ANUM);
                    temp[1] = (ib + 0.5) / float(BNUM);
                    temp[2] = (ic + 0.5) / float(CNUM);
                    for( int i=0; i<=2 ; ++i )
                    {
                        frac[i] = 0;
                        for( int j=0; j<=2 ; ++j )
                        {                        //frac = matmul(rotation,temp) + translation;
                            frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                        }
                        frac[i] = frac[i] + tran[i];
                    }
                    adjustFracCoorIntoAsu(frac);
                    int ia0 = floor(frac[0] * ANUM);
                    int ib0 = floor(frac[1] * BNUM);
                    int ic0 = floor(frac[2] * CNUM);
                    protMaskAsu[ia0][ib0][ic0] = protMaskAsuTemp[ia][ib][ic];
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void accumulateDensAsu()
{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsuSum[ia][ib][ic] += densAsu[ia][ib][ic];
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void computeAvgDensAsu(int num)
{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsu[ia][ib][ic] = densAsuSum[ia][ib][ic] / float(num);
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void computeDensAsuAvg(int num)
{
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsuAvg[ia][ib][ic] = densAsuAvg[ia][ib][ic] * (float(num-1) / float(num)) 
                        + densAsu[ia][ib][ic] * ( 1.0 / float(num));
                }
            }
        }
    }
    return;
}

/******************************************************************************/

void translateOrigDensAsuMinDiff()
{
    float rota[3][3], tran[3];

    //find the origin to minimize density-difference
    int iorig0;
    float densDiff, densDiffMin;

    for(int iorig = 0; iorig <= numOrigTran - 1 ; ++iorig) {
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = origTran[iorig][i][j]; //load rotation matrix
            }
            tran[i] = origTran[iorig][i][3]; //load translation matrix
        }
        for( int IALOOPASU ) {
            float x = (ia + 0.5) / float(ANUM);
            for( int IBLOOPASU ) {
                float y = (ib + 0.5) / float(BNUM);
                for( int ICLOOPASU ) {
                    float z = (ic + 0.5) / float(CNUM);
                    if( ASUCOND ) {
                        float frac[3], temp[3];
                        temp[0] = (ia + 0.5) / float(ANUM);
                        temp[1] = (ib + 0.5) / float(BNUM);
                        temp[2] = (ic + 0.5) / float(CNUM);
                        for( int i=0; i<=2 ; ++i )
                        {
                            frac[i] = 0;
                            for( int j=0; j<=2 ; ++j )
                            {                        //frac = matmul(rotation,temp) + translation;
                                frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                            }
                            frac[i] = frac[i] + tran[i];
                        }
                        adjustFracCoorIntoAsu(frac);
                        int ia0 = floor(frac[0] * ANUM);
                        int ib0 = floor(frac[1] * BNUM);
                        int ic0 = floor(frac[2] * CNUM);
                        densAsuTemp[ia0][ib0][ic0] = densAsu[ia][ib][ic];
                    }
                }
            }
        }
        densDiff = 0;
        for( int IALOOPASU ) {
            float x = (ia + 0.5) / float(ANUM);
            for( int IBLOOPASU ) {
                float y = (ib + 0.5) / float(BNUM);
                for( int ICLOOPASU ) {
                    float z = (ic + 0.5) / float(CNUM);
                    if( ASUCOND ) {
                        densDiff += abs(densAsuAvg[ia][ib][ic] - densAsuTemp[ia][ib][ic]);
                    }
                }
            }
        }
        if(iorig == 0) {
            iorig0 = 0;
            densDiffMin = densDiff;
        } else if(densDiff < densDiffMin) {
            iorig0 = iorig;
            densDiffMin = densDiff;
        }
    }

    cout << "Compare to the best run, relative origin is "<< iorig0 << endl;

    //after find the proper origin, translate origin
    for(int i = 0; i <= 2; ++i)
    {
        for(int j = 0; j <= 2; ++j) 
        {
            rota[i][j] = origTran[iorig0][i][j]; //load rotation matrix
        }
        tran[i] = origTran[iorig0][i][3]; //load translation matrix
    }
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    float frac[3], temp[3];
                    temp[0] = (ia + 0.5) / float(ANUM);
                    temp[1] = (ib + 0.5) / float(BNUM);
                    temp[2] = (ic + 0.5) / float(CNUM);
                    for( int i=0; i<=2 ; ++i )
                    {
                        frac[i] = 0;
                        for( int j=0; j<=2 ; ++j )
                        {                        //frac = matmul(rotation,temp) + translation;
                            frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                        }
                        frac[i] = frac[i] + tran[i];
                    }
                    adjustFracCoorIntoAsu(frac);
                    int ia0 = floor(frac[0] * ANUM);
                    int ib0 = floor(frac[1] * BNUM);
                    int ic0 = floor(frac[2] * CNUM);
                    densAsuTemp[ia0][ib0][ic0] = densAsu[ia][ib][ic];
                }
            }
        }
    }

    //update densAsu
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsu[ia][ib][ic] = densAsuTemp[ia][ib][ic];
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void translateOrigProtMaskAsuMinDiff()
{
    float rota[3][3], tran[3];

    //find the origin to minimize density-difference
    int iorig0;
    float maskDiff, maskDiffMin;

    for(int iorig = 0; iorig <= numOrigTran - 1 ; ++iorig) {
        for(int i = 0; i <= 2; ++i)
        {
            for(int j = 0; j <= 2; ++j) 
            {
                rota[i][j] = origTran[iorig][i][j]; //load rotation matrix
            }
            tran[i] = origTran[iorig][i][3]; //load translation matrix
        }
        for( int IALOOPASU ) {
            float x = (ia + 0.5) / float(ANUM);
            for( int IBLOOPASU ) {
                float y = (ib + 0.5) / float(BNUM);
                for( int ICLOOPASU ) {
                    float z = (ic + 0.5) / float(CNUM);
                    if( ASUCOND ) {
                        float frac[3], temp[3];
                        temp[0] = (ia + 0.5) / float(ANUM);
                        temp[1] = (ib + 0.5) / float(BNUM);
                        temp[2] = (ic + 0.5) / float(CNUM);
                        for( int i=0; i<=2 ; ++i )
                        {
                            frac[i] = 0;
                            for( int j=0; j<=2 ; ++j )
                            {                        //frac = matmul(rotation,temp) + translation;
                                frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                            }
                            frac[i] = frac[i] + tran[i];
                        }
                        adjustFracCoorIntoAsu(frac);
                        int ia0 = floor(frac[0] * ANUM);
                        int ib0 = floor(frac[1] * BNUM);
                        int ic0 = floor(frac[2] * CNUM);
                        protMaskAsuTemp[ia0][ib0][ic0] = protMaskAsu[ia][ib][ic];
                    }
                }
            }
        }
        maskDiff = 0;
        for( int IALOOPASU ) {
            float x = (ia + 0.5) / float(ANUM);
            for( int IBLOOPASU ) {
                float y = (ib + 0.5) / float(BNUM);
                for( int ICLOOPASU ) {
                    float z = (ic + 0.5) / float(CNUM);
                    if( ASUCOND ) {
                        maskDiff += abs(protMaskAsuAvg[ia][ib][ic] - protMaskAsuTemp[ia][ib][ic]);
                    }
                }
            }
        }
        if(iorig == 0) {
            iorig0 = 0;
            maskDiffMin = maskDiff;
        } else if(maskDiff < maskDiffMin) {
            iorig0 = iorig;
            maskDiffMin = maskDiff;
        }
    }

    cout << "Compare to the best run, relative origin is "<< iorig0 << endl;

    //after find the proper origin, translate origin
    for(int i = 0; i <= 2; ++i)
    {
        for(int j = 0; j <= 2; ++j) 
        {
            rota[i][j] = origTran[iorig0][i][j]; //load rotation matrix
        }
        tran[i] = origTran[iorig0][i][3]; //load translation matrix
    }
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    float frac[3], temp[3];
                    temp[0] = (ia + 0.5) / float(ANUM);
                    temp[1] = (ib + 0.5) / float(BNUM);
                    temp[2] = (ic + 0.5) / float(CNUM);
                    for( int i=0; i<=2 ; ++i )
                    {
                        frac[i] = 0;
                        for( int j=0; j<=2 ; ++j )
                        {                        //frac = matmul(rotation,temp) + translation;
                            frac[i] = frac[i] + rota[i][j] * temp[j]; //rotation
                        }
                        frac[i] = frac[i] + tran[i];
                    }
                    adjustFracCoorIntoAsu(frac);
                    int ia0 = floor(frac[0] * ANUM);
                    int ib0 = floor(frac[1] * BNUM);
                    int ic0 = floor(frac[2] * CNUM);
                    densAsuTemp[ia0][ib0][ic0] = densAsu[ia][ib][ic];
                }
            }
        }
    }

    //update densAsu
    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    densAsu[ia][ib][ic] = densAsuTemp[ia][ib][ic];
                }
            }
        }
    }

    return;
}

/******************************************************************************/

void cutAndOutputSingMoleMask(float atomCoorArray[][3], string fileNameString)
{

    ofstream outputFile( fileNameString.c_str() );
    if (!outputFile.is_open()) {
        cout << "Can't open the file " << fileNameString << endl;
        exit(EXIT_FAILURE);
    }
    outputFile << writeCryst1() << endl;

    //find a qubic region to hold a single molecule
    float orth[3],frac[3];
    for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom) {
        orth[0] = atomCoorArray[iatom][0];
        orth[1] = atomCoorArray[iatom][1];
        orth[2] = atomCoorArray[iatom][2];
        convertOrthToFrac(orth, frac);
        int ia = floor(frac[0] * ANUM);
        int ib = floor(frac[1] * BNUM);
        int ic = floor(frac[2] * CNUM);
        if(iatom == 0){
            anum0 = ia;
            bnum0 = ib;
            cnum0 = ic;
            anum1 = ia;
            bnum1 = ib;
            cnum1 = ic;
            continue;
        }
        if(ia < anum0) anum0 = ia;
        if(ib < bnum0) bnum0 = ib;
        if(ic < cnum0) cnum0 = ic;
        if(ia > anum1) anum1 = ia;
        if(ib > bnum1) bnum1 = ib;
        if(ic > cnum1) cnum1 = ic;
    }
    //settle ncs mask into a little bit loose qubic region
    anum0 = anum0 - 3;
    anum1 = anum1 + 3;
    bnum0 = bnum0 - 3;
    bnum1 = bnum1 + 3;
    cnum0 = cnum0 - 3;
    cnum1 = cnum1 + 3;

    float rota[3][3], tran[3];
    float orth0[3], orth2[3], frac1[3], frac2[3];
    float dist, distMin;
    int igrid = 0;

    for( int IALOOPASU ) {
        float x = (ia + 0.5) / float(ANUM);
        for( int IBLOOPASU ) {
            float y = (ib + 0.5) / float(BNUM);
            for( int ICLOOPASU ) {
                float z = (ic + 0.5) / float(CNUM);
                if( ASUCOND ) {
                    if(protMaskAsu[ia][ib][ic] >= 1){
                        distMin = 1000;
                        frac[0] = (ia + 0.5) / float(ANUM);
                        frac[1] = (ib + 0.5) / float(BNUM);
                        frac[2] = (ic + 0.5) / float(CNUM);
                        for(int isymm = 0; isymm <= numSymmOper - 1; ++isymm) {
                            for(int i = 0; i <= 2; ++i){
                                for(int j = 0; j <= 2; ++j){
                                    rota[i][j] = symmOper[isymm][i][j]; //load rotation matrix
                                } 
                                tran[i] = symmOper[isymm][i][3]; //load rotation matrix
                            }
                            for( int i=0; i<=2 ; ++i ){
                                frac1[i] = 0;
                                for( int j=0; j<=2 ; ++j ){
                                    frac1[i] = frac1[i] + rota[i][j] * frac[j]; //rotation
                                }
                                frac1[i] = frac1[i] + tran[i]; //translation
                            }
                            adjustFracCoorIntoUnitCell(frac1);  
                            //27 adjacent cells, 9 cells on top, 9 cells in the middle, and 9 cells at bottom
                            for(int icell = 0; icell <= numAdjaCell - 1; ++icell){
                                for(int i = 0; i <= 2; ++i){
                                    frac2[i] = frac1[i] + adjaCell[icell][i];
                                }
                                int ia0 = floor(frac2[0] * ANUM);
                                int ib0 = floor(frac2[1] * BNUM);
                                int ic0 = floor(frac2[2] * CNUM);
                                if(ia0 >= anum0 && ia0 <= anum1 && ib0 >= bnum0 && 
                                    ib0 <= bnum1 && ic0 >= cnum0 && ic0 <= cnum1){
                                    convertFracToOrth(frac2, orth2);
                                    for(int iatom = 0; iatom <= numAtomAsu - 1; ++iatom){
                                        orth[0] = atomCoorArray[iatom][0];
                                        orth[1] = atomCoorArray[iatom][1];
                                        orth[2] = atomCoorArray[iatom][2];
                                        dist = distOrth(orth, orth2);
                                        if(dist < distMin) {
                                            distMin = dist;
                                            orth0[0] = orth2[0];
                                            orth0[1] = orth2[1];
                                            orth0[2] = orth2[2];
                                        }
                                    }
                                }
                            }
                        }
                        igrid++;
                        if(igrid >= 9999) igrid = 0; //too many atoms
                        orth0[0] = orth0[0] + random() % 10 * 0.01; //make surface in pymol
                        orth0[1] = orth0[1] + random() % 10 * 0.01;
                        orth0[2] = orth0[2] + random() % 10 * 0.01;
                        float tempFact = 0;
                        float occupancy = 1.0;
                        outputFile << "ATOM  " << setfill(' ') << setw(5) << igrid
                        << ' ' << " C   ARG A   1    " 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth0[0]
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth0[1] 
                        << setfill(' ') << setw(8) << fixed << setprecision(3) << orth0[2] 
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << occupancy
                        << setfill(' ') << setw(6) << fixed << setprecision(2) << tempFact
                        << string(11, ' ') << 'C' << endl;
                    }
                }
            }
        }
    }
    outputFile.close();
    return;
}

/******************************************************************************/
