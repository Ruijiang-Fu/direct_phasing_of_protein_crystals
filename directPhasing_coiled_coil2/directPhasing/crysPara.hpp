#define PDB_CODE_REFE_HIST 	PDB_CODE
#define PDB_CODE "6c4z"

//cell parameters
#define CRYST_A 172.279
#define CRYST_B 172.279
#define CRYST_C 153.758
#define CRYST_ALPHA ( 90.00 / 180.00 * PI )
#define CRYST_BETA  ( 90.00 / 180.00 * PI )
#define CRYST_GAMMA ( 120.00 / 180.00 * PI )

//grid size
#define ANUM 172 //mod(ANUM,2)==0
#define BNUM 172 //mod(BNUM,2)==0
#define CNUM 156 //mod(CNUM,6)==0 for P6122

using namespace std;

const float resoCutoff = 3.3;
const float resoCutoffLow = 15.0;

const float solvContIni = 0.73; //sfcheck 75.3% pdb matthews 81.22%
const float solvContFina = 0.73;
float solvCont = solvContIni;

//weighted observed data
const float sigmWeigObsIni = 1.5;
float sigmWeigObs = sigmWeigObsIni;

//weighted average density
const float sigmWeigAvgIni = 4.0;
const float sigmWeigAvgFina = 3.0;
float sigmWeigAvg = sigmWeigAvgIni;
char updateSpeedWeigAvgDens[] = "normally"; 
//"normally" 100%; "slowly" 50% + 50%; "very slowly" 5% + 95%; "extremely slowly" 1% + 99%

const int numIter = 10000;

const float goodProtMaskThreshold = 0.90; //protMaskMatch > threshold
const float convThreshold = 70.0; //meanPhaseError < threshold

#define NCSANUM 180
#define NCSBNUM 180
#define NCSCNUM 180

//sigma of weighting function
const float sigmWeigNcs = 5.0;

//sigma of weighting function for updating centMassOrth
const float sigmWeigCentMass = 15.0;

const int numNcsOper = 2;
float ncsOper[numNcsOper][3][4] = 
{  	{ { 1.0, 0.0, 0.0, 0.0 }, 
	  { 0.0, 1.0, 0.0, 0.0 }, 
	  { 0.0, 0.0, 1.0, 0.0 } }, 
	{ { -0.2094,  0.8364, -0.5066,  51.1571},
	  { -0.6932, -0.4923, -0.5263, 141.5947},
	  { -0.6896,  0.2410,  0.6829,  42.6693} } };

const int numCentMass = 1;
const float centMassOrthTrue[numCentMass][3] = {22.242, 0.011, -1.513}; //center of molecule
float centMassOrth[numCentMass][3] = {22.242, 0.011, -1.513};

const float pointAOnNcsAxisTrue[3] = {22.242, 0.011, -1.513};
const float vectOnNcsAxisTrue[3] = {0,0,1}; //{-0.016, 0.036, 27.957};
float pointAOnNcsAxis[3] = {22.242, 0.011, -1.513};
float pointBOnNcsAxis[3] = {};
float vectOnNcsAxis[3] = {0,0,1};;

//self-rotation function at section kappa: lateral angle psi; azimuthal angle phi.
const float kappa = 180.0 / 180.0 * 3.1415926;
const float lateAnglPsi = 26.0 / 180.0 * 3.1415926;
const float azimAnglPhi = -175.0 / 180.0 * 3.1415926; //+-5 deg, and +-175 deg
//const float vectOnNcsAxisIni[3] = {-0.43670, -0.03821, 0.89879};

const int numNcsAxisCand = 101;//candNcsAxis + initNcsAxis

const float ncsAxisDeviDistMax = 10.0;	//10 angstroms
const float ncsAxisDeviAnglMax = 10.0 / 180.0 * 3.1415926; //10 degrees

const float rotaNcsAngl = 180.0 / 180.0 * 3.1415926;
const bool tranNcsFlag = false;
const float tranNcsDist = 0;

#include "../header/commPara.hpp"
#include "../header/P6122.hpp"

