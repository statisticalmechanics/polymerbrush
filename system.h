/******************************************
 * headfile containing most global variables
 ******************************************/
#include <stdio.h>

#define MAX_NumberOfParticles 80000
#define MAX_Nb 80000
#define MAX_NumberOfNeighbors 500 // max number of neighboring particles in the Verlet list
#define MAX_drBins 10000 // for g(r)
//#define MAX_Dimension 3
#define MAX_NumberOfChains 200
#define MAX_NumberOfMonomers 100


#define SQR(x) ((x)*(x))
#define CUBIC(x) ((x)*(x)*(x))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define MIN(x,y) ((x)<(y) ? (x) : (y))

extern int JobIndex;
extern double randomseed;
extern double runtime;
extern int WallSwitch,FrozenSwitch,TorsionSwitch,BendSwitch,COMSwitch;

extern int particlefrozen[MAX_NumberOfParticles]; // 0: not frozen 1: frozen
extern double atomsize[100],atommass[100],atomepsilon[100];
extern double bondstiff[100],bondlength[100];
extern double bondanglestiff[100],bondangle[100];
extern double torsionstiff1[100],torsionstiff2[100],torsionphase[100];
extern int atomtype[100];
extern int NAtomPerMonomer,NMonomerPerChain;
extern int NAtomPerChain,NChains,NChainX,NChainY;// 
extern int NSprings[MAX_NumberOfParticles];
extern int NSpringList[MAX_NumberOfParticles][5];
extern double SpringStiff[MAX_NumberOfParticles][5],SpringLength[MAX_NumberOfParticles][5];
extern int NBondAnglesA[MAX_NumberOfParticles],BondAngleListA[MAX_NumberOfParticles][5][2];
extern int NBondAnglesB[MAX_NumberOfParticles],BondAngleListB[MAX_NumberOfParticles][5][2];
extern double BondAngleStiffA[MAX_NumberOfParticles][5],BondAngleA[MAX_NumberOfParticles][5];
extern double BondAngleStiffB[MAX_NumberOfParticles][5],BondAngleB[MAX_NumberOfParticles][5];
extern int NBackBones,BackBoneList[200][100];
extern double TorsionStiff1[100],TorsionStiff2[100],TorsionPhase[100];
extern double RGyration[1000],dRend[1000];
extern double PolymerHeight[1000];
extern double COSij[500][100];

extern int rInitialType,vInitialType,ParticleInitialType;
extern int PackType; // SC, BCC or FCC
extern int PotentialType; //1: shifted-force L-J
extern int EnsembleType; // 

extern int NumberOfParticles; //N
extern int Nf; // degrees of freedom 3N-3 if velocity center of mass fixed
extern int NumberOfLatticeSites;
//extern int NA,NB;
//extern double fA,fB;// fraction

extern int NumberOfSteps;
extern int NumberOfInitialSteps;
extern int SampleMultiplier,MovieMultiplier;

extern double rc,rminLJ; // cutoff distance of pair potential u(r)=0 for r > rcutoff
//extern double rcA,rcB,rcAB;
//extern double ucA,ucB,ucAB; // u(r=rc) at cutoff
//extern double ducdrA,ducdrB,ducdrAB; // du(r=rc)/dr at cutoff
extern double rv; // radius of Verlet neighbor shell

/**************************************************/
extern double dt; // time increment
/**************************************************/

/**************************************************/
extern double kB; // Boltzmann constant, set to 1
extern double T,T0,Tf; // temperature
extern double Tinstant; // temperature
extern double Kinstant,Upotential,Virial,HyperVirial,PdotF,RdotPdotX,PdotP,Enthalpy;
extern double dT; // temperature increment
extern double beta; // 1/(kB*T)
extern double rho,rhoA,rhoB; // number density rho = N/V
extern double packingfraction; // phi = pi/6 rho in 3D
extern double V;  // Volume
//extern double L,L1,L2,L3; //simulation box length V = L^3
extern double L; // min LXLYLZ
extern double LX,LX1,LX2,LX3; //simulation box length V = LX LY LZ
extern double LY,LY1,LY2,LY3; //simulation box length V = LX LY LZ
extern double LZ,LZ1,LZ2,LZ3; //simulation box length V = LX LY LZ
extern double AspectRatioX,AspectRatioY,AspectRatioZ;
extern double latticeconst; // lattice constant
extern double P; // pressure
extern double Pinstant;
extern int BCType[3]; // boundary condition
extern int MinimumImage[3]; // Minimum Image distance switch
/**************************************************/

extern double sigma[MAX_NumberOfParticles]; // hardcore diameter  sigma_i
extern double epsilon[MAX_NumberOfParticles]; // attraction well depth epsilon_i
extern double mass[MAX_NumberOfParticles]; // particle mass  m_i
extern double identity[MAX_NumberOfParticles]; // particle identity 1:A 2:B
//extern double epsilonA,epsilonB,epsilonAB; // attraction well depth for binary mixture
//extern double sigmaA,sigmaB,sigmaAB;
//extern double massA,massB;
extern double g[MAX_drBins];//,gAA[MAX_drBins],gBB[MAX_drBins],gAB[MAX_drBins];
extern double dradial;
extern int drBins;

typedef struct
{
	double x;
	double y;
	double z;
} VECTOR;

extern VECTOR position[MAX_NumberOfParticles]; // t
extern VECTOR position0[MAX_NumberOfParticles]; // position since last Verlet list update
extern VECTOR position_old[MAX_NumberOfParticles]; // t-dt
extern VECTOR position_new[MAX_NumberOfParticles]; // t+dt
extern VECTOR velocity[MAX_NumberOfParticles];
extern VECTOR velocity_old[MAX_NumberOfParticles];
extern VECTOR velocity_new[MAX_NumberOfParticles];
extern VECTOR force[MAX_NumberOfParticles];
extern VECTOR force_old[MAX_NumberOfParticles];
extern VECTOR force_new[MAX_NumberOfParticles];

extern VECTOR atomposition[MAX_NumberOfParticles]; // atom position in a single monomer
extern int NBackboneBonds,NPhenylBonds,NRings;
extern int BackboneBondIndexA[MAX_Nb],BackboneBondIndexB[MAX_Nb]; // A-B
extern int PhenylBondIndexA[MAX_Nb],PhenylBondIndexB[MAX_Nb]; // A-B
extern int RingIndexA[MAX_Nb],RingIndexB[MAX_Nb],RingIndexC[MAX_Nb]; // A-B-C

/************ NPT Gear-Predictor-Corrector**************************/
extern double RX[MAX_NumberOfParticles],RY[MAX_NumberOfParticles],RZ[MAX_NumberOfParticles]; // position not put particle back to box
//extern double RX1[MAX_NumberOfParticles],RY1[MAX_NumberOfParticles],RZ1[MAX_NumberOfParticles]; //
extern double RX2[MAX_NumberOfParticles],RY2[MAX_NumberOfParticles],RZ2[MAX_NumberOfParticles]; //
extern double RX3[MAX_NumberOfParticles],RY3[MAX_NumberOfParticles],RZ3[MAX_NumberOfParticles]; //

extern double PX[MAX_NumberOfParticles],PY[MAX_NumberOfParticles],PZ[MAX_NumberOfParticles]; //
extern double PX1[MAX_NumberOfParticles],PY1[MAX_NumberOfParticles],PZ1[MAX_NumberOfParticles]; //
extern double PX2[MAX_NumberOfParticles],PY2[MAX_NumberOfParticles],PZ2[MAX_NumberOfParticles]; //
extern double PX3[MAX_NumberOfParticles],PY3[MAX_NumberOfParticles],PZ3[MAX_NumberOfParticles]; //

//extern double FX[MAX_NumberOfParticles],FY[MAX_NumberOfParticles],FZ[MAX_NumberOfParticles]; //
//extern double BOX,BOX1,BOX2,BOX3;

extern double chi,xi,chipxi; // chi and chi+xi
/******************************************************************/

extern int VerletCheckIndex; // 0: no need to update  1: need to update
extern int point[MAX_NumberOfParticles]; // point to the index number of where the verlet list of particle i starts from
extern int nlist; // verlet list index
extern int VerletList[MAX_NumberOfParticles*MAX_NumberOfNeighbors]; // verlet list stores all the neighbors of all particles


void ReadInput(void);
void Initialization(int job);
double BoxMuller(double mm, double ss);

double Distance(int i,int j); // return rij^2
double SigmaIJ(int i,int j);
double EpsilonIJ(int i,int j);
double Potential(int i,int j,double r); // u(rij)
double Potential_dr(int i,int j,double r); // du(rij)/dr
double Potential_dr2(int i,int j,double r); // d2u(rij)/dr2
double SpringPotential(int i,int nb,double r); //
double SpringPotential_dr(int i,int nb,double r); //
double BondAnglePotential(int i,int j,int k,double stiff,double equiangle);
VECTOR BondAngleForceA(int i,int j,int k,double stiff,double equiangle);
VECTOR BondAngleForceB(int i,int j,int k,double stiff,double equiangle); // angle i-j-k
double TorsionPotential(int nc, int nb); // input backbone C index 0,1,...NBackBones-1 on nc the chain
VECTOR TorsionForce(int nc, int nb);
double TorsionPotentialPhenyl(int nm, int na,int ni); // input monomer index nm and atom index na and overall atom index ni
VECTOR TorsionForcePhenyl(int nm,int na,int ni);
void PBC(double *xx,double box);

void VerletCheck(void);
void Force(void);
void SpringForce(void); 
void Kinetic(void);
void RadialDis(void);
void Writemovie(FILE *FilePtr);

/**************************************/
// integration algorithms in different ensembles
void VelocityVerlet(int part); // NVE
/**************************************/

void ScaleVelocity(void); //  need to use after Kinetic()
void ScaleVolume(void); //  need to use after Kinetic()
void ScaleNPT(void); //  need to use after Kinetic()

void PredictNVT(double DT);
void CorrectNVT(double DT,double XI);
void PredictNPH(double DT);
void CorrectNPH(double DT,double CHI);
void PredictNPT(double DT);
void CorrectNPT(double DT,double CHI,double CHIPXI);

void Store0(void); // initialization
void COMcheck(FILE *fp);
void Sample(void);
double P2COS(int type);
void RadiusGyration(void);
void BondCorrelation(void);
