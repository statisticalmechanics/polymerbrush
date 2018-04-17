#include "system.h"

int JobIndex;
double randomseed;
double runtime;
int WallSwitch,FrozenSwitch,TorsionSwitch,BendSwitch,COMSwitch;

int particlefrozen[MAX_NumberOfParticles]; // 0: not frozen 1: frozen
double atomsize[100],atommass[100],atomepsilon[100];
double bondstiff[100],bondlength[100];
double bondanglestiff[100],bondangle[100];
double torsionstiff1[100],torsionstiff2[100],torsionphase[100];
int atomtype[100];
int NAtomPerMonomer,NMonomerPerChain;
int NAtomPerChain,NChains,NChainX,NChainY;// 
int NSprings[MAX_NumberOfParticles];
int NSpringList[MAX_NumberOfParticles][5];
double SpringStiff[MAX_NumberOfParticles][5],SpringLength[MAX_NumberOfParticles][5];
int NBondAnglesA[MAX_NumberOfParticles],BondAngleListA[MAX_NumberOfParticles][5][2];
int NBondAnglesB[MAX_NumberOfParticles],BondAngleListB[MAX_NumberOfParticles][5][2];
double BondAngleStiffA[MAX_NumberOfParticles][5],BondAngleA[MAX_NumberOfParticles][5];
double BondAngleStiffB[MAX_NumberOfParticles][5],BondAngleB[MAX_NumberOfParticles][5];
int NBackBones,BackBoneList[200][100];
double TorsionStiff1[100],TorsionStiff2[100],TorsionPhase[100];
double RGyration[1000],dRend[1000];
double PolymerHeight[1000];
double COSij[500][100];

int rInitialType,vInitialType,ParticleInitialType;
int PackType; // SC, BCC or FCC
int PotentialType; //1: shifted-force L-J
int EnsembleType; // 

int NumberOfParticles; //N
int Nf; // degrees of freedom 3N-3 if velocity center of mass fixed
int NumberOfLatticeSites;
//int NA,NB;
//double fA,fB;// fraction

int NumberOfSteps;
int NumberOfInitialSteps;
int SampleMultiplier,MovieMultiplier;

double rc,rminLJ; // cutoff distance of pair potential u(r)=0 for r > rcutoff
//double rcA,rcB,rcAB;
//double ucA,ucB,ucAB; // u(r=rc) at cutoff
//double ducdrA,ducdrB,ducdrAB; // du(r=rc)/dr at cutoff
double rv; // radius of Verlet neighbor shell

/**************************************************/
double dt; // time increment
/**************************************************/

/**************************************************/
double kB; // Boltzmann constant, set to 1
double T,T0,Tf; // temperature
double Tinstant; // temperature
double Kinstant,Upotential,Virial,HyperVirial,PdotF,RdotPdotX,PdotP,Enthalpy;
double dT; // temperature increment
double beta; // 1/(kB*T)
double rho,rhoA,rhoB; // number density rho = N/V
double packingfraction; // phi = pi/6 rho in 3D
double V;  // Volume
//extern double L,L1,L2,L3; //simulation box length V = L^3
double L; // min LXLYLZ
double LX,LX1,LX2,LX3; //simulation box length V = LX LY LZ
double LY,LY1,LY2,LY3; //simulation box length V = LX LY LZ
double LZ,LZ1,LZ2,LZ3; //simulation box length V = LX LY LZ
double AspectRatioX,AspectRatioY,AspectRatioZ;
double latticeconst; // lattice constant
double P; // pressure
double Pinstant;
int BCType[3]; // boundary condition
int MinimumImage[3]; // Minimum Image distance switch
/**************************************************/

double sigma[MAX_NumberOfParticles]; // hardcore diameter  sigma_i
double epsilon[MAX_NumberOfParticles]; // attraction well depth epsilon_i
double mass[MAX_NumberOfParticles]; // particle mass  m_i
double identity[MAX_NumberOfParticles]; // particle identity 1:A 2:B
//double epsilonA,epsilonB,epsilonAB; // attraction well depth for binary mixture
//double sigmaA,sigmaB,sigmaAB;
//double massA,massB;
double g[MAX_drBins];//,gAA[MAX_drBins],gBB[MAX_drBins],gAB[MAX_drBins];
double dradial;
int drBins;

VECTOR position[MAX_NumberOfParticles]; // t
VECTOR position0[MAX_NumberOfParticles]; // position since last Verlet list update
VECTOR position_old[MAX_NumberOfParticles]; // t-dt
VECTOR position_new[MAX_NumberOfParticles]; // t+dt
VECTOR velocity[MAX_NumberOfParticles];
VECTOR velocity_old[MAX_NumberOfParticles];
VECTOR velocity_new[MAX_NumberOfParticles];
VECTOR force[MAX_NumberOfParticles];
VECTOR force_old[MAX_NumberOfParticles];
VECTOR force_new[MAX_NumberOfParticles];

VECTOR atomposition[MAX_NumberOfParticles]; // atom position in a single monomer
int NBackboneBonds,NPhenylBonds,NRings;
int BackboneBondIndexA[MAX_Nb],BackboneBondIndexB[MAX_Nb]; // A-B
int PhenylBondIndexA[MAX_Nb],PhenylBondIndexB[MAX_Nb]; // A-B
int RingIndexA[MAX_Nb],RingIndexB[MAX_Nb],RingIndexC[MAX_Nb]; // A-B

/************ NPT Gear-Predictor-Corrector**************************/
double RX[MAX_NumberOfParticles],RY[MAX_NumberOfParticles],RZ[MAX_NumberOfParticles]; // position not put particle back to box
//extern double RX1[MAX_NumberOfParticles],RY1[MAX_NumberOfParticles],RZ1[MAX_NumberOfParticles]; //
double RX2[MAX_NumberOfParticles],RY2[MAX_NumberOfParticles],RZ2[MAX_NumberOfParticles]; //
double RX3[MAX_NumberOfParticles],RY3[MAX_NumberOfParticles],RZ3[MAX_NumberOfParticles]; //

double PX[MAX_NumberOfParticles],PY[MAX_NumberOfParticles],PZ[MAX_NumberOfParticles]; //
double PX1[MAX_NumberOfParticles],PY1[MAX_NumberOfParticles],PZ1[MAX_NumberOfParticles]; //
double PX2[MAX_NumberOfParticles],PY2[MAX_NumberOfParticles],PZ2[MAX_NumberOfParticles]; //
double PX3[MAX_NumberOfParticles],PY3[MAX_NumberOfParticles],PZ3[MAX_NumberOfParticles]; //

//extern double FX[MAX_NumberOfParticles],FY[MAX_NumberOfParticles],FZ[MAX_NumberOfParticles]; //
//extern double BOX,BOX1,BOX2,BOX3;

double chi,xi,chipxi; // chi and chi+xi
/******************************************************************/

int VerletCheckIndex; // 0: no need to update  1: need to update
int point[MAX_NumberOfParticles]; // point to the index number of where the verlet list of particle i starts from
int nlist; // verlet list index
int VerletList[MAX_NumberOfParticles*MAX_NumberOfNeighbors]; // verlet list stores all the neighbors of all particles
