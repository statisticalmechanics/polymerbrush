/*********************************************************************
 * input simulation parameters from file "input"
 *********************************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void ReadInput(void)
{
 int i;
 FILE *fp;
 
 fp=fopen("input","r"); 
 fscanf(fp,"%d %d %d",&NChainX,&NChainY,&NMonomerPerChain);
 fscanf(fp,"%d",&EnsembleType);
 fscanf(fp,"%d %d",&rInitialType,&vInitialType);
 fscanf(fp,"%d",&PotentialType);
 fscanf(fp,"%lf %lf %lf",&T0,&Tf,&P);
// fscanf(fp,"%lf",&V);
 fscanf(fp,"%lf",&latticeconst);
 fscanf(fp,"%lf %lf",&rc,&rv);
 fscanf(fp,"%d %d",&NumberOfSteps,&NumberOfInitialSteps);
 fscanf(fp,"%d %d",&SampleMultiplier,&MovieMultiplier);
 fscanf(fp,"%d",&drBins);
 fscanf(fp,"%lf",&dt);
 fscanf(fp,"%lf",&randomseed);
 fscanf(fp,"%lf %lf %lf",&AspectRatioX,&AspectRatioY,&AspectRatioZ);
 fscanf(fp,"%d %d %d",&BCType[0],&BCType[1],&BCType[2]);
 fscanf(fp,"%d %d %d",&MinimumImage[0],&MinimumImage[1],&MinimumImage[2]);
 fscanf(fp,"%d %d %d %d %d",&WallSwitch,&FrozenSwitch,&BendSwitch,&TorsionSwitch,&COMSwitch);
 fclose(fp);
 
 //fscanf(fp,"%d %d %d",&NumberOfParticles,&NA,&NB);
 //fscanf(fp,"%lf %lf",&massA,&massB);
 //fscanf(fp,"%lf %lf %lf",&sigmaA,&sigmaB,&sigmaAB);
 //fscanf(fp,"%lf %lf %lf",&epsilonA,&epsilonB,&epsilonAB);
 //fscanf(fp,"%d %d %d %d",&rInitialType,&PackType,&NumberOfLatticeSites,&vInitialType);
 
 atomsize[0] = 1.0;
 NAtomPerMonomer = 8;
 NAtomPerChain = NMonomerPerChain*NAtomPerMonomer;
 NChains = NChainX*NChainY;
 NumberOfParticles = NAtomPerChain*NChains;
 printf("number of atom per monomer = %d\n",NAtomPerMonomer);
 printf("number of atom per chain = %d\n",NAtomPerChain);
 printf("number of chains = %d\n",NChains);
 printf("number of atom = %d\n",NumberOfParticles);


 rminLJ = pow(2.0,1.0/6.0);

 // input parameters are all in reduced units
 // rc and rv are multiplier of sigmaA

 //fA = 1.*NA/NumberOfParticles;
 //fB = 1.*NB/NumberOfParticles;

// Nf = 3*NumberOfParticles-3;
 Nf = 3*NumberOfParticles-(NChains*3);
 printf("degree of freedom Nf : %d \n",Nf);
 

 //InitializeRandomNumberGenerator(time(0l)); //time(0L)
 //InitializeRandomNumberGenerator(randomseed); //time(0L)

 //sigmaAB = (sigmaA+sigmaB)/2.;

 //rc = rc*sigmaA;
 //rv = rv*sigmaA;


 printf("\n");
 if(PotentialType == 0) printf("Lennard-Jones potential\n");
 if(PotentialType == 1) printf("shifted Lennard-Jones potential u(rc) = 0\n");
 if(PotentialType == 2) printf("shifted-force Lennard-Jones potential f(rc) = 0\n");
 printf("\n");

 if(EnsembleType == 0) printf("NVE ensemble\n");
 if(EnsembleType == 1) printf("isokinetic ensemble by velocity rescaling\n");
 if(EnsembleType == 2) printf("NVT ensemble via constraint\n");
 if(EnsembleType == 3) printf("NPH ensemble via constraint\n");
 if(EnsembleType == 4) printf("NPT ensemble via constraint\n");
 if(EnsembleType == 5) printf("NVT ensemble via constraint then switch to NVE\n");
 printf("\n");


 T = T0; // initial temperature
 kB = 1.0;
 beta = 1.0 / T / kB;

 //T = T + dT*JobIndex;
 //rho = rho + dT*JobIndex;

 //L = pow(V,1.0/3.0); 
 LX = pow(V/(AspectRatioX*AspectRatioY*AspectRatioZ),1.0/3.0)*AspectRatioX; 
 LY = pow(V/(AspectRatioX*AspectRatioY*AspectRatioZ),1.0/3.0)*AspectRatioY; 
 LZ = pow(V/(AspectRatioX*AspectRatioY*AspectRatioZ),1.0/3.0)*AspectRatioZ; 

 LX = latticeconst*NChainX;
 LY = latticeconst*NChainY;
 LZ = 2.0*NMonomerPerChain+2.0; 
 L = MIN(MIN(LX,LY),LZ); 
 V = LX*LY*LZ;
 rho = NumberOfParticles/V*CUBIC(atomsize[0]);
 dradial = L/2./drBins; // PBC is used, L/2

 //packingfraction = M_PI/6.*rho*CUBIC(sigmaA);
 //packingfraction = M_PI/6.*rho*(fA*CUBIC(sigmaA)+fB*CUBIC(sigmaB));

 //V = NumberOfParticles / rho;
 //L = pow(V,1.0/Dimension); 

 //dradial = L/2.0/drBins;


 printf("\n");
 printf("Number of MD steps: %d\n",NumberOfSteps);
 printf("Number of equilibrium MD steps: %d\n",NumberOfInitialSteps);
 printf("time increment dt = %lf\n",dt);
 printf("total time = %lf\n",dt*NumberOfSteps);
 printf("Sample multiplier (sampling frequency): %d\n",SampleMultiplier);
 printf("Movie multiplier (draw snapshot frequency): %d\n",MovieMultiplier);
 printf("\n");

 //printf("Dimension: %d\n",Dimension);

/*
 if(RadialDisSwitch == 1)
 {
  printf("Radial Distribution Function g(r) is calculated. \n");
  printf("%d bins of width %lf in histogram of g(r)\n",drBins,dradial);
 }
*/

 printf("\n");
 printf("T0 = %lf\tTf = %lf\n",T0,Tf);
 printf("T = %lf\tbeta = %lf\n",T,beta);
 printf("rho = %lf\n",rho);
 printf("external pressure P = %lf\n",P);
 printf("packing fraction = %lf\n",packingfraction);
 printf("V = %lf\n",V);
 printf("LX LY LZ = %lf %lf %lf\n",LX,LY,LZ);
 printf("rcutoff = %lf\n",rc);
 printf("rv = %lf\n",rv);
 
 printf("\n");
 //printf("ucA = %lf\tucB = %lf\tucAB = %lf\n",ucA,ucB,ucAB);
 printf("\n");


 return;
}
