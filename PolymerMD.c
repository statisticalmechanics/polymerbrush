/***************************************************
* 3D Molecular Dynamics (MD) simulation of
* polymers
* with continuous potential
* Lennard-Jones (LJ) 	            0
* shifted Lennard-Jones (LJ)        1
* shifted-force Lennard-Jones (LJ)  2
* NVT Gaussian Thermostat
* NPH Gaussian Barostat
* Kai Zhang, Yale University, 2014
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "system.h"
#include "ran_uniform.h"

//int main(void)
int main(int argc, char *argv[])
{
 int step;
 int i,j,k;
 char filename[20];
 FILE *fp;
 FILE *fpmovie;
 FILE *fpvelocity;
 FILE *fpsample;
 FILE *fpr; // position
 FILE *fpv; // velocity
 FILE *fpf; // force

 double Sum_K,Sum_U,Sum_E,Sum_W,Sum_T,Sum_V,Sum_rho,Sum_P,Sum_H;
 double P2backbone,Sum_P2backbone,P2phenyl,Sum_P2phenyl,P2ring,Sum_P2ring;
 double P2ringy,P2ringx;
 double Sum_P2ringy,Sum_P2ringx;
 double Sum_RGyration[1000],Sum_dRend[1000],AvgRGyration,AvgdRend;
 double Sum_PolymerHeight[1000],AvgPolymerHeight;
 double Sum_COSij[500][100],AvgCOSij[100];
 double Count;

 sscanf(argv[1],"%d",&JobIndex);
 
 printf("**************** 3D Molecular Dynamics simulation ********************");
 printf("\n");

 ReadInput();

 if(randomseed == 0.)  randomseed = (double)(JobIndex);
 printf("randomseed = %lf\n",randomseed);
 InitializeRandomNumberGenerator(randomseed); //time(0L)

 Initialization(JobIndex); // particle positions and velocities
 fflush(stdout);
 
 //sprintf(filename,"movie_%d.xyz",JobIndex);
 //fpmovie=fopen(filename,"w");
 // Writemovie(fpmovie);
 //fclose(fpmovie);
 //exit(1);

 Store0(); // initialize r=RX, p=PX, r'=RX1
 VerletCheckIndex = 1; // make list in the beginning
 Force(); // Verlet list is generated here
 Store0(); // initialize r''=RX2,r'''=RX3, PX, p'=f=PX1
 
 Kinetic();
 Pinstant = rho*Tinstant + Virial/V;
 fflush(stdout);
   
 if(EnsembleType == 3 || EnsembleType == 4) // NPH or NPT scale volume
 ScaleNPT();

 printf("initial Pinstant = %lf initial Volume = %lf\n",Pinstant,V);

 Sum_K = 0.;
 Sum_U = 0.;
 Sum_E = 0.;
 Sum_W = 0.;
 Sum_T = 0.;
 Sum_V = 0.;
 Sum_P = 0.;
 Sum_H = 0.;
 Sum_rho = 0.;
 Sum_P2backbone = 0.0;
 Sum_P2phenyl = 0.0;
 Sum_P2ring = 0.0;
 Sum_P2ringx = 0.0;
 Sum_P2ringy= 0.0;
   for(j=0;j<NChains;j++) 
   {
    Sum_RGyration[j] = 0.0;
    Sum_dRend[j] = 0.0;
    Sum_PolymerHeight[j] = 0.0;
     for(k=1;k<=NMonomerPerChain*2-2;k++)
      Sum_COSij[j][k] = 0.0;
   }
 Count = 0.;
 
 for(i=0;i<drBins;i++)
  g[i] = 0.;


 sprintf(filename,"movie_%d.xyz",JobIndex);
 fpmovie=fopen(filename,"w");
 sprintf(filename,"sample_%d.dat",JobIndex);
 fpsample=fopen(filename,"w");

// begin MD loop
 for(step=0;step<NumberOfSteps;step++)
 {
  runtime = step*dt;
  if(step%MovieMultiplier==0) Writemovie(fpmovie);

  if(step%SampleMultiplier==0)
  {
  P2backbone = P2COS(1);
  P2phenyl= P2COS(2);
  P2ring= P2COS(3);
  P2ringx= P2COS(4);
  P2ringy= P2COS(5);
  RadiusGyration();
  BondCorrelation();
  printf("%lf\n",runtime);
  //printf("%lf %lf %lf %lf\n",runtime,P2ring,P2ringx,P2ringy);
  fprintf(fpsample,"t = %lf K = %lf U = %lf E = %lf T = %lf P = %lf V = %lf H = %lf chi = %lf xi = %lf chi+xi = %lf P2b = %lf P2s = %lf Rg0 = %lf dRend0 = %lf z0 = %lf\n", \
  step*dt,Kinstant/NumberOfParticles,Upotential/NumberOfParticles,(Kinstant+Upotential)/NumberOfParticles,Tinstant,Pinstant,V,(Kinstant+Upotential+Pinstant*V)/NumberOfParticles,chi,xi,chipxi,P2backbone,P2phenyl,RGyration[0],dRend[0],PolymerHeight[0]);

 /********************/
  if(step>=NumberOfInitialSteps)
  {
   Sum_K += Kinstant;
   Sum_U += Upotential;
   Sum_E += Kinstant + Upotential;
   Sum_H += Kinstant + Upotential + Pinstant*V;
   Sum_W += Virial;
   Sum_T += Tinstant;
   Sum_V += V;
   Sum_rho += rho;
   Sum_P += Pinstant;
   Sum_P2backbone += P2backbone;
   Sum_P2phenyl+= P2phenyl;
   Sum_P2ring+= P2ring;
   Sum_P2ringx+= P2ringx;
   Sum_P2ringy+= P2ringy;

   for(j=0;j<NChains;j++) 
   {
    Sum_RGyration[j] += RGyration[j];
    Sum_dRend[j] += dRend[j];
    Sum_PolymerHeight[j] += PolymerHeight[j];
     for(k=1;k<=NMonomerPerChain*2-2;k++)
      Sum_COSij[j][k] += COSij[j][k];
   }

   Count += 1.;
   RadialDis(); // g(r)
  }//endif after equilibration
 /********************/
  }//endif every ? steps

  /*********************************/
  if(EnsembleType == 0) // NVE
  {
   VelocityVerlet(1); // r(t+dt) and v' using f(t)
   VerletCheck(); // check if Verlet list needs to be updated
   Force(); // f(t+dt)
   VelocityVerlet(2); // v(t+dt) using f(t+dt)
   Kinetic();
   Pinstant = rho*Tinstant + Virial/V;
  }//endif NVE
  else if(EnsembleType == 1) // isokinetic
  {
   VelocityVerlet(1); // r(t+dt) and v' using f(t)
   VerletCheck(); // check if Verlet list needs to be updated
   Force(); // f(t+dt)
   VelocityVerlet(2); // v(t+dt) using f(t+dt)
   Kinetic();
   //ScaleVelocity();
   if(step%SampleMultiplier==0) ScaleVelocity();
   Pinstant = rho*Tinstant + Virial/V;
  }//endif isokinetic
  else if(EnsembleType == 2)// NVT with constraint method
  {
   PredictNVT(dt);
   VerletCheck();
   Force();
   Kinetic();
   Pinstant = rho*Tinstant + Virial/V;
   xi = PdotF/PdotP;
   CorrectNVT(dt,xi);
   //if(step < NumberOfInitialSteps || step%SampleMultiplier==0) ScaleVelocity(); //avoid long time drifting
   if(step < NumberOfInitialSteps/4 || step%SampleMultiplier==0) ScaleVelocity(); //avoid long time drifting
  }// endif NVT
  else if(EnsembleType == 3)// NPH with constraint method
  {
   PredictNPH(dt);
   VerletCheck();
   Force(); // PdotF and RdotPdotX
   Kinetic(); // PdotP Tinstant Twalltop Twallbottom
   V = LX*LY*LZ;
   rho = NumberOfParticles/V*CUBIC(atomsize[0]);
   Pinstant = rho*Tinstant + Virial/V;
   chi = (2.*PdotF - RdotPdotX)/(9.*HyperVirial+9.*Pinstant*V+2.*PdotP);
   CorrectNPH(dt,chi);
   if(step < NumberOfInitialSteps || step%SampleMultiplier==0) ScaleVolume(); //avoid long time drifting
  }//endif NPH
  else if(EnsembleType == 4) // NPT with constrait method
  {
   PredictNPT(dt);
   VerletCheck(); // check if Verlet list needs to be updated
   Force(); // us in box position[i] to evaluate rij
   Kinetic(); // get PdotP
   V = LX*LY*LZ;
   rho = NumberOfParticles/V*CUBIC(atomsize[0]);
   Pinstant = rho*Tinstant + Virial/V;
   chi =  -RdotPdotX/(9.*HyperVirial+9.*Pinstant*V);
   chipxi = PdotF/PdotP; // chipxi = PdotF/(2.*Kinstant); mass not 1 not work
   CorrectNPT(dt,chi,chipxi);
   if(step < NumberOfInitialSteps || step%SampleMultiplier==0) ScaleNPT(); //avoid long time drifting
  }//endif NPT
  else if(EnsembleType == 5)// NVT with constraint method then Switch to NVE
  {
   PredictNVT(dt);
   VerletCheck();
   Force();
   Kinetic();
   Pinstant = rho*Tinstant + Virial/V;
   xi = PdotF/PdotP;
   CorrectNVT(dt,xi);
   //if(step < NumberOfInitialSteps || step%SampleMultiplier==0) ScaleVelocity(); //avoid long time drifting
   if(step < NumberOfInitialSteps/5 || step%SampleMultiplier==0) ScaleVelocity(); //avoid long time drifting
   //if(step >= NumberOfInitialSteps/5) EnsembleType = 0 ; //Switch to NVE
   if(step >= NumberOfInitialSteps/5) EnsembleType = 1 ; //Switch to NVK
  }// endif NVT -> NVE
  /*********************************/

  fflush(stdout);
 }//end MD loop step
 fclose(fpsample);
 
 Writemovie(fpmovie);
 fclose(fpmovie);

 /***********************************************************/
 sprintf(filename,"position1_%d.dat",JobIndex);
 fpr=fopen(filename,"w");
 sprintf(filename,"velocity1_%d.dat",JobIndex);
 fpv=fopen(filename,"w");
 sprintf(filename,"force1_%d.dat",JobIndex);
 fpf=fopen(filename,"w");
 for(i=0;i<NumberOfParticles;i++)
 {
  fprintf(fpr,"%lf %lf %lf\n",position[i].x,position[i].y,position[i].z);
  fprintf(fpv,"%lf %lf %lf\n",velocity[i].x,velocity[i].y,velocity[i].z);
  fprintf(fpf,"%lf %lf %lf\n",force[i].x,force[i].y,force[i].z);
 }
 fclose(fpr);
 fclose(fpv);
 fclose(fpf);
 /***********************************************************/
  
/*******************/
 Sum_K /= Count;
 Sum_U /= Count;
 Sum_E /= Count;
 Sum_H /= Count;
 Sum_W /= Count;
 Sum_T /= Count;
 Sum_V /= Count;
 Sum_rho /= Count;
 Sum_P /= Count;
 Sum_P2backbone /= Count;
 Sum_P2phenyl /= Count;
 Sum_P2ring /= Count;
 Sum_P2ringx /= Count;
 Sum_P2ringy /= Count;
 for(j=0;j<NChains;j++) 
 {
  Sum_RGyration[j] /= Count;
  Sum_dRend[j] /= Count;
  Sum_PolymerHeight[j] /= Count;
     for(k=1;k<=NMonomerPerChain*2-2;k++)
      Sum_COSij[j][k] /= Count;
 }
 AvgRGyration = 0.0; AvgdRend = 0.0;AvgPolymerHeight = 0.0;
 for(j=0;j<NChains;j++) 
 {
  AvgRGyration += Sum_RGyration[j]; 
  AvgdRend += Sum_dRend[j]; 
  AvgPolymerHeight += Sum_PolymerHeight[j]; 
   for(k=1;k<=NMonomerPerChain*2-2;k++)
    AvgCOSij[k] += Sum_COSij[j][k];
 }
 AvgRGyration /= NChains; AvgdRend /= NChains; AvgPolymerHeight /= NChains;
 for(k=1;k<=NMonomerPerChain*2-2;k++) AvgCOSij[k] /= NChains;
/*******************/

 printf("\n");
 printf("sampling number = %lf\n",Count);

 printf("\n");

 printf("T = %lf\trho = %lf\tP = %lf\t<K>/N = %lf\t<U>/N = %lf\t<E>/N = %lf\t<H>/N = %lf\t<T> = %lf\t<V> = %lf\t<P> = %lf <rho> = %lf\n", \
 T, rho,P,Sum_K/NumberOfParticles,Sum_U/NumberOfParticles,Sum_E/NumberOfParticles,Sum_H/NumberOfParticles,Sum_T,Sum_V,Sum_P,Sum_rho);

 printf("bond order P2 backbone %lf phenyl %lf ring %lf\n",Sum_P2backbone,Sum_P2phenyl,Sum_P2ring);
 printf("P2 ring x = %lf P2 ring y = %lf\n",Sum_P2ringx,Sum_P2ringy);
 
 printf("radius gyration Rg %lf end to end distance dRend %lf polymer height %lf \n",AvgRGyration,AvgdRend,AvgPolymerHeight);
 
 sprintf(filename,"cosij_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(k=1;k<=NMonomerPerChain*2-2;k++) 
  fprintf(fp,"j = %d cosij = %lf\n",k,AvgCOSij[k]);
 fclose(fp);
 
 sprintf(filename,"gr_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  g[i] /= Count;
  g[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial)*Sum_rho;
  fprintf(fp,"r = %lf\tg(r) = %lf\n",(i+0.5)*dradial,g[i]/NumberOfParticles);
 }
 fclose(fp);
 
/********************************************************************/
 // check velocity
  sprintf(filename,"vmaxwell1_%d.dat",JobIndex);
  fp=fopen(filename,"w");
  COMcheck(fp);
  fclose(fp);
/********************************************************************/

 printf("****************************** the end *******************************");


 return 0;
}
