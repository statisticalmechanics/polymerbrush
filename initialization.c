/***************************************************
* initialize particle parameters
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void Initialization(int job)
{
 char AtomID;
 double bondlength0;
 char filename[20];
 int i,j,k,m;
 int nx,ny,nz;
 double scale;
 double particlepercell;
 //double Kinstant,scale,Tinstant;
 double Mcom;
 VECTOR Vcom;

 int indexshift;
 VECTOR rjshift,rishift;
 double monomerheight,dlx,dly;

 FILE *fp;

 bondlength0 = 1.0;
 monomerheight = 2.0*bondlength0;
 dlx = LX/NChainX; dly = LY/NChainY;
 //NAtomPerMonomer = 8;
 atomposition[0].x = 0.0; atomposition[0].y = 0.0; atomposition[0].z = 0.0;
 atomposition[1].x = 0.0; atomposition[1].y = 0.0; atomposition[1].z = bondlength0;
 atomposition[2].x = bondlength0; atomposition[2].y = 0.0; atomposition[2].z = bondlength0;
 atomposition[3].x = 1.5*bondlength0; atomposition[3].y = -sqrt(3.0)/2.0*bondlength0; atomposition[3].z = bondlength0;
 atomposition[4].x = 2.5*bondlength0; atomposition[4].y =-sqrt(3.0)/2.0*bondlength0 ; atomposition[4].z = bondlength0;
 atomposition[5].x = 3.0*bondlength0; atomposition[5].y = 0.0; atomposition[5].z = bondlength0;
 atomposition[6].x = 2.5*bondlength0; atomposition[6].y = sqrt(3.0)/2.0*bondlength0; atomposition[6].z = bondlength0;
 atomposition[7].x = 1.5*bondlength0; atomposition[7].y = sqrt(3.0)/2.0*bondlength0; atomposition[7].z = bondlength0;

 //van der Waals ~ 1 kcal/mol, considering coordination number, epsilon ~ 0.1 kcal /mol sigma ~ 1.8 A
 atomsize[0] = 1.0; atomtype[0] = 1; atommass[0] = 1.0; atomepsilon[0] = 1.0; // C
 atomsize[1] = 1.0; atomtype[1] = 1; atommass[1] = 1.0; atomepsilon[1] = 1.0;
 atomsize[2] = 1.0; atomtype[2] = 2; atommass[2] = 1.0; atomepsilon[2] = 1.0;
 atomsize[3] = 1.0; atomtype[3] = 2; atommass[3] = 1.0; atomepsilon[3] = 1.0;
 atomsize[4] = 1.0; atomtype[4] = 2; atommass[4] = 1.0; atomepsilon[4] = 1.0;
 atomsize[5] = 1.0; atomtype[5] = 2; atommass[5] = 1.0; atomepsilon[5] = 1.0;
 atomsize[6] = 1.0; atomtype[6] = 2; atommass[6] = 1.0; atomepsilon[6] = 1.0;
 atomsize[7] = 1.0; atomtype[7] = 2; atommass[7] = 1.0; atomepsilon[7] = 1.0;

 // bond stretch v = 1/2 k (l-l0)^2, l0~sigma  C-C k ~ 300 kcal /mol / A^2 = 300*1.8^2 = 972 kcal/mol /sigma^2 
 bondlength[0] = 1.0;bondstiff[0] = 9000;
 bondlength[1] = 1.0;bondstiff[1] = 9000;
 bondlength[2] = 1.0;bondstiff[2] = 18000;
 bondlength[3] = 1.0;bondstiff[3] = 18000;
 bondlength[4] = 1.0;bondstiff[4] = 18000;
 bondlength[5] = 1.0;bondstiff[5] = 18000;
 bondlength[6] = 1.0;bondstiff[6] = 18000;
 bondlength[7] = 1.0;bondstiff[7] = 18000;
 bondlength[8] = 1.0;bondstiff[8] = 9000;
 
 // 12 bond angles v = 1/2 k (theta-theta0)^2, k ~ 0.01 kcal/mol / deg^2 (1 rad  = 57.3 deg) = 32.8 kcal /mol /rad^2
 bondangle[0] = 109.5/180.0*M_PI;bondanglestiff[0] = 400.0;
 bondangle[1] = 2.0*M_PI/3.0;bondanglestiff[1] = 18000.0; //600.0;
 bondangle[2] = 2.0*M_PI/3.0;bondanglestiff[2] = 18000.0; //600.0;
 bondangle[3] = 2.0*M_PI/3.0;bondanglestiff[3] = 18000.0; //600.0;
 bondangle[4] = 2.0*M_PI/3.0;bondanglestiff[4] = 18000.0; //600.0;
 bondangle[5] = 2.0*M_PI/3.0;bondanglestiff[5] = 18000.0; //600.0;
 bondangle[6] = 2.0*M_PI/3.0;bondanglestiff[6] = 18000.0; //600.0;
 bondangle[7] = 2.0*M_PI/3.0;bondanglestiff[7] = 18000.0; //600.0;
 bondangle[8] = 2.0*M_PI/3.0;bondanglestiff[8] = 18000.0; //600.0;
 bondangle[9] = 109.5/180.0*M_PI;bondanglestiff[9] =  400.0;
 bondangle[10] = 109.5/180.0*M_PI;bondanglestiff[10] = 400.0;
 bondangle[11] = 109.5/180.0*M_PI;bondanglestiff[11] = 400.0;
 
 // torsion dihedral, v = c1 cos(phi-phase)^3 - c2 cos(phi-phase) + c1-c2, c ~ 10 kcal /mol
 torsionstiff1[0] = 10.0; torsionstiff2[0] = 6.0; torsionphase[0] = M_PI; //0.0; 
 //torsionstiff1[0] = 0.0; torsionstiff2[0] = 0.0; torsionphase[0] = M_PI; //0.0; 
 torsionstiff1[1] = 10.0; torsionstiff2[1] = 6.0; torsionphase[1] = M_PI; // M_PI; 
 //torsionstiff1[1] = 0.0; torsionstiff2[1] = 0.0; torsionphase[1] = M_PI; // M_PI; 
/*
 fp = fopen("monomerchem","r");
 for(i=0;i<NAtomPerMonomer;i++)
 {
 fscanf(fp,"%s\n",&AtomID);
 if(AtomID == 'C') {atomsize[i] = 1.0; atomtype[i]=1; }
 }
 fclose(fp);
*/

 NBackboneBonds=0;
 NPhenylBonds=0;
 NRings=0;
 for(i=0;i<NChains;i++)
 {
  rishift.x = i%NChainX*dlx ;rishift.y = (i-i%NChainX)/NChainX*dly ; rishift.z = 0.0; 
  NBackBones = 0;
  for(j=0;j<NMonomerPerChain;j++)
  {
   rjshift.x = 0.0;rjshift.y = 0.0; rjshift.z = monomerheight*j; 
   indexshift = j*NAtomPerMonomer+i*NMonomerPerChain*NAtomPerMonomer;
   for(k=0;k<NAtomPerMonomer;k++)
   {
    m = k+indexshift;
    
    position[m].x = atomposition[k].x+rjshift.x+rishift.x;
    position[m].y = atomposition[k].y+rjshift.y+rishift.y;
    position[m].z = atomposition[k].z+rjshift.z+rishift.z;
    identity[m] = atomtype[k];
    sigma[m] = atomsize[k];
    mass[m] = atommass[k];
    epsilon[m] = atomepsilon[k];

    if(k==0)
    {
    BackboneBondIndexA[NBackboneBonds] = m;
    BackboneBondIndexB[NBackboneBonds] = m+1;
    NBackboneBonds++;
    }
    
    if(k==1)
    {
    PhenylBondIndexA[NPhenylBonds] = m;
    PhenylBondIndexB[NPhenylBonds] = m+1;
    NPhenylBonds++;
    
    RingIndexA[NRings] = m+1;
    RingIndexB[NRings] = m+3;
    RingIndexC[NRings] = m+5;
    NRings++;

    if(j!=NMonomerPerChain-1)
    {
    BackboneBondIndexA[NBackboneBonds] = m;
    BackboneBondIndexB[NBackboneBonds] = m+7;
    NBackboneBonds++;
    }
    }

   /*******************set torsion*********************/
    if(k==0)
    {
     BackBoneList[i][NBackBones] = m; 
     TorsionStiff1[NBackBones] = torsionstiff1[0];
     TorsionStiff2[NBackBones] = torsionstiff2[0];
     TorsionPhase[NBackBones] = torsionphase[0];
     NBackBones++;
    }
    if(k==1)
    {
     BackBoneList[i][NBackBones] = m; 
     TorsionStiff1[NBackBones] = torsionstiff1[1];
     TorsionStiff2[NBackBones] = torsionstiff2[1];
     TorsionPhase[NBackBones] = torsionphase[1];
     NBackBones++;
    }
   /*******************end set torsion*********************/

   /*******************set bond length*********************/
    if(k==0){
    if(j!=0)
    {
    NSprings[m]=2;
    NSpringList[m][0]=m-7;
    NSpringList[m][1]=m+1;
    SpringLength[m][0]=bondlength[8]; SpringStiff[m][0]=bondstiff[8];
    SpringLength[m][1]=bondlength[0]; SpringStiff[m][1]=bondstiff[0];
    }
    else
    {
    NSprings[m]=1;
    NSpringList[m][0]=m+1;
    SpringLength[m][0]=bondlength[0]; SpringStiff[m][0]=bondstiff[0];
    }
    }
    
    if(k==1){
    if(j!=NMonomerPerChain-1)
    {
    NSprings[m]=3;
    NSpringList[m][0]=m-1;
    NSpringList[m][1]=m+1;
    NSpringList[m][2]=m+7;
    SpringLength[m][0]=bondlength[0]; SpringStiff[m][0]=bondstiff[0];
    SpringLength[m][1]=bondlength[1]; SpringStiff[m][1]=bondstiff[1];
    SpringLength[m][2]=bondlength[8]; SpringStiff[m][2]=bondstiff[8];
    }
    else
    {
    NSprings[m]=2;
    NSpringList[m][0]=m-1;
    NSpringList[m][1]=m+1;
    SpringLength[m][0]=bondlength[0]; SpringStiff[m][0]=bondstiff[0];
    SpringLength[m][1]=bondlength[1]; SpringStiff[m][1]=bondstiff[1];
    }
    }
    
    if(k==2){
    NSprings[m]=3;
    NSpringList[m][0]=m-1;
    NSpringList[m][1]=m+1;
    NSpringList[m][2]=m+5;
    SpringLength[m][0]=bondlength[1]; SpringStiff[m][0]=bondstiff[1];
    SpringLength[m][1]=bondlength[2]; SpringStiff[m][1]=bondstiff[2];
    SpringLength[m][2]=bondlength[7]; SpringStiff[m][2]=bondstiff[7];
    }
    
    if(k==3 || k==4 || k==5 || k==6){
    NSprings[m]=2;
    NSpringList[m][0]=m-1;
    NSpringList[m][1]=m+1;
    SpringLength[m][0]=bondlength[2]; SpringStiff[m][0]=bondstiff[2];
    SpringLength[m][1]=bondlength[3]; SpringStiff[m][1]=bondstiff[3];
    }
    
    if(k==7){
    NSprings[m]=2;
    NSpringList[m][0]=m-1;
    NSpringList[m][1]=m-5;
    SpringLength[m][0]=bondlength[2]; SpringStiff[m][0]=bondstiff[2];
    SpringLength[m][1]=bondlength[3]; SpringStiff[m][1]=bondstiff[3];
    }
   /*******************end set bond length*********************/
   
   /*******************set bond angle A*********************/
    if(k==0){
    if(j!=0)
    {
    NBondAnglesA[m]=1;
    BondAngleListA[m][0][0]=m-7;BondAngleListA[m][0][1]=m+1; //order not matter
    BondAngleStiffA[m][0]=bondanglestiff[11];BondAngleA[m][0]=bondangle[11];
    }
    else
    NBondAnglesA[m]=0;
    }
    
    if(k==1){
    if(j!=NMonomerPerChain-1)
    {
    NBondAnglesA[m]=3;
    BondAngleListA[m][0][0]=m-1;BondAngleListA[m][0][1]=m+1;
    BondAngleStiffA[m][0]=bondanglestiff[0];BondAngleA[m][0]=bondangle[0];
    BondAngleListA[m][1][0]=m-1;BondAngleListA[m][1][1]=m+7;
    BondAngleStiffA[m][1]=bondanglestiff[10];BondAngleA[m][1]=bondangle[10];
    BondAngleListA[m][2][0]=m+7;BondAngleListA[m][2][1]=m+1;
    BondAngleStiffA[m][2]=bondanglestiff[9];BondAngleA[m][2]=bondangle[9];
    }
    else
    {
    NBondAnglesA[m]=1;
    BondAngleListA[m][0][0]=m-1;BondAngleListA[m][0][1]=m+1;
    BondAngleStiffA[m][0]=bondanglestiff[0];BondAngleA[m][0]=bondangle[0];
    }
    }
   
    if(k==2){
    NBondAnglesA[m]=3;
    BondAngleListA[m][0][0]=m-1;BondAngleListA[m][0][1]=m+1;
    BondAngleStiffA[m][0]=bondanglestiff[1];BondAngleA[m][0]=bondangle[1];
    BondAngleListA[m][1][0]=m-1;BondAngleListA[m][1][1]=m+5;
    BondAngleStiffA[m][1]=bondanglestiff[8];BondAngleA[m][1]=bondangle[8];
    BondAngleListA[m][2][0]=m+5;BondAngleListA[m][2][1]=m+1;
    BondAngleStiffA[m][2]=bondanglestiff[7];BondAngleA[m][2]=bondangle[7];
    }
    
    if(k==3 || k==4 || k==5 || k==6){
    NBondAnglesA[m]=1;
    BondAngleListA[m][0][0]=m-1;BondAngleListA[m][0][1]=m+1;
    BondAngleStiffA[m][0]=bondanglestiff[2];BondAngleA[m][0]=bondangle[2]; // =3 4 5 
    }
   
    if(k==7){
    NBondAnglesA[m]=1;
    BondAngleListA[m][0][0]=m-1;BondAngleListA[m][0][1]=m-5;
    BondAngleStiffA[m][0]=bondanglestiff[6];BondAngleA[m][0]=bondangle[6]; // =3 4 5 
    }
   /*******************end set bond angle A*********************/
   
   /*******************set bond angle B*********************/
    if(k==0){
    if(j==0)
    {
    NBondAnglesB[m]=2;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m+2; // order matters
    BondAngleListB[m][1][0]=m+1;BondAngleListB[m][1][1]=m+8; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[0];BondAngleB[m][0]=bondangle[0]; 
    BondAngleStiffB[m][1]=bondanglestiff[10];BondAngleB[m][1]=bondangle[10]; 
    }
    else if(j==NMonomerPerChain-1)
    {
    NBondAnglesB[m]=3;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m+2; // order matters
    BondAngleListB[m][1][0]=m-7;BondAngleListB[m][1][1]=m-6; // order matters
    BondAngleListB[m][2][0]=m-7;BondAngleListB[m][2][1]=m-8; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[0];BondAngleB[m][0]=bondangle[0]; 
    BondAngleStiffB[m][1]=bondanglestiff[9];BondAngleB[m][1]=bondangle[9]; 
    BondAngleStiffB[m][2]=bondanglestiff[10];BondAngleB[m][2]=bondangle[10]; 
    } 
    else
    {
    NBondAnglesB[m]=4;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m+2; // order matters
    BondAngleListB[m][1][0]=m+1;BondAngleListB[m][1][1]=m+8; // order matters
    BondAngleListB[m][2][0]=m-7;BondAngleListB[m][2][1]=m-6; // order matters
    BondAngleListB[m][3][0]=m-7;BondAngleListB[m][3][1]=m-8; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[0];BondAngleB[m][0]=bondangle[0]; 
    BondAngleStiffB[m][1]=bondanglestiff[10];BondAngleB[m][1]=bondangle[10]; 
    BondAngleStiffB[m][2]=bondanglestiff[9];BondAngleB[m][2]=bondangle[9]; 
    BondAngleStiffB[m][3]=bondanglestiff[10];BondAngleB[m][3]=bondangle[10]; 
    }
    }
  
    if(k==1){
    if(j==0)
    {
    NBondAnglesB[m]=3;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m+2; // order matters
    BondAngleListB[m][1][0]=m+1;BondAngleListB[m][1][1]=m+6; // order matters
    BondAngleListB[m][2][0]=m+7;BondAngleListB[m][2][1]=m+8; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[1];BondAngleB[m][0]=bondangle[1]; 
    BondAngleStiffB[m][1]=bondanglestiff[8];BondAngleB[m][1]=bondangle[8]; 
    BondAngleStiffB[m][2]=bondanglestiff[11];BondAngleB[m][2]=bondangle[11]; 
    }
    else if(j==NMonomerPerChain-1)
    {
    NBondAnglesB[m]=3;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m+2; // order matters
    BondAngleListB[m][1][0]=m+1;BondAngleListB[m][1][1]=m+6; // order matters
    BondAngleListB[m][2][0]=m-1;BondAngleListB[m][2][1]=m-8; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[1];BondAngleB[m][0]=bondangle[1]; 
    BondAngleStiffB[m][1]=bondanglestiff[8];BondAngleB[m][1]=bondangle[8]; 
    BondAngleStiffB[m][2]=bondanglestiff[11];BondAngleB[m][2]=bondangle[11]; 
    }
    else
    {
    NBondAnglesB[m]=4;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m+2; // order matters
    BondAngleListB[m][1][0]=m+1;BondAngleListB[m][1][1]=m+6; // order matters
    BondAngleListB[m][2][0]=m-1;BondAngleListB[m][2][1]=m-8; // order matters
    BondAngleListB[m][3][0]=m+7;BondAngleListB[m][3][1]=m+8; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[1];BondAngleB[m][0]=bondangle[1]; 
    BondAngleStiffB[m][1]=bondanglestiff[8];BondAngleB[m][1]=bondangle[8]; 
    BondAngleStiffB[m][2]=bondanglestiff[11];BondAngleB[m][2]=bondangle[11]; 
    BondAngleStiffB[m][3]=bondanglestiff[11];BondAngleB[m][3]=bondangle[11]; 
    }
    }

    if(k==2){
    if(j!=NMonomerPerChain-1)
    {
    NBondAnglesB[m]=4;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m+2; // order matters
    BondAngleListB[m][1][0]=m+5;BondAngleListB[m][1][1]=m+4; // order matters
    BondAngleListB[m][2][0]=m-1;BondAngleListB[m][2][1]=m-2; // order matters
    BondAngleListB[m][3][0]=m-1;BondAngleListB[m][3][1]=m+6; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[2];BondAngleB[m][0]=bondangle[2]; 
    BondAngleStiffB[m][1]=bondanglestiff[6];BondAngleB[m][1]=bondangle[6]; 
    BondAngleStiffB[m][2]=bondanglestiff[0];BondAngleB[m][2]=bondangle[0]; 
    BondAngleStiffB[m][3]=bondanglestiff[9];BondAngleB[m][3]=bondangle[9]; 
    }
    else
    {
    NBondAnglesB[m]=3;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m+2; // order matters
    BondAngleListB[m][1][0]=m+5;BondAngleListB[m][1][1]=m+4; // order matters
    BondAngleListB[m][2][0]=m-1;BondAngleListB[m][2][1]=m-2; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[2];BondAngleB[m][0]=bondangle[2]; 
    BondAngleStiffB[m][1]=bondanglestiff[6];BondAngleB[m][1]=bondangle[6]; 
    BondAngleStiffB[m][2]=bondanglestiff[0];BondAngleB[m][2]=bondangle[0]; 
    }
    }
    
    if(k==3){
    NBondAnglesB[m]=3;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m+2; // order matters
    BondAngleListB[m][1][0]=m-1;BondAngleListB[m][1][1]=m-2; // order matters
    BondAngleListB[m][2][0]=m-1;BondAngleListB[m][2][1]=m+4; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[3];BondAngleB[m][0]=bondangle[3]; 
    BondAngleStiffB[m][1]=bondanglestiff[1];BondAngleB[m][1]=bondangle[1]; 
    BondAngleStiffB[m][2]=bondanglestiff[7];BondAngleB[m][2]=bondangle[7]; 
    }
    
    if(k==4 || k==5){
    NBondAnglesB[m]=2;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m+2; // order matters
    BondAngleListB[m][1][0]=m-1;BondAngleListB[m][1][1]=m-2; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[4];BondAngleB[m][0]=bondangle[4]; // 5
    BondAngleStiffB[m][1]=bondanglestiff[2];BondAngleB[m][1]=bondangle[2]; // 3
    }

    if(k==6){
    NBondAnglesB[m]=2;
    BondAngleListB[m][0][0]=m+1;BondAngleListB[m][0][1]=m-4; // order matters
    BondAngleListB[m][1][0]=m-1;BondAngleListB[m][1][1]=m-2; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[6];BondAngleB[m][0]=bondangle[6]; 
    BondAngleStiffB[m][1]=bondanglestiff[4];BondAngleB[m][1]=bondangle[4]; 
    }
    
    if(k==7){
    NBondAnglesB[m]=3;
    BondAngleListB[m][0][0]=m-5;BondAngleListB[m][0][1]=m-4; // order matters
    BondAngleListB[m][1][0]=m-5;BondAngleListB[m][1][1]=m-6; // order matters
    BondAngleListB[m][2][0]=m-1;BondAngleListB[m][2][1]=m-2; // order matters
    BondAngleStiffB[m][0]=bondanglestiff[7];BondAngleB[m][0]=bondangle[7]; 
    BondAngleStiffB[m][1]=bondanglestiff[8];BondAngleB[m][1]=bondangle[8]; 
    BondAngleStiffB[m][2]=bondanglestiff[5];BondAngleB[m][2]=bondangle[5]; 
    }

   /*******************end set bond angle B*********************/
   } // end loop k atoms
  }//end loop j monomers
 }// end loop i polymers

   printf("%d backbond carbons per chain\n",NBackBones);
   
   printf("%d P2 backbone bonds (should be %d)\n",NBackboneBonds,(NBackBones-1)*NChains);
   printf("%d P2 phenyl bonds (should be %d)\n",NPhenylBonds,(NMonomerPerChain)*NChains);
   printf("%d P2 rings (should be %d)\n",NRings,(NMonomerPerChain)*NChains);
/*
 for(i=0;i<NChains;i++)
 {
   printf("%d: ",i);
  for(j=0;j<NBackBones;j++)
  {
   //printf("%d ",BackBoneList[i][j]);
   printf("%lf ",TorsionPhase[j]);
  }
   printf("\n");
 }
*/

/*
for(i=0;i<NumberOfParticles;i++)
{
  printf("i = %d",i);
 for(k=0;k<NSprings[i];k++)
   {
    j = NSpringList[i][k];
    if(Distance(i,j) > SQR(SpringLength[i][k])) printf("%lf > %lf ",sqrt(Distance(i,j)),SpringLength[i][k]);
   }
    printf("\n");
}
*/

for(i=0;i<NumberOfParticles;i++)
{
 particlefrozen[i] = 0.;
 if(FrozenSwitch == 1) if(i%NAtomPerChain==0) particlefrozen[i]=1;
}



/****************  position **********/
 if(rInitialType == 0) // read from file
 {
  sprintf(filename,"position_%d",job);
  fp = fopen(filename,"r");
  for(i=0;i<NumberOfParticles;i++)
  fscanf(fp,"%lf\t%lf\t%lf\n",&position[i].x,&position[i].y,&position[i].z);
  fclose(fp);
 }

  for(i=0;i<NumberOfParticles;i++) 
  {
   position_old[i].x = position[i].x - velocity[i].x*dt;
   position_old[i].y = position[i].y - velocity[i].y*dt;
   position_old[i].z = position[i].z - velocity[i].z*dt;
    
    //PBC
   if(BCType[0]==0) PBC(&(position_old[i].x),LX);
   if(BCType[1]==0) PBC(&(position_old[i].y),LY);
   if(BCType[2]==0) PBC(&(position_old[i].z),LZ);
  }
/******************************* end position***************************************/

  
/********************************** velocity ********************************************/

if(vInitialType == 0) // read from file
{
 // fp = fopen("velocity","r");
  sprintf(filename,"velocity_%d",job);
  fp = fopen(filename,"r");
  for(i=0;i<NumberOfParticles;i++)
  fscanf(fp,"%lf\t%lf\t%lf\n",&velocity[i].x,&velocity[i].y,&velocity[i].z);
  fclose(fp);
}
else if(vInitialType == 1) // random velocity
{
  /********************random velocity**************************************/
  Mcom = 0.;
  Vcom.x = 0.;
  Vcom.y = 0.;
  Vcom.z = 0.;

  for(i=0;i<NumberOfParticles;i++)
  {
   velocity[i].x = BoxMuller(0.,1.);
   velocity[i].y = BoxMuller(0.,1.);
   velocity[i].z = BoxMuller(0.,1.);

   Vcom.x += mass[i]*velocity[i].x;
   Vcom.y += mass[i]*velocity[i].y;
   Vcom.z += mass[i]*velocity[i].z;

   Mcom += mass[i];
  }

   Vcom.x /= Mcom;
   Vcom.y /= Mcom;
   Vcom.z /= Mcom;

  Kinstant = 0.;
  for(i=0;i<NumberOfParticles;i++) 
  {
   velocity[i].x -= Vcom.x;
   velocity[i].y -= Vcom.y;
   velocity[i].z -= Vcom.z;

   Kinstant += mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));
  }
 
   Kinstant *= 0.5; // instantenous kinetic energy

   // 0.5*kT*Nf = K = 0.5* sum_mv^2, Nf = 3N-3
   Tinstant = 2.0*Kinstant/Nf/kB;

   scale = sqrt(T/Tinstant);

   for(i=0;i<NumberOfParticles;i++) 
   {
    velocity[i].x *= scale;
    velocity[i].y *= scale;
    velocity[i].z *= scale;
   }
}//end if maxwell velocity
 
// check velocity
  sprintf(filename,"vmaxwell0_%d.dat",JobIndex);
  fp=fopen(filename,"w");
  COMcheck(fp);
  fclose(fp);
 /****************** end velocity ***************************************/
 
 fp=fopen("ur.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tu(r) = %lf\n",i*0.003,Potential(1,1,i*0.003));
 fclose(fp);
 
 fp=fopen("fr.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tf(r) = %lf\n",i*0.003,-Potential_dr(1,1,i*0.003));
 fclose(fp);
 
 fp=fopen("d2udr2.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\td2u(r) = %lf\n",i*0.003,Potential_dr2(1,1,i*0.003));
 fclose(fp);


 fp=fopen("uspring.dat","w"); 
 for(i=0;i<1000;i++)
  fprintf(fp,"r = %lf\tu(r) = %lf\n",i*0.001,SpringPotential(0,0,i*0.001*bondlength[0]));
 fclose(fp);
 
 fp=fopen("fspring.dat","w"); 
 for(i=0;i<1000;i++)
  fprintf(fp,"r = %lf\tu(r) = %lf\n",i*0.001,SpringPotential_dr(0,0,i*0.001*bondlength[0]));
 fclose(fp);

 RX[NumberOfParticles]=1.0;RY[NumberOfParticles]=0.0;RZ[NumberOfParticles]=1.0;
 RX[NumberOfParticles+1]=0.0;RY[NumberOfParticles+1]=0.0;RZ[NumberOfParticles+1]=1.0;
 RX[NumberOfParticles+2]=0.0;RY[NumberOfParticles+2]=0.0;RZ[NumberOfParticles+2]=0.0;
 BackBoneList[NChains][0]=NumberOfParticles; 
 BackBoneList[NChains][1]=NumberOfParticles+1; 
 BackBoneList[NChains][2]=NumberOfParticles+2; 
 BackBoneList[NChains][3]=NumberOfParticles+3; 
 fp=fopen("utorsion.dat","w"); 
 for(i=0;i<1000;i++)
  { 
  RX[NumberOfParticles+3]=cos(2.0*M_PI*i/1000.0);RY[NumberOfParticles+3]=sin(2.0*M_PI*i/1000.0);RZ[NumberOfParticles+3]=0.0;
  fprintf(fp,"theta = %lf\tu(theta) = %lf\n",i/1000.0*360.0-180.0,TorsionPotential(NChains,1));
  }
 fclose(fp);

return;
}

