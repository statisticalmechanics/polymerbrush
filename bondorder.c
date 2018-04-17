/***************************************************
* sampling bondorder 
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

double P2COS(int type)
{
 double p2;
 double rx,ry,rz;
 double r2;
 double costheta2;
 int i;

 double a1,a2,a3,b1,b2,b3;

 p2 = 0.0;
 if(type == 1) // backbone
 {
 for(i=0;i<NBackboneBonds;i++)
 {
  rx = RX[BackboneBondIndexB[i]] - RX[BackboneBondIndexA[i]];
  ry = RY[BackboneBondIndexB[i]] - RY[BackboneBondIndexA[i]];
  rz = RZ[BackboneBondIndexB[i]] - RZ[BackboneBondIndexA[i]];
  r2 = SQR(rx)+SQR(ry)+SQR(rz);
  costheta2 = SQR(rz)/r2;
  p2 += 1.5*costheta2 - 0.5;
 }
 p2 /= NBackboneBonds;
 }
 if(type == 2) // phenyl
 {
 for(i=0;i<NPhenylBonds;i++)
 {
  rx = RX[PhenylBondIndexB[i]] - RX[PhenylBondIndexA[i]];
  ry = RY[PhenylBondIndexB[i]] - RY[PhenylBondIndexA[i]];
  rz = RZ[PhenylBondIndexB[i]] - RZ[PhenylBondIndexA[i]];
  r2 = SQR(rx)+SQR(ry)+SQR(rz);
  costheta2 = SQR(rz)/r2;
  p2 += 1.5*costheta2 - 0.5;
 }
 // printf("side chain (%lf %lf %lf) ",rx,ry,rz);
 p2 /= NPhenylBonds;
 }
 if(type == 3) // ring
 {
 for(i=0;i<NRings;i++)
 {
  a1=RX[RingIndexB[i]]-RX[RingIndexA[i]];
  a2=RY[RingIndexB[i]]-RY[RingIndexA[i]];
  a3=RZ[RingIndexB[i]]-RZ[RingIndexA[i]];
  b1=RX[RingIndexC[i]]-RX[RingIndexA[i]];
  b2=RY[RingIndexC[i]]-RY[RingIndexA[i]];
  b3=RZ[RingIndexC[i]]-RZ[RingIndexA[i]];
  rx = a2*b3-a3*b2;
  ry = a3*b1-a1*b3;
  rz = a1*b2-a2*b1;
  r2 = SQR(rx)+SQR(ry)+SQR(rz);
  costheta2 = SQR(rz)/r2;
  p2 += 1.5*costheta2 - 0.5;
 }
 // printf("a (%lf %lf %lf) ",a1,a2,a3);
 // printf("b (%lf %lf %lf) ",b1,b2,b3);
 // printf("ring (%lf %lf %lf)\n",rx,ry,rz);
 p2 /= NRings;
 }
 if(type == 4) // ring wrt x
 {
 for(i=0;i<NRings;i++)
 {
  a1=RX[RingIndexB[i]]-RX[RingIndexA[i]];
  a2=RY[RingIndexB[i]]-RY[RingIndexA[i]];
  a3=RZ[RingIndexB[i]]-RZ[RingIndexA[i]];
  b1=RX[RingIndexC[i]]-RX[RingIndexA[i]];
  b2=RY[RingIndexC[i]]-RY[RingIndexA[i]];
  b3=RZ[RingIndexC[i]]-RZ[RingIndexA[i]];
  rx = a2*b3-a3*b2;
  ry = a3*b1-a1*b3;
  rz = a1*b2-a2*b1;
  r2 = SQR(rx)+SQR(ry)+SQR(rz);
  costheta2 = SQR(rx)/r2;
  p2 += 1.5*costheta2 - 0.5;
 }
 p2 /= NRings;
 }
 if(type == 5) // ring wrt y
 {
 for(i=0;i<NRings;i++)
 {
  a1=RX[RingIndexB[i]]-RX[RingIndexA[i]];
  a2=RY[RingIndexB[i]]-RY[RingIndexA[i]];
  a3=RZ[RingIndexB[i]]-RZ[RingIndexA[i]];
  b1=RX[RingIndexC[i]]-RX[RingIndexA[i]];
  b2=RY[RingIndexC[i]]-RY[RingIndexA[i]];
  b3=RZ[RingIndexC[i]]-RZ[RingIndexA[i]];
  rx = a2*b3-a3*b2;
  ry = a3*b1-a1*b3;
  rz = a1*b2-a2*b1;
  r2 = SQR(rx)+SQR(ry)+SQR(rz);
  costheta2 = SQR(ry)/r2;
  p2 += 1.5*costheta2 - 0.5;
 }
 p2 /= NRings;
 }

 return p2;
}

/***************************************************
* sampling radius of gyration Rg and end-to-end distance
****************************************************/

void RadiusGyration(void)
{
 int j,k,indexshift,i,i0,in;
 double rmx,rmy,rmz;
 
 for(j=0;j<NChains;j++)
 {
  indexshift = j*NAtomPerChain;
  i0 = indexshift; in = indexshift+NAtomPerChain-7;
  dRend[j] = sqrt(SQR(RX[in]-RX[i0])+SQR(RY[in]-RY[i0])+SQR(RZ[in]-RZ[i0])); 

  rmx=0.0;rmy=0.0;rmz=0.0;
  for(k=0;k<NAtomPerChain;k++)
  {
 //  if(k%NAtomPerMonomer==0)
   {
   i = k+indexshift;
   rmx += RX[i];rmy += RY[i];rmz += RZ[i];
   }
  }
  rmx /= NAtomPerChain;rmy /= NAtomPerChain;rmz /= NAtomPerChain;
 // rmx /= NMonomerPerChain;rmy /= NMonomerPerChain;rmz /= NMonomerPerChain;

  RGyration[j] = 0.0;
  PolymerHeight[j] = 0.0;
  for(k=0;k<NAtomPerChain;k++)
  {
 //  if(k%NAtomPerMonomer==0)
   {
   i = k+indexshift;
   RGyration[j] += SQR(RX[i]-rmx)+SQR(RY[i]-rmy)+SQR(RZ[i]-rmz);
   }
   if(RZ[i]>PolymerHeight[j]) PolymerHeight[j] = RZ[i];
  }
   RGyration[j] /= NAtomPerChain; RGyration[j] = sqrt(RGyration[j]);
 //  RGyration[j] /= NMonomerPerChain; RGyration[j] = sqrt(RGyration[j]);
 }

 return;
}

void BondCorrelation(void)
{
 int i,j,k;
 int i0,j0;
 double li,lj,lilj;
 double lix,liy,liz;
 double ljx,ljy,ljz;

 for(i=0;i<NChains;i++)
 {
  i0 = i*NAtomPerChain;
  lix = RX[i0+1]-RX[i0];liy = RY[i0+1]-RY[i0];liz = RZ[i0+1]-RZ[i0];
  li = sqrt(SQR(lix)+SQR(liy)+SQR(liz));
  
  for(j=1;j<=NMonomerPerChain*2-2;j++)
  {
   if(j%2==0){i0=i*NAtomPerChain+(j/2)*NAtomPerMonomer;j0=i0+1;}
   else {i0=i*NAtomPerChain+(j/2)*NAtomPerMonomer+1;j0=i0+7;}
   //printf("%d %d %d %d\n",i,j,i0,j0);
   ljx = RX[j0]-RX[i0];ljy = RY[j0]-RY[i0];ljz = RZ[j0]-RZ[i0];
   lj = sqrt(SQR(ljx)+SQR(ljy)+SQR(ljz));
   lilj = lix*ljx+liy*ljy+liz*ljz;
   COSij[i][j] = lilj/li/lj;
  }
 }


 return;
}
