/********************************************
 * NVT ensemble constraint method
 * Gear Predictor-Corrector
 * Predict Step
 ********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void PredictNVT(double DT)
{
 int i;
 double c1,c2,c3;

 c1 = DT;
 c2 = c1*DT/2.0;
 c3 = c2*DT/3.0;

 for(i=0;i<NumberOfParticles;i++)
 if(particlefrozen[i]!=1)
 {
  // original r
  //RX[i]  = RX[i]  + c1*RX1[i] + c2*RX2[i] + c3*RX3[i];
 // RY[i]  = RY[i]  + c1*RY1[i] + c2*RY2[i] + c3*RY3[i];
 // RZ[i]  = RZ[i]  + c1*RZ1[i] + c2*RZ2[i] + c3*RZ3[i];
  RX[i]  = RX[i]  + c1*velocity[i].x + c2*RX2[i] + c3*RX3[i];
  RY[i]  = RY[i]  + c1*velocity[i].y + c2*RY2[i] + c3*RY3[i];
  RZ[i]  = RZ[i]  + c1*velocity[i].z + c2*RZ2[i] + c3*RZ3[i];
  // back to box r
  position[i].x  =  RX[i];
  position[i].y  =  RY[i];
  position[i].z  =  RZ[i];
  if(BCType[0]==0) PBC(&(position[i].x),LX);
  if(BCType[1]==0) PBC(&(position[i].y),LY);
  if(BCType[2]==0) PBC(&(position[i].z),LZ);
  
//  RX1[i] = RX1[i] + c1*RX2[i] + c2*RX3[i];
//  RY1[i] = RY1[i] + c1*RY2[i] + c2*RY3[i];
//  RZ1[i] = RZ1[i] + c1*RZ2[i] + c2*RZ3[i];
  // v = r'
  velocity[i].x = velocity[i].x + c1*RX2[i] + c2*RX3[i];
  velocity[i].y = velocity[i].y + c1*RY2[i] + c2*RY3[i];
  velocity[i].z = velocity[i].z + c1*RZ2[i] + c2*RZ3[i];
  // v' = r''
  RX2[i] = RX2[i] + c1*RX3[i];
  RY2[i] = RY2[i] + c1*RY3[i];
  RZ2[i] = RZ2[i] + c1*RZ3[i];
  //momentum p
  PX[i]  = PX[i]  + c1*PX1[i] + c2*PX2[i] + c3*PX3[i];
  PY[i]  = PY[i]  + c1*PY1[i] + c2*PY2[i] + c3*PY3[i];
  PZ[i]  = PZ[i]  + c1*PZ1[i] + c2*PZ2[i] + c3*PZ3[i];
  // p'
  PX1[i] = PX1[i] + c1*PX2[i] + c2*PX3[i];
  PY1[i] = PY1[i] + c1*PY2[i] + c2*PY3[i];
  PZ1[i] = PZ1[i] + c1*PZ2[i] + c2*PZ3[i];
  // p''
  PX2[i] = PX2[i] + c1*PX3[i];
  PY2[i] = PY2[i] + c1*PY3[i];
  PZ2[i] = PZ2[i] + c1*PZ3[i];
 }
 return;
}

/********************************************
 * NVT ensemble
 * Gear Predictor-Corrector
 * Correct Step
 * 4-value, 1st order ODE
 ********************************************/

void CorrectNVT(double DT,double XI)
{
 int i;
 double c1,c2,c3;
 double coeff0,coeff2,coeff3;
 double corrx, corry,corrz,corpx,corpy,corpz;
 double RX1I,RY1I,RZ1I,PX1I,PY1I,PZ1I;
 double GEAR0,GEAR1,GEAR2,GEAR3;

 GEAR0 = 3.0/8.0;
 GEAR1 = 1.0; // not used
 GEAR2 = 3.0/4.0;
 GEAR3 = 1.0/6.0;
 
 c1 = DT;
 c2 = c1*DT/2.0;
 c3 = c2*DT/3.0;

 coeff0 = GEAR0*c1;
 //coeff1 = GEAR1*c1/c1 = 1;
 coeff2 = GEAR2*c1/c2;
 coeff3 = GEAR3*c1/c3;

 for(i=0;i<NumberOfParticles;i++)
  if(particlefrozen[i]!=1)
 {
  // r' = p/m
  RX1I = PX[i]/mass[i];
  RY1I = PY[i]/mass[i];
  RZ1I = PZ[i]/mass[i];
// corrx = RX1I - RX1[i];
// corry = RY1I - RY1[i];
// corrz = RZ1I - RZ1[i];
  corrx = RX1I - velocity[i].x;
  corry = RY1I - velocity[i].y;
  corrz = RZ1I - velocity[i].z;

  // r
  RX[i] = RX[i] + coeff0*corrx;
  RY[i] = RY[i] + coeff0*corry;
  RZ[i] = RZ[i] + coeff0*corrz;

  position[i].x = RX[i];
  position[i].y = RY[i];
  position[i].z = RZ[i];
  if(BCType[0]==0) PBC(&(position[i].x),LX);
  if(BCType[1]==0) PBC(&(position[i].y),LY);
  if(BCType[2]==0) PBC(&(position[i].z),LZ);

//  RX1[i] = RX1I; 
//  RY1[i] = RY1I; 
//  RZ1[i] = RZ1I; 
  velocity[i].x = RX1I; 
  velocity[i].y = RY1I; 
  velocity[i].z = RZ1I; 
  RX2[i] = RX2[i] + coeff2*corrx;
  RY2[i] = RY2[i] + coeff2*corry;
  RZ2[i] = RZ2[i] + coeff2*corrz;
  RX3[i] = RX3[i] + coeff3*corrx;
  RY3[i] = RY3[i] + coeff3*corry;
  RZ3[i] = RZ3[i] + coeff3*corrz;
  
  //p' = f - xi*p
  PX1I = force[i].x - XI*PX[i];
  PY1I = force[i].y - XI*PY[i];
  PZ1I = force[i].z - XI*PZ[i];
  corpx = PX1I - PX1[i];
  corpy = PY1I - PY1[i];
  corpz = PZ1I - PZ1[i];

  PX[i] = PX[i] + coeff0*corpx;
  PY[i] = PY[i] + coeff0*corpy;
  PZ[i] = PZ[i] + coeff0*corpz;
  PX1[i] = PX1I; 
  PY1[i] = PY1I; 
  PZ1[i] = PZ1I; 
  PX2[i] = PX2[i] + coeff2*corpx;
  PY2[i] = PY2[i] + coeff2*corpy;
  PZ2[i] = PZ2[i] + coeff2*corpz;
  PX3[i] = PX3[i] + coeff3*corpx;
  PY3[i] = PY3[i] + coeff3*corpy;
  PZ3[i] = PZ3[i] + coeff3*corpz;
 }

 return;
}
