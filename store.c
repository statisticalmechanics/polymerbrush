/********************************************
 * for Gear predictor-corrector algorithm
 * initialize RX, PX etc.
 * store position, velocity , force
 ********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void Store0(void) // initialization
{
 int i;

 for(i=0;i<NumberOfParticles;i++)
 {
  //r
  RX[i] = position[i].x;
  RY[i] = position[i].y;
  RZ[i] = position[i].z;
  //r'=v
//  RX1[i] = velocity[i].x;
//  RY1[i] = velocity[i].y;
//  RZ1[i] = velocity[i].z;
  // r'' = f/m
  RX2[i] = force[i].x/mass[i];
  RY2[i] = force[i].y/mass[i];
  RZ2[i] = force[i].z/mass[i];
  // r'''
  RX3[i] = 0.;
  RY3[i] = 0.;
  RZ3[i] = 0.;
  // p = mv
  PX[i]  = mass[i]*velocity[i].x;
  PY[i]  = mass[i]*velocity[i].y;
  PZ[i]  = mass[i]*velocity[i].z;
  // p' = f
  PX1[i] = force[i].x;
  PY1[i] = force[i].y;
  PZ1[i] = force[i].z;
  // p''
  PX2[i] = 0.;
  PY2[i] = 0.;
  PZ2[i] = 0.;
  // = p'''
  PX3[i] = 0.;
  PY3[i] = 0.;
  PZ3[i] = 0.;
 }

 //L1 = 0.;
 //L2 = 0.;
 //L3 = 0.;
 //L'
 LX1 = 0.;
 LY1 = 0.;
 LZ1 = 0.;
 //L''
 LX2 = 0.;
 LY2 = 0.;
 LZ2 = 0.;
 //L'''
 LX3 = 0.;
 LY3 = 0.;
 LZ3 = 0.;

 return;
}
