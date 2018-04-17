/***************************************************
* Scale velocity to make isokinetic
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

void ScaleVelocity(void) //  need to use after Kinetic()
{
 int i;
 double scale;

 /**************************************************
 Kinstant and Tinstant need to be done by Kinetic()
 ***************************************************/
 scale = sqrt(T/Tinstant);

 for(i=0;i<NumberOfParticles;i++) 
 {
  //p
  PX[i] *= scale;
  PY[i] *= scale;
  PZ[i] *= scale;
  //r'=v
  velocity[i].x *= scale;
  velocity[i].y *= scale;
  velocity[i].z *= scale;
 }
 Tinstant = T;
 Kinstant = Nf*kB/2.0*Tinstant;
 return;
}


void ScaleVolume(void) //  need to use after Kinetic()
{
 int i;
 double scale;

 /**************************************************
 Kinstant and Tinstant need to be done by Kinetic()
 ***************************************************/
 Pinstant = rho*Tinstant + Virial/V;

 scale = 1.0 + (Pinstant - P)/(Pinstant + HyperVirial/V)/3.0; 
 
 for(i=0;i<NumberOfParticles;i++) 
 {
  RX[i] *= scale;
  RY[i] *= scale;
  RZ[i] *= scale;
  position[i].x *= scale;
  position[i].y *= scale;
  position[i].z *= scale;
 }
 LX *= scale;
 LY *= scale;
 LZ *= scale;
 V = LX*LY*LZ;
 rho = NumberOfParticles/V*CUBIC(atomsize[0]);

 VerletCheck(); // check if Verlet list needs to be updated
 Force();
 Pinstant = rho*Tinstant + Virial/V;

 return;
}


void ScaleNPT(void) //  need to use after Kinetic()
{
 int i;
 double scale;

 /**************************************************
 Kinstant and Tinstant need to be done by Kinetic()
 ***************************************************/
 scale = sqrt(T/Tinstant);

 for(i=0;i<NumberOfParticles;i++) 
 {
  PX[i] *= scale;
  PY[i] *= scale;
  PZ[i] *= scale;
  // need to rescale r'?
 // velocity[i].x *= scale;
 // velocity[i].y *= scale;
 // velocity[i].z *= scale;
 }

 Tinstant = T;
 Kinstant = Nf*kB/2.0*Tinstant;

 Pinstant = rho*Tinstant + Virial/V;

 scale = 1.0 + (Pinstant - P)/(Pinstant + HyperVirial/V)/3.0; 
 
 for(i=0;i<NumberOfParticles;i++) 
 {
  RX[i] *= scale;
  RY[i] *= scale;
  RZ[i] *= scale;
  position[i].x *= scale;
  position[i].y *= scale;
  position[i].z *= scale;
 }
 LX *= scale;
 LY *= scale;
 LZ *= scale;
 V = LX*LY*LZ;
 rho = NumberOfParticles/V*CUBIC(atomsize[0]);

 VerletCheck(); // check if Verlet list needs to be updated
 Force();
 Pinstant = rho*Tinstant + Virial/V;

 return;
}


