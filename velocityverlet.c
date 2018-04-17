/***************************************
 * integrate equation of motion
 * via velocity Verlet algorithm
 * in two parts 1 and 2
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void VelocityVerlet(int part)
{
  int i;
  
  if(part == 1) // velocity Verlet step 1
  {
  for(i=0;i<NumberOfParticles;i++)
  if(particlefrozen[i]!=1)
  {
  //r(t+dt) = r(t) + v(t)*dt + f/2m *dt^2
   //position[i].x = position[i].x + velocity[i].x*dt + SQR(dt)/2./mass[i]*force[i].x;
   //position[i].y = position[i].y + velocity[i].y*dt + SQR(dt)/2./mass[i]*force[i].y;
   //position[i].z = position[i].z + velocity[i].z*dt + SQR(dt)/2./mass[i]*force[i].z;
   RX[i] += velocity[i].x*dt + SQR(dt)/2./mass[i]*force[i].x;
   RY[i] += velocity[i].y*dt + SQR(dt)/2./mass[i]*force[i].y;
   RZ[i] += velocity[i].z*dt + SQR(dt)/2./mass[i]*force[i].z;
   
   position[i].x = RX[i];
   position[i].y = RY[i];
   position[i].z = RZ[i];
   //PBC
   if(BCType[0]==0) PBC(&(position[i].x),LX);
   if(BCType[1]==0) PBC(&(position[i].y),LY);
   if(BCType[2]==0) PBC(&(position[i].z),LZ);

  // v' = v(t) + f(t)/2m * dt
   velocity[i].x = velocity[i].x + dt/2./mass[i]*force[i].x;
   velocity[i].y = velocity[i].y + dt/2./mass[i]*force[i].y;
   velocity[i].z = velocity[i].z + dt/2./mass[i]*force[i].z;
  }//end loop i
  } //endif part 1
  else if(part == 2) // velocity Verlet step 2, force needs to be updated
  {
  for(i=0;i<NumberOfParticles;i++)
  if(particlefrozen[i]!=1)
  {
  // v(t+dt) = v' + f(t+dt)/2m * dt
   velocity[i].x = velocity[i].x + dt/2./mass[i]*force[i].x;
   velocity[i].y = velocity[i].y + dt/2./mass[i]*force[i].y;
   velocity[i].z = velocity[i].z + dt/2./mass[i]*force[i].z;
  }//end loop i
  }

 return;
}
