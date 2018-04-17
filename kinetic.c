/***************************************************
* sampling
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

void Kinetic(void)
{
 int i;
 int j;
 double Mcom; // total mass
 VECTOR Vcom,Pcom;

 PdotP = 0.;
 Kinstant = 0.;
 


 if(COMSwitch == 1)
 {
  Mcom = 0.;
  Vcom.x = 0.; Vcom.y = 0.; Vcom.z = 0.;
  Pcom.x = 0.; Pcom.y = 0.; Pcom.z = 0.;
  for(i=0;i<NumberOfParticles;i++) 
  {
   Mcom += mass[i];
   Vcom.x += mass[i]*velocity[i].x;
   Vcom.y += mass[i]*velocity[i].y;
   Vcom.z += mass[i]*velocity[i].z;
   Pcom.x += PX[i];Pcom.y += PY[i];Pcom.z += PZ[i];
  }
   Vcom.x /= Mcom; Vcom.y /= Mcom;Vcom.z /= Mcom;
  for(i=0;i<NumberOfParticles;i++) 
  {
   velocity[i].x -=Vcom.x;
   velocity[i].y -=Vcom.y;
   velocity[i].z -=Vcom.z;
   PX[i] -= Pcom.x;
   PY[i] -= Pcom.y;
   PZ[i] -= Pcom.z;
  }
 }



 for(i=0;i<NumberOfParticles;i++) 
 {
  if(EnsembleType == 0 || EnsembleType == 1) // NVE or NVK using velocity Verlet
  Kinstant += mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));
  else // NVT or NPH or NPT using constraint
  {
  Kinstant += (SQR(PX[i]) + SQR(PY[i]) + SQR(PZ[i]))/mass[i];
  PdotP += (SQR(PX[i]) + SQR(PY[i]) + SQR(PZ[i]))/mass[i];
  }
 }// end loop i and if particle not frozen
 
 Kinstant *= 0.5; // instantenous kinetic energy
 // 0.5*kT*Nf = K = 0.5* sum_mv^2, Nf = 3N-3
 // Tinstant = 2.0*Kinstant/(3.*NumberOfParticles-3.)/kB;
 Tinstant = 2.0*Kinstant/Nf/kB;
 
 return;
}

//check center of mass velocity
void COMcheck(FILE *fp)
{
 int i;
 double Mcom; // total mass
 VECTOR Vcom;

  Mcom = 0.;
  Vcom.x = 0.;
  Vcom.y = 0.;
  Vcom.z = 0.;
  Kinstant = 0.;
  for(i=0;i<NumberOfParticles;i++) 
  {
   Mcom += mass[i];
   Vcom.x += mass[i]*velocity[i].x;
   Vcom.y += mass[i]*velocity[i].y;
   Vcom.z += mass[i]*velocity[i].z;
   Kinstant += mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));

   fprintf(fp,"%lf\n",sqrt(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z)));
  }

   Vcom.x /= Mcom;
   Vcom.y /= Mcom;
   Vcom.z /= Mcom;
   Kinstant *= 0.5;
   Tinstant = 2.0*Kinstant/Nf/kB;
   printf("center of mass velosity (%lf, %lf, %lf) and instantenous Tinstant = %lf\n",Vcom.x,Vcom.y,Vcom.z,Tinstant);
/********************************************************************/

 return;
}
