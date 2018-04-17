/***************************************
 * calculate of force fx,fy,fz of 
 * each particle i at corrent position t
 * calculate potential energy Upotential
 * and virial = 1/d <sum f*r> as well
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void Force(void) 
{
  int i,j,k,nb;
  int jbegin,jend,jlist;
  double rvij,rvij2; // Verlet list radius
  double rcij,rcij2; // rcutoff
  double sigmaij,sigmaij2;
  double rij,rij2,rij6; // avoid sqrt() if possible
  double xij,yij,zij;
 // double dudr; // du/dr
  double fij;
  double fxij,fyij,fzij;
  double pxij,pyij,pzij;
  double wij,uij,Xij,rijpij,pijfij;

  double Ubond,Ubondangle,Utorsion;
  VECTOR fangleA,fangleB;
  VECTOR ftorsion;
   
  Upotential = 0.;
  Virial = 0.;
  HyperVirial = 0.;
  PdotF = 0.;
  RdotPdotX = 0.;


  for(i=0;i<NumberOfParticles;i++)
  {
    force[i].x = 0.;
    force[i].y = 0.;
    force[i].z = 0.;
  }
  
  if(WallSwitch == 1)
  {
  for(i=0;i<NumberOfParticles;i++)
  {
    zij = position[i].z- (-rminLJ);
    if(zij < rminLJ) fzij = 48.0*1.0/zij*(pow(1.0/zij,12.0)-pow(1.0/zij,6.0)/2.0);
    else fzij = 0.;

    force[i].z += fzij;
  }
  }

  Ubond = 0.0;
  for(i=0;i<NumberOfParticles;i++)
  {
   for(nb=0;nb<NSprings[i];nb++)
   {
    j = NSpringList[i][nb];

    //should not use minimum image for bond 
    xij = position[i].x - position[j].x;
    yij = position[i].y - position[j].y;
    zij = position[i].z - position[j].z;
 //   if(MinimumImage[0]==0) xij = xij - LX*round(xij/LX);
 //   if(MinimumImage[1]==0) yij = yij - LY*round(yij/LY);
 //   if(MinimumImage[2]==0) zij = zij - LZ*round(zij/LZ);
    rij2 = SQR(xij)+SQR(yij)+SQR(zij);

 //   rij2 = Distance(i,j);
    rij = sqrt(rij2); // avoid sqrt if possible
    uij =  SpringPotential(i,nb,rij); Ubond += uij;
    fij =  -SpringPotential_dr(i,nb,rij);// printf("%lf\n",fij);
    fxij = fij*xij/rij;
    fyij = fij*yij/rij;
    fzij = fij*zij/rij;

    force[i].x += fxij;
    force[i].y += fyij;
    force[i].z += fzij;
   }
  }
  Ubond /= 2.0; // double acount
  
  Ubondangle = 0.0;
  if(BendSwitch == 1)
  {
  for(i=0;i<NumberOfParticles;i++)
  {
   for(nb=0;nb<NBondAnglesA[i];nb++)
   {
    j=BondAngleListA[i][nb][0];k=BondAngleListA[i][nb][1];
    Ubondangle += BondAnglePotential(i,j,k,BondAngleStiffA[i][nb],BondAngleA[i][nb]);
    fangleA = BondAngleForceA(i,j,k,BondAngleStiffA[i][nb],BondAngleA[i][nb]);
    force[i].x += fangleA.x;
    force[i].y += fangleA.y;
    force[i].z += fangleA.z;
   }
   
   for(nb=0;nb<NBondAnglesB[i];nb++)
   {
    j=BondAngleListB[i][nb][0];k=BondAngleListB[i][nb][1];
    fangleB = BondAngleForceB(i,j,k,BondAngleStiffB[i][nb],BondAngleB[i][nb]);
    force[i].x += fangleB.x;
    force[i].y += fangleB.y;
    force[i].z += fangleB.z;
   }
  }
  }
  
  
 Utorsion = 0.0;
if(runtime>10.0 && TorsionSwitch == 1)
{
  //torsion of backbone
  for(j=0;j<NChains;j++)
  {
   for(nb=0;nb<NBackBones;nb++)
   {
    i = BackBoneList[j][nb]; 
    uij = TorsionPotential(j,nb); 
    Utorsion += uij; 
    ftorsion = TorsionForce(j,nb);
    force[i].x += ftorsion.x;
    force[i].y += ftorsion.y;
    force[i].z += ftorsion.z;
   }
  }

 //torsion of side chain
  for(i=0;i<NumberOfParticles;i++)
  {
    //j = i/NAtomPerChain;// which chain
    j = (i%NAtomPerChain)/NAtomPerMonomer; // which monomer
    k = (i%NAtomPerChain)%NAtomPerMonomer; // which atom in the monomer
    uij = TorsionPotentialPhenyl(j,k,i); 
    Utorsion += uij; 
    ftorsion = TorsionForcePhenyl(j,k,i);
    force[i].x += ftorsion.x;
    force[i].y += ftorsion.y;
    force[i].z += ftorsion.z;
  }
}
   // use to calculate Xi in NVT ???
  for(i=0;i<NumberOfParticles;i++)
   PdotF += (PX[i]*force[i].x + PY[i]*force[i].y + PZ[i]*force[i].z)/mass[i];

  // need to update Verlet list, Verlet list is (re)generated
  if(VerletCheckIndex == 1)
  {
   for(i=0;i<NumberOfParticles;i++)
   {
    position0[i].x = RX[i];
    position0[i].y = RY[i];
    position0[i].z = RZ[i];
   }// end store old position loop

   nlist = -1; //Verlet list index   
   for(i=0;i<NumberOfParticles-1;i++)
   {
    point[i] = nlist + 1; //starting index of particle i in the list

    for(j=i+1;j<NumberOfParticles;j++)
    {
     sigmaij = SigmaIJ(i,j);
     rvij = rv*sigmaij;
     rvij2 = SQR(rvij);
     rcij = rc*sigmaij;
     rcij2 = SQR(rcij);

     xij = position[i].x - position[j].x;
     yij = position[i].y - position[j].y;
     zij = position[i].z - position[j].z;
     if(MinimumImage[0]==0) xij = xij - LX*round(xij/LX);
     if(MinimumImage[1]==0) yij = yij - LY*round(yij/LY);
     if(MinimumImage[2]==0) zij = zij - LZ*round(zij/LZ);
     rij2 = SQR(xij)+SQR(yij)+SQR(zij);

  //   rij2 = Distance(i,j); // not working for accuracy of position[i].x - position[y].x causing Uspring diverge

     if(rij2 < rvij2) // if within Verlet neighbor cell
     {
      nlist++;
      VerletList[nlist] = j;

      if(rij2 < rcij2) // try avoid sqrt for r>rc
      {
       rij = sqrt(rij2); // avoid sqrt if possible
       uij =  Potential(i,j,rij);
       fij =  -Potential_dr(i,j,rij);
       fxij = fij*xij/rij;
       fyij = fij*yij/rij;
       fzij = fij*zij/rij;

       force[i].x += fxij;
       force[i].y += fyij;
       force[i].z += fzij;

       force[j].x += -fxij;
       force[j].y += -fyij;
       force[j].z += -fzij;

   //    pxij = mass[i]*velocity[i].x - mass[j]*velocity[j].x;
   //    pyij = mass[i]*velocity[i].y - mass[j]*velocity[j].y;
   //    pzij = mass[i]*velocity[i].z - mass[j]*velocity[j].z;
       pxij = PX[i] - PX[j];
       pyij = PY[i] - PY[j];
       pzij = PZ[i] - PZ[j];

       rijpij = xij*pxij + yij*pyij + zij * pzij;
       pijfij = pxij*fxij + pyij*fyij + pzij*fzij; 
       wij = fij*rij;
       Xij = -wij + rij2*Potential_dr2(i,j,rij);

       Upotential += uij;
       Virial += wij;
       HyperVirial += Xij;
       PdotF += pijfij/mass[i]; // mass?
       RdotPdotX += rijpij*Xij/rij2/mass[i]; 
      } //endif rij < rcutoff
     } // endif rij < rv
    }//end loop j
   }//end loop i
   point[NumberOfParticles-1] = nlist +1;
 
  }// endif 1
  // no need to update, Verlet list is used
  else if (VerletCheckIndex == 0)
  {
   for(i=0;i<NumberOfParticles-1;i++)
   {
    jbegin = point[i];
    jend = point[i+1]-1;
    if(jbegin<=jend) // if list is not empty
    {
     for(jlist=jbegin;jlist<=jend;jlist++)
     {
      j = VerletList[jlist];

      xij = position[i].x - position[j].x;
      yij = position[i].y - position[j].y;
      zij = position[i].z - position[j].z;
     if(MinimumImage[0]==0) xij = xij - LX*round(xij/LX);
     if(MinimumImage[1]==0) yij = yij - LY*round(yij/LY);
     if(MinimumImage[2]==0) zij = zij - LZ*round(zij/LZ);
      rij2 = SQR(xij)+SQR(yij)+SQR(zij);
 //    rij2 = Distance(i,j);

      sigmaij = SigmaIJ(i,j);
      rcij = rc*sigmaij;
      rcij2 = SQR(rcij);
      if(rij2 < rcij2) // try avoid sqrt for r>rc
      {
       rij = sqrt(rij2); // avoid sqrt if possible
       uij =  Potential(i,j,rij);
       fij =  -Potential_dr(i,j,rij);
       fxij = fij*xij/rij;
       fyij = fij*yij/rij;
       fzij = fij*zij/rij;

       force[i].x += fxij;
       force[i].y += fyij;
       force[i].z += fzij;

       force[j].x += -fxij;
       force[j].y += -fyij;
       force[j].z += -fzij;
       
    //   pxij = mass[i]*velocity[i].x - mass[j]*velocity[j].x;
    //   pyij = mass[i]*velocity[i].y - mass[j]*velocity[j].y;
    //   pzij = mass[i]*velocity[i].z - mass[j]*velocity[j].z;
       pxij = PX[i] - PX[j];
       pyij = PY[i] - PY[j];
       pzij = PZ[i] - PZ[j];
       
       rijpij = xij*pxij + yij*pyij + zij * pzij;
       pijfij = pxij*fxij + pyij*fyij + pzij*fzij; 
       wij = fij*rij;
       Xij = -wij + rij2*Potential_dr2(i,j,rij);

       Upotential += uij;
       Virial += wij;
       HyperVirial += Xij;
       PdotF += pijfij/mass[i];// mass?
       RdotPdotX += rijpij*Xij/rij2/mass[i]; 
      } //endif rij < rcutoff
     }// end loop jlist
    }// endif list is not empty
   } //end loop i
  }// endif 0

 Upotential += Ubond;
 Upotential += Ubondangle;
 Upotential += Utorsion;

 Virial = Virial/3.0; // dimension
 HyperVirial = HyperVirial/9.0; // dimension^2

 return;

} // end function Force
