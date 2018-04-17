/***************************************************
* return distance r, u(r), du/dr, u_tail, w_tail
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void PBC(double *xx,double box) // put back to box
{
  /*
  if(*xx < 0.0)
   *xx += box;
  else if(*xx >= box)
   *xx -= box;
  */
  
  if(*xx < 0.0)
  {
  do
   *xx += box;
  while(*xx < 0.0);
  }
  else if(*xx >= box)
  {
  do
   *xx -= box; 
  while(*xx >= box);
  }

 return;
}

double Distance(int i,int j) // return rij^2
{
  double r2,dx,dy,dz;
   
  // miminum image convention under periodic boundary condition
  dx = position[i].x-position[j].x;
  dy = position[i].y-position[j].y;
  dz = position[i].z-position[j].z;
 
  if(MinimumImage[0]==0) dx = dx - LX*round(dx/LX);
  if(MinimumImage[1]==0) dy = dy - LY*round(dy/LY);
  if(MinimumImage[2]==0) dz = dz - LZ*round(dz/LZ);
 
  r2 = SQR(dx)+SQR(dy)+SQR(dz); // reduced unit
  //r = sqrt(SQR(dx)+SQR(dy)+SQR(dz));

  return r2;
}

double EpsilonIJ(int i,int j)
{
 double epsilonij;

 epsilonij = (epsilon[i]+epsilon[j])/2.0;
 //epsilonij = sqrt(epsilon[i]*epsilon[j]);

 return epsilonij;

}

double SigmaIJ(int i,int j)
{
 double sigmaij;

 sigmaij = (sigma[i]+sigma[j])/2.0;

 return sigmaij;
 
}

double Potential(int i,int j,double r) // u(rij)
{
  double u;
  double sigmaij,epsilonij,rcij;
  double uc,ducdr;

  sigmaij = SigmaIJ(i,j);
  rcij = rc*sigmaij;
  epsilonij = EpsilonIJ(i,j);
  uc = 4.0*epsilonij*(pow(sigmaij/rcij,12.0)-pow(sigmaij/rcij,6.0));
  ducdr = -48.0*epsilonij/rcij*(pow(sigmaij/rcij,12.0) - pow(sigmaij/rcij,6.0)/2.0);

  if(PotentialType == 0) // L-J
  {
   if(r < rcij) u = 4.0*epsilonij*(pow(sigmaij/r,12.0)-pow(sigmaij/r,6.0));
   else u = 0.;
  } //endif L-J
  else if(PotentialType == 1) // shifted L-J
  {
   if(r < rcij) u = 4.0*epsilonij*(pow(sigmaij/r,12.0)-pow(sigmaij/r,6.0)) - uc;
   else u = 0.;
  } //endif shifted L-J
  else if(PotentialType == 2) // shifted-force L-J  1/r^12 - 1/r^6
  {
   if(r < rcij) u = 4.0*epsilonij*(pow(sigmaij/r,12.0)-pow(sigmaij/r,6.0)) - uc-ducdr*(r-rcij);
   else u = 0.;
  } // endif shifted-force L-J

  return u;
}

double Potential_dr(int i,int j,double r) // du(rij)/dr
{
  double dudr;
  double sigmaij,epsilonij,rcij;
  double uc,ducdr;

  sigmaij = SigmaIJ(i,j);
  rcij = rc*sigmaij;
  epsilonij = EpsilonIJ(i,j);
 // uc = 4.0*epsilonij*(pow(sigmaij/rcij,12.0)-pow(sigmaij/rcij,6.0));
  ducdr = -48.0*epsilonij/rcij*(pow(sigmaij/rcij,12.0) - pow(sigmaij/rcij,6.0)/2.0);
 
  if(PotentialType == 0 || PotentialType == 1) // L-J and shiftedLJ
  {
    if(r < rcij) dudr = -48*epsilonij/r*(pow(sigmaij/r,12.0)-pow(sigmaij/r,6.0)/2.0);
    else dudr = 0.;
  } //endif L-J
  else if(PotentialType == 2) // shifted-force L-J  1/r^12 - 1/r^6
  {
    if(r < rcij) dudr = -48*epsilonij/r*(pow(sigmaij/r,12.0)-pow(sigmaij/r,6.0)/2.0) - ducdr;
    else dudr = 0.;
  } // endif shifted-force L-J

  return dudr;
}

double Potential_dr2(int i,int j,double r) // d2u(rij)/dr2
{
 double dudr2; // d2u/dr2
  double sigmaij,epsilonij,rcij;
  double uc,ducdr;

  sigmaij = SigmaIJ(i,j);
  rcij = rc*sigmaij;
  epsilonij = EpsilonIJ(i,j);
  //uc = 4.0*epsilonij*(pow(sigmaij/rcij,12.0)-pow(sigmaij/rcij,6.0));
  //ducdr = -48.0*epsilonij/rcij*(pow(sigmaij/rcij,12.0) - pow(sigmaij/rcij,6.0)/2.0);
  
 if(PotentialType == 0 || PotentialType == 1 || PotentialType == 2) // L-J
 {
   if(r < rcij) dudr2 = 4.0*epsilonij/SQR(r)*(12.0*13.0*pow(sigmaij/r,12.0)-6.0*7.0*pow(sigmaij/r,6.0));
   else dudr2 = 0.;
 } //endif L-J

 return dudr2;
}

double SpringPotential(int i,int nb,double r) //
{
 double u; //
 double hij; // spring stiffness
 double qij; // spring length

 double tol;
 tol = 1.0E-10;
 //r = r-tol;
    
 //j = NSpringList[i][nb];
 hij =  SpringStiff[i][nb];
 qij =  SpringLength[i][nb];

// if(r>qij) {printf("%d %d %lf>%lf\n",i,NSpringList[i][nb],r,qij);};
// u = -0.5*hij*SQR(qij)*log(SQR(qij)-SQR(r));
 //u = hij*r/(1.0-SQR(r)/SQR(qij)); 
 u = 0.5*hij*SQR(r-qij); 
 

 return u;
}

double SpringPotential_dr(int i,int nb,double r) //
{
 double dudr; //
 double hij; // spring stiffness
 double qij; // spring length

 double pp,qq;
 
 double tol;
 tol = 1.0E-10;
 //r = r-tol;
    
 //j = NSpringList[i][nb];
 hij =  SpringStiff[i][nb];
 qij =  SpringLength[i][nb];

 //dudr = hij*r/(1.0-SQR(r)/SQR(qij)); 

// qq = SQR(r)/SQR(qij);
// pp = 1.0-qq;
 //dudr = hij/pp + 2.0*hij*qq/SQR(pp); 
 
dudr = hij*(r-qij); 

 return dudr;
}

double BondAnglePotential(int i,int j,int k,double stiff,double equiangle)
{
 double u;
 double dijx,dijy,dijz,dikx,diky,dikz;
 double costheta;
 double theta;

// VECTOR dij,dik;

// dij.x = position[j].x-position[i].x;
// dij.y = position[j].y-position[i].y;
// dij.z = position[j].z-position[i].z;

// dik.x = position[k].x-position[i].x;
// dik.y = position[k].y-position[i].y;
// dik.z = position[k].z-position[i].z;

 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i];
 dikx = RX[k]-RX[i]; diky = RY[k]-RY[i]; dikz = RZ[k]-RZ[i];

 costheta = (dijx*dikx+dijy*diky+dijz*dikz)/sqrt(SQR(dijx)+SQR(dijy)+SQR(dijz))/sqrt(SQR(dikx)+SQR(diky)+SQR(dikz));

 if(fabs(costheta)> 1.0) {printf("cos = %lf\n",costheta); exit(1);}
 theta=acos(costheta);


 u = 0.5*stiff*SQR(theta-equiangle);

 return u;
}

VECTOR BondAngleForceA(int i,int j,int k,double stiff,double equiangle) // angle j-i-k
{
 VECTOR f;
 double dudcos;
 double costheta;
 double theta;
 double dijx,dijy,dijz,dikx,diky,dikz;
 double Cjj,Ckk,Cjk;

 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i];
 dikx = RX[k]-RX[i]; diky = RY[k]-RY[i]; dikz = RZ[k]-RZ[i];

 Cjk = dijx*dikx+dijy*diky+dijz*dikz;
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);
 Ckk = SQR(dikx)+SQR(diky)+SQR(dikz);

 costheta = Cjk/sqrt(Cjj)/sqrt(Ckk);

 if(fabs(costheta)> 1.0) {printf("cos = %lf\n",costheta); exit(1);}
 theta=acos(costheta);

 dudcos = stiff*(theta-equiangle)/sin(theta);

 f.x = -Cjk/Ckk*dikx - Cjk/Cjj*dijx + dikx + dijx;
 f.y = -Cjk/Ckk*diky - Cjk/Cjj*dijy + diky + dijy;
 f.z = -Cjk/Ckk*dikz - Cjk/Cjj*dijz + dikz + dijz;

 f.x = -f.x/sqrt(Cjj*Ckk)*dudcos;
 f.y = -f.y/sqrt(Cjj*Ckk)*dudcos;
 f.z = -f.z/sqrt(Cjj*Ckk)*dudcos;

 return f;
}

VECTOR BondAngleForceB(int i,int j,int k,double stiff,double equiangle) // angle i-j-k
{
 VECTOR f;
 double dudcos;
 double costheta;
 double theta;
 double dijx,dijy,dijz,djkx,djky,djkz;
 double Cjj,Ckk,Cjk;

 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i];
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j];

 Cjk = dijx*djkx+dijy*djky+dijz*djkz;
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);
 Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);

 costheta = Cjk/sqrt(Cjj)/sqrt(Ckk);

 if(fabs(costheta)> 1.0) {printf("cos = %lf\n",costheta); exit(1);}
 theta=acos(costheta);

 if(theta==0.0) theta=1.0;

 dudcos = stiff*(theta-(M_PI-equiangle))/sin(theta);

 f.x = Cjk/Cjj*dijx - djkx;
 f.y = Cjk/Cjj*dijy - djky;
 f.z = Cjk/Cjj*dijz - djkz;

 f.x = f.x/sqrt(Cjj*Ckk)*dudcos;
 f.y = f.y/sqrt(Cjj*Ckk)*dudcos;
 f.z = f.z/sqrt(Cjj*Ckk)*dudcos;

 return f;
}

double TorsionPotential(int nc, int nb) // input backbone C index 0,1,...NBackBones-1 on nc the chain
{
 int i,j,k,l;
 double c1,c2,phase;
 double u;
 double theta;
 double costheta;
 double dijx,dijy,dijz;
 double djkx,djky,djkz;
 double dklx,dkly,dklz;
 //double dijkx,dijky,dijkz;
 double djklx,djkly,djklz;
 double Cjj,Ckk,Cll,Cjk,Ckl,Cjl;
 double Djk,Dkl,Djl;

 if(nb>=1 && nb<=NBackBones-3)
 { 
 i=BackBoneList[nc][nb-1]; 
 j=BackBoneList[nc][nb]; 
 k=BackBoneList[nc][nb+1]; 
 l=BackBoneList[nc][nb+2]; 
 c1 = TorsionStiff1[nb];
 c2 = TorsionStiff2[nb];
 phase = TorsionPhase[nb];

 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i];
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j];
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k];

// dijkx = djky*dijz-djkz*dijy;dijky = -1.0*(djkx*dijz-djkz*dijx);dijkz = djkx*dijy-djky*dijx;
// djklx = dkly*djkz-dklz*djky;djkly = -1.0*(dklx*djkz-dklz*djkx);djklz = dklx*djky-dkly*djkx;
// costheta = -(dijkx*djklx+dijky*djkly+dijkz*djklz)/sqrt(SQR(dijkx)+SQR(dijky)+SQR(dijkz))/sqrt(SQR(djklx)+SQR(djkly)+SQR(djklz));
 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);
 Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);
 Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;
 Ckl = djkx*dklx+djky*dkly+djkz*dklz;
 Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);
 Dkl = Ckk*Cll-SQR(Ckl);
 Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 //u = c1*CUBIC(cos(theta-phase))-c2*cos(theta-phase)+c1-c2;
 u = c1*CUBIC(costheta)-c2*costheta+c1-c2;
 }
 else u=0.0;

 return u;
}

VECTOR TorsionForce(int nc, int nb)
{
 double fx,fy,fz;
 VECTOR fall;
 int i,j,k,l;
 double c1,c2,phase;
 double costheta;
 double dudcos;
 double dijx,dijy,dijz;
 double djkx,djky,djkz;
 double dklx,dkly,dklz;
 double dijkx,dijky,dijkz;
 //double djklx,djkly,djklz;
 double Cjj,Ckk,Cll,Cjk,Ckl,Cjl;
 double Djk,Dkl,Djl;

 fall.x=0.0;fall.y=0.0;fall.z=0.0;
 //  i  -  j  -  k  - l
 // a-3 - a-2 - a-1 - a
 if(nb>=3 && nb<=NBackBones-1)
 {
 i=BackBoneList[nc][nb-3]; 
 j=BackBoneList[nc][nb-2]; 
 k=BackBoneList[nc][nb-1]; 
 l=BackBoneList[nc][nb]; 
 c1 = TorsionStiff1[nb-2];
 c2 = TorsionStiff2[nb-2];
 phase = TorsionPhase[nb-2];

 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; // d(a-2)
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; // d(a-1)
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; // d(a)
 //dijkx = djky*dijz-djkz*dijy;dijky = -1.0*(djkx*dijz-djkz*dijx);dijkz = djkx*dijy-djky*dijx;
 //djklx = dkly*djkz-dklz*djky;djkly = -1.0*(dklx*djkz-dklz*djkx);djklz = dklx*djky-dkly*djkx;
 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);
 Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);
 Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;
 Ckl = djkx*dklx+djky*dkly+djkz*dklz;
 Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);
 Dkl = Ckk*Cll-SQR(Ckl);
 Djl = Ckl*Cjk-Cjl*Ckk;

 //costheta = -(dijkx*djklx+dijky*djkly+dijkz*djklz)/sqrt(SQR(dijkx)+SQR(dijky)+SQR(dijkz))/sqrt(SQR(djklx)+SQR(djkly)+SQR(djklz));
costheta = -Djl/sqrt(Dkl*Djk); 
costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);

 fx = Cjk*djkx - Ckk*dijx - Djl/Dkl*(Ckk*dklx - Ckl*djkx);
 fy = Cjk*djky - Ckk*dijy - Djl/Dkl*(Ckk*dkly - Ckl*djky);
 fz = Cjk*djkz - Ckk*dijz - Djl/Dkl*(Ckk*dklz - Ckl*djkz);
 
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 }
 
 //   i  -  j  -  k  - l
 //  a-2 - a-1 -  a - a+1
 if(nb>=2 && nb<=NBackBones-2)
 {
 i=BackBoneList[nc][nb-2]; 
 j=BackBoneList[nc][nb-1]; 
 k=BackBoneList[nc][nb]; 
 l=BackBoneList[nc][nb+1]; 
 c1 = TorsionStiff1[nb-1];
 c2 = TorsionStiff2[nb-1];
 phase = TorsionPhase[nb-1];

 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; // d(a-1)
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; // d(a)
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; // d(a+1)
// dijkx = djky*dijz-djkz*dijy;dijky = -1.0*(djkx*dijz-djkz*dijx);dijkz = djkx*dijy-djky*dijx;
// djklx = dkly*djkz-dklz*djky;djkly = -1.0*(dklx*djkz-dklz*djkx);djklz = dklx*djky-dkly*djkx;

 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);
 Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);
 Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;
 Ckl = djkx*dklx+djky*dkly+djkz*dklz;
 Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);
 Dkl = Ckk*Cll-SQR(Ckl);
 Djl = Ckl*Cjk-Cjl*Ckk;

 //costheta = -(dijkx*djklx+dijky*djkly+dijkz*djklz)/sqrt(SQR(dijkx)+SQR(dijky)+SQR(dijkz))/sqrt(SQR(djklx)+SQR(djkly)+SQR(djklz));
costheta = -Djl/sqrt(Dkl*Djk); 
costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
//dudcos = 3.0*c1*SQR(costheta)-c2;
dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);

 fx = Cjk*dklx - Cjk*djkx + Ckl*dijx +Ckk*dijx - 2.0*Cjl*djkx - Djl/Djk*(Cjj*djkx-Cjk*dijx) - Djl/Dkl*(Cll*djkx-Ckk*dklx-Ckl*dklx+Ckl*djkx);
 fy = Cjk*dkly - Cjk*djky + Ckl*dijy +Ckk*dijy - 2.0*Cjl*djky - Djl/Djk*(Cjj*djky-Cjk*dijy) - Djl/Dkl*(Cll*djky-Ckk*dkly-Ckl*dkly+Ckl*djky);
 fz = Cjk*dklz - Cjk*djkz + Ckl*dijz +Ckk*dijz - 2.0*Cjl*djkz - Djl/Djk*(Cjj*djkz-Cjk*dijz) - Djl/Dkl*(Cll*djkz-Ckk*dklz-Ckl*dklz+Ckl*djkz);
 
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;

 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 }
 
 //   i  -  j  -  k  - l
 //  a-1 -  a  - a+1 - a+2
 if(nb>=1 && nb<=NBackBones-3)
 {
 i=BackBoneList[nc][nb-1]; 
 j=BackBoneList[nc][nb]; 
 k=BackBoneList[nc][nb+1]; 
 l=BackBoneList[nc][nb+2]; 
 c1 = TorsionStiff1[nb];
 c2 = TorsionStiff2[nb];
 phase = TorsionPhase[nb];

 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; // d(a)
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; // d(a+1)
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; // d(a+2)
 //dijkx = djky*dijz-djkz*dijy;dijky = -1.0*(djkx*dijz-djkz*dijx);dijkz = djkx*dijy-djky*dijx;
 //djklx = dkly*djkz-dklz*djky;djkly = -1.0*(dklx*djkz-dklz*djkx);djklz = dklx*djky-dkly*djkx;
 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);
 Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);
 Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;
 Ckl = djkx*dklx+djky*dkly+djkz*dklz;
 Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);
 Dkl = Ckk*Cll-SQR(Ckl);
 Djl = Ckl*Cjk-Cjl*Ckk;

// costheta = -(dijkx*djklx+dijky*djkly+dijkz*djklz)/sqrt(SQR(dijkx)+SQR(dijky)+SQR(dijkz))/sqrt(SQR(djklx)+SQR(djkly)+SQR(djklz));
costheta = -Djl/sqrt(Dkl*Djk); 
costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
//dudcos = 3.0*c1*SQR(costheta)-c2;
dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);

 fx = -Cjk*dklx + Ckl*djkx - Ckl*dijx -Ckk*dklx + 2.0*Cjl*djkx - Djl/Dkl*(-Cll*djkx+Ckl*dklx) - Djl/Djk*(Ckk*dijx-Cjj*djkx-Cjk*djkx+Cjk*dijx);
 fy = -Cjk*dkly + Ckl*djky - Ckl*dijy -Ckk*dkly + 2.0*Cjl*djky - Djl/Dkl*(-Cll*djky+Ckl*dkly) - Djl/Djk*(Ckk*dijy-Cjj*djky-Cjk*djky+Cjk*dijy);
 fz = -Cjk*dklz + Ckl*djkz - Ckl*dijz -Ckk*dklz + 2.0*Cjl*djkz - Djl/Dkl*(-Cll*djkz+Ckl*dklz) - Djl/Djk*(Ckk*dijz-Cjj*djkz-Cjk*djkz+Cjk*dijz);
 
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;

 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 }

 //  i -  j  -  k  - l
 //  a - a+1 - a+2 - a+3
 if(nb>=0 && nb<=NBackBones-4)
 {
 i=BackBoneList[nc][nb]; 
 j=BackBoneList[nc][nb+1]; 
 k=BackBoneList[nc][nb+2]; 
 l=BackBoneList[nc][nb+3]; 
 c1 = TorsionStiff1[nb+1];
 c2 = TorsionStiff2[nb+1];
 phase = TorsionPhase[nb+1];

 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; // d(a+1)
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; // d(a+2)
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; // d(a+3)
 //dijkx = djky*dijz-djkz*dijy;dijky = -1.0*(djkx*dijz-djkz*dijx);dijkz = djkx*dijy-djky*dijx;
 //djklx = dkly*djkz-dklz*djky;djkly = -1.0*(dklx*djkz-dklz*djkx);djklz = dklx*djky-dkly*djkx;
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);
 Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);
 Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;
 Ckl = djkx*dklx+djky*dkly+djkz*dklz;
 Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);
 Dkl = Ckk*Cll-SQR(Ckl);
 Djl = Ckl*Cjk-Cjl*Ckk;

// costheta = -(dijkx*djklx+dijky*djkly+dijkz*djklz)/sqrt(SQR(dijkx)+SQR(dijky)+SQR(dijkz))/sqrt(SQR(djklx)+SQR(djkly)+SQR(djklz));
costheta = -Djl/sqrt(Dkl*Djk); 
costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
//dudcos = 3.0*c1*SQR(costheta)-c2;
dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);

 fx = -Ckl*djkx + Ckk*dklx - Djl/Djk*(-Ckk*dijx + Cjk*djkx);
 fy = -Ckl*djky + Ckk*dkly - Djl/Djk*(-Ckk*dijy + Cjk*djky);
 fz = -Ckl*djkz + Ckk*dklz - Djl/Djk*(-Ckk*dijz + Cjk*djkz);
 
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 }

 return fall;
}

double TorsionPotentialPhenyl(int nm, int na,int ni) // input monomer index nm and atom index na and overall atom index ni
{
 int i,j,k,l;
 double c1,c2,phase;
 double u;
 double theta;
 double costheta;
 double dijx,dijy,dijz;
 double djkx,djky,djkz;
 double dklx,dkly,dklz;
 //double dijkx,dijky,dijkz;
 double djklx,djkly,djklz;
 double Cjj,Ckk,Cll,Cjk,Ckl,Cjl;
 double Djk,Dkl,Djl;

 u = 0.0;

 if(na==0)
 {
  i=ni; j=ni+1; k=ni+2; l=ni+3; 
  c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
  dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i];
  djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j];
  dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k];
  Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
  Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
  Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
  costheta = -Djl/sqrt(Dkl*Djk); 
  costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
  u += c1*CUBIC(costheta)-c2*costheta+c1-c2; //u = c1*CUBIC(cos(theta-phase))-c2*cos(theta-phase)+c1-c2;
 
  i=ni; j=ni+1; k=ni+2; l=ni+7; 
  c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
  dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i];
  djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j];
  dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k];
  Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
  Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
  Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
  costheta = -Djl/sqrt(Dkl*Djk); 
  costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
  u += c1*CUBIC(costheta)-c2*costheta+c1-c2; //u = c1*CUBIC(cos(theta-phase))-c2*cos(theta-phase)+c1-c2;

  if(nm!=0)
  {
  i=ni; j=ni-7; k=ni-6; l=ni-5; 
  c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
  dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i];
  djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j];
  dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k];
  Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
  Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
  Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
  costheta = -Djl/sqrt(Dkl*Djk); 
  costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
  u += c1*CUBIC(costheta)-c2*costheta+c1-c2; //u = c1*CUBIC(cos(theta-phase))-c2*cos(theta-phase)+c1-c2;
  
  i=ni; j=ni-7; k=ni-6; l=ni-1; 
  c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
  dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i];
  djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j];
  dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k];
  Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
  Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
  Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
  costheta = -Djl/sqrt(Dkl*Djk); 
  costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
  u += c1*CUBIC(costheta)-c2*costheta+c1-c2; //u = c1*CUBIC(cos(theta-phase))-c2*cos(theta-phase)+c1-c2;
  }
 }

 if(na==1)
 {
 if(nm!=NMonomerPerChain-1)
 {

 }

 }
 
 //if(na==2)
 //if(na==3)
 //if(na==7)

 return u;
}


VECTOR TorsionForcePhenyl(int nm,int na,int ni)
{
 double fx,fy,fz;
 VECTOR fall;
 int i,j,k,l;
 double c1,c2,phase;
 double costheta;
 double dudcos;
 double dijx,dijy,dijz;
 double djkx,djky,djkz;
 double dklx,dkly,dklz;
 double dijkx,dijky,dijkz;
 //double djklx,djkly,djklz;
 double Cjj,Ckk,Cll,Cjk,Ckl,Cjl;
 double Djk,Dkl,Djl;

 fall.x=0.0;fall.y=0.0;fall.z=0.0;

 if(na==0)
 {
 i=ni;j=ni+1;k=ni+2;l=ni+3; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; // d(a+1)
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; // d(a+2)
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; // d(a+3)
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = -Ckl*djkx + Ckk*dklx - Djl/Djk*(-Ckk*dijx + Cjk*djkx);
 fy = -Ckl*djky + Ckk*dkly - Djl/Djk*(-Ckk*dijy + Cjk*djky);
 fz = -Ckl*djkz + Ckk*dklz - Djl/Djk*(-Ckk*dijz + Cjk*djkz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
  
 i=ni; j=ni+1; k=ni+2; l=ni+7; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; // d(a+1)
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; // d(a+2)
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; // d(a+3)
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = -Ckl*djkx + Ckk*dklx - Djl/Djk*(-Ckk*dijx + Cjk*djkx);
 fy = -Ckl*djky + Ckk*dkly - Djl/Djk*(-Ckk*dijy + Cjk*djky);
 fz = -Ckl*djkz + Ckk*dklz - Djl/Djk*(-Ckk*dijz + Cjk*djkz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 
 if(nm!=0)
 {
 i=ni; j=ni-7; k=ni-6; l=ni-5; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = -Ckl*djkx + Ckk*dklx - Djl/Djk*(-Ckk*dijx + Cjk*djkx);
 fy = -Ckl*djky + Ckk*dkly - Djl/Djk*(-Ckk*dijy + Cjk*djky);
 fz = -Ckl*djkz + Ckk*dklz - Djl/Djk*(-Ckk*dijz + Cjk*djkz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 
 i=ni; j=ni-7; k=ni-6; l=ni-1; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = -Ckl*djkx + Ckk*dklx - Djl/Djk*(-Ckk*dijx + Cjk*djkx);
 fy = -Ckl*djky + Ckk*dkly - Djl/Djk*(-Ckk*dijy + Cjk*djky);
 fz = -Ckl*djkz + Ckk*dklz - Djl/Djk*(-Ckk*dijz + Cjk*djkz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 }
 }//end if na==0

 if(na==1)
 {
 i=ni-1; j=ni; k=ni+1; l=ni+2; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = -Cjk*dklx + Ckl*djkx - Ckl*dijx -Ckk*dklx + 2.0*Cjl*djkx - Djl/Dkl*(-Cll*djkx+Ckl*dklx) - Djl/Djk*(Ckk*dijx-Cjj*djkx-Cjk*djkx+Cjk*dijx);
 fy = -Cjk*dkly + Ckl*djky - Ckl*dijy -Ckk*dkly + 2.0*Cjl*djky - Djl/Dkl*(-Cll*djky+Ckl*dkly) - Djl/Djk*(Ckk*dijy-Cjj*djky-Cjk*djky+Cjk*dijy);
 fz = -Cjk*dklz + Ckl*djkz - Ckl*dijz -Ckk*dklz + 2.0*Cjl*djkz - Djl/Dkl*(-Cll*djkz+Ckl*dklz) - Djl/Djk*(Ckk*dijz-Cjj*djkz-Cjk*djkz+Cjk*dijz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 
 i=ni-1; j=ni; k=ni+1; l=ni+6; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = -Cjk*dklx + Ckl*djkx - Ckl*dijx -Ckk*dklx + 2.0*Cjl*djkx - Djl/Dkl*(-Cll*djkx+Ckl*dklx) - Djl/Djk*(Ckk*dijx-Cjj*djkx-Cjk*djkx+Cjk*dijx);
 fy = -Cjk*dkly + Ckl*djky - Ckl*dijy -Ckk*dkly + 2.0*Cjl*djky - Djl/Dkl*(-Cll*djky+Ckl*dkly) - Djl/Djk*(Ckk*dijy-Cjj*djky-Cjk*djky+Cjk*dijy);
 fz = -Cjk*dklz + Ckl*djkz - Ckl*dijz -Ckk*dklz + 2.0*Cjl*djkz - Djl/Dkl*(-Cll*djkz+Ckl*dklz) - Djl/Djk*(Ckk*dijz-Cjj*djkz-Cjk*djkz+Cjk*dijz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 
 if(nm!=NMonomerPerChain-1)
 {
 i=ni+7; j=ni; k=ni+1; l=ni+2; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = -Cjk*dklx + Ckl*djkx - Ckl*dijx -Ckk*dklx + 2.0*Cjl*djkx - Djl/Dkl*(-Cll*djkx+Ckl*dklx) - Djl/Djk*(Ckk*dijx-Cjj*djkx-Cjk*djkx+Cjk*dijx);
 fy = -Cjk*dkly + Ckl*djky - Ckl*dijy -Ckk*dkly + 2.0*Cjl*djky - Djl/Dkl*(-Cll*djky+Ckl*dkly) - Djl/Djk*(Ckk*dijy-Cjj*djky-Cjk*djky+Cjk*dijy);
 fz = -Cjk*dklz + Ckl*djkz - Ckl*dijz -Ckk*dklz + 2.0*Cjl*djkz - Djl/Dkl*(-Cll*djkz+Ckl*dklz) - Djl/Djk*(Ckk*dijz-Cjj*djkz-Cjk*djkz+Cjk*dijz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;

 i=ni+7; j=ni; k=ni+1; l=ni+6; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = -Cjk*dklx + Ckl*djkx - Ckl*dijx -Ckk*dklx + 2.0*Cjl*djkx - Djl/Dkl*(-Cll*djkx+Ckl*dklx) - Djl/Djk*(Ckk*dijx-Cjj*djkx-Cjk*djkx+Cjk*dijx);
 fy = -Cjk*dkly + Ckl*djky - Ckl*dijy -Ckk*dkly + 2.0*Cjl*djky - Djl/Dkl*(-Cll*djky+Ckl*dkly) - Djl/Djk*(Ckk*dijy-Cjj*djky-Cjk*djky+Cjk*dijy);
 fz = -Cjk*dklz + Ckl*djkz - Ckl*dijz -Ckk*dklz + 2.0*Cjl*djkz - Djl/Dkl*(-Cll*djkz+Ckl*dklz) - Djl/Djk*(Ckk*dijz-Cjj*djkz-Cjk*djkz+Cjk*dijz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 }// end if not last monomer of the chain
 }// end if na == 1

 if(na==2)
 {
 i=ni-2; j=ni-1; k=ni; l=ni+1; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = Cjk*dklx - Cjk*djkx + Ckl*dijx +Ckk*dijx - 2.0*Cjl*djkx - Djl/Djk*(Cjj*djkx-Cjk*dijx) - Djl/Dkl*(Cll*djkx-Ckk*dklx-Ckl*dklx+Ckl*djkx);
 fy = Cjk*dkly - Cjk*djky + Ckl*dijy +Ckk*dijy - 2.0*Cjl*djky - Djl/Djk*(Cjj*djky-Cjk*dijy) - Djl/Dkl*(Cll*djky-Ckk*dkly-Ckl*dkly+Ckl*djky);
 fz = Cjk*dklz - Cjk*djkz + Ckl*dijz +Ckk*dijz - 2.0*Cjl*djkz - Djl/Djk*(Cjj*djkz-Cjk*dijz) - Djl/Dkl*(Cll*djkz-Ckk*dklz-Ckl*dklz+Ckl*djkz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 
 i=ni-2; j=ni-1; k=ni; l=ni+5; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = Cjk*dklx - Cjk*djkx + Ckl*dijx +Ckk*dijx - 2.0*Cjl*djkx - Djl/Djk*(Cjj*djkx-Cjk*dijx) - Djl/Dkl*(Cll*djkx-Ckk*dklx-Ckl*dklx+Ckl*djkx);
 fy = Cjk*dkly - Cjk*djky + Ckl*dijy +Ckk*dijy - 2.0*Cjl*djky - Djl/Djk*(Cjj*djky-Cjk*dijy) - Djl/Dkl*(Cll*djky-Ckk*dkly-Ckl*dkly+Ckl*djky);
 fz = Cjk*dklz - Cjk*djkz + Ckl*dijz +Ckk*dijz - 2.0*Cjl*djkz - Djl/Djk*(Cjj*djkz-Cjk*dijz) - Djl/Dkl*(Cll*djkz-Ckk*dklz-Ckl*dklz+Ckl*djkz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 
 if(nm!=NMonomerPerChain-1)
 {
 i=ni+6; j=ni-1; k=ni; l=ni+1; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = Cjk*dklx - Cjk*djkx + Ckl*dijx +Ckk*dijx - 2.0*Cjl*djkx - Djl/Djk*(Cjj*djkx-Cjk*dijx) - Djl/Dkl*(Cll*djkx-Ckk*dklx-Ckl*dklx+Ckl*djkx);
 fy = Cjk*dkly - Cjk*djky + Ckl*dijy +Ckk*dijy - 2.0*Cjl*djky - Djl/Djk*(Cjj*djky-Cjk*dijy) - Djl/Dkl*(Cll*djky-Ckk*dkly-Ckl*dkly+Ckl*djky);
 fz = Cjk*dklz - Cjk*djkz + Ckl*dijz +Ckk*dijz - 2.0*Cjl*djkz - Djl/Djk*(Cjj*djkz-Cjk*dijz) - Djl/Dkl*(Cll*djkz-Ckk*dklz-Ckl*dklz+Ckl*djkz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;

 i=ni+6; j=ni-1; k=ni; l=ni+5; 
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = Cjk*dklx - Cjk*djkx + Ckl*dijx +Ckk*dijx - 2.0*Cjl*djkx - Djl/Djk*(Cjj*djkx-Cjk*dijx) - Djl/Dkl*(Cll*djkx-Ckk*dklx-Ckl*dklx+Ckl*djkx);
 fy = Cjk*dkly - Cjk*djky + Ckl*dijy +Ckk*dijy - 2.0*Cjl*djky - Djl/Djk*(Cjj*djky-Cjk*dijy) - Djl/Dkl*(Cll*djky-Ckk*dkly-Ckl*dkly+Ckl*djky);
 fz = Cjk*dklz - Cjk*djkz + Ckl*dijz +Ckk*dijz - 2.0*Cjl*djkz - Djl/Djk*(Cjj*djkz-Cjk*dijz) - Djl/Dkl*(Cll*djkz-Ckk*dklz-Ckl*dklz+Ckl*djkz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 }// end if not last monomer of the chain
 }// end if na == 2
 
 if(na==3 || na==7)
 {
 if(na==3){i=ni-3; j=ni-2; k=ni-1; l=ni;} if(na==7){i=ni-7; j=ni-6; k=ni-5; l=ni;}
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = Cjk*djkx - Ckk*dijx - Djl/Dkl*(Ckk*dklx - Ckl*djkx);
 fy = Cjk*djky - Ckk*dijy - Djl/Dkl*(Ckk*dkly - Ckl*djky);
 fz = Cjk*djkz - Ckk*dijz - Djl/Dkl*(Ckk*dklz - Ckl*djkz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 
 if(nm!=NMonomerPerChain-1)
 {
 if(na==3){i=ni+5; j=ni-2; k=ni-1; l=ni;} if(na==7){i=ni+1; j=ni-6; k=ni-5; l=ni;}
 c1 = TorsionStiff1[0];c2 = TorsionStiff2[0];phase = TorsionPhase[0];
 dijx = RX[j]-RX[i]; dijy = RY[j]-RY[i]; dijz = RZ[j]-RZ[i]; 
 djkx = RX[k]-RX[j]; djky = RY[k]-RY[j]; djkz = RZ[k]-RZ[j]; 
 dklx = RX[l]-RX[k]; dkly = RY[l]-RY[k]; dklz = RZ[l]-RZ[k]; 
 Cjj = SQR(dijx)+SQR(dijy)+SQR(dijz);Ckk = SQR(djkx)+SQR(djky)+SQR(djkz);Cll = SQR(dklx)+SQR(dkly)+SQR(dklz);
 Cjk = dijx*djkx+dijy*djky+dijz*djkz;Ckl = djkx*dklx+djky*dkly+djkz*dklz;Cjl = dijx*dklx+dijy*dkly+dijz*dklz;
 Djk = Cjj*Ckk-SQR(Cjk);Dkl = Ckk*Cll-SQR(Ckl);Djl = Ckl*Cjk-Cjl*Ckk;
 costheta = -Djl/sqrt(Dkl*Djk); 
 costheta = costheta*cos(phase);//+sqrt(1.0-SQR(costheta))*sin(phase);
 dudcos = 3.0*c1*SQR(costheta)*cos(phase)-c2*cos(phase);//dudcos = 3.0*c1*SQR(costheta)-c2;
 fx = Cjk*djkx - Ckk*dijx - Djl/Dkl*(Ckk*dklx - Ckl*djkx);
 fy = Cjk*djky - Ckk*dijy - Djl/Dkl*(Ckk*dkly - Ckl*djky);
 fz = Cjk*djkz - Ckk*dijz - Djl/Dkl*(Ckk*dklz - Ckl*djkz);
 fx = fx/sqrt(Djk*Dkl)*dudcos;
 fy = fy/sqrt(Djk*Dkl)*dudcos;
 fz = fz/sqrt(Djk*Dkl)*dudcos;
 fall.x += fx;
 fall.y += fy;
 fall.z += fz;
 }// end if not last monomer of the chain
 }// end if na == 3 or 7
 
 
 return fall;
}

