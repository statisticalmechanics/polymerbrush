/***************************************
 * write movie snapshot in .xyz format
 *
 * number of particles
 * blank
 * name x y z
 * C 0.0000 0.0000 0.0000
 * S 0.0000 1.0000 0.2345
 *
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include "system.h"

void Writemovie(FILE *FilePtr)
{
  int i;
  
  fprintf(FilePtr,"%d\n",NumberOfParticles);
  fprintf(FilePtr,"%lf %lf %lf %lf %lf\n",LX,LY,LZ,Tinstant,Pinstant);

  for(i=0;i<NumberOfParticles;i++)
  {
    if(identity[i] == 1)	
     fprintf(FilePtr,"%s\t","C");
    else if(identity[i] == 2)	
     fprintf(FilePtr,"%s\t","N");
    else if(identity[i] == 3)	
     fprintf(FilePtr,"%s\t","O");
    else if(identity[i] == 4)	
     fprintf(FilePtr,"%s\t","S");
     
    fprintf(FilePtr,"%lf\t",position[i].x);
    fprintf(FilePtr,"%lf\t",position[i].y);
    fprintf(FilePtr,"%lf\n",position[i].z);
    
  // unboxed coordinates
  //  fprintf(FilePtr,"%lf\t",RX[i]);
  //  fprintf(FilePtr,"%lf\t",RY[i]);
  //  fprintf(FilePtr,"%lf\n",RZ[i]);
  }

}
