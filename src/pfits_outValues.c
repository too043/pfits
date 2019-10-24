//  Copyright (C) 2015,2016 George Hobbs
// This file is part of the pfits software package
//

/* pfits is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 * pfits is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 * You should have received a copy of the GNU General Public License 
 * along with pfits.  If not, see <http://www.gnu.org/licenses/>. 
*/

//  gcc -lm -o pfits_outValues pfits_outValues.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio



#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int i,j;
  int nFiles=0,nFiles_t=0;
  int subint;
  int s1,s2;
  float *out;
  float dedisp;
  int nTimeSamples,nFreqSamples;
  float deltaT;
  int delta;
  long nSamples;
  double dm=3490;
  // double dm=778;
  double fref;
  FILE *fout;
  int nc;
  
  // Allocate memory for these files
  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)  setFilename(argv[++i],dSet,debug);
      if (strcmp(argv[i],"-subint")==0) sscanf(argv[++i],"%d",&subint);
    }
  s1 = 105;
    s2 = 118;

    //s1 = 55;
    //  s2 = 59;
  
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);
  out = (float *)malloc(dSet->head->nchan*dSet->head->nsblk*8*(s2-s1+1));

  pfits_read1pol_float(out,0,dSet,s1,s2,2,&nSamples,&nTimeSamples,&nFreqSamples,debug,0);
  printf("Loaded %d %d %d %d (subint = %d)\n",nTimeSamples,nFreqSamples,dSet->head->nchan,dSet->head->nsblk,subint);

  for (i=0;i<96;i++)
    dSet->head->chanFreq[i] = 1374+144-288./96.*i-1.5;

  fref = (dSet->head->chanFreq[0]+dSet->head->chanFreq[dSet->head->nchan-1])/2.0;
  
  fout = fopen("rawData.dat","w");
  dm = 3523; //3523;
  printf("chanFreq = %g %g %g\n",dSet->head->chanFreq[0],dSet->head->chanFreq[dSet->head->nchan-1],fref);
    for (j=2560;j<2560+53248;j++)
  //  for (j=2048;j<2048+12288;j++)
    {
      dedisp=0;
      nc=0;
      for (i=0;i<96;i++)
	{	  
	  deltaT = 4.149e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[i]/1000.0,-2));
	  delta = -(int)(deltaT/dSet->head->tsamp);
	  if (j==2560)
	    printf("roll: %g %d %g\n",deltaT,delta,dSet->head->chanFreq[i]);
	  fprintf(fout,"%d %d %g\n",j,96-i,1-out[j*96+i]);
	  if ((delta+j) >= 0 && (delta+j) < nTimeSamples)
	    {
	      dedisp+=(1-out[(delta+j)*96+i]);
	      nc++;
	    }
	}
      fprintf(fout,"\n");
      printf("dedisp %d %d %g\n",j,nc,dedisp/(double)nc);
    }
  fclose(fout);

  printf("Complete\n");
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
  free(out);
  //  free(dSet);
}

