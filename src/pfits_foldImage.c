//
// Original
// 1181.813u 5.779s 19:47.78 99.9%	0+0k 520+6008000io 0pf+0w
//
// Not calculating phase properly
// 1181.813u 5.779s 14:04.78 99.9%	0+0k 520+6008000io 0pf+0w

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
// /usr/bin/gcc -lm -O5 -o pfits_foldImage pfits_foldImage.c pfits_loader.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"

#define MAX_PERIOD 1000

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  char fname[1024];
  int debug=0;
  float tsub;
  float requestTimeToProces;
  int i,j,k,ibin,p;
  int nbin;
  double *period;
  float **fold;
  unsigned long int **counter;
  int nchan;
  FILE *fout;
  long nS;
  int nT,nF;
  float *data;
  float time;
  int nPeriod=1000; //000;
  double tspan = 120;
  double mean[4096];
  unsigned long int meanCounter[4096];



  //  nbin = 1024;  
  period = (double *)malloc(sizeof(double)*MAX_PERIOD);
  //  for (i=0;i<nPeriod;i++)
  //     {
       //       period[i] = tspan/(i+1)/nbin;
       //       period[i] = 1.2930485876655;
  //       period[i] = 1+(double)i/1000.0; // 1.2930485876655;
       //       period[i] = 0.455021799+1e-4; //179960494;
  //     }
   //  period[0] = 1.;
   //  period[1] = 1.0001;
  //  period[0] = 1.292;
  //  period[i] = 1.2930485876655;

  
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);            
    }
  
  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));
  
  initialise(&dSet,debug);
  setFilename(fname,dSet,debug);
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  tsub = dSet->head->nsblk*dSet->head->tsamp;
  nchan = dSet->head->nchan;
  data = (float *)malloc(sizeof(float)*nchan*dSet->head->nsblk);
  printf("Subintegration time is %g s, nchan = %d\n",tsub,nchan);
  printf("Number of subtegrations = %d\n",dSet->head->nsub);
  fold = (float **)malloc(sizeof(float *)*nPeriod);
  counter = (unsigned long int **)malloc(sizeof(unsigned long int *)*nPeriod);
  for (i=0;i<nPeriod;i++)
     {
       fold[i] = (float *)calloc(nchan*nbin,sizeof(float));
       counter[i] = (unsigned long int *)calloc(nchan*nbin,sizeof(unsigned long int));
    }
  for (i=0;i<nchan;i++)
    {
      mean[i]=0.0;
      meanCounter[i]=0;
    }
  for (i=0;i<dSet->head->nsub;i++)
    {
      printf("Processing %d/%d\n",i,dSet->head->nsub);
      pfits_read1pol_float(data,0,dSet,i,i,1,&nS,&nT,&nF,debug,0);
      for (j=0;j<nT;j++)
	{
	  time = i*tsub + j*dSet->head->tsamp;
	  for (p=0;p<nPeriod;p++)
	    {
	      ibin = (int)(((time/period[p]) - (int)(time/period[p]))*nbin);
	      //	      if (p==1)
	      //		printf("%d %d %d %g %d %d bin\n",i*nT+j,p,ibin,period[p],nF,nT);
	      for (k=0;k<nF;k++)
		{
		  // Should also have a bin counter
		  fold[p][ibin*nchan+k] += data[j*nchan+k];
		  /*
		  counter[p][ibin*nchan+k] ++;
		  if (p==0) // Should check this!
		    {
		      mean[k] += data[j*nchan+k];
		      meanCounter[k]++;
		    }
		  */
		  //		  printf("At %d %d %g\n",ibin,ibin*nchan+k,fold[p][ibin*nchan+k]);
		}
	    }
	}	
    }

  printf("Finished loading data\n");
  for (p=0;p<nPeriod;p++)
    {
      sprintf(fname,"fold.txt.%f",period[p]);
      fout = fopen(fname,"w");
      for (i=0;i<nbin;i++)
	{
	  // Remove mean from each channel (should remove a baseline)
	  
	  for (j=0;j<nchan;j++)
	    {
	      fprintf(fout,"%d %d %g %d %d\n",i,j,fold[p][i*nchan+j]/(double)counter[p][i*nchan+j]-mean[j]/meanCounter[j],counter[p][i*nchan+j],meanCounter[j]);
	    }
	  fprintf(fout,"\n");
	}
      fclose(fout);
    }
  for (i=0;i<nPeriod;i++)
    free(fold[i]);
  free(fold);
  free(data);
  
  deallocateMemory(&dSet,debug);

}
