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

//gcc -lm -o pfits_getZeroDM pfits_getZeroDM.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lfftw3

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include "fftw3.h"

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,j,k;
  float *data;

  char outFile[128]="zeroDM.dat";
  FILE *fout;
  int subint;
  long nSamples;
  int nTime,nFreq;
  float sum1,sum2;
  long totCount;
  int sub0,sub1;
  int pol=0;
  int histogramTot[255];
  int i0,i1;
  float f0,f1;
  int setFreq=0;
  int sclOff=0;
  int tAv=1;
  int fft=0;
  long n,nc;

  
  fftw_complex *output;
  double mean=0;
  double *in_vals;
  fftw_plan transform_plan;
  
  i0 = -1;
  i1 = -1;

  sub0 = -1;
  sub1 = -1;

  
  for (i=0;i<255;i++)
    histogramTot[i]=0;
 
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outFile,argv[++i]);
      else if (strcasecmp(argv[i],"-tav")==0)
	sscanf(argv[++i],"%d",&tAv);
      else if (strcmp(argv[i],"-s1")==0)
	sscanf(argv[++i],"%d",&sub0);
      else if (strcmp(argv[i],"-s2")==0)
	sscanf(argv[++i],"%d",&sub1);
      else if (strcmp(argv[i],"-pol")==0)
	sscanf(argv[++i],"%d",&pol);
      else if (strcmp(argv[i],"-c1")==0)
	sscanf(argv[++i],"%d",&i0);
      else if (strcmp(argv[i],"-c2")==0)
	sscanf(argv[++i],"%d",&i1);
      else if (strcmp(argv[i],"-scale")==0)
	sclOff=1;
      else if (strcmp(argv[i],"-fft")==0)
	fft=1;
      else if (strcmp(argv[i],"-f0")==0)
	{sscanf(argv[++i],"%f",&f0); setFreq=1;}
      else if (strcmp(argv[i],"-f1")==0)
	{sscanf(argv[++i],"%f",&f1); setFreq=1;}
    }

  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  if (sub0 == -1)
    sub0 = 0;
  if (sub1 == -1)
    sub1 = dSet->head->nsub;
  
  if (setFreq==1)
    {
      i0=i1=-1;
      
      for (i=0;i<dSet->head->nchan-1;i++)
	{
	  if (i0==-1 && dSet->head->chanFreq[i] >= f0)
	    i0 = i;
	  if (i1==-1 && dSet->head->chanFreq[i] >= f1)
	    i1 = i;
	}
      if (i0==-1) i0=0;
      if (i1==-1) i1=dSet->head->nchan;
      printf("Have %d %d\n",i0,i1);
    }
  else
    {
      if (i0==-1) i0=0;
      if (i1==-1) i1=dSet->head->nchan;
    }
  data = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  
  totCount=0;

  fout = fopen(outFile,"w");

  if (fft==1)
    {
      in_vals = (double *)malloc(sizeof(double)*dSet->head->nsblk*dSet->head->nsub);
    }

  nc=0;
  for (subint=sub0;subint<sub1;subint++)
    {
      printf("Subint %d/%d\n",subint,sub1);
      pfits_read1pol_float(data,pol,dSet,subint,subint,1,&nSamples,&nTime,&nFreq,debug,sclOff);
      for (j=0;j<nTime;j+=tAv)
	{
	  sum1=0;
	  n=0;
	  for (k=j;k<j+tAv;k++)
	    {
	      if (k < nTime)
		{
		  for (i=0;i<dSet->head->nchan;i++)
		    {
		      if ((i >= i0 && i <= i1)) 
			{
			  sum1 += data[k*nFreq+i];
			  n++;
			}
		    }
		  totCount++;
		  
		}	      
	    }
	  if (fft==1)
	    {
	      mean+=sum1/(double)n;
	      in_vals[nc]=sum1/(double)n;
	    }
	  fprintf(fout,"%d %g %g\n",totCount,totCount*dSet->head->tsamp,sum1/(double)n);
	  nc++;
	}
    }
  fclose(fout);

  if (fft==1)
    {
      double freq;

      mean/=(double)nc;
      output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nc);
      printf("Allocation = %ld %d %d\n",(long)(dSet->head->nsblk*dSet->head->nsub),dSet->head->nsblk,dSet->head->nsub);
      printf("FFT n = %d\n",nc);
      for (i=0;i<nc;i++)
	in_vals[i]-=mean;
      printf("Removed mean\n");
      transform_plan = fftw_plan_dft_r2c_1d(nc,in_vals,output,FFTW_ESTIMATE);
      printf("Running fft\n");
      fftw_execute(transform_plan);
      printf("Complete running fft\n");
      fout = fopen("fft.dat","w");
      for (i=0;i<nc/2;i++)
	{
	  freq = (double)i/((sub1-sub0)*dSet->head->nsblk*dSet->head->tsamp);
	  fprintf(fout,"%g %g\n",freq,pow(output[i][0],2)+pow(output[i][1],2));
	}
      fclose(fout);
      free(in_vals);
      free(output);
    }
  
  free(data);
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
}
