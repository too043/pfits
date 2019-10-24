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

// gcc -lm -o pfits_extractImage pfits_extractImage.c pfits_setup.c pfits_loader.c T2toolkit.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lfftw3f

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include "T2toolkit.h"
#include <fftw3.h>

#define MAX_FRBS 1000
#define STRLEN 1024

typedef struct frbStruct {
  int   burstTimeSample;
  float burstDM;
  float burstWidth;
  float burstPercent;
} frbStruct;

int convertTo1Bits(unsigned char *inArray, unsigned char *outArray, unsigned int n);

int main(int argc,char *argv[])
{
  frbStruct *frb;
  dSetStruct *dSet;
  FILE *fin,*fout;
  int debug=0;
  int status=0;
  int i,ii,j,k,a;
  int setOut=0;
  int colnum_out;
  int colnum_datoffs;
  int colnum_datscl;
  char outFile[1024],oname[1024]; // Output filename
  fitsfile *outfptr;
  FILE *finBurst;
  float burstTime,burstDM;
  long sub;
  int bpos;
  int nTimeSamples;
  long nSamples;
  int nFrequencySamples;
  float *loadData;
  char *out;
  char *outputVals;
  float fch1,fch2;
  float flow,fhigh;
  float freq,dispersionDelay;
  double bt;
  long int pos;
  double t0,t1;
  unsigned int burstSamples;
  unsigned int changeVals[4][4];
  int ele;
  int totC;
  long randVal;
  float rnum;
  long int seed = TKsetSeed();
  int nFRB=0;
  int process=0;
  float tmax,tmin;
  long int ipos;
  int subint=0;
  int setOn = 0;
  int simRed=0;
  char frbFile[1024];
  int useFile=0;
  float *acf;
  float *redNoise;
  float minRed,maxRed;
  float modulationIndex = 1;
  float checkPercent;
  float acfScale = 0.2;
    // Initialise everything
  initialise(&dSet,debug);

  if (!(frb = (frbStruct *)malloc(sizeof(frbStruct)*MAX_FRBS)))
    {
      printf("Unable to allocate enough memory for %d FRBs\n",MAX_FRBS);
      exit(1);
    }

  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-nsub")==0)
	sscanf(argv[++i],"%d",&subint);
      else if (strcmp(argv[i],"-setOn")==0)
	sscanf(argv[++i],"%d",&setOn);
      else if (strcmp(argv[i],"-simRed")==0)
	{
	  sscanf(argv[++i],"%f",&acfScale);
	  sscanf(argv[++i],"%f",&modulationIndex);
	  simRed=1;
	}
      else if (strcmp(argv[i],"-frbFile")==0)
	{
	  strcpy(argv[++i],frbFile);
	  useFile=1;	  
	}
      else if (strcmp(argv[i],"-frb")==0)
	{
	  sscanf(argv[++i],"%d",&(frb[nFRB].burstTimeSample));
	  sscanf(argv[++i],"%f",&(frb[nFRB].burstDM));
	  sscanf(argv[++i],"%f",&(frb[nFRB].burstWidth));
	  sscanf(argv[++i],"%f",&(frb[nFRB].burstPercent));
	  nFRB++;
	}
      else if (strcmp(argv[i],"-o")==0)
	{
	  strcpy(outFile,argv[++i]);
	  setOut=1;
	}
    }

  if (useFile==1)
    {
      nFRB=0;
      fin = fopen(frbFile,"r");
      while (!feof(fin))
	{
	  if (fscanf(fin,"%d %f %f %f",&(frb[nFRB].burstTimeSample),&(frb[nFRB].burstDM),&(frb[nFRB].burstWidth),&(frb[nFRB].burstPercent))==4)
	    nFRB++;
	}
      fclose(fin);
    }

  
  if (setOut==0)
    errorStop("Must provide an output filename using -o\n",dSet,debug);

  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  fch1 = dSet->head->chanFreq[0];
  fch2 = dSet->head->chanFreq[dSet->head->nchan-1];
  if (fch1 < fch2)
    {
      flow = fch1;
      fhigh = fch2;
    }
  else
    {
      flow = fch2;
      fhigh = fch1;      
    }
  
  loadData = (float *)malloc(sizeof(float)*dSet->head->nchan*dSet->head->nsblk);
  outputVals = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nchan*dSet->head->nsblk);
  out = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nchan*dSet->head->nsblk);

  // Generate red noise model
  if (simRed==1)
    {
      fftwf_plan plan;
      fftwf_complex *spectrum;
      float acf;

      redNoise = (float *)malloc(sizeof(float)*dSet->head->nchan);
      
      //      model = setupRedNoiseModel(flow,fhigh,dSet->head->nchan,100,1e-6,-8./3.0,0);
      spectrum = (fftwf_complex*) fftwf_malloc((dSet->head->nchan)*sizeof(fftwf_complex));

      for (i=0;i<dSet->head->nchan;i++)
	{
	  acf = exp(-(i+1)*acfScale);
	  spectrum[i][0] = acf*TKgaussDev(&seed);
	  spectrum[i][1] = acf*TKgaussDev(&seed);
	}
      plan=fftwf_plan_dft_c2r_1d(dSet->head->nchan,spectrum,redNoise,FFTW_ESTIMATE);
      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
      
      for (i=0;i<dSet->head->nchan;i++)
	{	  
	  if (i==0)
	    minRed=maxRed=redNoise[i];
	  else
	    {
	      if (minRed > redNoise[i]) minRed = redNoise[i];
	      if (maxRed < redNoise[i]) maxRed = redNoise[i];	      
	    }
	  	      
	}
      free(spectrum);
    }
  

  printf("nbits = %d\n",dSet->head->nbits);
  printf("subint = %d\n",subint);
  pfits_read1pol_float(loadData,0,dSet,(float)subint,(float)subint,1,&nSamples,&nTimeSamples,&nFrequencySamples,debug);
  printf("Have loaded %d %d samples \n",nTimeSamples,nFrequencySamples);
  fout = fopen(outFile,"w");
  for (k=0;k<nFrequencySamples;k++)
    {     
      for (j=0;j<nTimeSamples;j++)
	{
	  if (nFRB > 0)
	    {
	      for (i=0;i<nFRB;i++)
		{
		  freq = (fch1 + k*dSet->head->chanbw)/1000.0;
		  dispersionDelay = 4.15e-3*frb[i].burstDM*(1.0/(pow(freq,2))-1.0/(pow(fch1/1000,2))); // Note different units for freq and fcentral
		  bt = frb[i].burstTimeSample*dSet->head->tsamp;
		  ipos = (long int)(((bt+dispersionDelay))/(double)dSet->head->tsamp);
		  //		  printf("Adding at ipos %d\n",(int)ipos);
		  if (j > ipos - frb[i].burstWidth && j < ipos + frb[i].burstWidth)
		    {

		      if (simRed==1)
			{
			  printf("red noise: %g %g %g %g %g\n",redNoise[k],((redNoise[k]-minRed)/(maxRed-minRed))*modulationIndex,minRed,maxRed,acfScale);
			  checkPercent = (frb[i].burstPercent/100.0+((redNoise[k]+minRed)/(maxRed-minRed))*modulationIndex);
			  if (checkPercent > 1) checkPercent = 1;
			  if (checkPercent < 0) checkPercent = 0;
			}
		      else
			checkPercent = frb[i].burstPercent/100.0;
		      
		      rnum = TKranDev(&seed);
		      if (rnum < checkPercent)
			loadData[j*nFrequencySamples+k] = setOn;
		    }
		}
	    }
	  if (j < nTimeSamples-1)
	    fprintf(fout,"%d,",(int)loadData[j*nFrequencySamples+k]);
	  else
	    fprintf(fout,"%d\n",(int)loadData[j*nFrequencySamples+k]);	      
	}
    }
  fclose(fout);
  

  
    // De-allocate the memory
  deallocateMemory(&dSet,debug);
  free(loadData);
  free(outputVals);
  free(out);
  free(frb);
  if (simRed==1)
    free(redNoise);
}



int convertTo1Bits(unsigned char *inArray, unsigned char *outArray, unsigned int n)
{

  int i,j;
  long nv=0;
  unsigned char tc;
  double bit_level=0;

  for (i=0;i<n/8;i++)
    {
      tc=0;
      for (j=0;j<8;j++)
	{
	  if (inArray[nv] > bit_level)
	    tc = tc | (1 << (7-j));
	  nv++;
	}
      outArray[i] = tc;
    }
  
  
  return 0;
}
