//
// Want to record the median and MAD for each frequency channel for each block
// Should replace bad samples with something drawn from the distribution given by the median and MAD
// Should store how many samples are being replaced
// CHECK SUPER BRIGHT FRB (or Vela pulsar)


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
//  gcc -lm -o pfits_singlePulseSearch pfits_singlePulseSearch.c pfits_setup.c pfits_loader.c T2toolkit.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -lpgplot -lfftw3 -O4

//gcc -lm -o pfits_plot pfits_plot.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include "fftw3.h"
#include "T2toolkit.h"


typedef struct candidateStruct {
  float t0;
  float dm;
  float sigma;
} candidateStruct;

void dedisperseAndSum(dSetStruct* dset,float *data,float dm,float fref,float *dedisp,int nTimeSamples,int nFreqSamples,int *nFinalTimeSamples,long int sampleCheck,char *wt);
void findZeroDMEvents(float *dedisp,int n,char *wt,int nFreq);
int compare(const void * a, const void * b);
void calculateMedianMad(dSetStruct *dSet,float *dataScls,float *dataSclsFlagged,int nTimeSamples,
			int nFreqSamples, int nsblk,float *medianChan,float *madChan);
void flagSigma(dSetStruct *dSet,float *dataSclsFlagged,int nTimeSamples,int nFreqSamples,int nsblk,float *medianChan,float *madChan);
void findCandidates(dSetStruct *dSet,float *dedisp,int nTimeSamples,float dm,int nsblk,candidateStruct *candidate,int *nCandidate);

#define MAX_CANDIDATES 10000

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,ii,j;
  char outFile[1024],oname[1024]; // Output filename
  int setOut=0;
  float dm=0;
  float fref=-1;
  int offScl = 1;
  long nSamples;
  double freq;
  int nTimeSamples,nFreqSamples,nFinalTimeSamples;
  long int sampleCheck=-1;
  float *dataNoScls,*dataScls,*dedisp,*bandpass,*medianChan,*madChan,*dataSclsFlagged;
  candidateStruct *candidate;
  int nCandidate=0;
  char *wt;
  FILE *fout;
  int prestoRun=1;
  
  candidate = (candidateStruct *)malloc(sizeof(candidateStruct)*MAX_CANDIDATES);
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-dm")==0)
	sscanf(argv[++i],"%f",&dm);
      else if (strcmp(argv[i],"-fref")==0)
	sscanf(argv[++i],"%f",&fref);
      else if (strcmp(argv[i],"-sampleCheck")==0)
	sscanf(argv[++i],"%ld",&sampleCheck);
      else if (strcasecmp(argv[i],"-noscales")==0)
	offScl=0;
      else if (strcmp(argv[i],"-o")==0)
	{
	  strcpy(outFile,argv[++i]);
	  setOut=1;
	}
    }

  if (setOut==0)
    errorStop("Must provide an output filename using -o\n",dSet,debug);

  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);


  // Allocate memory
  if (!(dataScls = (float *)malloc(sizeof(float)*8192*40*2048))) // CHECK SIZE PROPERLY
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }
  if (!(dataSclsFlagged = (float *)malloc(sizeof(float)*8192*40*2048))) // CHECK SIZE PROPERLY
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }
  if (!(dataNoScls = (float *)malloc(sizeof(float)*8192*40*2048))) // CHECK SIZE PROPERLY
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }
  wt = (char *)malloc(sizeof(char)*8192*40*2048);
  
  dedisp = (float *)malloc(sizeof(float)*8192*40);
  bandpass = (float *)malloc(sizeof(float)*2048); // FIX NCHAN

  medianChan = (float *)malloc(sizeof(float)*40*2048); // FIX NSUB AND NCHAN
  madChan = (float *)malloc(sizeof(float)*40*2048); // FIX NSUB AND NCHAN

  
  // Set the weights
  for (i=0;i<8192*40*2048;i++)
    wt[i] = 1;

  // Load in all the data

  if (prestoRun == 1)
    pfits_read1pol_float(dataNoScls,0,dSet,0,20,2,&nSamples,&nTimeSamples,&nFreqSamples,0,0);
  else
    {
      //  printf("Completed load without scales\n");
      printf("Loading with scales\n");
      pfits_read1pol_float(dataScls,0,dSet,0,20,2,&nSamples,&nTimeSamples,&nFreqSamples,0,1);
    }
  printf("Complete load\n");
  printf("Have loaded %d %d\n",nTimeSamples,nFreqSamples);


  if (prestoRun==0)
    {
      // Determine the median and MAD values per channel per subint
      // and flag based on these values
      calculateMedianMad(dSet,dataScls,dataSclsFlagged,nTimeSamples,nFreqSamples,8192,medianChan,madChan);
      
      // Now calculate the zeroDM timeseries
      // Now do a search for events that we would have picked up, but flag them instead
      flagSigma(dSet,dataSclsFlagged,nTimeSamples,nFreqSamples,8192,medianChan,madChan);
      
      
      printf("De-dispersing at zero DM\n");
      
      dm = 0;
      fref = 960;
      dedisperseAndSum(dSet,dataSclsFlagged,dm,fref,dedisp,nTimeSamples,nFreqSamples,&nFinalTimeSamples,sampleCheck,wt);
    }      
  else if (prestoRun==1)
    {
      // Should do some RFI flagging

      // Zero DM search
      dm = 0;
      fref = 960;
      for (dm = 0; dm < 520; dm+=10)
	{
	  dedisperseAndSum(dSet,dataNoScls,dm,fref,dedisp,nTimeSamples,nFreqSamples,&nFinalTimeSamples,sampleCheck,wt);
	  /*
	    fout = fopen("dedisp500_presto.dat","w");
	    for (i=0;i<nFinalTimeSamples;i++)
	    {
	    fprintf(fout,"%d %g\n",i,dedisp[i]);
	    }
	    fclose(fout);
	  */
	// This is not the standard PRESTO way of getting candidates  
	findCandidates(dSet,dedisp,nTimeSamples,dm,8192,candidate,&nCandidate);      
    }
    }
  printf("Number of candidates = %d\n",nCandidate);

  fout = fopen("candidates.txt","w");
  for (i=0;i<nCandidate;i++)
    {
      fprintf(fout,"%d %g %g %g\n",i,candidate[i].t0,candidate[i].dm,candidate[i].sigma);
    }
  fclose(fout);
  
  // Remove RFI
  


  /*

  // Bad channels
  // Look at the bandpass data set to identify frequency channels to remove
  fout = fopen("bandpass.dat","w");
  for (i=0;i<nFreqSamples;i++)
    {
      bandpass[i]=0;
      for (j=0;j<nTimeSamples;j++)
	bandpass[i] += dataScls[j*nFreqSamples+i];
      bandpass[i]/=(double)nTimeSamples;

      // Use known bad frequencies
      freq = dSet->head->chanFreq[i];
      if (freq > 1017 && freq < 1019) // Parkes airport DME
	wt[i] = 0;
      else if (freq > 1081 && freq < 1091) // Aliased signal
	wt[i] = 0;
      else if (freq > 1023 && freq < 1025) // Digitiser
	wt[i] = 0;
      fprintf(fout,"%d %g %g %d\n",i,freq,bandpass[i],wt[i]);
    }
  fclose(fout);
  */

  // Bad time samples
  // Deal with bandpass variations
  
  
  // Dedisperse


  /*
  printf("Finding zeroDM events\n");
  findZeroDMEvents(dedisp,nFinalTimeSamples,wt,nFreqSamples);
  */
  // PROBABLY SHOULD DIVIDE TIME SERIES BY MEDIAN VALUES POSSIBLY BOTH FREQUENCY CHANNELS AND ALSO NSBLKS

  /*
  printf("De-dispersing after flagging\n");
  dm = 0;
  //  dedisperseAndSum(dSet,dataNoScls,dm,fref,dedisp,nTimeSamples,nFreqSamples,&nFinalTimeSamples,sampleCheck,wt);
  dedisperseAndSum(dSet,dataScls,dm,fref,dedisp,nTimeSamples,nFreqSamples,&nFinalTimeSamples,sampleCheck,wt);

  //  dedisperseAndSum(dSet,dataNoScls,dm,fref,dedisp,nTimeSamples,nFreqSamples,&nFinalTimeSamples,sampleCheck,wt);
  */

  if (prestoRun==0)
    {
      printf("Complete dedispersion at %f\n",dm);
      fout = fopen("dedisp0.dat","w");
      for (i=0;i<nFinalTimeSamples;i++)
	{
	  /*      if (wt[i*nFreqSamples] == 0)
		  fprintf(fout,"\n");
		  else */
	  fprintf(fout,"%d %g\n",i,dedisp[i]);
	}
      fclose(fout);
      
      dm = 0;
      fref = 960;
      printf("De-dispersing at DM = 500cm-3pc\n");
      dedisperseAndSum(dSet,dataSclsFlagged,dm,fref,dedisp,nTimeSamples,nFreqSamples,&nFinalTimeSamples,sampleCheck,wt);
      fout = fopen("dedisp500.dat","w");
      for (i=0;i<nFinalTimeSamples;i++)
	{
	  /*      if (wt[i*nFreqSamples] == 0)
		  fprintf(fout,"\n");
		  else */
	  fprintf(fout,"%d %g\n",i,dedisp[i]);
	}
      fclose(fout);
      // Find candidates
      findCandidates(dSet,dedisp,nTimeSamples,dm,8192,candidate,&nCandidate);      
    }
  // De-allocate the memory
  deallocateMemory(&dSet,debug);

  free(dataScls);
  free(dataSclsFlagged);
  free(dataNoScls);
  free(wt);
  free(dedisp);
  free(bandpass);
  free(medianChan);
  free(madChan);
  free(candidate);
}

void findCandidates(dSetStruct *dSet,float *dedisp,int nTimeSamples,float dm,int nsblk,candidateStruct *candidate,int *nCandidate)
{
  int len = 8;
  int i,j,k,l,ii;
  float sum;
  float x,x2;
  float mean,sdev;
  FILE *fout;
  float sigma = 10;
  
  //  fout = fopen("candidates.txt","w");
  
  ii=0;
  for (i=0;i<nTimeSamples;i+=nsblk)
    {
      x = x2 = 0;
      for (j=i;j<i+nsblk;j++)
	{
	  x += dedisp[j];
	  x2 += pow(dedisp[j],2);
	}
      mean = x/(double)nsblk;
      sdev = sqrt(x2/(double)nsblk - pow(mean,2));
      printf("[%d] mean and sdev = %g/%g\n",ii,mean,sdev);
      for (j=i;j<i+nsblk-len;j+=len)
	{
	  sum=0;
	  for (k=j;k<j+len;k++)
	    sum+=(dedisp[k]-mean);
	  if (sum/sdev > sigma)
	    {
	      candidate[*nCandidate].dm = dm;
	      candidate[*nCandidate].t0 = (k+len/2.)*64e-6; // FIX HARDCODE
	      candidate[*nCandidate].sigma = sum/sdev;
	      (*nCandidate)++;
	    }
	  //	  fprintf(fout,"%d %g %g\n",k+len/2,sum,sdev);
	}
      ii++;
    }
  //  fclose(fout);
   

}

void flagSigma(dSetStruct *dSet,float *dataSclsFlagged,int nTimeSamples,int nFreqSamples,int nsblk,float *medianChan,float *madChan)
{
  int i,j,k,ii;
  float zerodm[nsblk];
  FILE *fout;
  double sx,sx2;
  long seed = TKsetSeed();
  double mean,sdev;
  double kval=3;
  
  ii=0;
  for (i=0;i<nTimeSamples;i+=nsblk)
    {
      printf(" .... form zero dm and flag subint %d\n",ii);
      sx=sx2=0;
      for (j=0;j<nsblk;j++)
	{
	  zerodm[j] = 0;
	  for (k=0;k<nFreqSamples;k++)
	    zerodm[j]+=dataSclsFlagged[(i+j)*nFreqSamples+k];
	  zerodm[j]/=(double)nFreqSamples;
	  sx += zerodm[j];
	  sx2 += pow(zerodm[j],2);
	  //	  vals[j] = dataScls[(i+j)*nFreqSamples+k];
	}
      mean = sx/(double)nsblk;
      sdev = sqrt(sx2/(double)nsblk - pow(mean,2));
      printf("mean and sdev = %g/%g\n",mean,sdev);
      for (j=0;j<nsblk;j++)
	{
	  if (zerodm[j] > kval * sdev + mean)
	    {
	      for (k=0;k<nFreqSamples;k++)
		dataSclsFlagged[(i+j)*nFreqSamples+k] =  (medianChan[ii*nFreqSamples+k] + TKgaussDev(&seed)*1.4826*madChan[ii*nFreqSamples+k])/medianChan[ii*nFreqSamples+k];
	    }
	}
      ii++;
    }
}

void calculateMedianMad(dSetStruct *dSet,float *dataScls,float *dataSclsFlagged,int nTimeSamples,
			int nFreqSamples, int nsblk,float *medianChan,float *madChan)
{
  int i,j,k,ii,i0,j0;
  float vals[nsblk];
  float deviation[nsblk];
  float ds;
  FILE *fout;
  double kval = 1.4826;
  long seed = TKsetSeed();
  int flagChannel=0;
  float freq;
  
  ii=0;
  printf("Calculating median and MAD values\n");
  //  fout = fopen("madData.txt","w");
  for (i=0;i<nTimeSamples;i+=nsblk)
    {
      printf(" .... subint %d\n",ii);
      for (k=0;k<nFreqSamples;k++)
	{
	  for (j=0;j<nsblk;j++)
	    vals[j] = dataScls[(i+j)*nFreqSamples+k];
	  qsort(vals,nsblk,sizeof(float),compare);
	  medianChan[ii*nFreqSamples+k] = vals[nsblk/2];
	  for (j=0;j<nsblk;j++)
	    deviation[j] = fabs(dataScls[(i+j)*nFreqSamples+k]-medianChan[ii*nFreqSamples+k]);
	  qsort(deviation,nsblk,sizeof(float),compare);
	  madChan[ii*nFreqSamples+k] = deviation[nsblk/2];
	  //	  fprintf(fout,"%d %d %g %g\n",ii,k,medianChan[ii*nFreqSamples+k],madChan[ii*nFreqSamples+k]);
	}
      //      fprintf(fout,"\n");

      // Choose to flag
      j0 = ii*nFreqSamples;

      for (j=0;j<nsblk;j++)
	{
	  i0 = (i+j)*nFreqSamples;
	  for (k=0;k<nFreqSamples;k++)
	    {	      
	      flagChannel=0;

	      freq = dSet->head->chanFreq[k];
	      if (freq > 1017 && freq < 1019) // Parkes airport DME
		flagChannel = 1;
	      else if (freq > 1081 && freq < 1091) // Aliased signal
		flagChannel = 1;
	      else if (freq > 1023 && freq < 1025) // Digitiser
		flagChannel = 1;
	      //	      if (k < 1024)
	      //	      if (k < 1024 && k > 512)
	      //	      if (k < 1024 && k > 900)
	      //	      if (k > 1024) // NO SIGNAL
	      //	      if (k > 800) // No signal
	      //	      if (k>512)
	      //	      if (k>2000) LOTS OF SIGNAL
	      //	      if (k > 1500) // WEAK
	      //	      if (k > 1700) // Moderately strong
	      //	      if (k > 1300)
	      //		flagChannel = 1;
	      
	      if (flagChannel == 1 || (dataScls[i0+k] > medianChan[j0+k] + 3*kval*madChan[j0+k])) // https://en.wikipedia.org/wiki/Median_absolute_deviation
		{
		  dataSclsFlagged[i0+k] = ((TKgaussDev(&seed)*kval*madChan[j0+k]/medianChan[j0+k])+1.0);
		}
	      else
		{
		  dataSclsFlagged[i0+k] = (dataScls[i0+k])/medianChan[j0+k];
		}
	    }
	}
      
      ii++;
    }
  //  fclose(fout);
  printf("Completed calculated median and MAD values\n");


}

void findZeroDMEvents(float *dedisp,int n,char *wt,int nFreq)
{
  double x,x2,mean,sdev;
  int i,j,k;
  int nsblk = 8192;
  float vals[nsblk];
  float deviation[nsblk];
  double median,mad;
  double kval = 1.4826;
  
  // Process in blocks of NSBLK
  for (i=0;i<n;i+=nsblk)
    {
      for (j=0;j<nsblk;j++)
	vals[j] = dedisp[i+j];
      printf("vals[0] = %g\n",vals[0]);
      qsort(vals,nsblk,sizeof(float),compare);
      median = vals[nsblk/2];
      for (j=0;j<nsblk;j++)
	deviation[j] = fabs(dedisp[i+j]-median);
      qsort(deviation,nsblk,sizeof(float),compare);
      mad = deviation[nsblk/2];
      for (j=0;j<nsblk;j++)
	{
	  if (dedisp[i+j] > median + 3*kval*mad) // https://en.wikipedia.org/wiki/Median_absolute_deviation
	    {
	      for (k=0;k<nFreq;k++)		
		wt[(i+j)*nFreq+k] = 0;
	    }
	}
      printf("vals[0] now = %g\n",vals[0]);
      printf("median = %g mad = %g\n",vals[nsblk/2],mad);
    }
	/*
	x = x2 = 0.0;
  printf("Processing %d time samples\n",n);
  for (i=0;i<n;i++)
    {
      x+=dedisp[i];
      x2+=pow(dedisp[i],2);
    }
  mean = x/(double)n;
  sdev = sqrt(x2/(double)n - pow(mean,2));
  printf("mean/sdev = %g/%g\n",mean,sdev); */
}

void dedisperseAndSum(dSetStruct* dSet,float *data,float dm,float fref,float *dedisp,int nTimeSamples,int nFreqSamples,int *nFinalTimeSamples,long int sampleCheck,char *wt)
{
  unsigned long int i,j,nc,i0;
  float tdiff;
  int idiff[nFreqSamples];
  FILE *fout;

  for (i=0;i<nFreqSamples;i++)
    {
      tdiff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[i]/1000.0,-2));
      idiff[i] = (int)(tdiff/dSet->head->tsamp);      
    }
  
  *nFinalTimeSamples = 0;
  for (i=0;i<nTimeSamples;i++)
    {
      dedisp[i] = 0;
      nc=0;
      //      if (i==sampleCheck)
      //	fout = fopen("sampleCheck.dat","w");
      for (j=0;j<nFreqSamples;j++)
	{
	  i0 = i - idiff[j];
	  if (i0 > 0 && i0 < nTimeSamples)
	    {
	      dedisp[i]+=(data[(i0)*nFreqSamples+j]); // *wt[i0*nFreqSamples+j]);
	      nc++;
	      //	      if (i==sampleCheck) fprintf(fout,"%d %d %d %g\n",i,j,idiff,data[(i-idiff)*nFreqSamples+j]);
	      //	      if (wt[(i-idiff)*nFreqSamples + j] > 0)
	      //		nc++;
	    }
	}
      //      if (i==sampleCheck)
      //	fclose(fout);
      (*nFinalTimeSamples)++;
      if (nc > 0)
      	dedisp[i]/=(double)nc;
      else dedisp[i] = 0;
    }
}

int compare(const void * a, const void * b) {

  if ( *(float*)a <  *(float*)b ) return -1;
  if ( *(float*)a == *(float*)b ) return 0;
  if ( *(float*)a >  *(float*)b ) return 1;
  
  // backup for function completeness
  return 0;
}
