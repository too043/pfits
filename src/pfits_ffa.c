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

//  /usr/bin/gcc -lm -o pfits_ffa pfits_ffa.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio
// /usr/bin/gcc -lm -o pfits_ffa_zwicky pfits_ffa.c pfits_setup.c pfits_loader.c -I/usr/lib/pgplot5 -L/usr/lib/pgplot5 -lcpgplot -lpgplot -L/pulsar/psr/software/stable/bullseye/lib -L/usr/lib/gcc/x86_64-linux-gnu/10/../lib -L/usr/lib/gcc/x86_64-linux-gnu/10/../../.. -lgfortran -lm -lquadmath -lcfitsio


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include <assert.h>

#define ZERO_PADDING 0
#define K 1.4826
#define TRUE 1
#define FALSE 0
#define MAX_FILTER_WIDTH 0.2
#define ARRAY_PADDING 2

typedef double ffadata;
typedef struct paddedArray {
  ffadata* dataarray;
  int nchan;
  int datasize;
  int fullsize;
  int scalefactor;
  int redflag;
  int window;
} paddedArray;


void ffa(float *data,  dSetStruct *dSet,int nsubLoaded,int lowperiod,int highperiod,int prelim_downsamples,int flip,int noFiles);
void mad(ffadata* array, int size);
ffadata postMadMatchedFilterMetric(ffadata* sourcearray, int startpos, int subsize);
paddedArray* createPaddedArray(int datasize, int fullsize,int nchan);
ffadata* getPaddedArrayDataArray(paddedArray* x);
void setRedFlag(paddedArray *x, int newflag);
void setWindow(paddedArray *x, int newwindow);
int getPaddedArrayDataSize(paddedArray* x);
void singleFFA(FILE* outputfile, FILE* profilefile, FILE* normprofilefile, paddedArray* sourcedata, int baseperiod, double (*metric)(ffadata*, int, int), int mfsize,int noFiles);
void massFFA(FILE* outputfile, FILE* profilefile, FILE* normprofilefile, paddedArray* sourcedata, int lowperiod, int highperiod, double (*metric)(ffadata*, int, int), int mfsize,  int prelim_ds, FILE* redfile, int PRESTO_flag, int timenorm_flag,int noFiles);
void deletePaddedArray(paddedArray* x);
int getPaddedArrayScaleFactor(paddedArray *x);
paddedArray* downsampleDataArray(paddedArray* sourcedata);
int getPaddedArrayFullSize(paddedArray* x);
ffadata* copyPaddedArrayDataArray(paddedArray* x);
ffadata generateZeroPadding();
ffadata resample(ffadata x, ffadata y);
void profiledump(FILE* profilefile, double period, int scalefactor, ffadata* sumarrays, int startpos, int subsize,int nchan);
int power2Resizer(int size, int period);
void slideAdd(ffadata* sourcearray, ffadata* resultarray, int sourcesubstartpos1, int sourcesubstartpos2, int resultsubstart, int subsize, int slide,int nchan);
int getRedFlag(paddedArray *x);
int getWindow(paddedArray *x);
ffadata add(ffadata x, ffadata y);
int compare(const void * a, const void * b);
void setPaddedArrayScaleFactor(paddedArray *x, int newfactor);

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,j,k;
  float *data;
  int noFiles=0;
    
  char outFile[128]="dumpValues";
  char oName[1024];
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
  float chbw;
  int setFreq=0;
  int sclOff=0;
  int tAv=1;
  long n,nc;
  int lowperiod = 1024;
  int highperiod = 1025;
  int prelim_downsamples = 0;
  int flip=0;
  
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
      else if (strcasecmp(argv[i],"-nofiles")==0)
	noFiles=1;
      else if (strcmp(argv[i],"-flip")==0)
	flip=1;
      else if (strcmp(argv[i],"-ds")==0)
	sscanf(argv[++i],"%d",&prelim_downsamples);
      else if (strcasecmp(argv[i],"-tav")==0)
	sscanf(argv[++i],"%d",&tAv);
      else if (strcmp(argv[i],"-lp")==0)
	sscanf(argv[++i],"%d",&lowperiod);
      else if (strcmp(argv[i],"-hp")==0)
	sscanf(argv[++i],"%d",&highperiod);
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
  chbw = dSet->head->chanFreq[1] - dSet->head->chanFreq[0];
  
  if (setFreq==1)
    {
      i0=i1=-1;
      
      for (i=0;i<dSet->head->nchan-1;i++)
	{
	  if (chbw > 0)
	    {
	      if (i0==-1 && dSet->head->chanFreq[i] >= f0)
		i0 = i;
	      if (i1==-1 && dSet->head->chanFreq[i] >= f1)
		i1 = i;
	    }
	  else
	    {
	      if (i0==-1 && dSet->head->chanFreq[i] <= f0)
		i0 = i;
	      if (i1==-1 && dSet->head->chanFreq[i] <= f1)
		i1 = i;
	    }
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

  // Load all into data
  data = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan*(sub1-sub0));
  
  totCount=0;


  //  fout = fopen(outFile,"wb");

  printf("Loading input data\n");
  for (subint=sub0;subint<sub1;subint++)
    pfits_read1pol_float(data+(subint-sub0)*dSet->head->nsblk*dSet->head->nchan,pol,dSet,subint,subint,1,&nSamples,&nTime,&nFreq,debug,sclOff);
  printf("Complete loading input data\n");
  ffa(data,dSet,sub1-sub0,lowperiod,highperiod,prelim_downsamples,flip,noFiles);
  //  fclose(fout);

  free(data);
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
}


void ffa(float *data,  dSetStruct *dSet,int nsubLoaded,int lowperiod,int highperiod,int prelim_downsamples,int flip,int noFiles)
{
  int i; // counter
  int datasize = 0;
  int nchan=dSet->head->nchan;
  unsigned long kk;
  ffadata* array;
  int loops = 1;
  int seedperiod = 500;
  int seedwidth = 15;
  FILE *outputfile = NULL;
  FILE *profilefile = NULL;
  FILE *normprofilefile = NULL;
  FILE *parrotfile = NULL;
  FILE *originalfile = NULL;
  FILE *originalderedfile = NULL;
  paddedArray* sourcedata = NULL;
  int mfsize = 0;
  int dered_flag = FALSE;
  double (*metric)(ffadata*, int, int);
  int dered_window = 1;
  int PRESTO_flag = FALSE;
  int timenorm_flag = FALSE;

  outputfile = fopen("ffa_output.dat","w");
  
  assert(lowperiod >= 2); // this is the smallest possible period size for the FFA - a period of 1 results in no array shifting
  assert(highperiod > lowperiod);
  assert(outputfile != NULL);
  assert(seedperiod > seedwidth);
  if (mfsize < 0) {
    printf("Matched filter size cannot be less than zero!\n");
    exit(0);
  }
  if (dered_window <= 0) {
    printf("De-reddening window must be greater than 0!\n");
    exit(0);
  }
  metric = postMadMatchedFilterMetric;
  
  datasize = dSet->head->nsblk*nsubLoaded;
  // GH HACK
  //    int paddedsize = (int)ceil(datasize*ARRAY_PADDING);
  int paddedsize = (int)ceil(datasize*1);
  // to assist with downsampling later, make padded size a multiple of two
  if ((paddedsize % 2) != 0) {
    paddedsize++;
  }
  printf("Creating padded array\n");
  sourcedata     = createPaddedArray(datasize, paddedsize,nchan);
  printf("Complete that bit\n");
  array          = getPaddedArrayDataArray(sourcedata);
  printf("Reading input file\n");
  for (kk=0;kk<datasize;kk++)
    {
      //	printf("Loaded %g\n",kk*100.0/datasize);
      //	fscanf(inputfile,"%f",&x[0]);
      for (i=0;i<nchan;i++)
	{
	  //	    array[i*datasize+kk] = (ffadata)x[i];
	  if (flip==0)
	    array[(unsigned long)kk*(unsigned long)nchan+(unsigned long)i] = (ffadata)data[kk*nchan+i]; // Could use something more efficint
	  else
	    array[(unsigned long)kk*(unsigned long)nchan+(unsigned long)i] = (ffadata)(-data[kk*nchan+i]+1);
	  //  array[(unsigned long)kk*(unsigned long)nchan+(unsigned long)i] = (ffadata)x[0];
	}
    }
  printf("Finished load. Now padding %ld %d\n",kk,paddedsize);
  //   while (kk < paddedsize)
  //      array[kk] = generateZeroPadding();
  printf("Finished reading input file\n");
  
  printf("AA\n");
  setRedFlag(sourcedata, dered_flag);
  printf("BB\n");
  setWindow(sourcedata, dered_window);
  printf("At this point\n");
  printf("Input file read successfully - %d samples in length.\n", getPaddedArrayDataSize(sourcedata));

  assert(getPaddedArrayDataSize(sourcedata) > highperiod);

  printf("Now scanning from a period of %d samples to %d original samples...\n", lowperiod, highperiod);
  
  // now ready to begin FFA
  profilefile = fopen("profileFile.ignore","w");
  massFFA(outputfile, profilefile, normprofilefile, sourcedata, lowperiod, highperiod, metric, mfsize, prelim_downsamples, originalderedfile, PRESTO_flag, timenorm_flag,noFiles);

  // file I/O should now be complete - close files
  fclose(outputfile);
  if (profilefile != NULL) {
    fclose(profilefile);
  }
  if (normprofilefile != NULL) {
    fclose(normprofilefile);
  }
  if (originalderedfile != NULL) {
    fclose(originalderedfile);
  }

  // cleanup
  deletePaddedArray(sourcedata);

  printf("\nFFA complete.\n");
 
  //  return 0;
  

}


ffadata postMadMatchedFilterMetric(ffadata* sourcearray, int startpos, int subsize) {

  assert(subsize > 0);
  assert(sourcearray != NULL);

  // metric needs to scan array using successively larger matched filters up to some set limit
  // can use the same principle as Kondratiev - 20%? 25%? Use nearest power of two?
  // Run with 20% to match Kondratiev, choosing next highest power of 2 above the 20% width as the max filter size

  int n_layers = (int)ceil(log2(subsize * MAX_FILTER_WIDTH));

  double max_SNR = 0; // stores the maximum SNR detection from this profile across all matched filters

  // loop counters
  int ii;
  int jj;

  // create a copy of the array
  ffadata* normarray = (ffadata*)malloc(sizeof(ffadata) * subsize);
  for (ii = 0; ii < subsize; ii++) {
    normarray[ii] = sourcearray[ii + startpos];
  }
 
  // normalise it using MAD 
  mad(normarray, subsize);

  // as a baseline, first perform a scan for a matched filter size of 1
  for (ii = 0; ii < subsize; ii++) {    
    if (normarray[ii] > max_SNR) {
      max_SNR = normarray[ii];
    }
  }

  // now begin the metric proper - begin scanning through successive matched filters of power 2
  // build copy and temp arrays to use for optimised pointer swapping
  ffadata* resultarray = (ffadata*)malloc(sizeof(ffadata) * subsize);
  ffadata* temparray;

  for (ii = 0; ii < n_layers; ii++) {
    int shift = 1 << ii;

    // convolve the matched filter and scan for max SNR simultaneously
    for (jj = 0; jj < subsize; jj++) {
      resultarray[jj] = normarray[jj] + normarray[(jj + shift + subsize)%subsize];
      
      if (resultarray[jj]/sqrt(shift * 2) > max_SNR) {
	max_SNR = resultarray[jj]/sqrt(shift * 2);
      }
    }
    
    // convolved array is complete
    // now swap arrays for next loop
    
    temparray = normarray;
    normarray = resultarray;
    resultarray = temparray;
  }

  // cleanup
  free(resultarray);
  free(normarray);

  // max SNR should by now have been isolated
  return max_SNR;

}

// Allocates the memory for a padded array struct and returns a pointer. Internal values are uninitialised.
paddedArray* createPaddedArray(int datasize, int fullsize,int nchan) {

  paddedArray* x = (paddedArray*)malloc(sizeof(paddedArray));
  assert(x != NULL);

  x->dataarray = (ffadata*)malloc(sizeof(ffadata)*fullsize*nchan);
  assert(x->dataarray != NULL);
  x->nchan = nchan;
  x->datasize = datasize;
  x->fullsize = fullsize;
  x->scalefactor = 1;
  x->redflag = FALSE;
  x->window = 1;

  return x;

}

ffadata* getPaddedArrayDataArray(paddedArray* x) {

  assert(x != NULL);

  return x->dataarray;

}

void setRedFlag(paddedArray *x, int newflag) {

  assert(x != NULL);

  if (newflag) {
    x->redflag = TRUE;
  } else {
    x->redflag = FALSE;
  }

  return;

}


void setWindow(paddedArray *x, int newwindow) {

  assert(x != NULL);
  assert(newwindow >= 0);

  x->window = newwindow;

  return;

}

int getPaddedArrayDataSize(paddedArray* x) {

  assert(x != NULL);

  return x->datasize;

}

void massFFA(FILE* outputfile, FILE* profilefile, FILE* normprofilefile, paddedArray* sourcedata, int lowperiod, int highperiod, double (*metric)(ffadata*, int, int), int mfsize,  int prelim_ds, FILE* redfile, int PRESTO_flag, int timenorm_flag,int noFiles) {

  // UPDATE - THIS SCRIPT MUST REFER ANY DE-REDDENING AND RESULTANT DOWNSAMPLING BACK TO THE ORIGINAL SOURCEDATA ARRAY FOR COMPUTATIONAL CORRECTNESS
  // DOUBLE UPDATE 15/04/2016 - THE DOWNSAMPLING FUNCTION NO LONGER INCLUDES AUTOMATIC DE-REDDENING
  // IT IS NOW APPLIED ONCE AFTER THE FULL DOWNSAMPLING LOOP IS COMPLETE

  // validity checks
  printf("Check 1\n");
  assert(outputfile != NULL);
  printf("Check 2\n");
  assert(sourcedata != NULL);
  printf("Check 3\n");
  if (fmod(lowperiod, getPaddedArrayScaleFactor(sourcedata)) != 0) {
    printf("ERROR: Low search period (%d) must be integer multiple of initial downsampling factor (%d).\nPlease correct and try again.\n", lowperiod, getPaddedArrayScaleFactor(sourcedata));
    exit(EXIT_FAILURE);
  }

  // This function needs to generate all the base periods to test, re-scale the array for each base period, and then execute the individual FFA algorithm for each base period
  // Function also needs to oversee downsampling by a factor of 2 each time the lowperiod is doubled.

  int i = lowperiod;
  paddedArray* workingdata = sourcedata;
  paddedArray* tempdata;
  int loopscalefactor = 0; //used for controlling the downsampling during FFA operation

  fprintf(outputfile, "# Period (original samples) | Downsample factor | Period (downsampled samples) | Metric\n");

  printf("Check 4\n");
  while (i < highperiod) {

    if (i == lowperiod*((int)pow(2, loopscalefactor))) {
      // this means that we have either just started the scan, or have reached a downsampling point
      printf("\nDownsampling initiated: i = %d\n", i);
      
      // the process will first involve de-reddening the datarray (but only if de-reddening is selected
      // clean up the current workingdata - it is no longer needed
      if (workingdata != sourcedata) {
	deletePaddedArray(workingdata);
      }

      // reset to sourcedata
      workingdata = sourcedata;

      // array has been de-reddened, now downsample according to the correct number of loops
      int jj = 0;
      printf("prelim_ds = %d loopscalefactor = %d\n",prelim_ds,loopscalefactor);
      while (jj < prelim_ds + loopscalefactor) {
	printf("Doing downsampling\n");
	tempdata = downsampleDataArray(workingdata);
	printf("Done downsampling\n");
	// we need to clean up the old workingdata struct
	if (workingdata != sourcedata) {
	  deletePaddedArray(workingdata);
	}
	workingdata = tempdata;

	jj++;
	printf("Downsample loop %d completed.\n", jj);
      }
      
      printf("Downsample factor: %d\n", prelim_ds + loopscalefactor);


      // modify the internal scale factor
      loopscalefactor++;

    }

    // downsampling, if neccessary, is now complete
    
    
    // we now need to pass the relevant parameters to the singleFFA function
    // baseperiod must be calculated to match with the current downsampling
    printf("\nCalling single FFA search for a period of %d original samples...\n", i);
    singleFFA(outputfile, profilefile, normprofilefile, workingdata, i/getPaddedArrayScaleFactor(workingdata), metric, mfsize,noFiles);

    // increment i according to the scale factor
    i = i + getPaddedArrayScaleFactor(workingdata);

  }

  // final clean up
  if ((workingdata != sourcedata) && (workingdata != NULL)) {
      deletePaddedArray(workingdata);
  }

  return;
}

void singleFFA(FILE* outputfile, FILE* profilefile, FILE* normprofilefile, paddedArray* sourcedata, int baseperiod, double (*metric)(ffadata*, int, int), int mfsize,int noFiles) {
  int nchan = sourcedata->nchan;
  
  // basic validity checks
  assert(outputfile != NULL);
  assert(sourcedata != NULL);

  printf("Entered singleFFA with baseperiod of %d samples and a scalefactor of %d...\n", baseperiod, getPaddedArrayScaleFactor(sourcedata));
  printf("Data size = %d %d\n",getPaddedArrayDataSize(sourcedata),baseperiod);
  // need the size of the array to use based on N/n = 2^x
  int size = power2Resizer(getPaddedArrayDataSize(sourcedata), baseperiod);
  // double check that the array size is still within memory limits 
  printf("Out with size = %d\n",size);
  assert(size <= getPaddedArrayFullSize(sourcedata));
  printf("Array size rescaled from %d to %d (%.1f%% change).\n", getPaddedArrayDataSize(sourcedata), size, abs(getPaddedArrayDataSize(sourcedata) - size)*100/(float)(getPaddedArrayDataSize(sourcedata)));

  // initialise counters
  int i, j, k;

  // setup variables controlling the scale of the FFA
  int branches = (int)size/baseperiod;
  int addition_iterations = (int)log2(branches);
  double period_increment = (double)1/((double)(branches - 1));

  // build this many new arrays matching the original size to store the cumulative additions - the first should be the source array
  ffadata* sumarrays[addition_iterations + 1];

  for (i = 0; i <= addition_iterations; i++) {
    if (i == 0) {
      sumarrays[i] = copyPaddedArrayDataArray(sourcedata);
    } else {
      sumarrays[i] = (ffadata*)malloc(sizeof(ffadata)*size*nchan);
    }
  }

  // NEW SECTION - HANDLES ZERO PADDING ISSUE
  // If array has been padded out, then a branch of the sourcearray data will contain part data and part zeroes, causing baseline jumps and false detections
  // This row must be entirely set to zeroes
  for (i = 0; i < branches; i++) {
    if ((((i + 1) * baseperiod) > getPaddedArrayDataSize(sourcedata)) && ((i * baseperiod) < getPaddedArrayDataSize(sourcedata))) {
      // if we are a partial row of data, clean it out
      for (j = i*baseperiod ; j < (i+1)*baseperiod ; j++) {
	{
	  for (k=0;k<nchan;k++)
	    (sumarrays[0])[j*nchan+k] = generateZeroPadding();
	}
      }
    }
  }

  // for file output
  double period;

  // start counting through the addition steps
  for (i = 1; i <= addition_iterations; i++) {

    ffadata* startarray = sumarrays[i-1];
    ffadata* endarray   = sumarrays[i];

    // a segment represents the self-contained module of array elements that are adding together at each addition step
    int segmentsize = (int)pow(2, i);
    int segments = (int)branches/segmentsize;

    for (j = 0; j < segments; j++) {

      for (k = 0; k < segmentsize; k++) {
	int slide = (int)ceil((float)k/2);
	int sourcecellpos1 = ((int)floor((float)k/2) + j*segmentsize)*baseperiod;

	// we have now honed in on the result cell, and have enough information to select the source cells to use in the addition and the slide amount
	// add sub array cells
	slideAdd(startarray, endarray, sourcecellpos1, sourcecellpos1 + baseperiod*segmentsize/2, (k + j*segmentsize)*baseperiod, baseperiod, slide,nchan);

	// if this is the last iteration of the FFA additions, we can output the metric now without having to re-scan the loop
	if (i == addition_iterations) {
	  period = k * period_increment + baseperiod; // this is the tested period in units of (downsampled) samples

	  // normalise the profile for post-MAD
	  //postMadProfileNormaliser(sumarrays[i], (k + j*segmentsize)*baseperiod, baseperiod, (int)ceil(sourcedata->datasize/((double)baseperiod)));
	  //postMadProfileNormaliser(sumarrays[i], (k + j*segmentsize)*baseperiod, baseperiod, branches);
	  
	  fprintf(outputfile, "%.10f %d %.10f %.10f\n", period*getPaddedArrayScaleFactor(sourcedata), getPaddedArrayScaleFactor(sourcedata), period, metric(sumarrays[i], (k + j*segmentsize)*baseperiod, baseperiod));
	  if (noFiles==0)
	    profiledump(profilefile, period*getPaddedArrayScaleFactor(sourcedata), getPaddedArrayScaleFactor(sourcedata), sumarrays[i], (k + j*segmentsize)*baseperiod, baseperiod,sourcedata->nchan);

	}
      }

    }

  }

  // individual FFA execution should now be complete

  // free memory - but don't delete the original array that is part of the paddedArray struct

  for (i = 0; i <= addition_iterations; i++) {
    free(sumarrays[i]);
  }
  return;
}

// leans up the DataArray once processing is complete
void deletePaddedArray(paddedArray* x) {

  assert(x != NULL);

  // delete memory-allocated contents
  free(x->dataarray);

  // delete struct itself
  free(x);

  // all memory cleared

  return;
}

void mad(ffadata* array, int size) {

  //printf("Entered MAD normalisation function...\n");

  assert(array != NULL);
  int i;

  // STEP 1 - Get the median of the array
  
  // copy array for sorting
  ffadata* sortedarray = (ffadata*)malloc(sizeof(ffadata)*size);
  for (i = 0; i < size; i++) {
    sortedarray[i] = array[i];
    
    // TEST SECTION
    //printf("%d %f\n", i, sortedarray[i]); 

  }

  //printf("Array copied...\n");

  // sort array
  qsort(sortedarray, size, sizeof(ffadata), compare);

  // TEST SECTION
  /*printf("Array sorted...displaying contents\n");
  for (i = 0; i < size; i++) {
     printf("%d %f\n", i, sortedarray[i]); 
     }*/

  // get median
  int median_pos = floor((double)size/((double)2));
  ffadata median = sortedarray[median_pos];

  //printf("Median obtained = %f...\n", median);
  
  // STEP 2 - Remove median and calculate the absolute value of the deviances

  ffadata* deviances = (ffadata*)malloc(sizeof(ffadata)*size);
  for (i = 0; i < size; i++) {
    array[i] = array[i] - median;
    deviances[i] = fabs(array[i]);
    //printf("Deviance %d: %f\n", i, array[i]); 
    //printf("Abs. Deviance %d: %f\n", i, deviances[i]); 
  }

  //printf("Median removed...\n");

  // STEP 3 - Get the median of the deviances (MAD)
  qsort(deviances, size, sizeof(ffadata), compare);
  ffadata median_deviance = deviances[median_pos];

  //printf("Median deviance obtained = %f...\n", median_deviance);

  // STEP 4 - Divide all elements by MAD * K
  for (i = 0; i < size; i++) {
    array[i] = array[i]/(median_deviance * K);
  }

  free(sortedarray);
  free(deviances);

  //printf("MAD normalisation complete.\n");
  return;
}

int getPaddedArrayScaleFactor(paddedArray *x) {

  assert(x != NULL);

  return(x->scalefactor);

}

paddedArray* downsampleDataArray(paddedArray* sourcedata) {
  int nchan = sourcedata->nchan;
  printf("In downsampling nchan = %d\n",nchan);
  assert(sourcedata != NULL);
  printf("In downsampling 2\n");
  // NOTE: Issue of whether to use noisy padding or zero padding is still undecided.
  // Code relating to noisy generation has been left inside this function but has been commented out.
  // Future functionality may include ability to nominate noise generation method via input variable.

  // Determine the sizes of the elements in the downsampled paddedArray based on the size of the elements in sourcedata
  // Need both datasizes to be mulitiples of 2
  int sourcedatasize = getPaddedArrayDataSize(sourcedata);
  int sourcefullsize = getPaddedArrayFullSize(sourcedata);

  if (sourcedatasize % 2 != 0) {
    sourcedatasize = sourcedatasize + 1;
  }
  if (sourcefullsize % 2 != 0) {
    sourcefullsize = sourcefullsize + 1;
  }
  
  // this may cause the sizes to be slightly larger than they were originally
  // if sourcedatasize > datasize, it should just read into the empty padding elements
  // if sourcefullsize > fullsize, then there should be a check in the downsampling loop that if an element does not exist in the original array, it is filled in as padding

  // Both sizes are now divisible by two - allocate new sizes
  int downdatasize = sourcedatasize/2;
  int downfullsize = sourcefullsize/2;
  
  // check that the padding ratio still holds, if not, correct this
  printf("In downsampling 3\n");
  while (downfullsize < downdatasize * ARRAY_PADDING) {
    downfullsize = downfullsize + 2;
  }
  
  // sizes configured - initialise the paddedArray
  printf("In downsampling 4\n");
  paddedArray* downdata = createPaddedArray(downdatasize, downfullsize,nchan);
  setRedFlag(downdata, getRedFlag(sourcedata));
  setWindow(downdata, getWindow(sourcedata));
  printf("In downsampling 5\n");
  ffadata* downdataarray = getPaddedArrayDataArray(downdata);
  printf("In downsampling 6\n");
  // now perform downsampling using abstracted "average" function to combine values 
  int i, index1, index2, element1, element2,k;
  
  for (k=0;k<nchan;k++)
    {
      for (i = 0; i < downdata->fullsize; i++) {
	index1 = i*2;
	index2 = index1 + 1;
	
	if (index1 > sourcedata->fullsize) {
	  element1 = generateZeroPadding();
	  element2 = generateZeroPadding();
	} else if (index2 > sourcedata->fullsize) {
	  element1 = sourcedata->dataarray[index1*nchan+k];
	  element2 = generateZeroPadding();
	} else {
	  element1 = sourcedata->dataarray[index1*nchan+k];
	  element2 = sourcedata->dataarray[index2*nchan+k];
	}
	
	// now seed new array with the averaged value
	downdataarray[i*nchan+k] = resample(element1, element2);
      }
    }
  printf("In downsampling 7\n");
  // downsampling complete - modify scalefactor
  setPaddedArrayScaleFactor(downdata, getPaddedArrayScaleFactor(sourcedata) * 2);

  // preliminary downsampling is now complete - can now execute the de-reddening
  /*if (getRedFlag(downdata)) {

    printf("De-reddening as part of downsampling...\n");
    
    paddedArray* outdata = dereddenDataArray(downdata);

    // dereddening now complete - replace downdata with outdata and cleanup
    deletePaddedArray(downdata);
    downdata = outdata;

    }
  */

  return downdata;

}

int getPaddedArrayFullSize(paddedArray* x) {

  assert(x != NULL);

  return x->fullsize;

}

ffadata* copyPaddedArrayDataArray(paddedArray* x) {

  assert(x != NULL);
  int fullsize = getPaddedArrayFullSize(x);
  int nchan = x->nchan;
  
  ffadata* copyarray = (ffadata*)malloc(sizeof(ffadata)*fullsize*nchan);
  ffadata* sourcearray = getPaddedArrayDataArray(x);
  assert(sourcearray != NULL);

  int i;
  for (i = 0; i < fullsize*nchan; i++) {
    copyarray[i] = sourcearray[i];
  }

  return copyarray;

}

ffadata generateZeroPadding() {

  return ZERO_PADDING;

}

ffadata resample(ffadata x, ffadata y) {

  // basic implementation - just return the sum
  return (x + y);

}

// prints out the full profiles produced by an FFA folding sequence to specified filestream
// Format will be "TrialPeriod(%.10f) ScaleFactor(%d) Bin1(%d) Bin2(%d) etc..."
void profiledump(FILE* profilefile, double period, int scalefactor, ffadata* sumarrays, int startpos, int subsize,int nchan)
{

  char fname[1024];
  FILE *fout;
  double mean[nchan];
  sprintf(fname,"george_ffa_out_%.6f_%d_%d_%d.txt",period,scalefactor,startpos,subsize);
  fout = fopen(fname,"w");
  assert(profilefile != NULL);
  assert(sumarrays != NULL);

  // output the first part of the line to the file
  //  fprintf(profilefile, "%.10f %d", period, scalefactor);

  // scan through the array and output datapoints in the profile to a continuous line of outout, separated by space delimiters
  int i,j;
  for (i=0;i<nchan;i++)
    mean[i] = 0;
  for (i = 0; i < subsize; i++)
    {
      for (j=0;j<nchan;j++)
	{
	  mean[j] += sumarrays[(i+startpos)*nchan+j];
	}
    }

    for (i = 0; i < subsize; i++)
    {
      for (j=0;j<nchan;j++)
	{

	  //	  fprintf(fout, "%d %d %g\n", i,j,sumarrays[(i+startpos)*nchan+j]);
	  fprintf(fout, "%d %d %g\n", i,j,sumarrays[(i+startpos)*nchan+j]-mean[j]/(double)subsize);
	  
	}
      fprintf(fout,"\n");
      //    fprintf(profilefile, " %d", (int)sourcearray[i+startpos]);
    } 
  
  fclose(fout);
  return;
}

int power2Resizer(int size, int period) {

  // work out the nearest ratios of size/period = 2^x that side above the existing size/period

  double actualpower = log2((double)size/(double)period);
  int highpower2 = (int)pow(2, ceil(actualpower));
  
  int highsize = highpower2 * period;
  printf("resizer with %g %d %d\n",actualpower,highpower2,highsize);


  return highsize;
}

void slideAdd(ffadata* sourcearray, ffadata* resultarray, int sourcesubstartpos1, int sourcesubstartpos2, int resultsubstart, int subsize, int slide,int nchan) {

  assert(sourcearray != NULL);
  assert(resultarray != NULL);

  // scroll through the elements of the result sub array, and add in the required elements of the source sub arrays, with the appropriate slide
  int i, j,k;

  for (k=0;k<nchan;k++)
    {
      for (i = 0; i < subsize; i++) {
	
	// result element is i
	// the position in the first source sub array will just be i
	// the position in the second source sub array will be i + slide, modulated by the subsize
	j = (i + slide) % subsize;
	
	int sourcesubpos1 = i + sourcesubstartpos1;
	int sourcesubpos2 = j + sourcesubstartpos2;
	
	// now add
	resultarray[(i + resultsubstart)*nchan+k] = add(sourcearray[sourcesubpos1*nchan+k], sourcearray[sourcesubpos2*nchan+k]);
	
      }
    }

  return;
}

int getRedFlag(paddedArray *x) {

  assert(x != NULL);

  return(x->redflag);

}

int getWindow(paddedArray *x) {

  assert(x != NULL);

  return(x->window);

}

ffadata add(ffadata x, ffadata y) {

  ffadata result = x + y;
  
  return result;
}

int compare(const void * a, const void * b) {

  if ( *(ffadata*)a <  *(ffadata*)b ) return -1;
  if ( *(ffadata*)a == *(ffadata*)b ) return 0;
  if ( *(ffadata*)a >  *(ffadata*)b ) return 1;
  
  // backup for function completeness
  return 0;
}

void setPaddedArrayScaleFactor(paddedArray *x, int newfactor) {
  
  assert(x != NULL);

  x->scalefactor = newfactor;

  return;

}
