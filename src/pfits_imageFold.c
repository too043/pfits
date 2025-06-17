#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cpgplot.h>
#include "pfits.h"

// gcc -O4 -lm -o pfits_imageFold pfits_imageFold.c -lcpgplot -lcfitsio pfits_loader.c pfits_setup.c


int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  char fname[1024];
  int debug=0;
  int nchan;
  int nsblk;
  float tsamp;
  float *dataVals;
  long int nSamples;
  int nTimeSamples;
  int nFreqSamples;
  int nsub;
  int scls=0;
  int i,ii,j,k;
  int nbin = 128;
  int ibin;
  int ntrial = 20;
  unsigned int **image;
  float *plotImage;
  double period = 0.089389421366021; // 0.016052722944505;
  double tval;
  float tr[6];
  float minVal,maxVal;

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

  nchan = dSet->head->nchan;
  nsblk = dSet->head->nsblk;
  tsamp = dSet->head->tsamp;
  nsub  = dSet->head->nsub;
  image = (unsigned int **)malloc(sizeof(unsigned int *)*ntrial);
  plotImage = (float *)malloc(sizeof(float)*nbin*nchan);
  for (i=0;i<ntrial;i++)
    {
      image[i] = (unsigned int *)calloc(sizeof(unsigned int),nbin*nchan);
    }
  
  if (!(dataVals = (float *)malloc(sizeof(float)*nchan*nsblk)))
    {
      printf("ERROR: unable to allocate memory\n");
      printf("Trying to allocate %d channels and %d time samples\n",nchan,nsblk);
      exit(1);
    }
    for (j=0;j<nsub;j++)
  //  for (j=0;j<10;j++)
    {
      printf("Reading subintegration %d/%d\n",j,nsub);

      // Should just read the lot into memory as bytes ...
      pfits_read1pol_float(dataVals,0,dSet,j,j,1,&nSamples,&nTimeSamples,&nFreqSamples,debug,scls);
      for (i=0;i<nTimeSamples;i++)
	{
	  tval = j*nsblk*tsamp+i*tsamp;
	  for (ii=0;ii<ntrial;ii++)
	    {
	      ibin = (int)((tval/(period+ii*1e-5)-(int)(tval/(period+ii*1e-5)))*nbin);
	      for (k=0;k<nFreqSamples;k++)
		{
		  image[ii][ibin*nchan+k]+=(int)(dataVals[i*nchan+k]*2+5);
		}
	    }
	}
    }
  printf("Now searching for min/max\n");
  minVal = 1e30; maxVal=-1e30;
  for (i=0;i<nbin;i++)
    {
      for (k=0;k<nFreqSamples;k++)
	{
	  if (minVal > image[0][i*nchan+k]) minVal = image[0][i*nchan+k];
	  if (maxVal < image[0][i*nchan+k]) maxVal = image[0][i*nchan+k];
	}
    }
  printf("min/max = %g/%g\n",minVal,maxVal);
  if (1==1)
    {
      tr[0] = 0; tr[1] = 0; tr[2] = 1;
      tr[3] = 0; tr[4] = 1; tr[5] = 0;
      cpgbeg(0,"1/xs",1,1);
      for (ii=0;ii<ntrial;ii++)
	{
	  for (i=0;i<nchan*nbin;i++)
	    plotImage[i] = image[ii][i];
	  cpgenv(0,nbin,0,nchan,0,1);
	  cpggray(plotImage,nchan,nbin,1,nchan,1,nbin,minVal,maxVal,tr);
	}
      cpgend();
    }
  free(dataVals);
  free(plotImage);
  for (i=0;i<ntrial;i++)
    free(image[i]);

  free(image);
  free(dSet);
}
