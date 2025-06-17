 // /usr/bin/gcc -lm -o pfits_inspect2ghz pfits_inspect2ghz.c pfits_loader.c pfits_setup.c T2toolkit.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -O3
 


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include "T2toolkit.h"
#include <cpgplot.h>

int main(int argc,char *argv[])
{
  long i,j,k;
  char fname[1024];
  dSetStruct *dSet;
  int debug=0;
  int nchan,nsub;
  float chbw;
  float f0;
  float flow,fhigh;
  int useChannel[1024];
  int nUseChannel;
  int nsblk;
  int bluetoothChannel = 22;
  float *inspectSeries;
  float *smoothSeries;
  float *loadData;
  long nsamples;
  int nTimeSamples;
  int nFreqSamples;
  float statistic;
  int boxWidth;
  float x,x2,*mean,*rms;
  float threshold=6;
  int position=0;
  long i0,j0;
  long seed = TKsetSeed();
  long n_notFlag=0;
  long n_flag=0;
  float pred[20];
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
    }
  printf("Processing file %s\n",fname);

  initialise(&dSet,debug);
  setFilename(fname,dSet,debug);
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  nchan = dSet->head->nchan;
  nsblk = dSet->head->nsblk;
  nsub = dSet->head->nsub;
  chbw = dSet->head->chanFreq[1]-dSet->head->chanFreq[0];
  f0 = dSet->head->chanFreq[0];
  printf("Nchan = %d, nsub = %d, f0 = %g, chanbw = %g, nsblk = %d\n",nchan,nsub,f0,chbw,nsblk);


  
  // Find Bluetooth LE channel of interest
  if (bluetoothChannel==22)
    {
      flow = 2449;
      fhigh = 2451;
    }
  else
    {
      printf("Don't know this channel\n");
      exit(1);
    }

  nUseChannel=0;
  for (i=0;i<nchan;i++)
    {
      if (dSet->head->chanFreq[i] >= flow && dSet->head->chanFreq[i] <= fhigh) // Check edges here
	{
	  useChannel[nUseChannel++] = i;
	}
    }
  printf("Using %d channels\n",nUseChannel);

  
  // Make a time series for this Bluetooth LE channel
  inspectSeries = (float *)malloc(sizeof(float)*nsub*nsblk);
  smoothSeries = (float *)malloc(sizeof(float)*nsub*nsblk);
  mean = (float *)malloc(sizeof(float)*nsub*nsblk);
  rms = (float *)malloc(sizeof(float)*nsub*nsblk);
  loadData = (float *)malloc(sizeof(float)*nsblk*nchan);
  
  for (i=0;i<nsub;i++)
    {
      //      pfits_read1pol_float(loadData,0,dSet,i,i,1,&nsamples,&nTimeSamples,&nFreqSamples,debug,1);
      pfits_read1pol_float(loadData,0,dSet,i,i,1,&nsamples,&nTimeSamples,&nFreqSamples,debug,0);
      for (j=0;j<nsblk;j++)
	{


	  inspectSeries[i*nsblk + j] = 0.0;
	  for (k=0;k<nUseChannel;k++)
	    inspectSeries[i*nsblk + j] += loadData[j*nFreqSamples + useChannel[k]];
	  inspectSeries[i*nsblk + j] /= (float)nUseChannel;
	  printf("output series: %d %g\n",i*nsblk+j,inspectSeries[i*nsblk+j]);
	}	    
    }
  
  
  
  // Search using a box car for the signal of interest
  // We assume the signals are rare and so can get the RMS from a slightly earlier bit of data
  //    boxWidth=150;
    boxWidth=50;
  for (i=boxWidth;i<nsub*nsblk-boxWidth;i++)
    {
      //    statistic=0;
      smoothSeries[i]=inspectSeries[i];
      x = inspectSeries[i-boxWidth];
      x2 = pow(inspectSeries[i-boxWidth],2);
      for (j=1;j<boxWidth;j++)
	{
	  smoothSeries[i]+=inspectSeries[i+j];
	  x+= inspectSeries[i-boxWidth+j];
	  x2+= pow(inspectSeries[i-boxWidth+j],2);
	}
      rms[i] = sqrt(1./(double)boxWidth * x2 - pow(1.0/boxWidth*x,2));
      mean[i] = x/boxWidth;
      smoothSeries[i] =  smoothSeries[i]/(double)boxWidth; //(smoothSeries[i]/(double)boxWidth-mean[i])/rms[i];
      printf("smooth: %d %g %g\n",i,smoothSeries[i],inspectSeries[i]);
	  //	  printf("Adding: %d %d %g %g %g\n",i,j,inspectSeries[i+j],inspectSeries[i+j]-126.,statistic);
	  //	  statistic += (inspectSeries[i+j]-126.)*1;
    }
  // Now identify significant peaks
  position=0;
  for (i=boxWidth;i<nsub*nsblk-boxWidth;i++)
    {
      if (smoothSeries[i] > threshold && position==0)
	{
	  position=1;
	  i0=i;
	}
      else if (smoothSeries[i] < threshold && position==1)
	{
	  float maxVal;
	  long bestPos;
	  position=2;
	  j0=i;
	  // Find the peak
	  maxVal = smoothSeries[i0];
	  bestPos = i0;
	  for (j=i0;j<j0;j++)
	    {
	      if (smoothSeries[j] > maxVal)
		{
		  maxVal = smoothSeries[j];
		  bestPos = j;
		}
	    }
	  printf("GOT MATCH: %d %d %d\n",i0,j0,bestPos);
	  // Now replace the original data
	  for (j=bestPos;j<bestPos+boxWidth;j++)
	    {
	      inspectSeries[j] = mean[i0] + TKgaussDev(&seed)*rms[i0];
	      n_flag++;
	    }
	  i = bestPos+boxWidth;
	}
      if (position==0)
	n_notFlag++;
      if (position==2)
	position=0;
    }
  //      if (statistic > threshold)
//	{
//	  printf("result: %d %g\n",i,statistic);
//	  for (j=0;j<boxWidth+boxWidth*0.5;j++) // FIX ME HERE
//	    inspectSeries[i+j] = 124+TKgaussDev(&seed)*5;
//	}
//    }
  for (i=0;i<nsub*nsblk;i++)
    printf("finalResult: %d %g\n",i,inspectSeries[i]);
  printf("Not flagged = %d, flaged = %d\n",n_notFlag,n_flag);
  
  free(inspectSeries);
  free(smoothSeries);
  free(mean);
  free(rms);
  free(loadData);
}
