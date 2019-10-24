// gcc -lm -o pfits_pulseExtraction pfits_pulseExtraction.c -lcpgplot -lcfitsio pfits_loader.c pfits_setup.c
//
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  char fname[1024];
  int debug=0;
  int i,j,k;
  int n;
  long npts=0;
  float *data1sub,*data;
  long nsamples;
  int nTimeSamples;
  int nFrequencySamples;
  int offScl=0;
  int sub=0;
  float *dedispTime;
  FILE *fout;
  int ch1 = 600; // 450;
  int ch2 = 800;
  //  float dm = 88.373; // For J1745-3040
  float dm = 16.1356; // For J1946+1805
  float fref,toff;
  int sampleOff;
  int bs=4;
  float final;
  
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
    }
  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));
  initialise(&dSet,debug);
  setFilename(fname,dSet,debug);

  pfitsOpenFile(dSet,debug);
  printf("Loading header\n");
  pfitsLoadHeader(dSet,debug);

  data1sub = (float *)malloc(sizeof(float)*dSet->head->nchan*dSet->head->nsblk);
  data = (float *)malloc(sizeof(float)*(ch2-ch1)*dSet->head->nsblk*dSet->head->nsub);
  dedispTime = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nsub);
  for (i=0;i<dSet->head->nsblk*dSet->head->nsub;i++)
    dedispTime[i]=0.0;
  
  for (i=0;i<dSet->head->nsub;i++) 
    {
      printf("Loading subint %d out of %d\n",i,dSet->head->nsub);
      pfits_read1pol_float(data1sub,0,dSet,i,i,1,&nsamples,&nTimeSamples,&nFrequencySamples,debug,offScl);

      for (j=0;j<dSet->head->nsblk;j++)
	{
	  for (k=ch1;k<ch2;k++)
	    data[i*dSet->head->nsblk*(ch2-ch1)+j*(ch2-ch1)+(k-ch1)]
	      = data1sub[j*dSet->head->nchan+k];
	}
    }
  printf("Completed loading\n");

  fref = (dSet->head->chanFreq[ch1] + dSet->head->chanFreq[ch2])/2.;
  
  npts=0;
  for (i=0;i<dSet->head->nsblk*dSet->head->nsub;i++)
    {
      n=0;
      for (j=0;j<ch2-ch1;j++)
	{
	  toff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[j+ch1]/1000.0,-2));
	  sampleOff = (int)(toff/dSet->head->tsamp);
	  if (i-sampleOff > 0)
	    {
	      dedispTime[i]+=data[(i-sampleOff)*(ch2-ch1)+j];
	      n++;
	    }
	}
      dedispTime[i]/=(float)n;
      npts++;
    }
  printf("Finished %d\n",(int)npts);

  fout = fopen("dedisp.dat","w");
  for (i=0;i<npts;i+=bs)
    {
      final =0.0;
      for (j=0;j<bs;j++)
	final+=dedispTime[i+j];
        fprintf(fout,"%d %g\n",i,final/(double)bs);
    }
  fclose(fout);
  
  cpgend();
  

  //  pfitsCloseFile(dSet,debug);
  free(dSet);
  free(data);
  free(data1sub);
  free(dedispTime);
}




