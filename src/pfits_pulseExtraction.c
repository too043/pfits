// gcc -lm -o pfits_pulseExtraction pfits_pulseExtraction.c -lcpgplot -lcfitsio pfits_loader.c pfits_setup.c

// gcc -lm -o pfits_pulseExtraction pfits_pulseExtraction.c T2toolkit.c -lcpgplot -lcfitsio pfits_loader.c pfits_setup.c -L/u/hob044/software/new_c/pfits/libsndfile-1.0.25/src/ -I/u/hob044/software/new_c/pfits/libsndfile-1.0.25/src/ -lsndfile

//
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include <sndfile.h>
#include "T2toolkit.h"

#define AMPLITUDE (1.0 * 0x7F000000)

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
  //  float dm = 16.1356; // For J1946+1805
  float dm = 67.97; // For J0835-4510
  float fref,toff;
  int sampleOff;
  int bs=4;
  float final;

  SNDFILE *file;
  SF_INFO sfinfo;
  int *buffer;
  float *convolve;
  int sample_count;
  int ival;
  float sample_rate = 44100;
  float min,max;
  float amp;
  long seed = TKsetSeed();
  long np,nc;
  float convolveW = 20;
  float nConvolve = 1000;

  //  nConvolve=0;
  
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
  convolve = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nsub);

  for (i=0;i<dSet->head->nsblk*dSet->head->nsub;i++)
    dedispTime[i]=0.0;

  dSet->head->nsub = 400;
  
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

  sample_count = 200000;
  
  buffer = (int *)malloc(sizeof(int)*sample_count);
  memset(&sfinfo,0,sizeof(sfinfo));
  sfinfo.samplerate = sample_rate;
  sfinfo.frames = sample_count;
  sfinfo.channels = 1; // For mono
  sfinfo.format = (SF_FORMAT_WAV | SF_FORMAT_PCM_24);
  if (!(file = sf_open("pulsar.wav",SFM_WRITE,&sfinfo)))
    {
      printf("Unable to open pulsar.wav file\n");
      exit(1);
    }
  
    


  
  
  fout = fopen("dedisp.dat","w");
  for (i=0;i<npts;i+=bs)
    {
      final =0.0;
      for (j=0;j<bs;j++)
	final+=dedispTime[i+j];
      fprintf(fout,"%d %g\n",i,final/(double)bs);
    }
  fclose(fout);
  printf("Making sound file\n");

  // Convolving
  printf("Convolving\n");
  for (i=0;i<sample_count;i++)
    {
      convolve[i] = 0.0;
      nc=0;
      if (nConvolve==0)
	convolve[i] = dedispTime[i];
      else
	{
	  for (j=-nConvolve/2.;j<nConvolve/2.;j++)
	    {
	      if (i+j > 0 && i+j < npts)
		{
		  convolve[i] += dedispTime[i+j]*exp(-pow(j,2)/2./pow(convolveW,2));
		  nc++;
		}
	    }
	  convolve[i]/=(double)nc;
	}
    }
  printf("Complete convolving\n");
  
  min=max=convolve[0];
  for (i=0;i<sample_count;i++)
    {   
      if (min > convolve[i]) min = convolve[i];
      if (max < convolve[i]) max = convolve[i];
    }

  fout = fopen("sound.dat","w");
  for (i=0;i<sample_count;i++)
    {
      ival = (int)((double)i/(double)sample_rate/dSet->head->tsamp);
      //      buffer[i] = (int)(AMPLITUDE*2*((dedispTime[ival]-min)/(max-min)-0.5));
      amp = (AMPLITUDE*(2*(convolve[ival]-min)/(max-min)));
      if (amp > AMPLITUDE) amp = AMPLITUDE;
      buffer[i] = (int)(TKgaussDev(&seed)*(amp)/2.);
      //      buffer[i] = (int)amp;
      printf("Have %d %g %g %g\n",buffer[i],convolve[ival],min,max);
      fprintf(fout,"%d %d %g %g %g %d %g\n",i,ival,convolve[i],convolve[ival],dedispTime[ival],buffer[i],amp);
    }
  fclose(fout);
  if (sf_write_int(file,buffer,sfinfo.channels*sample_count) != sfinfo.channels*sample_count)
    puts(sf_strerror(file));
  
  sf_close(file);
  //  pfitsCloseFile(dSet,debug);
  free(convolve);
  free(buffer);
  free(dSet);
  free(data);
  free(data1sub);
  free(dedispTime);
}




