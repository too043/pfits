//  gcc -lm -o pfits_read_1bitExtraction pfits_read_1bitExtraction.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -lpgplot

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

#define FILE_LEN 128

int main(int argc,char *argv[])
{
  FILE *fin;
  int i,j;
  float t0_plot,t1_plot;
  char  fname[1024];
  char  inName[FILE_LEN];
  int   nchan;
  float *chanFreq;
  float tsamp;
  int   stt_imjd;
  float stt_smjd;
  float stt_offs;
  float t0,t1;
  unsigned char *cVals;
  float *vals;
  long nSamples;
  int nbits,npol;
  int samplesperbyte;
  float tr[6];
  
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
    }

  fin = fopen(fname,"rb");
  fread(inName,sizeof(char),FILE_LEN,fin);
  fread(&nchan,sizeof(int),1,fin);
  //  printf("inName = %s, nchan = %d\n",inName,nchan);
  chanFreq = (float *)malloc(sizeof(float)*nchan);
  fread(chanFreq,sizeof(int),nchan,fin);
  fread(&tsamp,sizeof(float),1,fin);
  fread(&stt_imjd,sizeof(int),1,fin);
  fread(&stt_smjd,sizeof(float),1,fin);
  fread(&stt_offs,sizeof(float),1,fin);
  fread(&t0,sizeof(float),1,fin);
  fread(&t1,sizeof(float),1,fin);
  fread(&nbits,sizeof(int),1,fin);
  fread(&npol,sizeof(int),1,fin);

  samplesperbyte = 8/nbits;
  nSamples = (t1-t0)/tsamp;
  printf("Number of samples = %d %g %g %g\n",nSamples,t1,t0,tsamp);
  cVals = (unsigned char *)malloc(sizeof(unsigned char)*nSamples*nchan/samplesperbyte);
  fread(cVals,sizeof(unsigned char),nSamples*nchan/samplesperbyte,fin);
  printf("Read cVals\n");

  vals = (float *)malloc(sizeof(float)*nSamples*nchan);
  printf("Allocated vals\n");
  for (i=0;i<nSamples;i++)
    {
      pfits_bytesToFloats(samplesperbyte,nchan,cVals+i*nchan/samplesperbyte,vals+i*nchan);
    }
  fclose(fin);
  printf("Making plot\n");
  tr[0] = 0;
  tr[1] = 0;
  tr[2] = 1;
  tr[3] = 0;
  tr[4] = 1;
  tr[5] = 0;

  cpgbeg(0,"/xs",1,1);
  cpgenv(0,nSamples,0,nchan,0,1);
  cpggray(vals,nchan,nSamples,1,nchan,1,nSamples,0,3,tr);
  cpgend();

  free(chanFreq);  
  free(cVals);
  free(vals);
}
