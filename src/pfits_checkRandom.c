// gcc -lm -o pfits_checkRandom pfits_checkRandom.c -lcpgplot -lcfitsio pfits_loader.c pfits_setup.c
//
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

//void drawColourMap(dSetStruct *dSet,int pol);
void checkRandom(dSetStruct *dSet,int nsub,int dk);

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  char fname[1024];
  int debug=0;
  int i;
  int nsub=0;
  int dk=512; 
  
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-nsub")==0)
	sscanf(argv[++i],"%d",&nsub);
      else if (strcmp(argv[i],"-nsamp")==0)
	sscanf(argv[++i],"%d",&dk);
      
    }
  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));
  initialise(&dSet,debug);
  setFilename(fname,dSet,debug);

  pfitsOpenFile(dSet,debug);
  printf("Loading header\n");
  pfitsLoadHeader(dSet,debug);
  checkRandom(dSet,nsub,dk);
  
  cpgend();
  

  //  pfitsCloseFile(dSet,debug);
  free(dSet);

}


void checkRandom(dSetStruct *dSet,int nsub,int dk)
{
  float *val;
  int nchan = dSet->head->nchan;
  int nsblk = dSet->head->nsblk;
  int i,j,k,kk;
  long nSamples;
  int nTimeSamples,nFreqSamples;
  int debug=0;
  char fname[1024];
  FILE *fout;
  long t=0;
  int n1=0;

  printf("dk = %d\n",dk);
  
  val = (float *)malloc(sizeof(float)*nchan*nsblk);
  for (i=nsub;i<nsub+1;i++)
    {
      pfits_read1pol_float(val,0,dSet,i,i,1,&nSamples,&nTimeSamples,&nFreqSamples,debug,0);
      printf("Loaded %d %d %d\n",i,nTimeSamples,nFreqSamples);
      for (k=0;k<nTimeSamples;k+=dk)
	{
	  sprintf(fname,"random.%05ld.dat",t);
	  fout = fopen(fname,"w");
	  n1=0;
	  for (j=0;j<dk;j++)
	    {
	      for (kk=0;kk<nFreqSamples;kk++)
		{
		  if ((int)val[(k+j)*nFreqSamples+kk] == 1)
		    n1++;
		  //		  printf("Trying: %d\n",(k+j)*nFreqSamples+kk);
		  fprintf(fout,"%d",(int)val[(k+j)*nFreqSamples+kk]);
		  //		  printf("Written\n");
		}
	    }
	  printf("%s k = %d n1 = %d %g\n",fname,k,n1,(double)n1/(nFreqSamples*dk));
	  //	  fprintf(fout,"\n");

	  fclose(fout);
	  t++;
	}
    }
  free(val);
  
}


/*
void drawColourMap(dSetStruct *dSet,int pol)
{
  float tr[6];
  int nchan = dSet->head->nchan;
  int nbin = dSet->head->nbin;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  char title[128];
  int i,j;
  int debug=0;
  int colnum;
  int subint=0;
  int status=0;
  short int nVal=0;
  short int *sVals;
  int initflag=0;
  float minV,maxV;
  float minV_scl,maxV_scl;
  float nfVal = 0;
  float *plotArr;
  
  float mx,my,mx2,my2;
  char key;
  int plot=1;
  int plotExpected=-1;
  int plotY=2;
  int zap;

  int timeOff=0;
  
  char  plotStr[1024];
  float plot_minX,plot_minY,plot_maxX,plot_maxY;
  int nsblk=8192;
  long int nSamples;
  int nTimeSamples;
  int nFreqSamples;
  long nplot=0;
  int dt = 64; // This is the offset in samples
  long tmax = 5*8192;
  long plotPts=2048;
  int s0=420;
  int s1=430;
  printf("nchan = %d\n",nchan);
  
  plotArr = (float *)malloc(sizeof(float)*nchan*nsblk*((s1-s0)+1));
  pfits_read1pol_float(plotArr,0,dSet,s0,s1,2,&nSamples,&nTimeSamples,&nFreqSamples,debug,0);
  printf("Loaded %d time samples\n",nTimeSamples);
  //  pfits_read1pol_float(plotArr[j],0,dSet[j],subintCount,subintCount,1,&nSamples,&nTimeSamples,&nFreqSamples,debug,0);
  tr[0] = 0;  tr[1] = 0;  tr[2] = dSet->head->tsamp;
  tr[3] = dSet->head->chanFreq[0];  tr[4] = dSet->head->chanbw;  tr[5] = 0;

  for (timeOff=0;timeOff<tmax;timeOff+=dt)
    {
      printf("Timeoff = %d\n",timeOff);
      sprintf(plotStr,"plot%5.5d.png/png",nplot);
      
      cpgbeg(0,plotStr,1,1);
      cpgask(0);
      cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
      
      //      plot_minX=timeOff*3;
      //      plot_maxX=timeOff*3+128;
      plot_minX=timeOff*dSet->head->tsamp;
      plot_maxX=(timeOff+plotPts)*dSet->head->tsamp;

      plot_minY=dSet->head->chanFreq[0];
      //      plot_maxY=96;
      plot_maxY=dSet->head->chanFreq[nchan-1];
      
      cpgenv(plot_minX,plot_maxX,plot_minY,plot_maxY,0,1);
      cpglab("Time (seconds)","Frequency (channel)","");
      //      cpgimag(plotArr,nchan,nsblk,1,nchan,1,nsblk,0,1,tr);
      //      cpgimag(plotArr,nchan,nTimeSamples,1,nchan,1,nTimeSamples,-2.5,2.5,tr);
      cpggray(plotArr,nchan,nTimeSamples,1,nchan,1,nTimeSamples,0,1,tr);
      nplot++;
    }
  free(plotArr);
}
*/
