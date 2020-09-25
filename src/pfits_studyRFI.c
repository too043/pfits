// gcc -lm -o pfits_studyRFI pfits_studyRFI.c -lcpgplot -lcfitsio pfits_loader.c pfits_setup.c
//
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

void drawColourMap(dSetStruct *dSet,int pol,float fmin,float fmax,float tmin,float tmax,int plotTransmitters,char *grDev);
void showTransmitter(float fc,float bw,char *label,float t0,float t1);
void makeZeroDM(dSetStruct *dSet,int pol,float fmin,float fmax,float tmin,float tmax,int scls,char *outName);

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  char fname[1024];
  int debug=0;
  int i;
  int colourMap=0;
  int zeroDM =0;
  float fmin = 704;
  float fmax = 4032;
  float tmin = 0;
  float tmax = 0.5;
  int plotTransmitters=0;
  char grDev[128]="/xs";
  char outName[1024]="DM0.dat";
  int scls=1;
  
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-cm")==0)
	colourMap=1;
      else if (strcmp(argv[i],"-no_scales")==0)
	scls=0;
      else if (strcasecmp(argv[i],"-zerodm")==0)
	zeroDM=1;
      else if (strcmp(argv[i],"-f0")==0)
	sscanf(argv[++i],"%f",&fmin);
      else if (strcmp(argv[i],"-f1")==0)
	sscanf(argv[++i],"%f",&fmax);
      else if (strcmp(argv[i],"-t0")==0)
	sscanf(argv[++i],"%f",&tmin);
      else if (strcmp(argv[i],"-t1")==0)
	sscanf(argv[++i],"%f",&tmax);
      else if (strcmp(argv[i],"-transmitters")==0)
	plotTransmitters=1;
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outName,argv[++i]);
      else if (strcmp(argv[i],"-g")==0)
	strcpy(grDev,argv[++i]);
    }
  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));
  initialise(&dSet,debug);
  setFilename(fname,dSet,debug);

  pfitsOpenFile(dSet,debug);
  printf("Loading header\n");
  pfitsLoadHeader(dSet,debug);
  printf("Header loaded\n");
  if (colourMap==1)
    drawColourMap(dSet,0,fmin,fmax,tmin,tmax,plotTransmitters,grDev);

  // Make zero DM results
  if (zeroDM==1)
    makeZeroDM(dSet,0,fmin,fmax,tmin,tmax,scls,outName);

  
  cpgend();
  

  //  pfitsCloseFile(dSet,debug);
  free(dSet);

}

void makeZeroDM(dSetStruct *dSet,int pol,float fmin,float fmax,float tmin,float tmax,int scls,char *outName)
{
  int nchan = dSet->head->nchan;
  int nbin  = dSet->head->nbin;
  int nsblk = dSet->head->nsblk;
  int i,j;  
  float *dataVals;
  double sumVal;
  FILE *fout;
  long int nSamples;
  int nTimeSamples;
  int nFreqSamples;
  int debug=0;
  double freq;
  int nc=0;
  int sub=0;
  double tsamp = dSet->head->tsamp;
  
  nTimeSamples = (int)((tmax-tmin)/tsamp+0.5);
  if (!(dataVals = (float *)malloc(sizeof(float)*nchan*nTimeSamples)))
    {
      printf("ERROR: unable to allocate memory\n");
      printf("Trying to allocate %d channels and %d time samples\n",nchan,nTimeSamples);
      exit(1);
    }
  printf("Nchan = %d\n",nchan);
  printf("nTSamples = %d\n",nTimeSamples);
  pfits_read1pol_float(dataVals,0,dSet,tmin,tmax,3,&nSamples,&nTimeSamples,&nFreqSamples,debug,scls);
  printf("nTSamples now = %d\n",nTimeSamples);

  fout = fopen(outName,"w");
            
  for (i=0;i<nTimeSamples;i++)
    {
      sumVal=0;
      nc=0;
      for (j=0;j<nFreqSamples;j++)
	{
	  freq = dSet->head->chanFreq[j];
	  //	      if (freq > 2435 && freq < 2440)
	  //	      if (freq > 835 && freq < 840)
	  if (freq > fmin && freq < fmax)
	    {
	      sumVal += dataVals[i*nFreqSamples+j];
	      nc++;
	    }
	}
      fprintf(fout,"%g %g\n",sub*nsblk*dSet->head->tsamp+i*dSet->head->tsamp,sumVal/(float)nc);
    }
  free(dataVals);
  fclose(fout);
}


void drawColourMap(dSetStruct *dSet,int pol,float fmin,float fmax,float tmin,float tmax,int plotTransmitters,char *grDev)
{
  float tr[6];
  int nchan = dSet->head->nchan;
  int nbin = dSet->head->nbin;
  int nsblk = dSet->head->nsblk;
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
  float fx[2],fy[2];
  
  float mx,my,mx2,my2;
  char key;
  int plot=1;
  int plotExpected=-1;
  int plotY=2;
  int zap;

  int timeOff=0;
  
  char  plotStr[1024];
  float plot_minX,plot_minY,plot_maxX,plot_maxY;
  long int nSamples;
  int nTimeSamples;
  int nFreqSamples;
  double tsamp = dSet->head->tsamp;
  
  nTimeSamples = (int)((tmax-tmin)/tsamp+0.5);
  plotArr = (float *)malloc(sizeof(float)*nchan*nTimeSamples);
  printf("nTSamples = %d\n",nTimeSamples);
  pfits_read1pol_float(plotArr,0,dSet,tmin,tmax,3,&nSamples,&nTimeSamples,&nFreqSamples,debug,0);
  printf("nTSamples now = %d\n",nTimeSamples);
  //  for (i=0;i<nSamples;i++)
  //    printf("%g\n",plotArr[i]);

  //  pfits_read1pol_float(plotArr[j],0,dSet[j],subintCount,subintCount,1,&nSamples,&nTimeSamples,&nFreqSamples,debug,0);
  tr[0] = tmin;  tr[1] = 0;  tr[2] = dSet->head->tsamp;
  tr[3] = dSet->head->chanFreq[0];  tr[4] = dSet->head->chanbw;  tr[5] = 0;

  {
    printf("Timeoff = %d\n",timeOff);
    
    cpgbeg(0,grDev,1,1);
    cpgsci(1.2);
    cpgscf(2);
    cpgslw(2);

    cpgask(0);
    cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
    
    plot_minX=tmin;
    plot_maxX=tmax;
    
    plot_minY=fmin;
    plot_maxY=fmax;
    
    cpgenv(plot_minX,plot_maxX,plot_minY,plot_maxY,0,1);
    cpglab("Time (seconds)","Frequency (MHz)","");
    //    cpgimag(plotArr,nchan,nTimeSamples,1,nchan,1,nTimeSamples,0,3,tr);
    cpgimag(plotArr,nchan,nTimeSamples,1,nchan,1,nTimeSamples,0,255,tr);
    //    cpggray(plotArr,nchan,nsblk,1,nchan,1,nsblk,0,3,tr);
  }

  if (plotTransmitters==1)
    {
      showTransmitter(763,10,"Optus",tmin,tmax);
      showTransmitter(778,20,"Telstra",tmin,tmax);
	showTransmitter(847.8,0.4,"Parkes radio",tmin,tmax);
	showTransmitter(849.5,0.2,"Police",tmin,tmax);
	showTransmitter(872.5,5,"Vodafone",tmin,tmax);
	showTransmitter(882.5,14.9,"Telstra",tmin,tmax);
	showTransmitter(905,3,"Amateur radio/walkie talkie",tmin,tmax);
	showTransmitter(930.5,0.2,"Police",tmin,tmax);
	showTransmitter(947.6,8.4,"Optus",tmin,tmax);
	showTransmitter(956.2,5,"Vodafone",tmin,tmax);
	showTransmitter(1018,0.5,"PKS airport",tmin,tmax);
	showTransmitter(1575.42,32.736,"Beidou",tmin,tmax);
	showTransmitter(1176.45,20.46,"Beidou",tmin,tmax);
	showTransmitter(1268.52,20.46,"Beidou",tmin,tmax);
	showTransmitter(1542,34,"Inmarsat",tmin,tmax);
	showTransmitter(1643.25,33.5,"Inmarsat",tmin,tmax);
	showTransmitter(1575,32,"Galileo",tmin,tmax);
	showTransmitter(1280,40,"Galileo",tmin,tmax);
	showTransmitter(1189,50,"Galileo",tmin,tmax);
	showTransmitter(1575.42,15.456,"GPS",tmin,tmax);
	showTransmitter(1227.6,11,"GPS",tmin,tmax);
	showTransmitter(1176.45,12.5,"GPS",tmin,tmax);
	showTransmitter(1626.25,0.5,"Iridium",tmin,tmax);
	showTransmitter(1542.5,35,"Thurya",tmin,tmax);
	showTransmitter(1597.21875,16.3125,"Glonass",tmin,tmax);
	showTransmitter(1245.78125,5.6875,"Glonass",tmin,tmax);
	showTransmitter(1815,20,"Telstra",tmin,tmax);
	showTransmitter(1835,20,"Telstra",tmin,tmax);
	showTransmitter(1857.5,25,"Optus",tmin,tmax);
	showTransmitter(2115,10,"Vodafone",tmin,tmax);
	showTransmitter(2167.5,4.8,"Vodafone",tmin,tmax);
	showTransmitter(2122.5,5,"Telstra",tmin,tmax);
	showTransmitter(2127.5,5,"Telstra",tmin,tmax);
	showTransmitter(2142.5,5,"Optus",tmin,tmax);
	showTransmitter(2147.5,5,"Optus",tmin,tmax);
	showTransmitter(2152.5,5,"Optus",tmin,tmax);
	showTransmitter(2160.0,9.9,"Telstra",tmin,tmax);
	showTransmitter(2312,19.2,"NBN",tmin,tmax);
	showTransmitter(2331.2,19.2,"NBN",tmin,tmax);
	showTransmitter(2350.4,19.2,"NBN",tmin,tmax);
	showTransmitter(2369.6,19.2,"NBN",tmin,tmax);
	showTransmitter(2386.7,15,"NBN",tmin,tmax);

	showTransmitter(2650.0,40,"Telstra",tmin,tmax);
	showTransmitter(2680.0,20,"Optus",tmin,tmax);

	showTransmitter(3455.0,19.9,"NBN",tmin,tmax);
	showTransmitter(3560.0,19.9,"NBN",tmin,tmax);
	
	showTransmitter(921.5,13,"WiFi 802.11",tmin,tmax);
	showTransmitter(2442.0,82,"WiFi",tmin,tmax);

	showTransmitter(1090,50./1000.,"ADS-B",tmin,tmax);
    }

  free(plotArr);
}

void showTransmitter(float fc,float bw,char *label,float t0,float t1)
{
  float fx[2],fy[2];

  fx[0] = t0; fx[1] = t1+(t1-t0)*0.02;
  fy[0] = fy[1] = fc-bw/2.;
  cpgslw(4);
  cpgsci(4);
  cpgline(2,fx,fy);
  fy[0] = fy[1] = fc+bw/2.;
  cpgline(2,fx,fy);
  //  cpgsci(5);
  cpgslw(2);
  cpgptxt(fx[1],fc,90,0.5,label);
  cpgsci(1);
  cpgslw(1);

}
