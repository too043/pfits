// gcc -lm -o pfits_makeMovie pfits_makeMovie.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -lpgplot

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>


void drawColourMap(dSetStruct *dSet,int pol,int outTxt);

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  char fname[1024];
  int debug=0;
  int outTxt=0;
  int i;

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-txt")==0) // THIS DOESN'T DO ANYTHING YET
	outTxt=1;
    }
  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));
  initialise(&dSet,debug);
  setFilename(fname,dSet,debug);

  pfitsOpenFile(dSet,debug);
  printf("Loading header\n");
  pfitsLoadHeader(dSet,debug);
  drawColourMap(dSet,0,outTxt);
  if (outTxt==0)
    cpgend();
  

  //  pfitsCloseFile(dSet,debug);
  free(dSet);

}




void drawColourMap(dSetStruct *dSet,int pol,int outTxt)
{
  float tr[6];
  int nchan = dSet->head->nchan;
  int nbin = dSet->head->nbin;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  char title[128];
  float *plotArr;
  float *plotArr_scl;
  float *plotArr_meanScl;
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
  float *dat_offs;
  float *dat_scl;
  float *mean_scl;
  float nfVal = 0;
  float *expectedX,*expectedY;
  
  float mx,my,mx2,my2;
  char key;
  int plot=1;
  int plotExpected=-1;
  int plotY=2;
  int zap;

  int freqOff=0;
  
  char  plotStr[1024];
  float plot_minX,plot_minY,plot_maxX,plot_maxY;
  
  printf("nchan = %d, nbin = %d\n",nchan,nbin);
  
  plotArr = (float *)malloc(sizeof(float)*nbin*nchan);
  plotArr_scl = (float *)malloc(sizeof(float)*nbin*nchan);
  expectedX = (float *)malloc(sizeof(float)*nchan);
  expectedY = (float *)malloc(sizeof(float)*nchan);
  dat_offs = (float *)malloc(sizeof(float)*nchan);
  dat_scl = (float *)malloc(sizeof(float)*nchan);
  mean_scl = (float *)malloc(sizeof(float)*nchan);
  
  sVals = (short int *)malloc(sizeof(short int)*nbin*nchan);

  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum,&status);
  fits_read_col(dSet->fp,TFLOAT,colnum,subint+1,1,nchan,&nfVal,dat_offs,&initflag,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum,&status);
  fits_read_col(dSet->fp,TFLOAT,colnum,subint+1,1,nchan,&nfVal,dat_scl,&initflag,&status);
  
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum,&status);

  fits_read_col(dSet->fp,TSHORT,colnum,subint+1,1,nchan*nbin,&nVal,sVals,&initflag,&status);

  
  for (i=0;i<nchan;i++)
    {
      for (j=0;j<nbin;j++)
	{
	  plotArr[i*nbin+j] = (float)sVals[i*nbin+j];
	  plotArr_scl[i*nbin+j] = plotArr[i*nbin+j]*dat_scl[i]+dat_offs[i];
	  if (i==0 && j==0)
	    {
	      minV = maxV = plotArr[i*nbin+j];
	      minV_scl = maxV_scl = plotArr_scl[i*nbin+j];
	    }
	  else
	    {
	      if (minV > plotArr[i*nbin+j]) minV = plotArr[i*nbin+j];
	      if (maxV < plotArr[i*nbin+j]) maxV = plotArr[i*nbin+j];
	      if (minV_scl > plotArr_scl[i*nbin+j]) minV_scl = plotArr_scl[i*nbin+j];
	      if (maxV_scl < plotArr_scl[i*nbin+j]) maxV_scl = plotArr_scl[i*nbin+j];
	    }
	}
    }
  printf("min/max = %g %g\n",minV_scl,maxV_scl);
  tr[0] = 0;  tr[1] = 64e-6;  tr[2] = 0;
  tr[3] = 4032;  tr[4] = 0;  tr[5] = -(4032.0-704.0)/dSet->head->nchan;

  for (freqOff=0;freqOff<3200;freqOff++)
    {
      sprintf(plotStr,"plot%5.5d.png/png",freqOff);
      
      cpgbeg(0,plotStr,1,1);
      cpgask(0);
      cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
      
      plot_minX=0;
      plot_minY=704;
      plot_maxX=0.25;
      plot_maxY=800;
      
      cpgenv(plot_minX,plot_maxX,704+freqOff,704+freqOff+128,0,1);
      cpglab("Time (seconds)","Frequency (MHz)","");
      cpgimag(plotArr,nbin,nchan,1,nbin,1,nchan,minV,maxV,tr);
      
    }
  
  free(plotArr);
  free(plotArr_scl);
  free(sVals);
  free(dat_offs);
  free(dat_scl);
  free(mean_scl);
  free(expectedX);
  free(expectedY);
}
