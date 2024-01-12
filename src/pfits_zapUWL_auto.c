// Code to zap pulse profiles from the UWL receiver
//
//gcc -lm -o pfits_zapUWL_auto pfits_zapUWL_auto.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -O3

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include <stdint.h>
#include "fitsio.h"

#define VNUM "v1.1"

void help();
void get_dm_tbin(dSetStruct *dSet,double *dm,double *tbin);
void saveFile(float *datWts, int nsub, int nchan,dSetStruct *dSet);
void autoZapTransmitters(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL);
void autoZapSatellites(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL);
void autoZapWiFi(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll);
void autoZapHandsets(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL);
void autoZapAircraft(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL);
void autoZapDigitisers(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL);
void autoZapUnexplained(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL);
void autoZapBandEdge(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,float bandEdge);

void help()
{
  printf("\npfits_zapUWL_auto    %s\n",VNUM);
  printf("PFITS version:       %s\n",SOFTWARE_VER);
  printf("Author:              George Hobbs\n");
  printf("Software to automatically zap known RFI in the UWL band\n");

  printf("\nCommand line arguments:\n\n");
  printf("-h                   This help\n");
  printf("-f                   PSRFITS filename to zap\n");
  printf("-tbin                TBIN\n");
  printf("-zapAll              Remove all known RFI\n");
  printf("-keepSL              Do not remove known spectral lines (use in conjunction with -zapAll)\n");

  printf("\nExample 1: brute force zap of a pulsar observation\n\n");
  printf("pfits_zapUWLauto -f <psrfits_file> -zapAll\n\n");

  printf("\nExample 2: zap but preserve known spectral lines\n\n");
  printf("pfits_zapUWLauto -f <psrfits_file> -zapAll -keepSL\n\n");
  
  exit(1);
}


int main(int argc,char *argv[])
{
  float tsub;
  float chanbw;
  double tbin;
  double dm;
  
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,j,i0;

  short int *sval,*svalPsr;
  short int n_val =0;
  float n_fval=0;
  float *datScl;
  float *datOffs;
  float *datWts;
  float *fChan,fref;
  int initflag = 0;
  long nchan,nsub,nbin,npol;
  int colnum_data;
  int colnum_datOffs;
  int colnum_datScl;
  int colnum_datWts;
  long ii,jj,kk;
  float f0,f1;
  int polNum = 0;
  int zapAll=0;
  int keepSL=0;
  float bandEdge = 10;

  // Searching for saturation
  float **spectrum;
  float *specMean1,*specMean2,*specMean3;
  float *sort1,*sort2,*sort3;
  float mean0,mean1,mean2,mean3;
  float val,freq;
  float low1,low2,low3;
  int onlyScl=1;
  int n1,n2,n3;
  float sigma1=2; // Should calculate these threshold values better
  float sigma2=5;
  float sigma3=5;

  if (argc==1)
    help();

  // Initialise everything
  initialise(&dSet,debug);

  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-h")==0)
	help();
      else if (strcmp(argv[i],"-f")==0)
	{setFilename(argv[++i],dSet,debug);}
      else if (strcmp(argv[i],"-tbin")==0)
	sscanf(argv[++i],"%lf",&tbin);
      else if (strcasecmp(argv[i],"-zapAll")==0)
	zapAll=1;
      else if (strcasecmp(argv[i],"-keepSL")==0)
	keepSL=1;
    }
       
  
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  nchan = dSet->head->nchan;
  nbin  = dSet->head->nbin;
  nsub  = dSet->head->nsub;
  npol  = dSet->head->npol;
  chanbw = dSet->head->chanFreq[1]-dSet->head->chanFreq[0];
  
  get_dm_tbin(dSet,&dm,&tbin);
  
  printf("Loaded header\n");
  printf("Number of channels = %d\n",nchan);
  printf("Number of bins = %d\n",nbin);
  printf("Number of sub-integrations = %d\n",nsub);
  printf("Number of polarisations = %d\n",npol);
  printf("Dispersion measure = %g\n",dm);
  
  fChan    = (float *)malloc(sizeof(float)*nchan);
  sval     = (short int *)calloc(nchan*nbin,sizeof(short int));
  spectrum = (float **)malloc(sizeof(float *)*nsub);
  specMean1 = (float *)malloc(sizeof(float)*nsub);
  specMean2 = (float *)malloc(sizeof(float)*nsub);
  specMean3 = (float *)malloc(sizeof(float)*nsub);
  for (i=0;i<nsub;i++)
    spectrum[i] = (float *)malloc(sizeof(float)*nchan);

  datScl  = (float *)malloc(sizeof(float)*nchan*4);
  datOffs = (float *)malloc(sizeof(float)*nchan*4);
  datWts  = (float *)malloc(sizeof(float)*nchan*nsub);
  
  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_data,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);

  for (ii=0;ii<nsub;ii++)
    {      
      fits_read_col(dSet->fp,TFLOAT,colnum_datWts,ii+1,1,nchan,&n_fval,datWts+ii*nchan,&initflag,&status);
      // Automatic zapping
      autoZapTransmitters(dSet,datWts+ii*nchan,nchan,chanbw,zapAll,keepSL);
      autoZapDigitisers(dSet,datWts+ii*nchan,nchan,chanbw,zapAll,keepSL);
      autoZapAircraft(dSet,datWts+ii*nchan,nchan,chanbw,zapAll,keepSL);
      autoZapSatellites(dSet,datWts+ii*nchan,nchan,chanbw,zapAll,keepSL);
      autoZapWiFi(dSet,datWts+ii*nchan,nchan,chanbw,zapAll);
      autoZapUnexplained(dSet,datWts+ii*nchan,nchan,chanbw,zapAll,keepSL);
      autoZapHandsets(dSet,datWts+ii*nchan,nchan,chanbw,zapAll,keepSL);
      autoZapBandEdge(dSet,datWts+ii*nchan,nchan,chanbw,zapAll,bandEdge);

      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,ii+1,1,nchan*npol,&n_fval,datScl,&initflag,&status);
      mean0=mean1=mean2=mean3 = 0.0;
      n1=n2=n3=0;
      for (i=0;i<nchan;i++)
	{
	  freq = dSet->head->chanFreq[i];
	  val = datScl[i]+datScl[i+nchan];
	  spectrum[ii][i] = val;
	  mean0+=val;
	  //	    printf("wts = %g\n",datWts[i]);
	  if (datWts[ii*nchan+i] > 0)
	    {
	      if (freq > 703 && freq < 1344)
		{mean1+=val;n1++;}
	      else if (freq >= 1344 && freq <2368)
		{mean2+=val;n2++;}
	      else
		{mean3+=val;n3++;}
	    }
	  //	    printf("%.5f %g\n",dSet->head->chanFreq[i],spectrum);
	  //	    printf("%.5f %g %g\n",dSet->head->chanFreq[i],datScl[i],
	  //datScl[i+nchan],datScl[i+2*nchan],datScl[i+3*nchan]);
	}
      mean0/=(float)nchan;
      mean1/=(float)n1;
      mean2/=(float)n2;
      mean3/=(float)n3;
      specMean1[ii] = mean1;
      specMean2[ii] = mean2;
      specMean3[ii] = mean3;
      if (ii==i0)
	{
	  low1=mean1;
	  low2=mean2;
	  low3=mean3;
	}
      else
	{
	  if (mean1 < low1) low1 = mean1;
	  if (mean2 < low2) low2 = mean2;
	  if (mean3 < low3) low3 = mean3;
	}
      printf("DAT_SCL %d %g %g %g %g\n",ii,mean0,mean1,mean2,mean3);
    }
  
  // Should calculate the mean and sdev of the lowest few points (knowing that RFI is always positive)
  // in the specMean values

  
  // Now check for band saturation
  // Should check first subint carefully -- perhaps start from 1 not 0
  for (ii=0;ii<nsub;ii++)
    {
      if (specMean1[ii] > sigma1*low1)
	{
	  printf("Likely saturation in subband %d in RF band 1 (%g, %g)\n",ii,specMean1[ii],low1);
	  for (i=0;i<nchan;i++)
	    {
	      freq = dSet->head->chanFreq[i];
	      if (freq > 703 && freq < 1344)
		datWts[ii*nchan+i]=0;
	    }
	}
      if (specMean2[ii] > sigma2*low2)
	{
	  printf("Likely saturation in subband %d in RF band 2 (%g, %g)\n",ii,specMean2[ii],low2);
	  for (i=0;i<nchan;i++)
	    {
	      freq = dSet->head->chanFreq[i];
	      if (freq >= 1344 && freq <2368)
		datWts[ii*nchan+i]=0;
	    }
	}
      if (specMean3[ii] > sigma3*low3)
	{
	  printf("Likely saturation in subband %d in RF band 3 (%g, %g)\n",ii,specMean3[ii],low3);
	  for (i=0;i<nchan;i++)
	    {
	      freq = dSet->head->chanFreq[i];
	      if (freq > 2368)
		datWts[ii*nchan+i]=0;
	    }
	}
    }
  

  
  // Save the file
  printf("Flagging complete: saving output data file\n");
  saveFile(datWts,nsub,nchan,dSet);
  printf("Complete\n");
  
  // De-allocate the memory
  //  pfitsCloseFile(dSet,debug);
  deallocateMemory(&dSet,debug);

  //  free(dVals);
  //  free(sdevArr);
  free(sval);
  free(fChan);
  free(datScl);
  free(datOffs);
  free(datWts);
  for (i=0;i<nsub;i++)
    free(spectrum[i]);
  free(spectrum);
  free(specMean1); free(specMean2); free(specMean3);

}

void autoZapBandEdge(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,float bandEdge)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;

  //
  // Fixed mobile transmission towers
  // We see persistent emission from these at all telescope pointing angles
  //
  for (i=0;i<27;i++)
    {
      f0[nT] = 704+128*i-bandEdge/2.; f1[nT++] = 704+128*i+bandEdge/2.; 
    }
  
  printf("Number of band edges being removed = %d\n",nT);
  

  for (i=0;i<nchan;i++)
    {
      freq = dSet->head->chanFreq[i];
      for (j=0;j<nT;j++)
	{
	  // Should probably add/subtract one channel bandwidth here to ensure at least
	  // one channel is zapped
	  if (freq+chanbw/2. >= f0[j] && freq-chanbw/2. <= f1[j])
	    {
	      datWts[i] = 0;		  	      
	      nzap++;
	    }
	}
    }
  printf(" ... have flagged %d channels\n",nzap);
  
}

void autoZapTransmitters(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;

  //
  // Fixed mobile transmission towers
  // We see persistent emission from these at all telescope pointing angles
  // NOTE: we can check against INSPECTA's spectral lines list
  //
  f0[nT] = 758; f1[nT++] = 768;       //  Optus  
  f0[nT] = 768; f1[nT++] = 788;       //  Telstra
  f0[nT] = 869.95; f1[nT++] = 875.05; //  Vodafone
  f0[nT] = 875.05; f1[nT++] = 889.95; //  Telstra
  f0[nT] = 915; f1[nT++] = 928;       //  Agricultural 920 MHz RFI
  f0[nT] = 943.4; f1[nT++] = 951.8;   //  Optus
  f0[nT] = 953.7; f1[nT++] = 958.7;   //  Vodafone
  f0[nT] = 1081.7; f1[nT++] = 1086.7; // Aliased signal
  if (keepSL==1)
    {
      // adjusted to account for spectral line at 1720.530
      f0[nT] = 1725; f1[nT++] = 1736.2;   // Aliased signal
    }
  else
    {
      f0[nT] = 1720; f1[nT++] = 1736.2;   // Aliased signal
    }
  f0[nT] = 1805; f1[nT++] = 1825;     //  Telstra
  if (keepSL==1)
    {
      // adjusted to account for spectral line at 1849.634
      f0[nT] = 1855; f1[nT++] = 1865;     //  Optus
    }
  else
    {
      f0[nT] = 1845; f1[nT++] = 1865;     //  Optus
    }
  f0[nT] = 1973.7; f1[nT++] = 1991.0; // Aliased signal
  f0[nT] = 2110; f1[nT++] = 2120;     //  Vodafone
  f0[nT] = 2140; f1[nT++] = 2145;     //  Optus
  f0[nT] = 2145; f1[nT++] = 2150;     //  Optus
  f0[nT] = 2164.9; f1[nT++] = 2170.1; //  Vodafone
  f0[nT] = 2302.05; f1[nT++] = 2321.95; //  NBN
  f0[nT] = 2322.05; f1[nT++] = 2341.95; //  NBN
  f0[nT] = 2342.05; f1[nT++] = 2361.95; //  NBN
  f0[nT] = 2362.05; f1[nT++] = 2381.95; //  NBN
  f0[nT] = 2487.00; f1[nT++] = 2496.00; //  NBN alias
  f0[nT] = 2670; f1[nT++] = 2690;     //  Optus
  f0[nT] = 3445.05; f1[nT++] = 3464.95; //  NBN
  f0[nT] = 3550.05; f1[nT++] = 3569.95; //  NBN

  if (zapAll==1)
    {
      f0[nT] = 804.4; f1[nT++] = 804.6;               //  NSW Police Force
      f0[nT] = 2152.5-2.5; f1[nT++] = 2152.5+2.5;     //  Telstra - not always on
      f0[nT] = 847.8-0.2; f1[nT++] = 847.8+0.2;     //  Radio broadcast Parkes
      f0[nT] = 849.5-0.1; f1[nT++] = 849.5+0.1;               //  NSW Police Force
      f0[nT] = 848.6-0.230/2.; f1[nT++] = 848.6+0.230/2.; //  Radio broadcast Mount Coonambro
      
      // Note see Licence number 1927906/1 in the ACMA database
      f0[nT] = 2127.5-2.5; f1[nT++] = 2127.5+2.5; // Parkes: "Station open to official correspondence exclusively"

      f0[nT] = 3575; f1[nT++] = 3640;               //  Telstra from Orange or Dubbo
    }


  
  printf("Number of fixed transmitters being removed = %d\n",nT);
  

  for (i=0;i<nchan;i++)
    {
      freq = dSet->head->chanFreq[i];
      for (j=0;j<nT;j++)
	{
	  // Should probably add/subtract one channel bandwidth here to ensure at least
	  // one channel is zapped
	  if (freq+chanbw/2. >= f0[j] && freq-chanbw/2. <= f1[j])
	    {
	      datWts[i] = 0;		  	      
	      nzap++;
	    }
	}
    }
  printf(" ... have flagged %d channels\n",nzap);
}

void autoZapDigitisers(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;

  //
  // Digitiser-related signals that we always see
  //
  f0[nT] = 1023; f1[nT++] = 1025;
  f0[nT] = 1919.9; f1[nT++] = 1920.1;
  f0[nT] = 3071.9; f1[nT++] = 3072.05;       
  
  printf("Number of digitiser-related signals being removed = %d\n",nT);
  
  for (i=0;i<nchan;i++)
    {
      freq = dSet->head->chanFreq[i];
      for (j=0;j<nT;j++)
	{
	  // Should probably add/subtract one channel bandwidth here to ensure at least
	  // one channel is zapped
	  if (freq+chanbw/2. >= f0[j] && freq-chanbw/2. <= f1[j])
	    {
	      datWts[i] = 0;		  	      
	      nzap++;
	    }
	}
    }
  printf(" ... have flagged %d channels\n",nzap);
}

void autoZapUnexplained(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;

  // We do not know what is causing this
  f0[nT] = 824.95; f1[nT++] = 825.05;

  f0[nT] = 1225; f1[nT++] = 1230; // This is potentially interesting - should look more into this one
  f0[nT] = 1399.9; f1[nT++] = 1400.15;
  f0[nT] = 1498.0; f1[nT++] = 1499.0;
  f0[nT] = 1499.8; f1[nT++] = 1500.2;
  f0[nT] = 1880; f1[nT++] = 1904;
  f0[nT] = 2032.0; f1[nT++] = 2033.0;
  f0[nT] = 2063.0; f1[nT++] = 2064.0;
  if (keepSL==0)
    {
      // spectral line at 2078.068 (CNCHO)
      f0[nT] = 2077.0; f1[nT++] = 2078.0;
      f0[nT] = 2079.0; f1[nT++] = 2080.0;
    }
  f0[nT] = 2093.0; f1[nT++] = 2094.0;
  f0[nT] = 2160.0; f1[nT++] = 2161;
  f0[nT] = 2191.0; f1[nT++] = 2192;
  f0[nT] = 2205.0; f1[nT++] = 2206.0;
  f0[nT] = 2207.0; f1[nT++] = 2208.0;
  f0[nT] = 2221.0; f1[nT++] = 2222.0;
  f0[nT] = 2226.3; f1[nT++] = 2226.7;
  printf("Number of unexplained signals being removed = %d\n",nT);
  
  for (i=0;i<nchan;i++)
    {
      freq = dSet->head->chanFreq[i];
      for (j=0;j<nT;j++)
	{
	  // Should probably add/subtract one channel bandwidth here to ensure at least
	  // one channel is zapped
	  if (freq+chanbw/2. >= f0[j] && freq-chanbw/2. <= f1[j])
	    {
	      datWts[i] = 0;		  	      
	      nzap++;
	    }
	}
    }
  printf(" ... have flagged %d channels\n",nzap);
}

void autoZapSatellites(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;

  //
  // Satellites that we always see
  //
  if (keepSL==0)
    {
      // spectral lines at 1624.518 and 1626.161 (^17OH)
      f0[nT] = 1618; f1[nT++] = 1626.5; // Iridium
    }
  if (zapAll==1)
    {
      f0[nT] = 1164; f1[nT++] = 1189;
      f0[nT] = 1189; f1[nT++] = 1214;
      f0[nT] = 1240; f1[nT++] = 1260;
      f0[nT] = 1260; f1[nT++] = 1300;
      if (keepSL==1)
        { 
	  // exclude spectral lines          
	  f0[nT] = 1525; f1[nT++] = 1533;    // Inmarsat
	  f0[nT] = 1544; f1[nT++] = 1566;    // Inmarsat
	  f0[nT] = 1574; f1[nT++] = 1580;    // Inmarsat
	  f0[nT] = 1588; f1[nT++] = 1606;    // Inmarsat
	  f0[nT] = 1616; f1[nT++] = 1620;    // Inmarsat
	  f0[nT] = 1630; f1[nT++] = 1633;    // Inmarsat
	  f0[nT] = 1643; f1[nT++] = 1646.5;  // Inmarsat
        }
      else
        {
          f0[nT] = 1525; f1[nT++] = 1646.5;  // Inmarsat - this is too wide
	}
    } 
  printf("Number of satellite signals being removed = %d\n",nT);
  
  for (i=0;i<nchan;i++)
    {
      freq = dSet->head->chanFreq[i];
      for (j=0;j<nT;j++)
	{
	  // Should probably add/subtract one channel bandwidth here to ensure at least
	  // one channel is zapped
	  if (freq+chanbw/2. >= f0[j] && freq-chanbw/2. <= f1[j])
	    {
	      datWts[i] = 0;		  	      
	      nzap++;
	    }
	}
    }
  printf(" ... have flagged %d channels\n",nzap);
}

void autoZapWiFi(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;

  if (zapAll==1)
    {
      f0[nT] = 2401; f1[nT++] = 2483; // Entire band
    }
  printf("Number of WiFi signals being removed = %d\n",nT);
  
  for (i=0;i<nchan;i++)
    {
      freq = dSet->head->chanFreq[i];
      for (j=0;j<nT;j++)
	{
	  // Should probably add/subtract one channel bandwidth here to ensure at least
	  // one channel is zapped
	  if (freq+chanbw/2. >= f0[j] && freq-chanbw/2. <= f1[j])
	    {
	      datWts[i] = 0;		  	      
	      nzap++;
	    }
	}
    }
  printf(" ... have flagged %d channels\n",nzap);
}

void autoZapHandsets(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;

  if (zapAll==1)
    {
      if (keepSL==1)
        {
          f0[nT] = 706;    f1[nT++] = 713;    // 4G Optus
          f0[nT] = 713;    f1[nT++] = 720;
          f0[nT] = 728;    f1[nT++] = 733;
	}
      else
        {
	  f0[nT] = 703;    f1[nT++] = 713;    // 4G Optus
          f0[nT] = 704.5;  f1[nT++] = 708;    // Alias
          f0[nT] = 713;    f1[nT++] = 733;
        }
      f0[nT] = 825.1;  f1[nT++] = 829.9; // Vodafone 4G
      if (keepSL==1)
        {
          f0[nT] = 838; f1[nT++] = 844.95; // Telstra 3G
        }
      else
        {
          f0[nT] = 830.05; f1[nT++] = 844.95; // Telstra 3G
        }
      f0[nT] = 847.6;  f1[nT++] = 848.0;
      f0[nT] = 898.4;  f1[nT++] = 906.4;  // Optus 3G
      f0[nT] = 906.8;  f1[nT++] = 915.0;  // Vodafone 3G
      f0[nT] = 953;    f1[nT++] = 960.1;  // Alias
      if (keepSL==1)
        {
          f0[nT] = 1710;   f1[nT++] = 1716;  // 4G Telstra
          f0[nT] = 1724;   f1[nT++] = 1725;  // 4G Telstra
        }
      else
        {
          f0[nT] = 1710;   f1[nT++] = 1725;  // 4G Telstra
	}
      f0[nT] = 1745;   f1[nT++] = 1755;  // 4G Optus
      f0[nT] = 2550;   f1[nT++] = 2570;  // 4G Optus
      
    }
  printf("Number of Handset signals being removed = %d\n",nT);
  
  for (i=0;i<nchan;i++)
    {
      freq = dSet->head->chanFreq[i];
      for (j=0;j<nT;j++)
	{
	  // Should probably add/subtract one channel bandwidth here to ensure at least
	  // one channel is zapped
	  if (freq+chanbw/2. >= f0[j] && freq-chanbw/2. <= f1[j])
	    {
	      datWts[i] = 0;		  	      
	      nzap++;
	    }
	}
    }
  printf(" ... have flagged %d channels\n",nzap);
}


void autoZapAircraft(dSetStruct *dSet,float *datWts,int nchan,float chanbw,int zapAll,int keepSL)
{
  int maxTransmitters = 128;
  float f0[maxTransmitters],f1[maxTransmitters];
  float freq;
  int nT=0;
  int i,j;
  int nzap=0;

  //
  // Aircraft signals that we always see
  //
  f0[nT] = 1017; f1[nT++] = 1019;       // Parkes ground response
  f0[nT] = 1029; f1[nT++] = 1031;       // Parkes ground response

  if (zapAll==1)
    {
      f0[nT] = 1026.8; f1[nT++] = 1027.2;       // Unknown DME signal
      f0[nT] = 1027.8; f1[nT++] = 1028.2;       // Unknown DME signal
      f0[nT] = 1032.8; f1[nT++] = 1033.2;       // Unknown DME signal
      f0[nT] = 1040.8; f1[nT++] = 1041.2;       // Strong, unknown DME signal
      f0[nT] = 1061.8; f1[nT++] = 1062.2;       // Sydney DME
      f0[nT] = 1067.8; f1[nT++] = 1068.2;       // Richmond DME
      f0[nT] = 1071.8; f1[nT++] = 1072.2;       // Wagga Wagga DME
      f0[nT] = 1079.2; f1[nT++] = 1080.2;       // Unknown DME
      f0[nT] = 1080.8; f1[nT++] = 1081.2;       // Parkes DME
      f0[nT] = 1081.8; f1[nT++] = 1082.2;       // Sydney RW07 DME
      f0[nT] = 1103.8; f1[nT++] = 1104.2;       // Unknown DME
      f0[nT] = 1102.8; f1[nT++] = 1103.2;       // Unknown DME
      f0[nT] = 1120.8; f1[nT++] = 1121.2;       // Wagga Wagga DME
      f0[nT] = 1134.8; f1[nT++] = 1135.2;       // Nowra DME
      f0[nT] = 1137.6; f1[nT++] = 1138.4;       // Canberra DME
      f0[nT] = 1149.8; f1[nT++] = 1150.2;       // CHECK IF THIS IS A DME
      f0[nT] = 1150.8; f1[nT++] = 1151.2;       // CHECK IF THIS IS A DME
      
    }
  printf("Number of aircraft signals being removed = %d\n",nT);
  

  for (i=0;i<nchan;i++)
    {
      freq = dSet->head->chanFreq[i];
      for (j=0;j<nT;j++)
	{
	  // Should probably add/subtract one channel bandwidth here to ensure at least
	  // one channel is zapped
	  if (freq+chanbw/2. >= f0[j] && freq-chanbw/2. <= f1[j])
	    {
	      datWts[i] = 0;		  	      
	      nzap++;
	    }
	}
    }
  printf(" ... have flagged %d channels\n",nzap);
}

void get_dm_tbin(dSetStruct *dSet,double *dm,double *tbin)
{
  int status=0;
  float fdm,ftbin;

  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_read_key(dSet->fp,TFLOAT,"DM",&fdm,NULL,&status);
  *dm = fdm;
  if (*tbin == 0)
    {
      fits_read_key(dSet->fp,TFLOAT,"TBIN",&ftbin,NULL,&status);
      *tbin = ftbin;
    }
  printf("Loaded DM = %g\n",*dm);
  printf("Loaded TBIN = %g\n",*tbin);
}

void saveFile(float *datWts, int nsub, int nchan, dSetStruct *dSet)
{
  char outname[1024];
  char temp[1024];
  fitsfile *outfptr;
  int status=0;
  int hdu=1;
  int colnum_datWts;
  long ii;

  sprintf(outname,"!%s.autoZap",dSet->fileName);
  fits_create_file(&outfptr, outname, &status);
  if (status) {fits_report_error(stderr, status); printf("NOT WRITING FILE\n"); status=0; return;}
  
  /* Copy every HDU until we get an error */
  while( !fits_movabs_hdu(dSet->fp, hdu++, NULL, &status) )
    {
      fits_copy_hdu(dSet->fp, outfptr, 0, &status);
    }
  /* Reset status after normal error */
  if (status == END_OF_FILE) status = 0;

  // Update the DAT_WTS column
  fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);
  for (ii=0;ii<nsub;ii++)
    fits_write_col(outfptr,TFLOAT,colnum_datWts,ii+1,1,nchan,datWts+ii*nchan,&status);
  
  fits_close_file(outfptr, &status);
  if (status) {fits_report_error(stderr, status); printf("FILE MAY BE CORRUPT\n"); status=0; return;}
  printf("Data file written (%s)\n",outname);
  
}
