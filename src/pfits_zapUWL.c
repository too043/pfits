// Code to zap pulse profiles from the UWL receiver
//
//gcc -lm -o pfits_zapUWL pfits_zapUWL.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -O3

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include <stdint.h>
#include "fitsio.h"


void help();
void get_dm_tbin(dSetStruct *dSet,double *dm,double *tbin);
void saveFile(float *datWts, int nsub, int nchan,dSetStruct *dSet);

void help()
{
}


int main(int argc,char *argv[])
{
  float tsub = 30; // Should read from file
  
  // Define a small clean band (MHz) for extracting cal on/off bins
  float cleanBand_start = 2650;
  float cleanBand_end   = 2700;

  // Cal on/off ratio in a given channel in a given subband must
  // be greater than this value
  float sn_calOnOffRatio = 2;                             

  // Zap band edge percent (+/- val/2)
  float bandEdgePercent = 5;

  float *sdevArr;
  int   nChanPsr;
  int   nBinPsr;
  int   nSubPsr;
  float *dVals;

  float min,max;
  
  int   bin_calOff2On;
  int   bin_calOn2Off;
  int   *calOnOrOff;

  float *calOnMean, *calOffMean;
  float *calOnSdev, *calOffSdev;
  int   nOn,nOff;
  
  dSetStruct *dSetCal;
  dSetStruct *dSetPsr;
  int debug=0;
  int status=0;
  int i,j;

  float *loadAll_pol1;
  float *loadAll_pol2;
  short int *sval,*svalPsr;
  short int n_val =0;
  float n_fval=0;
  float *datScl;
  float *datOffs;
  float *datWtsCal;
  float *datWtsPsr;
  float *datSclPsr;
  float *datOffsPsr;
  float *fChan,fref;
  int initflag = 0;
  long nchan,nsub,nbin,npol;
  int colnum_data;
  int colnum_datOffs;
  int colnum_datScl;
  int colnum_datWts;
  long ii,jj,kk;
  float f0,f1,chanbw;
  int polNum = 0;
  double dm=0;
  double tbin=0;
  double mean1,mean2,mean1_2,mean2_2;  
  float val;

  int allZero;
  
  double sum_on,sum2_on;
  double sum_off,sum2_off;
  double *cleanCalSum;
  float minClean, maxClean;

  FILE *fin;
  float tsys[4096];
  float tcal[4096];
  int   nTsysTcal=0;
  float temp;
  float tcal_p1,tcal_p2,tsys_p1,tsys_p2;

  int haveCal=0;
  int havePsr=0;
  
  // Initialise everything
  initialise(&dSetCal,debug);
  initialise(&dSetPsr,debug);
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-cal")==0)
	{setFilename(argv[++i],dSetCal,debug); haveCal=1;}
      if (strcmp(argv[i],"-psr")==0)
	{setFilename(argv[++i],dSetPsr,debug); havePsr=1;}
      else if (strcmp(argv[i],"-tbin")==0)
	sscanf(argv[++i],"%lf",&tbin);
    }
       

  fin = fopen("scal_ssys.dat","r");
  nTsysTcal=0;
  while (!feof(fin))
    {
      if (fscanf(fin,"%f %f %f %f %f %f %f %f %f %f",&temp,&tcal_p1,&temp,&tcal_p2,&temp,
		 &temp,&tsys_p1,&temp,&tsys_p2,&temp)==10)
	{
	  tcal[nTsysTcal] = (tcal_p1+tcal_p2)/2.;
	  tsys[nTsysTcal] = (tsys_p1+tsys_p2)/2.;
	  nTsysTcal++;
	}
    
    }
  fclose(fin);
  
  pfitsOpenFile(dSetCal,debug);
  pfitsLoadHeader(dSetCal,debug);
  if (havePsr==1)
    {
      pfitsOpenFile(dSetPsr,debug);
      pfitsLoadHeader(dSetPsr,debug);
    }
  // Now load the data to plot
  // Average in polarisation
  // Produce data files also time averaged and frequency averaged

  nchan = dSetCal->head->nchan;
  nbin  = dSetCal->head->nbin;
  nsub  = dSetCal->head->nsub;
  npol  = dSetCal->head->npol;

  get_dm_tbin(dSetCal,&dm,&tbin);
  
  printf("Loaded header\n");
  printf("Number of channels = %d\n",nchan);
  printf("Number of bins = %d\n",nbin);
  printf("Number of sub-integrations = %d\n",nsub);
  printf("Number of polarisations = %d\n",npol);
  printf("Dispersion measure = %g\n",dm);

  loadAll_pol1 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
  loadAll_pol2 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
  cleanCalSum = (double *)malloc(sizeof(double)*nbin);
  calOnOrOff  = (int *)malloc(sizeof(int)*nbin);
  calOnMean   = (float *)malloc(sizeof(float)*nchan*nsub);
  calOffMean   = (float *)malloc(sizeof(float)*nchan*nsub);
  calOnSdev   = (float *)malloc(sizeof(float)*nchan*nsub);
  calOffSdev   = (float *)malloc(sizeof(float)*nchan*nsub);


  for (ii=0;ii<nbin;ii++)
    cleanCalSum[ii] = 0.0;

  
  fChan    = (float *)malloc(sizeof(float)*nchan);
  sval     = (short int *)calloc(nchan*nbin,sizeof(short int));

     
  datScl = (float *)malloc(sizeof(float)*nchan*nsub);
  datOffs = (float *)malloc(sizeof(float)*nchan*nsub);
  datWtsCal = (float *)malloc(sizeof(float)*nchan*nsub);

  if (havePsr==1)
    {
      datWtsPsr = (float *)malloc(sizeof(float)*dSetPsr->head->nchan*dSetPsr->head->nsub);
      datSclPsr = (float *)malloc(sizeof(float)*dSetPsr->head->nchan*dSetPsr->head->nsub);
      datOffsPsr = (float *)malloc(sizeof(float)*dSetPsr->head->nchan*dSetPsr->head->nsub);
    }
  
  fits_movnam_hdu(dSetCal->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSetCal->fp,CASEINSEN,"DATA",&colnum_data,&status);
  fits_get_colnum(dSetCal->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
  fits_get_colnum(dSetCal->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
  fits_get_colnum(dSetCal->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);

  printf("Loading all data\n");
  for (ii=0;ii<nchan;ii++)
    fChan[ii] = dSetCal->head->chanFreq[ii];
  fref = dSetCal->head->freq;
    
  f0 = fChan[0];
  f1 = fChan[nchan-1];
  chanbw = fChan[1]-fChan[0];
  
  printf("Reference freq = %g\n",fref);
  
  for (ii=0;ii<nsub;ii++)
    {
      // Polarisation 1
      polNum=0;
      fits_read_col(dSetCal->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
      fits_read_col(dSetCal->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);
      fits_read_col(dSetCal->fp,TFLOAT,colnum_datWts,ii+1,1,nchan,&n_fval,datWtsCal+ii*nchan,&initflag,&status);
      fits_read_col(dSetCal->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);
	   
      for (kk=0;kk<nchan;kk++)
	{
	  mean1=0;
	  for (jj=0;jj<nbin;jj++)
	    {
 	      val = (sval[kk*nbin+jj]*datScl[ii*nchan+kk]+datOffs[ii*nchan+kk]);
	      loadAll_pol1[ii*nbin*nchan+kk*nbin+jj] = val;
	      mean1 += val;
	    }
	  for (jj=0;jj<nbin;jj++)
	    loadAll_pol1[ii*nbin*nchan+kk*nbin+jj] -= mean1/(double)nbin;
	}
      
      // Polarisation 2
      polNum=1;
      fits_read_col(dSetCal->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
      fits_read_col(dSetCal->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);
      // Don't need to read the weights as not polarisation dependent
      fits_read_col(dSetCal->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);

      for (kk=0;kk<nchan;kk++)
	{
	  mean2=0;
	  for (jj=0;jj<nbin;jj++)
	    {
 	      val = (sval[kk*nbin+jj]*datScl[ii*nchan+kk]+datOffs[ii*nchan+kk]);
	      loadAll_pol2[ii*nbin*nchan+kk*nbin+jj] = val;
	      mean2+=val;
	    }	      
	  for (jj=0;jj<nbin;jj++)
	    {
	      loadAll_pol2[ii*nbin*nchan+kk*nbin+jj] -= mean2/(double)nbin;
	      //	      printf("%d %d %d %g %g\n",ii,kk,jj,loadAll_pol1[ii*nbin*nchan+kk*nbin+jj],loadAll_pol2[ii*nbin*nchan+kk*nbin+jj]);
	    }
	}
    }

  // Obtain average cal from a clean band
  for (ii=0;ii<nsub;ii++)
    {
      for (kk=0;kk<nchan;kk++)
	{
	  if (fChan[kk] > cleanBand_start && fChan[kk] <= cleanBand_end)
	    {
	      for (jj=0;jj<nbin;jj++)
		{
		  cleanCalSum[jj]+=(loadAll_pol1[ii*nbin*nchan+kk*nbin+jj]+loadAll_pol2[ii*nbin*nchan+kk*nbin+jj]);
		}
	    }
	}
    }

  minClean = 1e99;
  maxClean = -1e99;
  for (jj=0;jj<nbin;jj++)
    {
      if (minClean > cleanCalSum[jj]) minClean = cleanCalSum[jj];
      if (maxClean < cleanCalSum[jj]) maxClean = cleanCalSum[jj];
    }
  for (jj=0;jj<nbin-1;jj++)
    {
      if (cleanCalSum[jj] < (minClean+maxClean)*0.5 &&
	  cleanCalSum[jj+1] > (minClean+maxClean)*0.5)
	bin_calOff2On = jj;

      if (cleanCalSum[jj] > (minClean+maxClean)*0.5 &&
	  cleanCalSum[jj+1] < (minClean+maxClean)*0.5)
	bin_calOn2Off = jj;
    }
  // Sort out if each bin is on (=1)/off (=0) or unknown (=-1)
  for (jj=0;jj<nbin;jj++)
    {
      if (bin_calOff2On < bin_calOn2Off)
	{
	  if (jj > bin_calOff2On+(bin_calOn2Off-bin_calOff2On)*0.1 && jj < bin_calOn2Off - (bin_calOn2Off-bin_calOff2On)*0.1)
	    calOnOrOff[jj] = 1;
	  else if (jj < bin_calOff2On-(bin_calOn2Off-bin_calOff2On)*0.1)
	    calOnOrOff[jj] = 0;
	  else if (jj > bin_calOn2Off+(bin_calOn2Off-bin_calOff2On)*0.1)
	    calOnOrOff[jj] = 0;
	  else
	    calOnOrOff[jj] = -1;
	}
      printf("cal: %d %g %d\n",jj,cleanCalSum[jj],calOnOrOff[jj]);

    }

  // Now process all the cal channels
  for (ii=0;ii<nsub;ii++)
  {
      for (kk=0;kk<nchan;kk++)
	{
	  // Check if WiFi is on
	  // Check from 2385 to 2395 as clean bands near WiFi
	  // and 2401 to 2495 as WiFi bands.
	  // Wipe those if needed
	  // GOT TO HERE

	  //
	  // Zap band edges because of critically sampled filterbank
	  //
	  for (jj=0;jj<26;jj++) // Number of subbands
	    {
	      if (fChan[kk] > 704+jj*128-bandEdgePercent/100.*128/2. &&
		  fChan[kk] < 704+jj*128+bandEdgePercent/100.*128/2.)
		datWtsCal[ii*nchan+kk]=0;
	    }


	  	 
	  sum_on = sum2_on = sum_off = sum2_off = 0.0;
	  nOn = 0;
	  nOff = 0;
	    for (jj=0;jj<nbin;jj++)
	      {
		val = loadAll_pol1[ii*nbin*nchan+kk*nbin+jj]+loadAll_pol2[ii*nbin*nchan+kk*nbin+jj];
 		if (calOnOrOff[jj] == 1)
		  {
		    sum_on += val;
		    sum2_on += val*val;
		    nOn++;
		  }
 		else if (calOnOrOff[jj] == 0)
		  {
		    sum_off += val;
		    sum2_off += val*val;
		    nOff++;
		  }
	      }
	    calOnMean[ii*nchan+kk] = sum_on/(double)nOn;
	    calOffMean[ii*nchan+kk] = sum_off/(double)nOff;
	    calOnSdev[ii*nchan+kk] = sqrt(1./(double)nOn*sum2_on - pow(sum_on/(double)nOn,2));
	    calOffSdev[ii*nchan+kk] = sqrt(1./(double)nOff*sum2_off - pow(sum_off/(double)nOff,2));
	    printf("%d %d %g %g %g %g calOnOff_mean_sdev\n",ii,kk,calOnMean[ii*nchan+kk],calOnSdev[ii*nchan+kk],calOffMean[ii*nchan+kk],calOffSdev[ii*nchan+kk]);
	}
    }

    for (jj=0;jj<nbin;jj++)
    {
      kk = 1300;
      ii = 0;
      val = loadAll_pol1[ii*nbin*nchan+kk*nbin+jj]+loadAll_pol2[ii*nbin*nchan+kk*nbin+jj];
      printf("%d %g test1300\n",jj,val);
    }

  
  // Clear out really bad sub-ints and channels
  for (ii=0;ii<nsub;ii++)
    {
      for (kk=0;kk<nchan;kk++)
	{
	  val = (calOnMean[ii*nchan+kk]-calOffMean[ii*nchan+kk])/calOffSdev[ii*nchan+kk];
	  // Zap if incorrect S/N (/2 = scaling factor)
	  if (tcal[kk] <0.1 || tsys[kk] < 0.1 || tcal[kk]/tsys[kk]*sqrt(tsub/nbin*1e6)/2 > val)
	    {
	      datWtsCal[ii*nchan+kk]=0;
	    }
	  
	  if (val < sn_calOnOffRatio)
	    {
	      datWtsCal[ii*nchan+kk]=0;
	    }
	  if (ii==2) printf("%d %d %g %g %g outResults\n",ii,kk,datWtsCal[ii*nchan+kk],val,tcal[kk]/tsys[kk]);
	}
    }
  
  saveFile(datWtsCal,nsub,nchan,dSetCal);

  // Now process the pulsar file
  if (havePsr==1)
    {
      fits_movnam_hdu(dSetPsr->fp,BINARY_TBL,"SUBINT",1,&status);
      fits_get_colnum(dSetPsr->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);
      fits_get_colnum(dSetPsr->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
      fits_get_colnum(dSetPsr->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
      fits_get_colnum(dSetPsr->fp,CASEINSEN,"DATA",&colnum_data,&status);
      
      for (ii=0;ii<dSetPsr->head->nsub;ii++)
	{
	  fits_read_col(dSetPsr->fp,TFLOAT,colnum_datWts,ii+1,1,dSetPsr->head->nchan,
			&n_fval,datWtsPsr+ii*dSetPsr->head->nchan,&initflag,&status);
	}
      
      // set weights to zero if all sub-bands are weighted to zero in the cal
      for (jj=0;jj<nchan;jj++)
	{
	  allZero=1;
	  for (kk=0;kk<nsub;kk++)
	    {
	      if (datWtsCal[kk*nchan+jj]!=0)
		{allZero=0; break;}
	    }
	  if (allZero==1)
	    {
	      for (ii=0;ii<dSetPsr->head->nchan;ii++)
		{
		  for (kk=0;kk<dSetPsr->head->nsub;kk++)
		    {
		      if (dSetPsr->head->chanFreq[ii] >= fChan[jj] && dSetPsr->head->chanFreq[ii] <= fChan[jj])
			datWtsPsr[kk*dSetPsr->head->nchan+ii] = 0;		  
		    }
		}
	    }
	}
      
      
      // For any channel that hasn't been zeroed check for a large outlier in the total signal in a given channel
      // Get an array of sdev as a function of frequency and subint -- FOR BRIGHT PULSARS ALLOW THE USER TO DEFINE THE BASELINE RANGE
      //
      
      nChanPsr = dSetPsr->head->nchan;
      nBinPsr  = dSetPsr->head->nbin;
      nSubPsr  = dSetPsr->head->nsub;
      
      sdevArr = (float *)malloc(sizeof(float)*nChanPsr*nSubPsr);
      dVals   = (float *)malloc(sizeof(float)*nBinPsr);
      svalPsr = (short int *)calloc(nChanPsr*nBinPsr,sizeof(short int));
      
      for (ii=0;ii<nSubPsr;ii++)
	{
	  // Polarisation 1
	  polNum=0;
	  fits_read_col(dSetPsr->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nChanPsr,nChanPsr,&n_fval,datOffsPsr+ii*nChanPsr,&initflag,&status);
	  fits_read_col(dSetPsr->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nChanPsr,nChanPsr,&n_fval,datSclPsr+ii*nChanPsr,&initflag,&status);
	  fits_read_col(dSetPsr->fp,TSHORT,colnum_data,ii+1,1+polNum*nChanPsr*nBinPsr,nChanPsr*nBinPsr,&n_val,svalPsr,&initflag,&status);
	  
	  for (kk=0;kk<nChanPsr;kk++)
	    {
	      mean1=0;
	      mean1_2=0;
	      for (jj=0;jj<nBinPsr;jj++)
		{
		  val = (svalPsr[kk*nBinPsr+jj]*datSclPsr[ii*nChanPsr+kk]+datOffsPsr[ii*nChanPsr+kk]);
		  mean1 += val;
		  mean1_2 += val*val;
		}
	      sdevArr[kk*nSubPsr+ii] = sqrt(1./(double)nBinPsr*mean1_2 - pow(mean1/(double)nBinPsr,2));
	    }
	  
	  // Polarisation 1
	  polNum=1;
	  fits_read_col(dSetPsr->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nChanPsr,nChanPsr,&n_fval,datOffsPsr+ii*nChanPsr,&initflag,&status);
	  fits_read_col(dSetPsr->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nChanPsr,nChanPsr,&n_fval,datSclPsr+ii*nChanPsr,&initflag,&status);
	  fits_read_col(dSetPsr->fp,TSHORT,colnum_data,ii+1,1+polNum*nChanPsr*nBinPsr,nChanPsr*nBinPsr,&n_val,svalPsr,&initflag,&status);
	  
	  for (kk=0;kk<nChanPsr;kk++)
	    {
	      mean2=0;
	      mean2_2=0;
	      max=-1e99;
	      for (jj=0;jj<nBinPsr;jj++)
		{
		  val = (svalPsr[kk*nBinPsr+jj]*datSclPsr[ii*nChanPsr+kk]+datOffsPsr[ii*nChanPsr+kk]);
		  if (max < val) max = val;
		  mean2 += val;
		  mean2_2 += val*val;
		  
		}
	      val = sdevArr[kk*nSubPsr+ii];
	      sdevArr[kk*nSubPsr+ii] = (sdevArr[kk*nSubPsr+ii] + sqrt(1./(double)nBinPsr*mean2_2 - pow(mean2/(double)nBinPsr,2)))/2.;
	      //	  if (kk==2710)
	      //	    printf("loaded gh: %d %d %d %g %g %g\n",ii,kk,jj,sdevArr[kk*nSubPsr+ii],sqrt(1./(double)nBinPsr*mean2_2 - pow(mean2/(double)nBinPsr,2)),val);
	      
	      //	  sdevArr[kk*nSubPsr+ii] = max;
	    }
	  
	}
      printf("Zapping channels in pulsar file\n");
      for (ii=0;ii<nChanPsr;ii++)
	{
	  allZero=1;
	  for (kk=0;kk<nSubPsr;kk++)
	    {
	      if (datWtsPsr[kk*nChanPsr+ii]!=0)
		{allZero=0; break;}
	    }
	  if (allZero==0)
	    {
	      min = 1e99;
	      max = -1e99;
	      for (kk=0;kk<dSetPsr->head->nsub;kk++)
		{
		  if (min > sdevArr[ii*nSubPsr+kk]) min = sdevArr[ii*nSubPsr+kk];
		  if (max < sdevArr[ii*nSubPsr+kk]) max = sdevArr[ii*nSubPsr+kk];	      
		}
	      // Should set this 2 properly
	      for (kk=0;kk<dSetPsr->head->nsub;kk++)
		{
		  //	      if (ii==2710)
		  //		  printf("test George: %d %g %g %g\n",kk,sdevArr[ii*nSubPsr+kk],min,max);
		  if (sdevArr[ii*nSubPsr+kk] > 2*min)
		    {
		      datWtsPsr[kk*nChanPsr+ii] = 0;
		      //		if (ii==405)
		      //		    printf("Flagging channel %d, subint %d (out of %d)\n",ii,kk,nSubPsr);
		    }
		}
	      
	    }
	}
      {
	FILE *fout;
	fout = fopen("sdev.dat","w");
	
	for (ii=0;ii<dSetPsr->head->nchan;ii++)
	  {
	    for (kk=0;kk<dSetPsr->head->nsub;kk++)
	      {
		if (datWtsPsr[kk*dSetPsr->head->nchan+ii] == 0)
		  fprintf(fout,"%d %d %g zapped\n",ii,kk,0.0);
		else
		  fprintf(fout,"%d %d %g\n",ii,kk,sdevArr[ii*nSubPsr+kk]);
	      }
	    fprintf(fout,"\n");
	  }
	fclose(fout);
      }
      
      
      printf("Unloading pulsar file\n");
      saveFile(datWtsPsr,dSetPsr->head->nsub,dSetPsr->head->nchan,dSetPsr);        
    }
  
  // De-allocate the memory
  //  pfitsCloseFile(dSet,debug);
  deallocateMemory(&dSetCal,debug);

  if (havePsr==1)
    {
      deallocateMemory(&dSetPsr,debug);
      free(datWtsPsr);
      free(datSclPsr);
      free(datOffsPsr);
      free(svalPsr);
    }
  
  free(dVals);
  free(sdevArr);

  free(calOnOrOff);
  free(loadAll_pol1);
  free(loadAll_pol2);
  free(sval);
  free(fChan);
  free(datScl);
  free(datOffs);
  free(datWtsCal);

  free(cleanCalSum);
  free(calOnMean); free(calOffMean); free(calOnSdev); free(calOffSdev);
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
