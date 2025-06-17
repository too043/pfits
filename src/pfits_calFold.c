//  Copyright (C) 2015,2016 George Hobbs
// This file is part of the pfits software package
//

/* pfits is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 * pfits is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 * You should have received a copy of the GNU General Public License 
 * along with pfits.  If not, see <http://www.gnu.org/licenses/>. 
*/
// Newton:
// gcc -lm -o pfits_calFold pfits_calFold.c pfits_loader.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include "fftw3.h"

void determineOn_Off(dSetStruct *dSet,float *result_p0,float *result_p1,long nSamples);
void measureCalibrationParams(dSetStruct *dSet,long sub0,long sub1,double baseline0,double baseline1,double on0,double on1);
double fortran_mod( double a, double p);


int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,ii,jj,kk,j,k;
  double on0,on1;
  double baseline0,baseline1;
  float *loadAll_pol1;
  float *loadAll_pol2;
  float *loadAll_pol3;
  float *loadAll_pol4;
  float *datScl;
  float *datOffs;
  float *datWts;
  short int *sval;
  int colnum_data,colnum_datOffs,colnum_datScl,colnum_datWts;
  long nchan,nbin,nsub,npol;   
  float *fChan;
  float n_fval=0;
  short int n_val =0;
  float val;
  int polNum;
  int initflag =0;

  float f1,f2;
  float calOn1=-1.0,calOn2;
  float calOff1=-1.0,calOff2;

  double calOn_aa,calOff_aa;
  double calOn_bb,calOff_bb;
  double calOn_rab,calOff_rab;
  double calOn_iab,calOff_iab;
  double aa,bb,rab,iab;
  int nCalOn,nCalOff;
  double avCalOnAA,avCalOnBB,avCalOffAA,avCalOffBB;
  
  double deltaG2,deltaPsi;
  double sI,sQ,sU,sV;

  float autoCal1=-1,autoCal2=-1;
  int calOnOff[4096];

  FILE *fout;
  FILE *fin;
  char outName[1024]="calFold.dat";
  char tcalName[1024]="NULL";
  float tcalFreq[3328];
  float tcalAA[3328];
  float tcalBB[3328];
  int nTcal=0;
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outName,argv[++i]);
      else if (strcmp(argv[i],"-tcal")==0)
	strcpy(tcalName,argv[++i]);
      else if (strcmp(argv[i],"-autoCal")==0)
	{
	  sscanf(argv[++i],"%f",&autoCal1);
	  sscanf(argv[++i],"%f",&autoCal2);
	}
      else if (strcmp(argv[i],"-calOn")==0) // Pulse Phase 
 	{
	  sscanf(argv[++i],"%f",&calOn1);
	  sscanf(argv[++i],"%f",&calOn2);
	}
      else if (strcmp(argv[i],"-calOff")==0)
	{
	  sscanf(argv[++i],"%f",&calOff1);
	  sscanf(argv[++i],"%f",&calOff2);
	}
    }

  if (strcmp(tcalName,"NULL")!=0)
    {
      fin = fopen(tcalName,"r");
      while (!feof(fin))
	if (fscanf(fin,"%f %f %f",&tcalFreq[nTcal],&tcalAA[nTcal],&tcalBB[nTcal])==3)
	  nTcal++;
      
      fclose(fin);
      printf("Loaded %d tcal measurements\n",nTcal);
    }

  if ((calOn1 == -1 || calOff1 ==-1) && (autoCal1 == -1 && autoCal2 == -1))
    {
      printf("ERROR: must set calOn and calOff (or -autoCal)\n");
      exit(1);
    }
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);
  
  nchan = dSet->head->nchan;
  nbin  = dSet->head->nbin;
  nsub  = dSet->head->nsub;
  npol  = dSet->head->npol;

  if (calOn1 != -1)
    {
      for (i=0;i<nbin;i++)
	{
	  if (i >= calOn1*nbin && i <= calOn2*nbin)
	    calOnOff[i] = 1; // On
	  else if (i >= calOff1*nbin && i <= calOff2*nbin)
	    calOnOff[i] = 2; // Off
	  else
	    calOnOff[i] = 0; // Not set
	}
    }

  
  loadAll_pol1 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
  loadAll_pol2 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
  loadAll_pol3 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
  loadAll_pol4 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);

  sval     = (short int *)calloc(nchan*nbin,sizeof(short int));

  datScl = (float *)malloc(sizeof(float)*nchan*nsub);
  datOffs = (float *)malloc(sizeof(float)*nchan*nsub);
  datWts = (float *)malloc(sizeof(float)*nchan*nsub);
  fChan    = (float *)malloc(sizeof(float)*nchan);

  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_data,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);

  for (i=0;i<nchan;i++)
    fChan[i] = dSet->head->chanFreq[i];

  for (i=0;i<nsub;i++)
    {
      // Polarisation 1
      polNum=0;
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,i+1,1+polNum*nchan,nchan,&n_fval,datOffs+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,i+1,1+polNum*nchan,nchan,&n_fval,datScl+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datWts,i+1,1,nchan,&n_fval,datWts+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TSHORT,colnum_data,i+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);
      for (kk=0;kk<nchan;kk++)
	{
	  for (jj=0;jj<nbin;jj++)
	    {
 	      val = (sval[kk*nbin+jj]*datScl[i*nchan+kk]+datOffs[i*nchan+kk]);
	      loadAll_pol1[i*nbin*nchan+kk*nbin+jj] = val;
	    }
	}

      // Polarisation 2
      polNum=1;
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,i+1,1+polNum*nchan,nchan,&n_fval,datOffs+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,i+1,1+polNum*nchan,nchan,&n_fval,datScl+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datWts,i+1,1,nchan,&n_fval,datWts+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TSHORT,colnum_data,i+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);
      for (kk=0;kk<nchan;kk++)
	{
	  for (jj=0;jj<nbin;jj++)
	    {
 	      val = (sval[kk*nbin+jj]*datScl[i*nchan+kk]+datOffs[i*nchan+kk]);
	      loadAll_pol2[i*nbin*nchan+kk*nbin+jj] = val;
	    }
	}

      // Polarisation 3
      polNum=2;
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,i+1,1+polNum*nchan,nchan,&n_fval,datOffs+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,i+1,1+polNum*nchan,nchan,&n_fval,datScl+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datWts,i+1,1,nchan,&n_fval,datWts+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TSHORT,colnum_data,i+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);
      for (kk=0;kk<nchan;kk++)
	{
	  for (jj=0;jj<nbin;jj++)
	    {
 	      val = (sval[kk*nbin+jj]*datScl[i*nchan+kk]+datOffs[i*nchan+kk]);
	      loadAll_pol3[i*nbin*nchan+kk*nbin+jj] = val;
	    }
	}

      // Polarisation 4
      polNum=3;
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,i+1,1+polNum*nchan,nchan,&n_fval,datOffs+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,i+1,1+polNum*nchan,nchan,&n_fval,datScl+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datWts,i+1,1,nchan,&n_fval,datWts+i*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TSHORT,colnum_data,i+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);
      for (kk=0;kk<nchan;kk++)
	{
	  for (jj=0;jj<nbin;jj++)
	    {
 	      val = (sval[kk*nbin+jj]*datScl[i*nchan+kk]+datOffs[i*nchan+kk]);
	      loadAll_pol4[i*nbin*nchan+kk*nbin+jj] = val;
	    }
	}         
    }
  printf("All data loaded\n");
  if (autoCal1 > -1) // Automatically determine the cal on and off 
    {
      double avCal[nbin];
      double transitionPeriod=0.05;
      double min,max;
      int n=0;
      int nTransit;
      int i1,i2,di1,di2;
      
      for (i=0;i<nbin;i++)
	avCal[i] = 0.0;
      for (k=0;k<nchan;k++)
	{
	  if (fChan[k] >= autoCal1 && fChan[k] <= autoCal2)
	    {
	      for (j=0;j<nsub;j++)
		{
		  for (i=0;i<nbin;i++)
		    avCal[i] += loadAll_pol1[j*nbin*nchan+k*nbin+i] + loadAll_pol2[j*nbin*nchan+k*nbin+i];
		}

	    }
	}
      min=max=avCal[0];
      for (i=0;i<nbin;i++)
	{
	  if (min > avCal[i]) min = avCal[i];
	  if (max < avCal[i]) max = avCal[i];
	}
      for (i=0;i<nbin;i++)
	{
	  if (avCal[i] > min + (max-min)/2.)
	    calOnOff[i] = 1; // On
	  else 
	    calOnOff[i] = 2; // Off
	}
      // Identify transition periods
      i1=i2=-1;
      for (i=0;i<nbin-1;i++)
	{
	  if (calOnOff[i] != calOnOff[i+1])
	    {
	      if (i1==-1) i1=i;
	      else i2 = i;
	    }	    
	}
      for (i=0;i<nbin;i++)
	{
	  di1 = fabs(i1-i);
	  if (di1 < 0) di1 += nbin;
	  if (di1 >= nbin) di1 -= nbin;
	  di2 = fabs(i2-i);
	  if (di2 < 0) di2 += nbin;
	  if (di2 >= nbin) di2 -= nbin;

	  if (fabs(di1) < (int)(transitionPeriod*nbin) || fabs(di2) < (int)(transitionPeriod * nbin))
	    calOnOff[i] = 0;
	  printf("autoCal: %d %g %d %d %d %d %d\n",i,avCal[i],calOnOff[i],i1,i2,di1,di2);	  
	  
	}
      
    }

  fout = fopen(outName,"w");
  for (k=0;k<nchan;k++)
    {
      aa=bb=rab=iab=0.0;
      avCalOnAA=avCalOnBB=avCalOffAA=avCalOffBB=0.0;
      for (j=0;j<nsub;j++)
	{
	  nCalOn=0;
	  nCalOff=0;
	  calOn_aa = calOn_bb = calOn_rab = calOn_iab = 0.0;
	  calOff_aa = calOff_bb = calOff_rab = calOff_iab = 0.0;
	  for (i=0;i<nbin;i++)
	    {
	      if (calOnOff[i] == 1)
		{
		  calOn_aa  += loadAll_pol1[j*nbin*nchan+k*nbin+i];
		  calOn_bb  += loadAll_pol2[j*nbin*nchan+k*nbin+i];
		  calOn_rab += loadAll_pol3[j*nbin*nchan+k*nbin+i];
		  calOn_iab += loadAll_pol4[j*nbin*nchan+k*nbin+i];
		  nCalOn++;
		}	      
	      if (calOnOff[i] == 2)
		{
		  calOff_aa  += loadAll_pol1[j*nbin*nchan+k*nbin+i];
		  calOff_bb  += loadAll_pol2[j*nbin*nchan+k*nbin+i];
		  calOff_rab += loadAll_pol3[j*nbin*nchan+k*nbin+i];
		  calOff_iab += loadAll_pol4[j*nbin*nchan+k*nbin+i];
		  nCalOff++;
		}	      
	      if (j==0 && k == 21395)
		{
		  printf("GEORGE: %d %g %g %g %g %d\n",i,loadAll_pol1[j*nbin*nchan+k*nbin+i],loadAll_pol2[j*nbin*nchan+k*nbin+i],loadAll_pol3[j*nbin*nchan+k*nbin+i],loadAll_pol4[j*nbin*nchan+k*nbin+i],calOnOff[i]);
		}
	    }
	  
	  calOn_aa/=(double)nCalOn;
	  calOn_bb/=(double)nCalOn;
	  calOn_rab/=(double)nCalOn;
	  calOn_iab/=(double)nCalOn;

	  calOff_aa/=(double)nCalOff;
	  calOff_bb/=(double)nCalOff;
	  calOff_rab/=(double)nCalOff;
	  calOff_iab/=(double)nCalOff;

	  if (k==21395)
	    printf("GEORGE2: %g %g %d %d\n",calOn_aa,calOff_aa,nCalOn,nCalOff);
	

	  avCalOnAA += calOn_aa;
	  avCalOnBB += calOn_bb;
	  avCalOffAA += calOff_aa;
	  avCalOffBB += calOff_bb;
	  
	  aa  += (calOn_aa-calOff_aa);
	  bb  += (calOn_bb-calOff_bb);
	  rab += (calOn_rab-calOff_rab);
	  iab += (calOn_iab-calOff_iab);
	}

      // NOT ALL SUBINTS ARE EQUAL -- BUT SAME FOR AA,BB,RAB AND IAB
      aa/=(float)nsub;
      bb/=(float)nsub;
      rab/=(float)nsub;
      iab/=(float)nsub;

      avCalOnAA/=(float)nsub;
      avCalOnBB/=(float)nsub;
      avCalOffAA/=(float)nsub;
      avCalOffBB/=(float)nsub;
      
      sI = aa+bb;
      sQ = aa-bb;
      sU = 2*rab;
      sV = 2*iab;
      
      deltaG2  = sQ/sI;
      deltaPsi = atan2(sV,sU);
      if (nTcal > 0)
	{
	  fprintf(fout,"%.4f %g %g %g %g %g %g\n",fChan[k],avCalOffAA/(avCalOnAA-avCalOffAA)*tcalAA[k],	avCalOffBB/(avCalOnBB-avCalOffBB)*tcalBB[k],avCalOffAA,avCalOffBB,avCalOnAA,avCalOnBB);
	}
      else	
	fprintf(fout,"%d %.4f %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",k,fChan[k],deltaG2,deltaPsi*180/M_PI,avCalOnAA,avCalOnBB,avCalOffAA,avCalOffBB,aa,bb,rab,iab,sI,sQ,sU,sV);

    }
  fclose(fout);
  free(fChan);
  free(datScl);
  free(datOffs);
  free(datWts);
  free(loadAll_pol1);
  free(loadAll_pol2);
  free(loadAll_pol3);
  free(loadAll_pol4);
  free(sval);
  
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
}

