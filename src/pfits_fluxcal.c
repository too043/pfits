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
// gcc -lm -o pfits_fluxcal pfits_fluxcal.c pfits_loader.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include "fftw3.h"

#define MAX_FILES 32
#define MAX_BINS 8192

int main(int argc,char *argv[])
{
  FILE *fout;
  char outName[1024];
    
  dSetStruct **dSet;
  int debug=0;
  int status=0;
  int i,ii,j;
  char fname[MAX_FILES][1024];
  float freqPolAv[MAX_BINS];
  double srcOn_calOn_pol1,srcOff_calOn_pol1,srcOn_calOff_pol1,srcOff_calOff_pol1;
  double srcOn_calOn_pol2,srcOff_calOn_pol2,srcOn_calOff_pol2,srcOff_calOff_pol2;
  double scal_pol1,scal_pol2;
  double ssys_pol1,ssys_pol2;
  double ssrc,freq;
  float weight;
  short int *sval;
  short int n_val =0;
  int polNum = 0;
  float n_fval=0;
  int nFiles=0;
  long nchan,nsub,nbin,npol;
  long jj,kk;
  float *datScl;
  float *datOffs;
  float *datWts;
  double sumWts=0;
  float val;
  float *loadAll_pol1;
  float *loadAll_pol2;
  float **calOnSpec_pol1;
  float **calOffSpec_pol1;
  float **calOnSpec_pol2;
  float **calOffSpec_pol2;
  float **calOnSpec_wt;
  float **calOffSpec_wt;
  int colnum_data;
  int colnum_datOffs;
  int colnum_datScl;
  int colnum_datWts;
  int initflag = 0;
  float calOnBin0,calOnBin1;
  float calOffBin0,calOffBin1;

  // Initialise everything

  // Read inputs
  for (i=1;i<argc;i++)
    strcpy(fname[nFiles++],argv[i]);

  dSet = (dSetStruct **)malloc(sizeof(dSetStruct *)*nFiles);

  for (i=0;i<nFiles;i++)
    {
      initialise(&(dSet[i]),debug);
      setFilename(fname[i],dSet[i],debug);
      pfitsOpenFile(dSet[i],debug);
      pfitsLoadHeader(dSet[i],debug);
      printf("Completed loading header information for %s\n",fname[i]);
      // Average signal in time
      nchan = dSet[i]->head->nchan;
      nbin  = dSet[i]->head->nbin;
      nsub  = dSet[i]->head->nsub;
      npol  = dSet[i]->head->npol;
      if (i==0) // Assuming that all observations have the same number of channels
	{	    
	    calOnSpec_pol1 = (float **)malloc(sizeof(float *)*nFiles);
	    calOffSpec_pol1 = (float **)malloc(sizeof(float *)*nFiles);
	    calOnSpec_pol2 = (float **)malloc(sizeof(float *)*nFiles);
	    calOffSpec_pol2 = (float **)malloc(sizeof(float *)*nFiles);
	    calOnSpec_wt = (float **)malloc(sizeof(float *)*nFiles);
	    calOffSpec_wt = (float **)malloc(sizeof(float *)*nFiles);
	    for (j=0;j<nFiles;j++)
	      {
		calOnSpec_pol1[j] = (float *)calloc(sizeof(float),nchan);
		calOffSpec_pol1[j] = (float *)calloc(sizeof(float),nchan);
		calOnSpec_pol2[j] = (float *)calloc(sizeof(float),nchan);
		calOffSpec_pol2[j] = (float *)calloc(sizeof(float),nchan);
		calOnSpec_wt[j] = (float *)calloc(sizeof(float),nchan);
		calOffSpec_wt[j] = (float *)calloc(sizeof(float),nchan);
	      }
	}
      
      loadAll_pol1 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
      loadAll_pol2 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
      
      datScl = (float *)malloc(sizeof(float)*nchan*nsub);
      datOffs = (float *)malloc(sizeof(float)*nchan*nsub);
      datWts = (float *)malloc(sizeof(float)*nchan*nsub);
      sval     = (short int *)calloc(nchan*nbin,sizeof(short int));
    
      fits_movnam_hdu(dSet[i]->fp,BINARY_TBL,"SUBINT",1,&status);
      fits_get_colnum(dSet[i]->fp,CASEINSEN,"DATA",&colnum_data,&status);
      fits_get_colnum(dSet[i]->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
      fits_get_colnum(dSet[i]->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
      fits_get_colnum(dSet[i]->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);
      for (ii=0;ii<nsub;ii++)
	{
	  // Polarisation 1
	  polNum=0;
	  fits_read_col(dSet[i]->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
	  fits_read_col(dSet[i]->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);
	  fits_read_col(dSet[i]->fp,TFLOAT,colnum_datWts,ii+1,1,nchan,&n_fval,datWts+ii*nchan,&initflag,&status);
	  fits_read_col(dSet[i]->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);
	  for (kk=0;kk<nchan;kk++)
	    {
	      for (jj=0;jj<nbin;jj++)
		{
		  val = (sval[kk*nbin+jj]*datScl[ii*nchan+kk]+datOffs[ii*nchan+kk]);
		  loadAll_pol1[ii*nbin*nchan+kk*nbin+jj] = val;
		}
	    }
	  // Polarisation 2
	  polNum=1;
	  fits_read_col(dSet[i]->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
	  fits_read_col(dSet[i]->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);
	  fits_read_col(dSet[i]->fp,TFLOAT,colnum_datWts,ii+1,1,nchan,&n_fval,datWts+ii*nchan,&initflag,&status);
	  fits_read_col(dSet[i]->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);
	  for (kk=0;kk<nchan;kk++)
	    {
	      for (jj=0;jj<nbin;jj++)
		{
		  val = (sval[kk*nbin+jj]*datScl[ii*nchan+kk]+datOffs[ii*nchan+kk]);
		  loadAll_pol2[ii*nbin*nchan+kk*nbin+jj] = val;
		}
	    }
	}
      // Determine cal on and cal off bins
      sprintf(outName,"%s.sumCal",fname[i]);
      fout = fopen(outName,"w");
      for (jj=0;jj<nbin;jj++)
	{
	  sumWts=0;
	  for (ii=0;ii<nsub;ii++)
	    {
	      for (kk=0;kk<nchan;kk++)
		{
		  // DO I NEED TO DIVIDE SOMEWHERE BY THE SUM OF THE WEIGHTS &******* CHECK THIS
		  weight = datWts[ii*nchan+kk]; // Should squrae??
		  if (ii==0)
		    freqPolAv[jj] = loadAll_pol1[ii*nbin*nchan+kk*nbin+jj]*weight; // + loadAll_pol2[ii*nbin*nchan+kk*nbin+jj]*datWts[ii*nchan+kk];
		  else
		    freqPolAv[jj] += (loadAll_pol1[ii*nbin*nchan+kk*nbin+jj]*weight); //+loadAll_pol2[ii*nbin*nchan+kk*nbin+jj]*datWts[ii*nchan+kk]);
		  sumWts += weight;
		  // I'm out by 1.125!! WHERE FROM?
		}
	    }
	  fprintf(fout,"%s %d %g %g\n",fname[i],jj,freqPolAv[jj]/sumWts,sumWts);
	}
      fclose(fout);
      
      calOnBin0  = 0;     // HARDCODING -- MUST FIX ***
      calOnBin1  = 400;
      calOffBin0 = 500;
      calOffBin1 = 900;

      // Now determine the cal on and cal off spectra
      for (jj=0;jj<nbin;jj++)
	{
	  sumWts=0;
	  for (ii=0;ii<nsub;ii++)
	    {
	      for (kk=0;kk<nchan;kk++)
		{
		  weight = datWts[ii*nchan+kk]; // Should squrae??
		  if (jj >= calOnBin0 && jj <= calOnBin1)
		    {
		      calOnSpec_pol1[i][kk] += loadAll_pol1[ii*nbin*nchan+kk*nbin+jj]*weight;
		      calOnSpec_pol2[i][kk] += loadAll_pol2[ii*nbin*nchan+kk*nbin+jj]*weight;
		      calOnSpec_wt[i][kk] += weight;
		    }
		  if (jj >= calOffBin0 && jj <= calOffBin1)
		    {
		      calOffSpec_pol1[i][kk] += loadAll_pol1[ii*nbin*nchan+kk*nbin+jj]*weight;
		      calOffSpec_pol2[i][kk] += loadAll_pol2[ii*nbin*nchan+kk*nbin+jj]*weight;
		      calOffSpec_wt[i][kk] += weight;
		    }		  
		}
	    }
	}
      sprintf(outName,"%s.spec",fname[i]);
      fout = fopen(outName,"w");

      for (kk=0;kk<nchan;kk++)
	{
	  fprintf(fout,"%s %d %g %g %g %g\n",fname[i],kk,calOnSpec_pol1[i][kk]/calOnSpec_wt[i][kk],calOffSpec_pol1[i][kk]/calOffSpec_wt[i][kk],calOnSpec_pol2[i][kk]/calOnSpec_wt[i][kk],calOffSpec_pol2[i][kk]/calOffSpec_wt[i][kk]);
	}
      fclose(fout);
      
      free(loadAll_pol1);
      free(loadAll_pol2);
      free(datScl);
      free(datOffs);
      free(datWts);
      free(sval);
      printf("Processing next file\n");
    }

  // Now form ON and OFF spectra
  fout = fopen("onOff_1.spec","w");
  for (kk=0;kk<nchan;kk++)
	{
	  freq = 704 + kk;
	  ssrc = pow(10,1.3499-1.0387*log10(freq/1000.)-0.3467*pow(log10(freq/1000.),2)+0.0861*pow(log10(freq/1000.),3));
	  
	  srcOn_calOn_pol1   = calOnSpec_pol1[1][kk]/calOnSpec_wt[1][kk];
	  srcOn_calOff_pol1  = calOffSpec_pol1[1][kk]/calOffSpec_wt[1][kk];
	  srcOff_calOn_pol1  = (calOnSpec_pol1[0][kk]/calOnSpec_wt[0][kk] + calOnSpec_pol1[2][kk]/calOnSpec_wt[2][kk])*0.5;
	  srcOff_calOff_pol1 = (calOffSpec_pol1[0][kk]/calOffSpec_wt[0][kk] + calOffSpec_pol1[2][kk]/calOffSpec_wt[2][kk])*0.5;
	  scal_pol1 = (srcOff_calOn_pol1 - srcOff_calOff_pol1)/(srcOff_calOff_pol1) * (srcOff_calOff_pol1/(srcOn_calOff_pol1-srcOff_calOff_pol1));

	  srcOn_calOn_pol2   = calOnSpec_pol2[1][kk]/calOnSpec_wt[1][kk];
	  srcOn_calOff_pol2  = calOffSpec_pol2[1][kk]/calOffSpec_wt[1][kk];
	  srcOff_calOn_pol2  = (calOnSpec_pol2[0][kk]/calOnSpec_wt[0][kk] + calOnSpec_pol2[2][kk]/calOnSpec_wt[2][kk])*0.5;
	  srcOff_calOff_pol2 = (calOffSpec_pol2[0][kk]/calOffSpec_wt[0][kk] + calOffSpec_pol2[2][kk]/calOffSpec_wt[2][kk])*0.5;
	  scal_pol2 = (srcOff_calOn_pol2 - srcOff_calOff_pol2)/(srcOff_calOff_pol2) * (srcOff_calOff_pol2/(srcOn_calOff_pol2-srcOff_calOff_pol2));

	  ssys_pol1 = (srcOff_calOff_pol1/(srcOn_calOff_pol1-srcOff_calOff_pol1));
	  ssys_pol2 = (srcOff_calOff_pol2/(srcOn_calOff_pol2-srcOff_calOff_pol2));

	  fprintf(fout,"%d %g %g %g %g %g %g %g %g %g %g %g %g\n",kk,srcOn_calOn_pol1,srcOn_calOff_pol1,srcOff_calOn_pol1,srcOff_calOff_pol1,scal_pol1*ssrc/2.,
		  srcOn_calOn_pol2,srcOn_calOff_pol2,srcOff_calOn_pol2,srcOff_calOff_pol2,scal_pol2*ssrc/2.,ssys_pol1*ssrc/2.,ssys_pol2*ssrc/2.);
	}
	    
  fclose(fout);

  fout = fopen("onOff_2.spec","w");
  for (kk=0;kk<nchan;kk++)
	{
	  freq = 704 + kk;
	  ssrc = pow(10,1.3499-1.0387*log10(freq/1000.)-0.3467*pow(log10(freq/1000.),2)+0.0861*pow(log10(freq/1000.),3));

	  srcOn_calOn_pol1  = calOnSpec_pol1[3][kk]/calOnSpec_wt[3][kk];
	  srcOn_calOff_pol1 = calOffSpec_pol1[3][kk]/calOffSpec_wt[3][kk];
	  srcOff_calOn_pol1 = (calOnSpec_pol1[4][kk]/calOnSpec_wt[4][kk] + calOnSpec_pol1[2][kk]/calOnSpec_wt[2][kk])*0.5;
	  srcOff_calOff_pol1 = (calOffSpec_pol1[4][kk]/calOffSpec_wt[4][kk] + calOffSpec_pol1[2][kk]/calOffSpec_wt[2][kk])*0.5;
	  scal_pol1 = (srcOff_calOn_pol1 - srcOff_calOff_pol1)/(srcOff_calOff_pol1) * (srcOff_calOff_pol1/(srcOn_calOff_pol1-srcOff_calOff_pol1));

	  srcOn_calOn_pol2  = calOnSpec_pol2[3][kk]/calOnSpec_wt[3][kk];
	  srcOn_calOff_pol2 = calOffSpec_pol2[3][kk]/calOffSpec_wt[3][kk];
	  srcOff_calOn_pol2 = (calOnSpec_pol2[4][kk]/calOnSpec_wt[4][kk] + calOnSpec_pol2[2][kk]/calOnSpec_wt[2][kk])*0.5;
	  srcOff_calOff_pol2 = (calOffSpec_pol2[4][kk]/calOffSpec_wt[4][kk] + calOffSpec_pol2[2][kk]/calOffSpec_wt[2][kk])*0.5;
	  scal_pol2 = (srcOff_calOn_pol2 - srcOff_calOff_pol2)/(srcOff_calOff_pol2) * (srcOff_calOff_pol2/(srcOn_calOff_pol2-srcOff_calOff_pol2));


	  ssys_pol1 = (srcOff_calOff_pol1/(srcOn_calOff_pol1-srcOff_calOff_pol1));
	  ssys_pol2 = (srcOff_calOff_pol2/(srcOn_calOff_pol2-srcOff_calOff_pol2));

	  fprintf(fout,"%d %g %g %g %g %g %g %g %g %g %g %g %g\n",kk,srcOn_calOn_pol1,srcOn_calOff_pol1,srcOff_calOn_pol1,srcOff_calOff_pol1,scal_pol1*ssrc/2.,
		  srcOn_calOn_pol2,srcOn_calOff_pol2,srcOff_calOn_pol2,srcOff_calOff_pol2,scal_pol2*ssrc/2.,ssys_pol1*ssrc/2.,ssys_pol2*ssrc/2.);
	}
	    
  fclose(fout);

  
  // De-allocate the memory
  printf("Deallocating\n");

  for (i=0;i<nFiles;i++)
    {
      deallocateMemory(&dSet[i],debug);
      free(calOnSpec_pol1[i]);
      free(calOffSpec_pol1[i]);
      free(calOnSpec_pol2[i]);
      free(calOffSpec_pol2[i]);
      free(calOnSpec_wt[i]);
      free(calOffSpec_wt[i]);
      
    }
  
	    
  free(calOnSpec_pol1);
  free(calOffSpec_pol1);
  free(calOnSpec_pol2);
  free(calOffSpec_pol2);
  free(calOnSpec_wt);
  free(calOffSpec_wt);

  free(dSet);
  //  printf("Free'ing results\n");

  //  free(result_p0);
  //  free(result_p1);
}

