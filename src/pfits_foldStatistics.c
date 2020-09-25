// Code to zap pulse profiles
//
//gcc -lm -o pfits_foldStatistics pfits_foldStatistics.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -O3

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include <stdint.h>
#include "fitsio.h"


void help()
{
}


int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,j;
  float f0,f1,chanbw;
  
  short int *load_pol1;
  short int *load_pol2;
  short int *load_pol3;
  short int *load_pol4;

  short int n_val =0;
  float n_fval=0;
  long nchan,nsub,nbin,npol;
  long ii,jj,kk;
  int colnum_data;
  int colnum_datOffs;
  int colnum_datScl;
  int colnum_datWts;
  int initflag = 0;
  int polNum = 0;
  float *datScl;
  float *datOffs;
  float *datWts;
  float *fChan,fref;

  int applyOffsScls=0;
  int applyWts=0;
  
  short int lowestVal_pol1,highestVal_pol1;
  
  
  // Initialise everything
  initialise(&dSet,debug);
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-scale")==0)
	applyOffsScls=1;
      else if (strcmp(argv[i],"-weight")==0)
	applyWts=1;
    }
       
 
  
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  // Now load the data to plot
  // Average in polarisation
  // Produce data files also time averaged and frequency averaged

  nchan = dSet->head->nchan;
  nbin  = dSet->head->nbin;
  nsub  = dSet->head->nsub;
  npol  = dSet->head->npol;

  
  printf("Loaded header\n");
  printf("Number of channels = %d\n",nchan);
  printf("Number of bins = %d\n",nbin);
  printf("Number of sub-integrations = %d\n",nsub);
  printf("Number of polarisations = %d\n",npol);
  printf("Size of short int = %d\n",sizeof(short int));

  
  load_pol1 = (short int *)malloc(sizeof(short int)*nchan*nbin*nsub);
  load_pol2 = (short int *)malloc(sizeof(short int)*nchan*nbin*nsub);
  load_pol3 = (short int *)malloc(sizeof(short int)*nchan*nbin*nsub);
  load_pol4 = (short int *)malloc(sizeof(short int)*nchan*nbin*nsub);

  fChan    = (float *)malloc(sizeof(float)*nchan);

  datScl = (float *)malloc(sizeof(float)*nchan*nsub);
  datOffs = (float *)malloc(sizeof(float)*nchan*nsub);
  datWts = (float *)malloc(sizeof(float)*nchan*nsub);
  
  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_data,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);

  printf("Loading all data\n");
  for (ii=0;ii<nchan;ii++)
    fChan[ii] = dSet->head->chanFreq[ii];
  fref = dSet->head->freq;
    
  f0 = fChan[0];
  f1 = fChan[nchan-1];
  chanbw = fChan[1]-fChan[0];
  
  printf("Reference freq = %g\n",fref);
  
  for (ii=0;ii<nsub;ii++)
    {
      // Polarisation 1
      polNum=0;
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datWts,ii+1,1,nchan,&n_fval,datWts+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,load_pol1,&initflag,&status);
	         
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);
      // Don't need to read the weights as not polarisation dependent
      polNum=1;
      fits_read_col(dSet->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,load_pol2,&initflag,&status);
      polNum=2;
      fits_read_col(dSet->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,load_pol3,&initflag,&status);
      polNum=3;
      fits_read_col(dSet->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,load_pol4,&initflag,&status);

      for (i=0;i<nchan;i++)
	{
	  for (j=0;j<nbin;j++)
	    {
	      if (applyOffsScls==1)
		{
		  load_pol1[i*nbin+j] = (load_pol1[i*nbin+j]*datScl[i+ii*nchan]+datOffs[i+ii*nchan]);
		}
	      if (j==0)
		{
		  lowestVal_pol1 = highestVal_pol1 = load_pol1[i*nbin+j];
		}
	      else
		{
		  if (lowestVal_pol1 > load_pol1[i*nbin+j]) lowestVal_pol1 = load_pol1[i*nbin+j];
		  if (highestVal_pol1 < load_pol1[i*nbin+j]) highestVal_pol1 = load_pol1[i*nbin+j];
		}
	      if (ii==40)
	      // Should check for zapped samples
		{
		  if (applyWts==1)
		    {
		      if (datWts[ii*nchan+i]!=0)
			printf("stats: %d %d\n",j,load_pol1[i*nbin+j]);
		    }
		  else
		    printf("stats: %d %d\n",j,load_pol1[i*nbin+j]);
		}
	    }
	  //	  printf("stats: %d %d %d %d %d\n",ii,i,highestVal_pol1-lowestVal_pol1,lowestVal_pol1,highestVal_pol1);
	}
      
    }

      






  // De-allocate the memory
  //  pfitsCloseFile(dSet,debug);
  deallocateMemory(&dSet,debug);


  free(load_pol1);
  free(load_pol2);
  free(load_pol3);
  free(load_pol4);
  free(fChan);
  free(datScl);
  free(datOffs);
  free(datWts);

}

