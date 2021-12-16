// Code to form a bandpass
//
//gcc -lm -o pfits_formBandpass pfits_formBandpass.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -O3

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include <stdint.h>
#include "fitsio.h"


void help();

void help()
{
}


int main(int argc,char *argv[])
{
  FILE *fout;
  float tsub;
  float chanbw;
  double tbin;
  double dm;
  
  dSetStruct *dSet;
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
  float *datWts;
  float *az;
  float *ze;
  float *timeVal;
  float *fChan,fref;
  int initflag = 0;
  long nchan,nsub,nbin,npol;
  int colnum_data;
  int colnum_datOffs;
  int colnum_datScl;
  int colnum_datWts;
  int colnum_telAz;
  int colnum_telZe;
  long ii,jj,kk,i0;
  int polNum = 0;
  int zapAll=0;
  float bandEdge = 10;
  float **spectrum;
  float *specMean1,*specMean2,*specMean3;
  float mean0,mean1,mean2,mean3,mean4;
  float val,freq;
  float val1,val2;
  float mean4_1,mean4_2;
  float low1,low2,low3;
  int onlyScl=1;
  int n1,n2,n3,n4;
  float sigma1=5;
  float sigma2=5;
  float sigma3=5;

  char outFile[1024]="NULL";
  int ignoreFirstSub=0;
  int ignoreLastSub=0;
  int useWts=0;
  float f1=-1,f2=-1;
  int out2pol=0;
  
  // Initialise everything
  initialise(&dSet,debug);

  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{setFilename(argv[++i],dSet,debug);}
      else if (strcmp(argv[i],"-ignoreFirst")==0)
	ignoreFirstSub=1;
      else if (strcmp(argv[i],"-ignoreLast")==0)
	ignoreLastSub=1;
      else if (strcmp(argv[i],"-useWts")==0)
	useWts=1;
      else if (strcmp(argv[i],"-out2pol")==0)
	out2pol=1;
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outFile,argv[++i]);
      else if (strcmp(argv[i],"-freq")==0)
	{
	  sscanf(argv[++i],"%f",&f1);
	  sscanf(argv[++i],"%f",&f2);
	}
    }
  if (strcmp(outFile,"NULL")==0)
    sprintf(outFile,"%s.bandpass",dSet->fileName);
         
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  nchan = dSet->head->nchan;
  nbin  = dSet->head->nbin;
  nsub  = dSet->head->nsub;
  npol  = dSet->head->npol;
  chanbw = dSet->head->chanFreq[1]-dSet->head->chanFreq[0];
  
  
  printf("Loaded header\n");
  printf("Number of channels = %d\n",nchan);
  printf("Number of bins = %d\n",nbin);
  printf("Number of sub-integrations = %d\n",nsub);
  printf("Number of polarisations = %d\n",npol);

  spectrum = (float **)malloc(sizeof(float *)*nsub);
  specMean1 = (float *)malloc(sizeof(float)*nsub);
  specMean2 = (float *)malloc(sizeof(float)*nsub);
  specMean3 = (float *)malloc(sizeof(float)*nsub);
  for (i=0;i<nsub;i++)
    spectrum[i] = (float *)malloc(sizeof(float)*nchan);
  
  loadAll_pol1 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
  loadAll_pol2 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
  
  fChan    = (float *)malloc(sizeof(float)*nchan);
  sval     = (short int *)calloc(nchan*nbin,sizeof(short int));

  datScl  = (float *)malloc(sizeof(float)*nchan*4); // *4 for npol
  datOffs = (float *)malloc(sizeof(float)*nchan*4);
  datWts  = (float *)malloc(sizeof(float)*nchan);

  timeVal  = (float *)malloc(sizeof(float)*nsub);
  az  = (float *)malloc(sizeof(float)*nsub);
  ze  = (float *)malloc(sizeof(float)*nsub);

  for (i=0;i<nsub;i++)
    timeVal[i] = dSet->head->stt_smjd;
   
  
  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_data,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"TEL_AZ",&colnum_telAz,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"TEL_ZEN",&colnum_telZe,&status);

  printf("Output file = %s\n",outFile);
  fout = fopen(outFile,"w");
  
  if (onlyScl==1) // Determine the bandpass just from the scale parameter
    {
      // Ignoring first subint for now -- usually short!
      if (ignoreFirstSub==1)
	i0=1;
      else
	i0=0;
      
      if (ignoreLastSub==1)
	nsub--;
	  
      for (ii=i0;ii<nsub;ii++)
      //      ii = 1;
	{
	  fits_read_col(dSet->fp,TFLOAT,colnum_datWts,ii+1,1,nchan,&n_fval,datWts,&initflag,&status);
	  fits_read_col(dSet->fp,TFLOAT,colnum_datScl,ii+1,1,nchan*npol,&n_fval,datScl,&initflag,&status);
	  fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,ii+1,1,nchan*npol,&n_fval,datOffs,&initflag,&status);
	  fits_read_col(dSet->fp,TFLOAT,colnum_telAz,ii+1,1,1,&n_fval,&az[ii],&initflag,&status);
	  fits_read_col(dSet->fp,TFLOAT,colnum_telZe,ii+1,1,1,&n_fval,&ze[ii],&initflag,&status);
	  mean0=mean1=mean2=mean3 = mean4= 0.0;
	  mean4_1=0;
	  mean4_2=0;
	  n1=n2=n3=n4=0;

	  for (i=0;i<nchan;i++)
	    {
	      freq = dSet->head->chanFreq[i];
	      val = datScl[i]+datScl[i+nchan];
	      val1 = datScl[i]+datOffs[i];
	      val2 = datScl[i+nchan]+datOffs[i+nchan];
	      //	      printf("Here with %g %g %d\n",datScl[i],datOffs[i],status);
	      spectrum[ii][i] = val;
	      //	      printf("val = %g\n",val);
	      mean0+=val;
	      //	    printf("wts = %g\n",datWts[i]);
	      if (useWts==0 || (useWts==1 && datWts[i] > 0))
		{
		  if (freq > 703 && freq < 1344)
		    {mean1+=val;n1++;}
		  else if (freq >= 1344 && freq <2368)
		    {mean2+=val;n2++;}
		  else
		    {mean3+=val;n3++;}
		  if (f1==-1 && f2 == -1)
		    {mean4+=val;n4++;}
		  else if (freq >= f1 && freq <= f2)
		    {
		      mean4+=val;n4++;
		      mean4_1+=val1;
		      mean4_2+=val2;
		    }
		}
	      if (out2pol==0)
		fprintf(fout,"%.5d %.5f %g %g\n",ii,dSet->head->chanFreq[i],spectrum[ii][i],datWts[i]);
	      else
		fprintf(fout,"%.5d %.5f %g %g %g\n",ii,dSet->head->chanFreq[i],val1,val2,datWts[i]);
	      //	    printf("%.5f %g %g\n",dSet->head->chanFreq[i],datScl[i],
	      //datScl[i+nchan],datScl[i+2*nchan],datScl[i+3*nchan]);
	    }
	  fprintf(fout,"\n");
	  if (out2pol==0)
	    printf("%s %.5d %g %g %g [AzZe] %g\n",dSet->fileName,ii,mean4/(double)n4,az[ii],ze[ii],timeVal[ii]);
	  else
	    printf("%s %.5d %g %g %g %g [AzZe] %g\n",dSet->fileName,ii,mean4_1/(double)n4,mean4_2/(double)n4,az[ii],ze[ii],timeVal[ii]);
	}

    }
  fclose(fout);
  
  /*
  //  for (ii=0;ii<nsub;ii++)
  ii=0;
  {
      // First polarisation
      polNum=0;
      
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);


      fits_read_col(dSet->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);

      for (i=0;i<nchan;i++)
	{
	  val=0.0;
	  for (j=0;j<nbin;j++)
	    {
	      val += sval[i*nbin+j]*datScl[ii*nchan+i]+datOffs[ii*nchan+i];
	    }
	  p1 = val;
	  printf("%.5f %g %g\n",dSet->head->chanFreq[i],p1,datScl[ii*nchan+i]);
	}
     
    }
  */
  
  // De-allocate the memory
  //  pfitsCloseFile(dSet,debug);
  deallocateMemory(&dSet,debug);

  //  free(dVals);
  //  free(sdevArr);
  free(az); free(ze);
  free(timeVal);
  free(loadAll_pol1);
  free(loadAll_pol2);
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

