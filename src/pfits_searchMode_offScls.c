// Code to form a bandpass
//
//gcc -lm -o pfits_searchMode_offScls pfits_searchMode_offScls.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -O3

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include <stdint.h>
#include "fitsio.h"


#define MAX_FILES 128

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

  char fname[MAX_FILES][128];
  int nFiles=0;
  
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,j;

  short int n_val =0;
  float n_fval=0;
  float *datScl;
  float *datOffs;
  float *datWts;
  float *fChan,fref;
  int initflag = 0;
  long nchan,nsub,nsblk,npol;
  float tsamp;
  int colnum_datOffs;
  int colnum_datScl;
  int colnum_datWts;
  long ii,jj,kk,i0;
  int polNum = 0;
  int zapAll=0;
  float bandEdge = 10;
  float meanScl,meanOff,mean;
  float valScl,valOff,freq;
  float val;
  float low1,low2,low3;
  int onlyScl=1;
  int nVal;
  float sigma1=5;
  float sigma2=5;
  float sigma3=5;

  float *plotValsX;
  float *plotValsY1;
  float *plotValsY2;
  int nPlot=0;
  
  char outFile[1024]="NULL";
  int ignoreFirstSub=0;
  int ignoreLastSub=0;
  int useWts=0;
  float f1=-1,f2=-1;
  int out=1;
  int imjd,imjd0;
  float smjd;
  
  // Initialise everything
  initialise(&dSet,debug);

  // Read inputs
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-ignoreFirst")==0)
	ignoreFirstSub=1;
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outFile,argv[++i]);
      else if (strcmp(argv[i],"-ignoreLast")==0)
	ignoreLastSub=1;
      else if (strcmp(argv[i],"-useWts")==0)
	useWts=1;
      else if (strcmp(argv[i],"-out")==0)
	sscanf(argv[++i],"%d",&out);
      else if (strcmp(argv[i],"-freq")==0)
	{
	  sscanf(argv[++i],"%f",&f1);
	  sscanf(argv[++i],"%f",&f2);
	}
      else
	{
	  strcpy(fname[nFiles++],argv[i]);
	}
    }

  printf("outfile = %s\n",outFile);
  if (strcmp(outFile,"NULL")==0)
    sprintf(outFile,"offScls.bandpass",dSet->fileName);

  if (!(fout = fopen(outFile,"w")))
    {
      printf("Unable to open file >%s<\n",outFile);
      exit(1);
    }
  

  
  for (jj=0;jj<nFiles;jj++)
    {
      printf("Processing: %s\n",fname[jj]);
      setFilename(fname[jj],dSet,debug);     

      pfitsOpenFile(dSet,debug);
      pfitsLoadHeader(dSet,debug);

      nchan = dSet->head->nchan;
      tsamp = dSet->head->tsamp;
      nsblk  = dSet->head->nsblk;
      nsub  = dSet->head->nsub;
      npol  = dSet->head->npol;
      chanbw = dSet->head->chanFreq[1]-dSet->head->chanFreq[0];
      imjd = dSet->head->stt_imjd;
      smjd = dSet->head->stt_smjd;

      if (jj==0)
	imjd0 = imjd;
      
      printf("Loaded header\n");
      printf("Number of channels = %d\n",nchan);
      printf("NSBLK = %d\n",nsblk);
      printf("Number of sub-integrations = %d\n",nsub);
      printf("Number of polarisations = %d\n",npol);

      fChan    = (float *)malloc(sizeof(float)*nchan);
      
      datScl  = (float *)malloc(sizeof(float)*nchan*npol);
      datOffs = (float *)malloc(sizeof(float)*nchan*npol);
      datWts  = (float *)malloc(sizeof(float)*nchan);

      fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);
      
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
	      meanScl = 0.0;
	      meanOff = 0.0;
	      mean = 0.0;
	      nVal=0;
	      for (i=0;i<nchan;i++)
		{
		  freq = dSet->head->chanFreq[i];
		  valScl = datScl[i]; //+datScl[i+nchan]; 
		  valOff = datOffs[i]; //+datScl[i+nchan];
		  val = valScl+valOff;
		  
		  //	      printf("val = %g\n",val);
		  //	    printf("wts = %g\n",datWts[i]);
		  if (useWts==0 || (useWts==1 && datWts[i] > 0))
		    {
		      if (freq >= f1 && freq <= f2)
			{
			  meanScl+=valScl;
			  meanOff+=valOff;
			  mean+=val;
			  nVal++;
			}
		    }
		  if (out==2 && (freq >= f1 && freq <= f2))
		    fprintf(fout,"%s %.5d %.5f %g %g %g\n",fname[jj],ii,dSet->head->chanFreq[i],valScl,valOff,datWts[i]);
		  //	    printf("%.5f %g %g\n",dSet->head->chanFreq[i],datScl[i],
		  //datScl[i+nchan],datScl[i+2*nchan],datScl[i+3*nchan]);
		}
	      if (out==2) fprintf(fout,"\n");
	      if (out==1) fprintf(fout,"%s %.5d %g %g %g %g\n",fname[jj],ii,(imjd-imjd0)*86400+smjd+ii*nsblk*tsamp,meanScl/(double)nVal,meanOff/(double)nVal,mean/(double)val);
	      
	    }
	  
	}
      fprintf(fout,"\n");
      printf("I\n");      
      // De-allocate the memory
      //  pfitsCloseFile(dSet,debug);
      
      //  free(dVals);
      //  free(sdevArr);
      
      free(fChan);
      free(datScl);
      free(datOffs);
      free(datWts);
            printf("J\n");
    }
  printf("Closing\n");
  fclose(fout);
  deallocateMemory(&dSet,debug);


  //  float *plotValsX;
  //  float *plotValsY1;
  //  float *plotValsY2;
  //  int nPlot=0;

  
}

