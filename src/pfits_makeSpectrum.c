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

//  gcc -lm -o pfits_makeSpectrum pfits_makeSpectrum.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -lpgplot

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int plot=0;
  int status=0;
  int i,j;
  float *data;
  int subint;
  long nSamples;
  int nTime,nFreq;
  float *sum,*sumSub,*max,*min;
  float *meanSub;
  long totCount;
  int sub0,sub1;
  int pol=0;
  int histogramTot[255];
  int applyScaling=1;
  int colnum_datoffs,colnum_datscl;
  float *datOffs,*datScl;
  float n_fval=0;
  int initflag = 0;
  FILE *fout;
  FILE *fout_ts;
  FILE *fout_samp;
  char grDev[128]="/xs";
    
  sub0 = 0;
  sub1 = 1;
  
  // Initialise everything
  initialise(&dSet,debug);

  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-s0")==0)
	sscanf(argv[++i],"%d",&sub0);
      else if (strcmp(argv[i],"-s1")==0)
	sscanf(argv[++i],"%d",&sub1);
      else if (strcmp(argv[i],"-pol")==0)
	sscanf(argv[++i],"%d",&pol);
      else if (strcmp(argv[i],"-noScale")==0)
	applyScaling=0;
      else if (strcmp(argv[i],"-g")==0)
	{
	  plot=1;
	  strcpy(grDev,argv[++i]);
	}
    }

  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datoffs,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datscl,&status);
  
  data = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  datOffs = (float *)malloc(sizeof(float)*dSet->head->nchan);
  datScl = (float *)malloc(sizeof(float)*dSet->head->nchan);

  sum = (float *)malloc(sizeof(float)*dSet->head->nchan);
  sumSub = (float *)malloc(sizeof(float)*dSet->head->nchan);
  max = (float *)malloc(sizeof(float)*dSet->head->nchan);
  min = (float *)malloc(sizeof(float)*dSet->head->nchan);
  meanSub = (float *)malloc(sizeof(float)*dSet->head->nchan*((sub1-sub0)+1));
  for (i=0;i<dSet->head->nchan;i++)
    {sum[i]=sumSub[i]=0.0; max[i]=0.0; min[i]=1e99;}
  totCount=0;


  fout_ts = fopen("timeSeries.dat","w");
  fout_samp = fopen("sample.dat","w");
  for (subint=sub0;subint<sub1;subint++)
    {
      printf("Subint %d\n",subint);
      // SHOULD HAVE *NPOL HERE -- ASSUME = 1
      fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datoffs,subint+1,1,dSet->head->nchan,&n_fval,datOffs,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datscl,subint+1,1,dSet->head->nchan,&n_fval,datScl,&initflag,&status);

      /*      for (i=0;i<dSet->head->nchan;i++)
	      {
	      printf("SCALING %d %g %g %d %d %d\n",i,datOffs[i],datScl[i],colnum_datoffs,colnum_datscl,status);
	      }
	      exit(1);
      */
      
      pfits_read1pol_float(data,pol,dSet,subint,subint,1,&nSamples,&nTime,&nFreq,debug,0);
      // Must apply SCL and OFFs parameters
      //
      for (j=0;j<nTime;j++)
	{
	  for (i=0;i<dSet->head->nchan;i++)
	    {
	      //     if (dSet->head->chanFreq[i] > 1134 && dSet->head->chanFreq[i] < 1136)
	      //	         if (dSet->head->chanFreq[i] > 1089 && dSet->head->chanFreq[i] < 1091)
	      //	      if (dSet->head->chanFreq[i] > 1017 && dSet->head->chanFreq[i] < 1019)
	      if (dSet->head->chanFreq[i] > 2349 && dSet->head->chanFreq[i] < 2351)
		{
		  fprintf(fout_ts,"%g %g %g ",subint*nTime*dSet->head->tsamp+j*dSet->head->tsamp,dSet->head->chanFreq[i],data[j*nFreq+i]*datScl[i]+datOffs[i]);
		}
	      if (subint==sub0 && (j == 64))
		{
		  fprintf(fout_samp,"%d %g %g %g\n",j,dSet->head->chanFreq[i],data[j*nFreq+i]*datScl[i]+datOffs[i],data[j*nFreq+i]);
		}
	      if (applyScaling==1)
		{
		  //		  printf("Scaling %g %g\n",datScl[i],datOffs[i]);
		  data[j*nFreq+i] = data[j*nFreq+i]*datScl[i]+datOffs[i];
		}
	      if (data[j*nFreq+i] > max[i])
		max[i] = data[j*nFreq+i];
	      if (data[j*nFreq+i] < min[i])
		min[i] = data[j*nFreq+i];
	      sum[i]+=data[j*nFreq+i];
	    }
	  totCount++;
	}
      if (subint==sub0)
	fprintf(fout_samp,"\n\n");
    }
  fclose(fout_ts);
  fclose(fout_samp);
  fout = fopen("output.dat","w");
  for (i=0;i<dSet->head->nchan;i++)
    {
      sum[i]/=(double)totCount;
      fprintf(fout,"%d %g %g %g %g spectrum\n",i,dSet->head->chanFreq[i],sum[i],max[i],min[i]);
      sum[i]=log10(sum[i]);
    }
  if (plot==1)
    {
      cpgbeg(0,grDev,1,1);
      cpgslw(2);
      cpgscf(2);
      cpgsvp(0.10,0.95,0.10,0.35);
      cpgswin(704,1344,2,9.6);
      cpgbox("ABCTSN",0,0,"ABCTSNL",0,0);
      cpglab("Frequency (MHz)","","");
      cpgline(dSet->head->nchan,dSet->head->chanFreq,sum);

      cpgsvp(0.10,0.95,0.40,0.65);
      //	cpgswin(1344,2368,miny,maxy);
      cpgswin(1344,2368,3,8);
      cpgbox("ABCTSN",0,0,"ABCTSNL",0,0);
      cpglab("","Signal strength","");
      cpgline(dSet->head->nchan,dSet->head->chanFreq,sum);
      
      cpgsvp(0.10,0.95,0.7,0.95);
      //	cpgswin(2368,4096,miny,maxy);
      cpgswin(2368,4096,3,8);
      cpgbox("ABCTNS",0,0,"ABCTSNL",0,0);
      cpgline(dSet->head->nchan,dSet->head->chanFreq,sum);
      cpgend();
    }

  
  fclose(fout);
  free(data);
  free(max);
  free(min);
  free(sum);
  free(sumSub);
  free(meanSub);
  free(datOffs);
  free(datScl);
  // De-allocate the memory
  deallocateMemory(&dSet,debug);

  
}
