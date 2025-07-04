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

//  gcc -lm -o pfits_frb pfits_frb.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -lpgplot
//gcc -lm -o pfits_frb pfits_frb.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

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
  int i,j,k;
  int nbin=1;
  float *plotArr;
  float *plotArr2;
  float *plotArrBin;
  int nTimeSamples;
  int nFreqSamples;
  long nSamples;
  float tr[6];
  float lineOffset=0;
  float lineOffset2=0;
  int   setLO2=0;
  int nchan;
  int nsblk;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float t1 = 733.74;
  float t2 = 733.90;
  float frbTime = 733.768224;
  float tdispCentre = 733.82;
  float tdispWidth = 0.02;
  float plotWidth;
  float dm=-1;
  long nalloc;
  float fx[4096],fy[4096];
  float miny,maxy;
  float fref=-1;
  float pos0,tdiff;
  int nDedispSamp;
  int idiff;
  float *sumSigX,*sumSigY;
  int j0;
  char dmStr[1024];
  char grDev[128]="/xs";
  char output[128];
  int outputSet=0;
  int publish=0;
  int colourScheme=0;
  FILE *fout;
  int csvOutput=0;
  char csvFile[128];
  int flip=0;
  int stats=0;

  float minVal = 0;
  float maxVal = 3;
  
  printf("Starting\n");
  
  // Allocate memory for these files
  //  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));

  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-dm")==0)
	sscanf(argv[++i],"%f",&dm);
      else if (strcmp(argv[i],"-csv")==0)
	{
	  csvOutput=1;
	  strcpy(csvFile,argv[++i]);
	}
      else if (strcmp(argv[i],"-flip")==0)
	flip=1;
      else if (strcmp(argv[i],"-frb")==0)
	sscanf(argv[++i],"%f",&frbTime);
      else if (strcmp(argv[i],"-w1")==0) // This is the width of the final plot
	sscanf(argv[++i],"%f",&tdispWidth);
      else if (strcmp(argv[i],"-w2")==0) // This is the data window that gets loaded in
	sscanf(argv[++i],"%f",&plotWidth);
      else if (strcmp(argv[i],"-nbin")==0)
	sscanf(argv[++i],"%d",&nbin);
      else if (strcmp(argv[i],"-lo")==0)
	sscanf(argv[++i],"%f",&lineOffset);
      else if (strcmp(argv[i],"-lo2")==0)
	{sscanf(argv[++i],"%f",&lineOffset2); setLO2=1;}
      else if (strcmp(argv[i],"-publish")==0)
	publish=1;
      else if (strcmp(argv[i],"-publish2")==0)
	publish=2;
      else if (strcmp(argv[i],"-colour")==0)
	sscanf(argv[++i],"%d",&colourScheme);
      else if (strcmp(argv[i],"-minVal")==0)
	sscanf(argv[++i],"%f",&minVal);
      else if (strcmp(argv[i],"-maxVal")==0)
	sscanf(argv[++i],"%f",&maxVal);
      else if (strcmp(argv[i],"-stats")==0)
	stats=1;
      else if (strcmp(argv[i],"-fref")==0)
	sscanf(argv[++i],"%f",&fref);
      else if (strcmp(argv[i],"-g")==0)
	strcpy(grDev,argv[++i]);
      else if (strcmp(argv[i],"-o")==0)
	{
	  outputSet=1;
	  strcpy(output,argv[++i]);
	}
    }
  sprintf(dmStr,"DM = %.2f cm\\u-3\\dpc",dm);
  
  // Open the files and load header information
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);
  printf("Nchan = %d\n",dSet->head->nchan);
  nchan = dSet->head->nchan;
  nsblk = dSet->head->nsblk;
  printf("nchan d = %d\n",(int)dSet->head->nchan);
  if (fref==-1)
    fref = (dSet->head->chanFreq[0]+dSet->head->chanFreq[dSet->head->nchan-1])/2.0;
  //  tdiff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[0]/1000.0,-2));
  tdispCentre = frbTime; // + tdiff;
  printf("Reference frequency = %g\n",fref);
  
  t1 = tdispCentre-plotWidth/2.0;
  t2 = tdispCentre+plotWidth/2.0;


  nalloc = (int)((float)((t2-t1)*nchan/dSet->head->tsamp) + 10000);
  plotArr = (float *)malloc(sizeof(float)*nalloc);
  plotArr2 = (float *)malloc(sizeof(float)*nalloc);
  if (!(plotArrBin = (float *)malloc(sizeof(float)*nalloc)))
    {
      printf("UNABLE TO ALLOCATE MEMORY\n");
      exit(1);
    }
  pfits_read1pol_float(plotArr,0,dSet,t1,t2,3,&nSamples,&nTimeSamples,&nFreqSamples,1,0);    
  if (flip==1)
    {
      for (i=0;i<nSamples;i++)
	{
	  if (plotArr[i] == 0) plotArr[i]=1;
	  else plotArr[i]=0;
	}
    }
  
  // Plot 1a: time versus frequency - no dedispersion
  //
  tr[0] = t1;
  tr[1] = 0;
  tr[2] = dSet->head->tsamp;
  tr[3] = dSet->head->chanFreq[0];  tr[4] = dSet->head->chanbw;  tr[5] = 0;
  cpgbeg(0,grDev,1,1);
  //  cpgsvp(0.1,0.7,0.55,0.90);
  if (publish>0)
    {
      cpgsch(1.2);
      cpgslw(2);
      cpgscf(2);
      nDedispSamp = (int)(tdispWidth/dSet->head->tsamp);
      
      printf("nDedispSamp = %d\n",nDedispSamp);
      printf("nchan = %d\n",(int)dSet->head->nchan);
      printf("nTimeSamples = %d\n",nTimeSamples);
      printf("nbin = %d\n",nbin);
      
      for (i=0;i<dSet->head->nchan;i++)
	{
	  for (j=0;j<nTimeSamples/nbin;j++)
	    {
	      //	      printf("j = %d -- adding at %d %d\n",j,j*dSet->head->nchan+i,nalloc);
	      plotArrBin[j*dSet->head->nchan+i]=0;
	      //	      printf("Before k\n");
	      for (k=0;k<nbin;k++)
		plotArrBin[j*dSet->head->nchan+i]+=plotArr[(nbin*j+k)*dSet->head->nchan+i];
	      plotArrBin[j*dSet->head->nchan+i]/=nbin;
	    }
	}

      nDedispSamp = (int)(tdispWidth/dSet->head->tsamp);
      //      printf("nDedispSamp = %d\n",nDedispSamp);
      sumSigX = (float *)malloc(sizeof(float)*nDedispSamp);
      sumSigY = (float *)malloc(sizeof(float)*nDedispSamp);
      //      printf("nbin = %d\n",nbin);
      for (j=0;j<nDedispSamp/nbin;j++)
	{
	  sumSigY[j]=0;
	  sumSigX[j]=tdispCentre-tdispWidth/2.+dSet->head->tsamp*j*nbin;
	}
      j0 = ((tdispCentre-tdispWidth/2.0-t1)/dSet->head->tsamp)/nbin;
      printf("Nchan = %d\n",dSet->head->nchan);
      for (i=0;i<dSet->head->nchan;i++)
	{
	  for (j=0;j<nDedispSamp/nbin;j++)
	    {
	      tdiff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[i]/1000.0,-2));
	      idiff = (int)(tdiff/dSet->head->tsamp/nbin);
	      //	      printf("Setting with i = %d idiff = %d j = %d j0 = %d\n",i,idiff,j,j0);
	      //	      printf("Array positions are %d and %d (%d %d %d %d %d)\n",j*dSet->head->nchan+i,
	      //		     (j0+j-idiff)*dSet->head->nchan+i,j0,j,idiff,dSet->head->nchan,i);
	      plotArr2[j*dSet->head->nchan+i] = plotArrBin[(j0+j-idiff)*dSet->head->nchan+i];
	      //	      printf("Have set plotArr2 - now sumSigY\n");
	      sumSigY[j]+=plotArr2[j*dSet->head->nchan+i];
	      //	      printf("Have set\n");
	    }	  
	}
      //      printf("POS 1\n");
      tr[0] = tdispCentre-tdispWidth/2.0;
      tr[1] = 0;
      tr[2] = dSet->head->tsamp*nbin;
      tr[3] = dSet->head->chanFreq[0];  tr[4] = dSet->head->chanbw;  tr[5] = 0;

      if (publish==1)
	cpgsvp(0.1,0.9,0.15,0.75);
      else if (publish==2)
	cpgsvp(0.1,0.9,0.45,0.75);
      cpgswin(tdispCentre-tdispWidth/2.,tdispCentre+tdispWidth/2.,tr[3]+10,tr[3]+nchan*tr[4]-10);
      if (publish==1){
	cpgbox("ABCTSN",0,0,"ABCTNS",0,0);
	cpglab("Time from start of observation (s)","Frequency (MHz)","");
      }
      else if (publish==2){
	cpgbox("ABCTS",0,0,"ABNCTS",0,0);
	cpglab("","Frequency (MHz)","");
      }
      else
	{
	  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
	  cpglab("Time from start of observation (s)","","");
	}
      
      printf("Using %d %d\n",nchan,nTimeSamples);
      if (colourScheme==0)
	{
	  cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
	  cpggray(plotArr2,nchan,nDedispSamp/nbin,1,nchan,1,nDedispSamp/nbin,minVal,maxVal,tr);
	}
      else if (colourScheme==1)
	{
	  cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
	  cpgimag(plotArr2,nchan,nDedispSamp/nbin,1,nchan,1,nDedispSamp/nbin,0,3,tr);
	}
      else if (colourScheme==2)
	{
	  heat_r[0] = 0.0; heat_r[1] = 0.1; heat_r[2] = 0.2; heat_r[3] =  0.3; heat_r[4] = 1.0;
	  heat_g[0] = 0.0; heat_g[1] = 0.1; heat_g[2] = 0.2; heat_g[3] =  0.3; heat_g[4] = 1.0;
	  heat_b[0] = 0.0; heat_b[1] = 0.1; heat_b[2] = 0.2; heat_b[3] =  0.3; heat_b[4] = 1.0;

	  cpgctab(heat_l,heat_r,heat_g,heat_b,5,-1.0,0.5);
	  cpgimag(plotArr2,nchan,nDedispSamp/nbin,1,nchan,1,nDedispSamp/nbin,0,3,tr);
	}
      else if (colourScheme==3)
	{
	  heat_r[0] = 0.0; heat_r[1] = 0.05; heat_r[2] = 0.1; heat_r[3] =  0.15; heat_r[4] = 1.0;
	  heat_g[0] = 0.0; heat_g[1] = 0.05; heat_g[2] = 0.1; heat_g[3] =  0.15; heat_g[4] = 1.0;
	  heat_b[0] = 0.0; heat_b[1] = 0.05; heat_b[2] = 0.1; heat_b[3] =  0.15; heat_b[4] = 1.0;

	  cpgctab(heat_l,heat_r,heat_g,heat_b,5,-1.0,0.5);
	  cpgimag(plotArr2,nchan,nDedispSamp/nbin,1,nchan,1,nDedispSamp/nbin,0,3,tr);
	}
      else if (colourScheme==4)
	{
	  heat_r[0] = 0.0; heat_r[1] = 0.1; heat_r[2] = 0.2; heat_r[3] =  0.3; heat_r[4] = 1.0;
	  heat_g[0] = 0.0; heat_g[1] = 0.1; heat_g[2] = 0.2; heat_g[3] =  0.3; heat_g[4] = 1.0;
	  heat_b[0] = 0.0; heat_b[1] = 0.1; heat_b[2] = 0.2; heat_b[3] =  0.3; heat_b[4] = 1.0;

	  cpgctab(heat_l,heat_r,heat_g,heat_b,5,-1.0,0.5);
	  cpgimag(plotArr2,nchan,nDedispSamp/nbin,1,nchan,1,nDedispSamp/nbin,0,5,tr);
	}
      if (csvOutput==1)
	{
	  FILE *fout;
	  fout = fopen(csvFile,"w");
	  for (i=0;i<nchan;i++)
	    {
	      for (j=0;j<nDedispSamp/nbin;j++)
		{
		  fprintf(fout,"%g",plotArr2[j*nchan+i]);
		  if (j < nDedispSamp/nbin-1)
		    fprintf(fout,", ");
		  else
		    fprintf(fout,"\n");
		}
	    }
	  fclose(fout);
	}

      if (stats==1)
	{
	  long npts=0;
	  long n1=0;
	  for (i=0;i<nchan;i++)
	    {
	      for (j=0;j<nDedispSamp/nbin;j++)
		{
		  printf("%d %d %g\n",i,j,plotArr2[j*nchan+i]);
		  npts++;
		  if (plotArr2[j*nchan+i] > 0.99 && plotArr2[j*nchan+i] < 1.01)
		    n1++;
		}
	    }
	  printf("Number of points = %d, number of 1s = %d\n",npts,n1);
	}

      
      for (i=0;i<nDedispSamp/nbin;i++)
	{
	  if (i==0)
	    {
	      miny = maxy = sumSigY[i];
	    }
	  else
	    {
	      if (miny > sumSigY[i]) miny = sumSigY[i];
	      if (maxy < sumSigY[i]) maxy = sumSigY[i];
	    }
	}
      //
      // Plotting the profile at the top
      cpgsvp(0.1,0.9,0.75,0.9);
      cpgswin(tdispCentre-tdispWidth/2.,tdispCentre+tdispWidth/2.,miny,maxy+(maxy-miny)*0.1);
      cpgbox("ABCTS",0,0,"ABCTS",0,0);
      //  cpglab("","",dmStr);
      cpgline(nDedispSamp/nbin,sumSigX,sumSigY);
      printf("Done\n");

      if (publish == 2)
	{
	  cpgsvp(0.1,0.9,0.15,0.44);	  
	  cpgswin(tdispCentre-tdispWidth/2.,tdispCentre+tdispWidth/2.,tr[3]+10,tr[3]+nchan*tr[4]-10);
	  cpgbox("ABCTNS",0,0,"ABNCTS",0,0);
	  cpglab("Time from start of observation (s)","Frequency (MHz)","");
	  for (i=0;i<dSet->head->nchan;i++)
	    {
	      for (j=0;j<nTimeSamples/nbin;j++)
		{
		  plotArrBin[j*dSet->head->nchan+i]=0;
		  for (k=0;k<nbin;k++)
		    plotArrBin[j*dSet->head->nchan+i]+=plotArr[(nbin*j+k)*dSet->head->nchan+i];
		  plotArrBin[j*dSet->head->nchan+i]/=nbin;
		}
	    }
	  tr[0] = t1;
	  tr[1] = 0;
	  tr[2] = dSet->head->tsamp*nbin;
	  tr[3] = dSet->head->chanFreq[0];  tr[4] = dSet->head->chanbw;  tr[5] = 0;
	  
	  cpgimag(plotArr,nchan,nTimeSamples,1,nchan,1,nTimeSamples,minVal,maxVal,tr);
	  //	  cpgimag(plotArrBin,nchan,nTimeSamples/2,1,nchan,1,nTimeSamples/2,minVal,maxVal,tr);
	}
    }
  else
    {
      cpgsvp(0.05,0.45,0.55,0.90);
      cpgswin(t1,t2,tr[3],tr[3]+nchan*tr[4]);
      cpgbox("ABCTS",0,0,"ABCTSN",0,0);
      cpglab("","Frequency (MHz)",dSet->fileName);
      if (colourScheme==0)
	{
	  cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
	}
      else if (colourScheme==1)
	{
	  cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
	}
      else if (colourScheme==2)
	{
	  heat_r[0] = 0.0; heat_r[1] = 0.1; heat_r[2] = 0.2; heat_r[3] =  0.3; heat_r[4] = 1.0;
	  heat_g[0] = 0.0; heat_g[1] = 0.1; heat_g[2] = 0.2; heat_g[3] =  0.3; heat_g[4] = 1.0;
	  heat_b[0] = 0.0; heat_b[1] = 0.1; heat_b[2] = 0.2; heat_b[3] =  0.3; heat_b[4] = 1.0;

	  cpgctab(heat_l,heat_r,heat_g,heat_b,5,-1.0,0.5);
	}
      else if (colourScheme==3)
	{
	  heat_r[0] = 0.0; heat_r[1] = 0.05; heat_r[2] = 0.1; heat_r[3] =  0.15; heat_r[4] = 1.0;
	  heat_g[0] = 0.0; heat_g[1] = 0.05; heat_g[2] = 0.1; heat_g[3] =  0.15; heat_g[4] = 1.0;
	  heat_b[0] = 0.0; heat_b[1] = 0.05; heat_b[2] = 0.1; heat_b[3] =  0.15; heat_b[4] = 1.0;

	  cpgctab(heat_l,heat_r,heat_g,heat_b,5,-1.0,0.5);
	}

      //      cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
      printf("Using %d %d\n",nchan,nTimeSamples);
      // Top left plot - shows the dispersed pulse on a 2-D freq-time plot
      cpgimag(plotArr,nchan,nTimeSamples,1,nchan,1,nTimeSamples,minVal,maxVal,tr);
      if (outputSet==1)
	{
	  fout = fopen(output,"w");
	  for (j=0;j<nTimeSamples;j++)
	    {
	      for (i=0;i<nchan;i++)
		fprintf(fout,"%d %d %g\n",i,j,plotArr[j*dSet->head->nchan+i]);
	      fprintf(fout,"\n");	 
	    }
	  fclose(fout);
	}
      
      
      // Overlay FRB analytic signal
      tdiff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[0]/1000.0,-2));
      tdispCentre = frbTime + tdiff;
      
      pos0 = tdispCentre; //(t1+t2)/2.0;
      for (i=0;i<dSet->head->nchan;i++)
	{
	  fx[i] = pos0 - 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[i]/1000.0,-2))+lineOffset;
	  fy[i] = dSet->head->chanFreq[i];
	}
      cpgsci(2); cpgslw(4);
      cpgline(dSet->head->nchan,fx,fy);
      if (setLO2 == 1)
	{
	  for (i=0;i<dSet->head->nchan;i++)
	    {
	      fx[i] = pos0 - 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[i]/1000.0,-2))+lineOffset2;
	      fy[i] = dSet->head->chanFreq[i];
	    }
	  cpgsci(2); cpgslw(4);
	  cpgline(dSet->head->nchan,fx,fy);      
	}
      
      cpgsci(1); cpgslw(1);
      
      
      // Plot 1b: time versus frequency - no dedispersion
      //
      for (i=0;i<dSet->head->nchan;i++)
	{
	  for (j=0;j<nTimeSamples/nbin;j++)
	    {
	      plotArrBin[j*dSet->head->nchan+i]=0;
	      for (k=0;k<nbin;k++)
		plotArrBin[j*dSet->head->nchan+i]+=plotArr[(nbin*j+k)*dSet->head->nchan+i];
	      plotArrBin[j*dSet->head->nchan+i]/=nbin;
	    }
	}
      
      
      tr[0] = t1;
      tr[1] = 0;
      tr[2] = dSet->head->tsamp*nbin;
      tr[3] = dSet->head->chanFreq[0];  tr[4] = dSet->head->chanbw;  tr[5] = 0;
      //  cpgsvp(0.1,0.7,0.15,0.54);

      // This is the bottom left plot
      cpgsvp(0.05,0.45,0.15,0.54);
      cpgswin(t1,t2,tr[3],tr[3]+nchan*tr[4]);
      cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
      cpglab("Time from start of observation (s)","Frequency (MHz)","");
      //      cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
      printf("Using %d %d\n",nchan,nTimeSamples/2);
      cpgimag(plotArrBin,nchan,nTimeSamples/nbin,1,nchan,1,nTimeSamples/nbin,minVal,maxVal,tr);
      // Overlay FRB analytic signal
      tdiff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[0]/1000.0,-2));
      tdispCentre = frbTime + tdiff;
      
      pos0 = tdispCentre; //(t1+t2)/2.0;
      for (i=0;i<dSet->head->nchan;i++)
	{
	  fx[i] = pos0 - 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[i]/1000.0,-2)) + lineOffset;
	  fy[i] = dSet->head->chanFreq[i];
	}
      cpgsci(2); cpgslw(4);
      cpgline(dSet->head->nchan,fx,fy);
      if (setLO2 == 1)
	{
	  for (i=0;i<dSet->head->nchan;i++)
	    {
	      fx[i] = pos0 - 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[i]/1000.0,-2))+lineOffset2;
	      fy[i] = dSet->head->chanFreq[i];
	    }
	  cpgsci(2); cpgslw(4);
	  cpgline(dSet->head->nchan,fx,fy);      
	}
      
      cpgsci(1); cpgslw(1);
      
      
      nDedispSamp = (int)(tdispWidth/dSet->head->tsamp);
      
      // Plot 2a: same, but with dedispersion
      sumSigX = (float *)malloc(sizeof(float)*nDedispSamp);
      sumSigY = (float *)malloc(sizeof(float)*nDedispSamp);
      
      //    tdispCentre = 733.82;
      printf("frbTime = %g %g\n",tdiff,tdispCentre);
    
      for (j=0;j<nDedispSamp;j++)
	{
	  sumSigY[j]=0;
	  sumSigX[j]=tdispCentre-tdispWidth/2.+dSet->head->tsamp*j;
	}
      j0 = (tdispCentre-tdispWidth/2.0-t1)/dSet->head->tsamp;
      printf("j0 = %d\n",j0);
      
      for (i=0;i<dSet->head->nchan;i++)
	{
	  for (j=0;j<nDedispSamp;j++)
	    {
	      tdiff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[i]/1000.0,-2));
	      idiff = (int)(tdiff/dSet->head->tsamp);
	      
	      //	  plotArr2[j+i*nDedispSamp] = plotArr[(j0+idiff)+i*nDedispSamp]; 
	      //	  plotArr2[j*dSet->head->nchan+i] = plotArr[(j0+idiff)*dSet->head->nchan+i];
	      plotArr2[j*dSet->head->nchan+i] = plotArr[(j0+j-idiff)*dSet->head->nchan+i];
	      sumSigY[j]+=plotArr2[j*dSet->head->nchan+i];
	    }
	}
      
      tr[0] = tdispCentre-tdispWidth/2.0;
      tr[1] = 0;
      tr[2] = dSet->head->tsamp;
      tr[3] = dSet->head->chanFreq[0];  tr[4] = dSet->head->chanbw;  tr[5] = 0;
      
      cpgsvp(0.52,0.75,0.15,0.75);
      
      cpgswin(tdispCentre-tdispWidth/2.,tdispCentre+tdispWidth/2.,tr[3],tr[3]+nchan*tr[4]);
      cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
      cpglab("Time from start of observation (s)","Frequency (MHz)","");
      //      cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
      printf("Using %d %d\n",nchan,nTimeSamples);
      cpgimag(plotArr2,nchan,nDedispSamp,1,nchan,1,nDedispSamp,minVal,maxVal,tr);
      
      
      for (i=0;i<nDedispSamp;i++)
	{
	  if (i==0)
	    {
	      miny = maxy = sumSigY[i];
	    }
	  else
	    {
	      if (miny > sumSigY[i]) miny = sumSigY[i];
	      if (maxy < sumSigY[i]) maxy = sumSigY[i];
	    }
	}
      
      cpgsvp(0.52,0.75,0.75,0.9);
      cpgswin(tdispCentre-tdispWidth/2.,tdispCentre+tdispWidth/2.,miny,maxy);
      cpgbox("ABCTS",0,0,"ABCTSN",0,0);
      cpglab("","",dmStr);
      cpgline(nDedispSamp,sumSigX,sumSigY);
  

      
      
      // Plot 2b: same, but with dedispersion
      sumSigX = (float *)malloc(sizeof(float)*nDedispSamp);
      sumSigY = (float *)malloc(sizeof(float)*nDedispSamp);
      
      for (j=0;j<nDedispSamp/nbin;j++)
	{
	  sumSigY[j]=0;
	  sumSigX[j]=tdispCentre-tdispWidth/2.+dSet->head->tsamp*j*nbin;
	}
      j0 = ((tdispCentre-tdispWidth/2.0-t1)/dSet->head->tsamp)/nbin;
      
      for (i=0;i<dSet->head->nchan;i++)
	{
	  for (j=0;j<nDedispSamp/nbin;j++)
	    {
	      tdiff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(dSet->head->chanFreq[i]/1000.0,-2));
	      idiff = (int)(tdiff/dSet->head->tsamp/nbin);
	  
	  plotArr2[j*dSet->head->nchan+i] = plotArrBin[(j0+j-idiff)*dSet->head->nchan+i];
	  sumSigY[j]+=plotArr2[j*dSet->head->nchan+i];
	}
    }

  tr[0] = tdispCentre-tdispWidth/2.0;
  tr[1] = 0;
  tr[2] = dSet->head->tsamp*nbin;
  tr[3] = dSet->head->chanFreq[0];  tr[4] = dSet->head->chanbw;  tr[5] = 0;

  cpgsvp(0.76,0.98,0.15,0.75);

  cpgswin(tdispCentre-tdispWidth/2.,tdispCentre+tdispWidth/2.,tr[3],tr[3]+nchan*tr[4]);
  cpgbox("ABCTSN",0,0,"ABCTS",0,0);
  cpglab("Time from start of observation (s)","","");
  //  cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
  printf("Using %d %d\n",nchan,nTimeSamples);
  cpgimag(plotArr2,nchan,nDedispSamp/nbin,1,nchan,1,nDedispSamp/nbin,minVal,maxVal,tr);


  for (i=0;i<nDedispSamp/nbin;i++)
    {
      if (i==0)
	{
	  miny = maxy = sumSigY[i];
	}
      else
	{
	  if (miny > sumSigY[i]) miny = sumSigY[i];
	  if (maxy < sumSigY[i]) maxy = sumSigY[i];
	}
    }
  
  cpgsvp(0.76,0.98,0.75,0.9);
  cpgswin(tdispCentre-tdispWidth/2.,tdispCentre+tdispWidth/2.,miny,maxy);
  cpgbox("ABCTS",0,0,"ABCTS",0,0);
  //  cpglab("","",dmStr);
  cpgline(nDedispSamp/nbin,sumSigX,sumSigY);
  
    }


  
  cpgend();
  

  // Do the plot
  // doPlot(dSet,nFiles,plot,debug);
 
  // Close the file
  //  pfitsCloseFile(dSet,debug);
  
  // De-allocate the memory
  deallocateMemory(&dSet,debug);  
  free(plotArr);
  free(plotArr2);
  free(plotArrBin);
  free(sumSigX);
  free(sumSigY);
}

