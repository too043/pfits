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

//  gcc -lm -o pfits_makeImages pfits_makeImages.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -lpgplot
//gcc -lm -o pfits_makeImages pfits_makeImages.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

#define MAX_BEAM 13

int main(int argc,char *argv[])
{
  dSetStruct **dSet;
  char str[1024];
  int debug=0;
  int i,j,k;
  int kk1,kk2;
 float **plotArr;
  int nchan;
  int nsblk;
  int nsub;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  int stats=0;
  long int nSamples;
  int nTimeSamples,nFreqSamples;
  float tr[6];
  char fname[1024];
  int csv=0;
  FILE *fout;
  int maxNsub = -1;
  int nbeam=1;
  float delta=0.01;
  int plotTimeSamples = 128;
  int curPos=0;
  int subintCount=0;
  int ringPosI=0;
  int beamNum=0;
  int screen=0;
  long imgNum=0;
  float **ringBuffer;
  float **plotBuffer;
  int  shiftSample=1;
  
  char temp[1024];
  
  printf("Starting\n");
  
  // Allocate memory for these files
  dSet = (dSetStruct **)malloc(sizeof(dSetStruct *)*MAX_BEAM);
  
  // Initialise everything
  for (i=0;i<MAX_BEAM;i++)    
    initialise(&(dSet[i]),debug);
  printf("Completed initialisation");

  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet[beamNum++],debug);
      else if (strcmp(argv[i],"-csv")==0)
	csv=1;
      else if (strcmp(argv[i],"-maxNsub")==0)
	sscanf(argv[++i],"%d",&maxNsub);
      else if (strcmp(argv[i],"-screen")==0)
	screen=1;
      else if (strcmp(argv[i],"-shiftSample")==0)
	sscanf(argv[++i],"%d",&shiftSample);
    }
  nbeam = beamNum;
  
  // Open the files and load header information
  for (i=0;i<beamNum;i++)
    {
      pfitsOpenFile(dSet[i],debug);
      pfitsLoadHeader(dSet[i],debug);
    }
  nchan = dSet[0]->head->nchan;
  nsblk = dSet[0]->head->nsblk;
  nsub  = dSet[0]->head->nsub;

  tr[0] = 0;  tr[1] = 0;  tr[2] = 1;
  tr[3] = 0;  tr[4] = 1;  tr[5] = 0;

  plotArr = (float **)malloc(sizeof(float *)*nbeam);
  ringBuffer = (float **)malloc(sizeof(float *)*nbeam);
  plotBuffer = (float **)malloc(sizeof(float *)*nbeam);
  for (i=0;i<nbeam;i++)
    {
      plotArr[i]=(float *)malloc(sizeof(float)*nsblk*nchan);
      ringBuffer[i]=(float *)malloc(sizeof(float)*plotTimeSamples*nchan);
      plotBuffer[i]=(float *)malloc(sizeof(float)*plotTimeSamples*nchan);
    }
  printf("Nchan = %d, nsblk = %d, nsub = %d\n",nchan,nsblk,nsub);
  //  subintCount=421;
  
  for (j=0;j<nbeam;j++)
    {
      pfits_read1pol_float(plotArr[j],0,dSet[j],subintCount,subintCount,1,&nSamples,&nTimeSamples,&nFreqSamples,debug,0);
      
      for (i=0;i<plotTimeSamples*nchan;i++)
	ringBuffer[j][i] = plotArr[j][curPos*nchan+i];
    }
  ringPosI = 0;
  
  //  while (curPos+plotTimeSamples <= nsblk*nsub)
  while (curPos+plotTimeSamples <= nsblk*1)
    {
      printf("Current position = %d, subint = %d\n",curPos,subintCount);
      for (j=0;j<nbeam;j++)
	{
	  memcpy(plotBuffer[j],ringBuffer[j]+ringPosI*nchan,(plotTimeSamples*nchan-ringPosI*nchan)*sizeof(float));
	  if (ringPosI > 0)
	    memcpy(plotBuffer[j]+(plotTimeSamples*nchan-ringPosI*nchan),ringBuffer[j],(ringPosI*nchan)*sizeof(float));
	}
      if (screen==1)
	cpgbeg(0,"1/xs",1,1);
      else
	{
	  sprintf(fname,"movieImg/image_%06d.png/png",imgNum);
	  cpgbeg(0,fname,1,1);	  
	}
      cpgsvp(0,0.15,0.9,1.0);
      cpgswin(0,1,0,1);
      sprintf(str,"%.3f sec",curPos*dSet[0]->head->tsamp);
      cpgtext(0.1,0.1,str);

      if (nbeam==1)
	{
	  cpgsvp(0.1,0.9,0.1,0.8); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTSN",0,0,"ABCTSN",0,0);  // Beam 1
	  cpggray(plotBuffer[0],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}
      else
	{
	  cpgsvp(0.4+delta,0.6-delta,0.4+delta,0.6-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 1
	  cpggray(plotBuffer[0],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}
      if (nbeam >= 2)
	{
	  cpgsvp(0.2+delta,0.4-delta,0.25+delta,0.45-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 2
	  cpggray(plotBuffer[1],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}
      if (nbeam >= 3)
	{
	  cpgsvp(0.4+delta,0.6-delta,0.15+delta,0.35-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 3
	  cpggray(plotBuffer[2],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}

      
      if (nbeam >= 4)
	{      
	  cpgsvp(0.6+delta,0.8-delta,0.25+delta,0.45-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 4
	  cpggray(plotBuffer[3],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}

      if (nbeam >= 5)
	{
	  cpgsvp(0.6+delta,0.8-delta,0.55+delta,0.75-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 5
	  cpggray(plotBuffer[4],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}

      if (nbeam >= 6)
	{
	  cpgsvp(0.4+delta,0.6-delta,0.65+delta,0.85-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 6
	  cpggray(plotBuffer[5],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}

      if (nbeam >= 7)
	{
	  cpgsvp(0.2+delta,0.4-delta,0.55+delta,0.75-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 7
	  cpggray(plotBuffer[6],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}

      if (nbeam >= 8)
	{
	  cpgsvp(0.0+delta,0.2-delta,0.4+delta,0.6-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 8
	  cpggray(plotBuffer[7],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}

      if (nbeam >= 9)
	{
	  cpgsvp(0.2+delta,0.4-delta,0.0+delta,0.2-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 9
	  cpggray(plotBuffer[8],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}

      if (nbeam >= 10)
	{
	  cpgsvp(0.6+delta,0.8-delta,0.0+delta,0.2-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 10
	  cpggray(plotBuffer[9],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}

      if (nbeam >= 11)
	{
	  cpgsvp(0.8+delta,1.0-delta,0.4+delta,0.6-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 11
	  cpggray(plotBuffer[10],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}

      if (nbeam >= 12)
	{
	  cpgsvp(0.6+delta,0.8-delta,0.8+delta,1.0-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 12
	  cpggray(plotBuffer[11],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}

      if (nbeam >= 13)
	{
	  cpgsvp(0.2+delta,0.4-delta,0.8+delta,1.0-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 13
	  cpggray(plotBuffer[12],nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	}
      //      scanf("%s",&temp);
      cpgend();

      for (k=0;k<shiftSample;k++)
	{
	  curPos+=1;
	  ringPosI++;
	  if (ringPosI == plotTimeSamples) ringPosI = 0;
	  if (curPos==nsblk)
	    {
	      subintCount++;
	      curPos=0;
	      for (i=0;i<nbeam;i++)
		pfits_read1pol_float(plotArr[i],0,dSet[0],subintCount,subintCount,1,&nSamples,&nTimeSamples,&nFreqSamples,debug,0);
	    }
	  
	  for (j=0;j<nbeam;j++)
	    {
	      for (i=0;i<nchan;i++)
		{
		  if ((ringPosI) == plotTimeSamples)
		    ringBuffer[j][(ringPosI)*nchan+i]=plotArr[j][(curPos)*nchan+i];
		}
	    }
	}
      imgNum+=1;
      //      exit(1);
    }

 

  /*  if (maxNsub < 0) maxNsub = nsub;
  for (i=0;i<maxNsub;i++)
    {
      pfits_read1pol_float(plotArr,0,dSet,i,i,1,&nSamples,&nTimeSamples,&nFreqSamples,debug);
      
      if (csv==1)
	{
	  sprintf(fname,"image_%05d.csv",i);
	  fout = fopen(fname,"w");
 	  for (k=0;k<nFreqSamples;k++)
	    {
	      for (j=0;j<nTimeSamples;j++)
		{
		  if (j < nTimeSamples-1)
		    fprintf(fout,"%d,",(int)plotArr[j*nFreqSamples+k]);
		  else
		    fprintf(fout,"%d\n",(int)plotArr[j*nFreqSamples+k]);	      
		}
	      //	      fprintf(fout,"\n");
	    }
	  fclose(fout);
	}
      else
	{
	  //	  sprintf(fname,"image_%05d.gif/gif",i);
	  sprintf(fname,"/xs");
	  cpgbeg(0,fname,1,1);
	  if (nbeam==1)
	    {
	      cpgenv(0,nTimeSamples,0,nFreqSamples,0,1);
	      cpggray(plotArr,nchan,nTimeSamples,1,nchan,1,nTimeSamples,0,1,tr);
	    }
	  else if (nbeam==13)
	    {
	      cpgsvp(0.4+delta,0.6-delta,0.4+delta,0.6-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 1
	      cpggray(ringBuffer,nchan,plotTimeSamples,1,nchan,1,plotTimeSamples,0,1,tr);
	      
	      cpgsvp(0.2+delta,0.4-delta,0.25+delta,0.45-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 2
	      cpgsvp(0.4+delta,0.6-delta,0.15+delta,0.35-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 3
	      cpgsvp(0.6+delta,0.8-delta,0.25+delta,0.45-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 4
	      cpgsvp(0.6+delta,0.8-delta,0.55+delta,0.75-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 5
	      cpgsvp(0.4+delta,0.6-delta,0.65+delta,0.85-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 6
	      cpgsvp(0.2+delta,0.4-delta,0.55+delta,0.75-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 7
	      cpgsvp(0.0+delta,0.2-delta,0.4+delta,0.6-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 8
	      cpgsvp(0.2+delta,0.4-delta,0.0+delta,0.2-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 9
	      cpgsvp(0.6+delta,0.8-delta,0.0+delta,0.2-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 10
	      cpgsvp(0.8+delta,1.0-delta,0.4+delta,0.6-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 11
	      cpgsvp(0.6+delta,0.8-delta,0.8+delta,1.0-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 12
	      cpgsvp(0.2+delta,0.4-delta,0.8+delta,1.0-delta); cpgswin(0,plotTimeSamples,0,nFreqSamples);  cpgbox("ABCTS",0,0,"ABCTS",0,0);  // Beam 13
	    }
	  cpgend();
	}
    }
  */
  for (i=0;i<nbeam;i++)
    {
      free(plotArr[i]);
      free(ringBuffer);
      free(plotBuffer);
    }
  free(plotArr);
  free(ringBuffer);
  free(plotBuffer);
}

