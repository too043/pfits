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


// /usr/bin/gcc -lm -o pfits_findSignal pfits_findSignal.c pfits_setup.c pfits_loader.c  -lcfitsio -lcpgplot

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

void checkAdd(unsigned long *storePos,unsigned long *n,unsigned long val,unsigned char *dataVals,int nchan,int chn,unsigned long max,unsigned long *treeStart,int base,unsigned long candBaseOldVal,unsigned long *candBaseNew);
void keepUnique(unsigned long *storePos,unsigned long *n);
int cmpfunc (const void * a, const void * b);
void makeCandidate(unsigned long startPos,unsigned long endPos,unsigned char *dataVals,int nchan,char *candName,unsigned long basePos,int nsblk,double tsamp);

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  unsigned long int i,j,k;
  long nSamp;
  int nTime,nFreq,nSub;
  int debug=0;
  unsigned char *dataVals;
  float *dataValsF;
  FILE *fout;
  long s1=0,s2=0;
  char fname[1024];
  long int nSamples;
  int nTimeSamples;
  int nFreqSamples;
  int scls=0;
  int width=3;
  unsigned long *storePos;
  unsigned long *storePos2;
  unsigned long *treeStart;
  unsigned long nStorePos=0;
  unsigned long nStorePos2=0;
  unsigned long *nStoreOld,*nStoreNew;
  unsigned long *candBase,*candBase2;
  unsigned long *candBaseOld,*candBaseNew;
  unsigned long *storePosOld,*storePosNew;
  int dirn=1;
  float *arr;
  float tr[6];
  float fx[2],fy[2];
  int cSubint[4096];
  unsigned long int sample0[4096];
  unsigned long int sample1[4096];
  int candSubint=0;
  int minChan=-1;
  int maxChan=-1;
  int plotSub=-1;
  int showCands=0;
  unsigned long max;
  float psx1 = -1;
  float psx2 = -1;
  float psy1 = -1;
  float psy2 = -1;
  int base=0;
  
  //  fout = fopen("sb421.txt","w");
  // Allocate memory for these files
  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcasecmp(argv[i],"-showcands")==0)
	showCands=1;
      else if (strcmp(argv[i],"-ps")==0)
	sscanf(argv[++i],"%d",&plotSub);
      else if (strcmp(argv[i],"-pszoomx")==0)
	{
	  sscanf(argv[++i],"%f",&psx1);
	  sscanf(argv[++i],"%f",&psx2);
	}
      else if (strcmp(argv[i],"-pszoomy")==0)
	{
	  sscanf(argv[++i],"%f",&psy1);
	  sscanf(argv[++i],"%f",&psy2);
	}
    }

  // Open the file and load header information
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);
  nSub = dSet->head->nsub;

  max = dSet->head->nsblk*dSet->head->nsub;
  
  dataValsF = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  arr = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  dataVals = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nsblk*dSet->head->nchan*dSet->head->nsub);
  storePos = (unsigned long *)malloc(sizeof(unsigned long)*dSet->head->nsblk*dSet->head->nchan*dSet->head->nsub); // Shouldn't need nchan here ??
  treeStart = (unsigned long *)malloc(sizeof(unsigned long)*dSet->head->nsblk*dSet->head->nchan*dSet->head->nsub); // Shouldn't need nchan here ??
  storePos2 = (unsigned long *)malloc(sizeof(unsigned long)*dSet->head->nsblk*dSet->head->nchan*dSet->head->nsub);
  candBase = (unsigned long *)malloc(sizeof(unsigned long)*dSet->head->nsblk*dSet->head->nchan*dSet->head->nsub); // Shouldn't need nchan here ??
  candBase2 = (unsigned long *)malloc(sizeof(unsigned long)*dSet->head->nsblk*dSet->head->nchan*dSet->head->nsub); // Shouldn't need nchan here ??


  for (j=0;j<nSub;j++)
    {
      pfits_read1pol_float(dataValsF,0,dSet,j,j,1,&nSamples,&nTimeSamples,&nFreqSamples,debug,scls);
      for (i=0;i<dSet->head->nsblk*dSet->head->nchan;i++)
	dataVals[j*dSet->head->nsblk*dSet->head->nchan+i] = 1-(unsigned char)dataValsF[i];
	//	dataVals[j*dSet->head->nsblk*dSet->head->nchan+i] = (unsigned char)dataValsF[i];
      if (j==plotSub)
	{
	  for (i=0;i<dSet->head->nchan;i++)
	    {
	      for (k=0;k<dSet->head->nsblk;k++)
		{
		  arr[k*dSet->head->nchan+i] = dataValsF[k*dSet->head->nchan+i];
		  //		  fprintf(fout,"%d %d %f\n",k,i,dataValsF[k*dSet->head->nchan+i]);
		}
	      //		  fprintf(fout,"\n");
	    }
	}
    }
  //  fclose(fout);
  printf("Loaded all data\n");
  nStorePos=0;
  for (i=0;i<((long int)nSub*(long int)dSet->head->nsblk)-3;i++) // FIX ME ON THE 3!
    {
      //printf("Have %d\n",(int)dataVals[i]);
      if (dataVals[i*dSet->head->nchan]==1 &&
	  dataVals[(i+1)*dSet->head->nchan]==1 &&
	  dataVals[(i+2)*dSet->head->nchan]==1)
	{
	  storePos[nStorePos] = i;
	  candBase[nStorePos++] = i;
	  storePos[nStorePos] = i+1;
	  candBase[nStorePos++] = i+1;
	  storePos[nStorePos] = i+2;
	  candBase[nStorePos++] = i+2;
	}
      
    }
  if (plotSub > -1)
    {
      cpgbeg(0,"/xs",1,1);
      //    cpgenv(0,4096,0,96,0,1);
      if (psx1==-1) psx1 = 0;
      if (psx2==-1) psx2 = dSet->head->nsblk;
      if (psy1==-1) psy1 = 0;
      if (psy2==-1) psy2 = dSet->head->nchan;
      cpgenv(psx1,psx2,psy1,psy2,0,1);
      //    cpgenv(2000,2700,0,96,0,1);
      //    cpgenv(2300,2600,40,80,0,1);
      //      cpgenv(2400,2500,0,3,0,1);
      tr[0] = 0; tr[1] = 0; tr[2] = 1; tr[3] = 0; tr[4] = 1; tr[5] = 0;
      printf("Making image\n");
      cpgimag(arr,dSet->head->nchan,dSet->head->nsblk,1,dSet->head->nchan,1,dSet->head->nsblk,0,1,tr);
      printf("Made image\n");

      cpgsci(3);
    
      //  fout = fopen("sb421.results","w");
      printf("%d Number of stored positions = %d\n",0,nStorePos);
      for (i=0;i<nStorePos;i++)
	{
	  if ((int)((float)storePos[i]/(float)dSet->head->nsblk) == plotSub)
	    {
	      //	  fprintf(fout,"%ld %d p1\n",(storePos[i]-plotSub*dSet->head->nsblk),0);
	      fx[0] = (float)(storePos[i]-plotSub*dSet->head->nsblk)+0.5+0.5;//0.5?
	      fy[0] = (float)0+0.5+0.5; // 0.0;
	      //	  printf("Plotting at %g %g\n",fx[0],fy[0]);
	      cpgpt(1,fx,fy,9);
	    }
	  //      if (storePos[i] > 1726812-50 && storePos[i] < 1726812+50)
	  //	printf("Positions = %ld %d\n",storePos[i],(int)((float)storePos[i]/(float)dSet->head->nsblk));
	}
    }
  if (dirn==1)
    {
      nStoreOld = &nStorePos;
      nStoreNew = &nStorePos2;
      storePosOld = storePos;
      storePosNew = storePos2;
      candBaseOld = candBase;
      candBaseNew = candBase2;
    }
  //  for (j=1;j<dSet->head->nchan;j++)
  if (minChan==-1)      minChan = 0;
  if (maxChan==-1)      maxChan=dSet->head->nchan;
  base=0;
  for (j=minChan+1;j<maxChan;j++) 
    {
      printf("j = %Ld\n",j);
      *nStoreNew=0;

      for (i=0;i<*nStoreOld;i++)
	{
	  if (dataVals[(storePosOld[i]+0L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]+1L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]+2L)*dSet->head->nchan+j]==1)
	    {
	      checkAdd(storePosNew,nStoreNew,storePosOld[i],dataVals,dSet->head->nchan,j,max,treeStart,base,candBaseOld[i],candBaseNew);
	    }
	  if (dataVals[(storePosOld[i]-1L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]+0L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]+1L)*dSet->head->nchan+j]==1)
	    {
	      checkAdd(storePosNew,nStoreNew,storePosOld[i]-1,dataVals,dSet->head->nchan,j,max,treeStart,base,candBaseOld[i],candBaseNew);
	      //	      printf("passing with i-1 = %Ld\n",i-1);
	      //A	      checkAdd(storePosNew,nStoreNew,storePosOld,i-1,max);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+1);
	    }	      
	  if (dataVals[(storePosOld[i]+1L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]+2L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]+3L)*dSet->head->nchan+j]==1)
	    {
	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+1,dataVals,dSet->head->nchan,j,max,treeStart,base,candBaseOld[i],candBaseNew);
	      //A	      checkAdd(storePosNew,nStoreNew,storePosOld,i+1,max);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+2);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+3);
	    }
	  if (dataVals[(storePosOld[i]-2L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]-1L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]+0L)*dSet->head->nchan+j]==1)
	    {
	      checkAdd(storePosNew,nStoreNew,storePosOld[i]-2,dataVals,dSet->head->nchan,j,max,treeStart,base,candBaseOld[i],candBaseNew);
	      //A	      checkAdd(storePosNew,nStoreNew,storePosOld,i-2,max);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]-1);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]);
	    }
	  if (dataVals[(storePosOld[i]+2L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]+3L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]+4L)*dSet->head->nchan+j]==1)
	    {
	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+2,dataVals,dSet->head->nchan,j,max,treeStart,base,candBaseOld[i],candBaseNew);
	      //    checkAdd(storePosNew,nStoreNew,storePosOld[i]+2);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+3);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+4);
	      } 
	  if (dataVals[(storePosOld[i]-3L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]-2L)*dSet->head->nchan+j]==1 &&
	      dataVals[(storePosOld[i]-1L)*dSet->head->nchan+j]==1)
	    {
	      checkAdd(storePosNew,nStoreNew,storePosOld[i]-3,dataVals,dSet->head->nchan,j,max,treeStart,base,candBaseOld[i],candBaseNew);
	      //A	      checkAdd(storePosNew,nStoreNew,storePosOld,i-3,max);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]-2);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]-1);
	    }
	  if  (dataVals[(storePosOld[i]+3L)*dSet->head->nchan+j]==1 &&
	       dataVals[(storePosOld[i]+4L)*dSet->head->nchan+j]==1 &&
	       dataVals[(storePosOld[i]+5L)*dSet->head->nchan+j]==1)
	    {
	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+3,dataVals,dSet->head->nchan,j,max,treeStart,base,candBaseOld[i],candBaseNew);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+3);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+4);
	      //	      checkAdd(storePosNew,nStoreNew,storePosOld[i]+5);
	      } 

	}
	  base=1;
      // Keep unique
      keepUnique(storePosNew,nStoreNew);
    
      printf("%d Number of stored positions = %d\n",j,(*nStoreNew));
      printf("Min/Max channels = %d %d\n",minChan,maxChan);
      //
      if (plotSub>-1)
	{
	  for (i=0;i<*nStoreNew;i++)
	    {
	      if (((int)((float)storePosNew[i]/(float)dSet->head->nsblk)) == plotSub)
		{
		  //		  fprintf(fout,"%ld %d p2\n",(storePosNew[i]-plotSub*dSet->head->nsblk),(int)j);
		  fx[0] = (float)(storePosNew[i]-plotSub*dSet->head->nsblk)+0.5+0.5; //Why second 0.5?
		  fy[0] = (float)j+0.5+0.5; // 0.0;
		  //		  printf("Plotting at %g %g\n",fx[0],fy[0]);
		  cpgpt(1,fx,fy,9);

		}
	      //	      if (storePosNew[i] > 1726812-50 && storePosNew[i] < 1726812+50)
	      if (*nStoreNew < 500)
		printf("Positions = %ld %d (%d, %d)\n",storePosNew[i],(int)((float)storePosNew[i]/(float)dSet->head->nsblk),j,(*nStoreNew));
	    }
	}

	dirn*=-1;
	if (dirn==1)
	{
	  nStoreOld = &nStorePos;
	  nStoreNew = &nStorePos2;
	  storePosOld = storePos;
	  storePosNew = storePos2;
	  candBaseOld = candBase;
	  candBaseNew = candBase2;
	}
      else
	{
	  nStoreOld = &nStorePos2;
	  nStoreNew = &nStorePos;
	  storePosOld = storePos2;
	  storePosNew = storePos;
	  candBaseOld = candBase2;
	  candBaseNew = candBase;
      
	}
    }
  cpgend();
  // Output final candidates
  candSubint=0;
  {
    int sb;
    int found;
    for (i=0;i<*nStoreOld;i++)
      {
	//	printf("Positions = %ld subint = %d\n",storePosOld[i],(int)((float)storePosOld[i]/(float)dSet->head->nsblk));
	sb = (int)((float)storePosOld[i]/(float)dSet->head->nsblk);
	found=-1;
	for (j=0;j<candSubint;j++)
	  {
	    if (cSubint[j] == sb) // Should actually look if candidate goes across subints
	      {
		found=j;
		break;
	      }
	  }
	if (found==-1)
	  {
	    cSubint[candSubint]=sb;
	    candSubint++;
	  }
      }
    for (i=0;i<candSubint;i++)
      printf("Candidate %d in subint %d\n",i,cSubint[i]);

    // Make images of the candidates
    // Need to store the starting point of each tree
    //
    if (showCands==1)
    {
      int gap=50; // Maximum size of gap before making a new image
      unsigned long startPos,endPos;
      char candName[1024];
      int candNumber=0;
      int startBase=0;
      startPos = storePosOld[0];
      startBase = candBaseOld[0];
      
      for (i=1;i<*nStoreOld;i++)
	{
	  //	printf("Positions = %ld subint = %d\n",storePosOld[i],(int)((float)storePosOld[i]/(float)dSet->head->nsblk));
	  if (storePosOld[i]-storePosOld[i-1] < gap)
	    endPos=storePosOld[i];
	  else
	    {
	      printf("Have %Ld %Ld\n",startPos,endPos);
	      sprintf(candName,"%s_%04d.png/png",dSet->fileName,candNumber);
	      makeCandidate(startPos,endPos,dataVals,dSet->head->nchan,candName,startBase,dSet->head->nsblk,dSet->head->tsamp);
	      candNumber++;
	      startPos=endPos=storePosOld[i];
	      startBase = candBaseOld[i];
	      printf("Setting startBase to be %Ld %Ld\n",startBase,i);
	    }
	}
      endPos=storePosOld[i-1];
      printf("Have %Ld %Ld\n",startPos,endPos);
      sprintf(candName,"%s_%04d.png/png",dSet->fileName,candNumber);
      makeCandidate(startPos,endPos,dataVals,dSet->head->nchan,candName,startBase,dSet->head->nsblk,dSet->head->tsamp);
      candNumber++;
    }
  }
  

  //  fclose(fout);
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
  free(dataVals);
  free(arr);
  free(dataValsF);
  free(storePos);
  free(treeStart);
  free(storePos2);
  free(candBase); free(candBase2);
}

void makeCandidate(unsigned long startPos,unsigned long endPos,unsigned char *dataVals,int nchan,char *candName,unsigned long basePos,int nsblk,double tsamp)
{
  float *arr;
  float tr[6];
  unsigned long nlen = (endPos-startPos)*32L;
  unsigned long width = endPos-startPos;
  unsigned long i,j;
  char title[1024];
  float fx[2],fy[2];

  arr = (float *)malloc(sizeof(float)*nlen*nchan);
  for (i=0;i<nlen;i++)
    {
      for (j=0;j<nchan;j++)
	arr[i*nchan+j] = (float)dataVals[(i+startPos-width*16)*nchan+j];
    }
  tr[0] = (startPos-width*16)*tsamp; tr[1] = 0; tr[2] = (float)tsamp; tr[3] = 0; tr[4] = 1; tr[5] = 0;
  printf("candName = >%s<\n",candName);

  //  if (strcmp(candName,"SMC021_00861.sf_18.ps/cps")==0)
    {
      cpgbeg(0,candName,1,1);
      printf("startPos = %Ld, endPos = %Ld\n",startPos,endPos);
      // BASEPOS NOT CORRECT BECAUSE OF THE QSORT
      printf("basePos = %Ld\n",basePos);
      printf("nlen = %Ld\n",nlen);
      printf("Width = %Ld\n",width);
      cpgenv(tr[0],tr[0]+nlen*tsamp,0,nchan,0,1);
      sprintf(title,"%s startSubint=%d",candName,(int)((startPos-width*16)/(double)nsblk));
      cpglab("Time from start (seconds)","Frequency channel",title);

      //      cpgimag(arr,nchan,nlen,1,nchan,1,nlen,0,1,tr);
      cpggray(arr,nchan,nlen,1,nchan,1,nlen,0,1,tr);
      // Draw a grid
      cpgsci(7);
      cpgslw(2);
      cpgsls(2);
      for (i=0;i<5;i++)
	{
	  fx[0]=fx[1] = tr[0]+nlen*tsamp/5.0*(i+1);
	  fy[0]=0; fy[1]=nchan;
	  cpgline(2,fx,fy);
	}
      for (i=0;i<5;i++)
	{
	  fx[0]=tr[0]; fx[1]= tr[0]+nlen*tsamp;
	  fy[0]=fy[1]=(i+1)*nchan/5.0; 
	  cpgline(2,fx,fy);
	}
      cpgslw(1);
      cpgsci(1);
      cpgsls(1);
      cpgend();
    }
  free(arr);
}

//	      checkAdd(storePosNew,nStoreNew,storePosOld[i],dataVals,dSet->head->nchan,j,max);
void checkAdd(unsigned long *storePos,unsigned long *n,unsigned long val,unsigned char *dataVals,int nchan,int chn,unsigned long max,unsigned long *treeStart,int base,unsigned long candBaseOldVal,unsigned long *candBaseNew)
{
  unsigned long int i;
  int maxLen = 20;
  int found=0;
  //  for (i=pos;i<*n;i++)

  //  	      checkAdd(storePosNew,nStoreNew,storePosOld[i],dataVals,dSet->head->nchan,j,max);
  i=0;
  do 
    {
      storePos[*n]=val+i;
      //      if (base==0)	
      //	candBaseOldVal = val;
      //      if (chn == 7)
      //	printf("Swapping with %Ld\n",candBaseOldVal);
      candBaseNew[*n] = candBaseOldVal;
      (*n)++;
      i++;
    } while (i < maxLen && val+i < max && dataVals[(val+i)*nchan+chn]==1);  // SHOULD MAKE MAXIMUM LENGTH MORE SENSIBLE
  //  	       dataVals[(storePosOld[i]+4L)*dSet->head->nchan+j]==1 &&
  //storePos[(*n)++]=val[pos];
}

void keepUnique(unsigned long *storePos,unsigned long *n)
{
  unsigned long int origN = *n;
  unsigned long int newN = 1;
  unsigned long int i,j;
  int keep=1;
  
  //  printf("Original n = %ld\n",origN);
  qsort(storePos, origN, sizeof(unsigned long int), cmpfunc);
  //  printf("Sorted\n");
  
  // SHOULD DO THIS WITH A SORT ** DO IT
  
  for (i=1;i<origN;i++)
    {
      if (storePos[newN-1]!=storePos[i])
	{
	  storePos[newN] = storePos[i];
	  newN++;
	}
    }
  //  printf("New n = %ld\n",newN);
  *n = newN; 
}

int cmpfunc (const void * a, const void * b) {
   return ( *(unsigned long int*)a - *(unsigned long int*)b );
}
