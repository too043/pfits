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

//  gcc -lm -o pfits_convertToBinary pfits_convertToBinary.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -lpgplot

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
  int i;
  int nFiles=0,nFiles_t=0;
  float sval=-1;
  float t1=-1,t2=-1;
  float s1=-1,s2=-1;
  int offScl=0;
  int nohead=0;
  int interact=1;
  int colnum;
  int status=0;
  int nchan;
  int samplesperbyte=1;
  unsigned char *cVals;
  int subint = 0;
  unsigned char nval = '0';
  int initflag=0;
  FILE *fout,*fout2;
  unsigned long int nsamp=0;
  unsigned long int totsamp=100;
  int format=1; // Binary
  int checkSaturate=0;
  
  // Allocate memory for these files
  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));
  initialise(&dSet,debug);
  // Count number of files
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)          setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-scale")==0) offScl=1;
      else if (strcmp(argv[i],"-nsamp")==0) sscanf(argv[++i],"%Ld",&totsamp);
      else if (strcmp(argv[i],"-binary")==0) format=1;
      else if (strcmp(argv[i],"-text")==0) format=2;
      else if (strcmp(argv[i],"-checkSaturate")==0) checkSaturate=1;
    }

  printf("Looking for saturation: %d\n",checkSaturate);
  // Open the files and load header information
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);
  
  //  plotArr_p0 = (float **)malloc(sizeof(float *)*nFiles);
  //  pfits_read1pol_float(plotArr_p0[i],0,dSet[i],plot[i]->t1,plot[i]->t2,plot[i]->rangeType,&(plot[i]->nSamples),&(plot[i]->nTimeSamples),&(plot[i]->nFrequencySamples),debug,offScl); 
  // Close the file
  //    pfitsCloseFile(dSet,debug);
  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum,&status);
  
  nchan = dSet->head->nchan;
  printf("Number of channels = %d\n",nchan);
  
  cVals = (unsigned char *)malloc(sizeof(unsigned char)*nchan/samplesperbyte); // Note:       reading one polarisation
  if (format==2)
    fout = fopen("rawData.dat","w");
  else if (format == 1)
    fout = fopen("rawData.dat","wb");
  printf("Number of samples to read = %d\n",totsamp);

  if (checkSaturate==1) fout2 = fopen("saturate.dat","w");
  
  for (nsamp = 0;nsamp < totsamp;nsamp++)
    {
      fits_read_col_byt(dSet->fp,colnum,
			subint+1,
			nsamp*nchan+1,nchan,nval,cVals,&initflag,&status);

      if (format==2)
	{
	  for (i=0;i<nchan;i++)
	    fprintf(fout,"%d %d %d\n",nsamp,i,cVals[i]);
	  fprintf(fout,"\n");
	}
      else if (format==1)
	fwrite(cVals,1,nchan,fout);
      if (checkSaturate==1)
	{
	  for (i=0;i<nchan;i++)
	    {
	      if (cVals[i] == 255)
		fprintf(fout2,"%d %d\n",nsamp,i);
	    }
	}
    }
  fclose(fout);
  if (checkSaturate==1) fclose(fout2);
  
  printf("De-allocating memory\n");
  deallocateMemory(&dSet,debug);  
  free(cVals);
  free(dSet);
}

