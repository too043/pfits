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

//  gcc -lm -o pfits_extract1bit pfits_extract1bit.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -lpgplot
//gcc -lm -o pfits_extract1bit pfits_extract1bit.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  char outfile[1024];
  int debug=0;
  int i,j,k;
  int nchan,nsblk;
  int npol;
  float t0,t1;
  float tsamp;
  FILE *fout;
  long s0,s1;
  int firstSamp;
  int lastSamp;
  int subint;
  int initflag=0;
  int i0,i1;
  int colnum;
  int status=0;
  int samplesperbyte;
  int nbits;
  unsigned char *cVals;
  unsigned char nval = '0';
  int polNum=0;
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outfile,argv[++i]);
      else if (strcmp(argv[i],"-t0")==0)
	sscanf(argv[++i],"%f",&t0);
      else if (strcmp(argv[i],"-t1")==0)
	sscanf(argv[++i],"%f",&t1);
    }
  
  // Open the files and load header information
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);
  nchan = dSet->head->nchan;
  nsblk = dSet->head->nsblk;
  tsamp = dSet->head->tsamp;
  nbits = dSet->head->nbits;
  npol = dSet->head->npol;
  
  samplesperbyte = 8/nbits;
    
  printf("Tsamp = %g\n",tsamp);
  
  s0 = (long)(int)(t0/(tsamp*nsblk));
  s1 = (long)(int)(t1/(tsamp*nsblk));
  firstSamp = (int)((t0-s0*tsamp*nsblk)/tsamp);
  lastSamp = (int)((t1-s1*tsamp*nsblk)/tsamp+0.5);
  printf("Loading from (%d,%d) to (%d,%d)\n",s0,firstSamp,s1,lastSamp);
  
  // Find the relevant bit of the file and dump to the output file
  fout = fopen(outfile,"wb");
  fwrite(dSet->fileName,sizeof(char),FILE_LEN,fout);
  fwrite(&nchan,sizeof(int),1,fout);
  fwrite(dSet->head->chanFreq,sizeof(int),nchan,fout);
  fwrite(&tsamp,sizeof(float),1,fout);
  fwrite(&(dSet->head->stt_imjd),sizeof(int),1,fout);
  fwrite(&(dSet->head->stt_smjd),sizeof(float),1,fout);
  fwrite(&(dSet->head->stt_offs),sizeof(float),1,fout);
  printf("Writing t0 and t1\n");
  fwrite(&t0,sizeof(float),1,fout);
  fwrite(&t1,sizeof(float),1,fout);

  fwrite(&(dSet->head->nbits),sizeof(int),1,fout);
  fwrite(&(dSet->head->npol),sizeof(int),1,fout);
  printf("Reading the file\n");
  
  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum,&status);
  cVals = (unsigned char *)malloc(sizeof(unsigned char)*nchan/samplesperbyte); // Note:       reading one polarisation
  
  for (subint=s0;subint<=s1;subint++)
    {
      printf("Processing subint %d\n",subint);
      if (subint==s0)
	i0 = firstSamp;
      else
	i0 = 0;
      if (subint==s1)
	i1 = lastSamp;
      else
	i1 = nsblk;
      for (i=i0;i<i1;i++)
	{
	  fits_read_col_byt(dSet->fp,colnum,
			    subint+1,
			    1+i*npol*nchan/samplesperbyte+polNum*nchan/samplesperbyte,
			    nchan/samplesperbyte,
			    nval,cVals,&initflag,&status);
	  fwrite(cVals,sizeof(unsigned char),nchan/samplesperbyte,fout);
	}
      
    }
  fclose(fout);
  
 
  // Close the file
  //  pfitsCloseFile(dSet,debug);
  
  // De-allocate the memory
  deallocateMemory(&dSet,debug);  
}

