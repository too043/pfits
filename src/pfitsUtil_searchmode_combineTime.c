// Code to combine two or more PSRFITS files in the time direction to
// produce a single output file.
//
// Usage:
// ./pfitsUtil_searchmode_combineTime -o <output name> file1 file2 file3 file4 file5
//
// This software will work with PSRFITS fold mode files and PSRFITS search mode files
//
// Version 1: 5rd May 2018, G. Hobbs
// Version 2: 5th Dec 2018, G. Hobbs - fixed OFFS_SUB table
//
// gcc -lm -o pfitsUtil_searchmode_combineTime pfitsUtil_searchmode_combineTime.c -lcfitsio
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "fitsio.h"
#include <stdint.h>

#define VERSION 1.0

int main(int argc,char *argv[])
{
  fitsfile *infptr,*outfptr;
  int status=0;
  char outname[1024];
  char usename[1024];
  char inname[1000][1024];
  int colnum = 0;
  int nIn=0;  
  int ii,jj;
  int i,j,k,p;
  long inrows,outrows;
  int hdu=1;
  int extraRows=0;
  long width;
  int outrowNum=0;
  unsigned char *buffer=0;
  long nrows = 0;
  double nulldouble = 0.0;
  int initflag = 0;
  double *offs_subs;
  double offs_sub_step;
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(outname,argv[++i]);
      else
	{
	  strcpy(inname[nIn],argv[i]);
	  nIn++;
	  printf("nIn = %d\n",nIn);
	}
    }
  sprintf(usename,"!%s",outname);
  
  printf("fitsTimeCombine version %d\n",(int)(VERSION));
  // Copy the first file to the output file
  /* Open the input file */
  printf("Loading %s\n",inname[0]);
  if ( !fits_open_file(&infptr, inname[0], READONLY, &status) )
    {
      /* Create the output file */
      if ( !fits_create_file(&outfptr, usename, &status) )
	{
	  /* Copy every HDU until we get an error */
	  while( !fits_movabs_hdu(infptr, hdu++, NULL, &status) )
	    {
	      fits_copy_hdu(infptr, outfptr, 0, &status);
	      //	      if (hdu == 6) break; // Don't go past SUBINT --- FIX THE HARDCODE HERE // Doesn't seem to save time (??)
	    }
	      /* Reset status after normal error */
	  if (status == END_OF_FILE) status = 0;	  
	}
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
      fits_read_key(infptr, TLONG, "NAXIS1", &width, NULL, &status);
      if (!(buffer = (unsigned char*)malloc(width)))
	{
	  printf("Memory allocation error\n");
	  exit(1);
	}
      fits_close_file(infptr, &status);
    }
  fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);

  // Count how many extra rows are required.
  // *** This bit could be skipped if we know we have 1 subint per file ***
  //
  fits_get_num_rows(outfptr,  &outrows,  &status);
  outrowNum = outrows;
  printf("outrowNum = %d\n",outrowNum);
  
  for (i=1;i<nIn;i++)
    {
      fits_open_file(&infptr,  inname[i], READONLY,  &status) ; // Doing a lot of file opening and closing here. Could improve this.
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
      fits_get_num_rows(infptr,  &inrows,  &status);
      extraRows+=inrows;
      fits_close_file(infptr, &status); 
    }
  // Insert this many extra rows
  printf("Inserting %d extra rows\n",extraRows);
  fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
  // Try moving this down and see if it speeds up  --- it does
  //  fits_insert_rows(outfptr,outrows,extraRows,&status);

  if (status) {fits_report_error(stderr, status); exit(1);}

  // Now join the files
  for (i=1;i<nIn;i++)
    {
      printf("Processing file %d %s\n",i+1,inname[i]);
      fits_open_file(&infptr,  inname[i], READONLY,  &status) ;
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
      fits_get_num_rows(infptr,  &inrows,  &status);
      fits_get_num_rows(outfptr,  &outrows,  &status);
      printf("Have %d %d %d\n",inrows,outrows,outrowNum);
      for (ii=1, jj=outrowNum+1;ii<=inrows;ii++,jj++)
	{
	  fits_insert_rows(outfptr,jj-1,1,&status);
	  fits_read_tblbytes(infptr,ii,1,width,buffer,&status);
	  fits_write_tblbytes(outfptr,jj,1,width,buffer,&status);
	}
      outrowNum=jj-1;
      printf("Completed data copy\n");
      fits_close_file(infptr, &status);
      if (status) {fits_report_error(stderr, status); exit(1);}
      } 

  
    /* if error occured, print out error message */
  if (status) fits_report_error(stderr, status);

  // Now sort out the OFFS_SUB parameters
  fits_get_num_rows(outfptr, &nrows, &status);
  fits_get_colnum (outfptr, CASEINSEN, "OFFS_SUB", &colnum, &status);
  offs_subs = (double *)malloc(sizeof(double)*nrows);
  fits_read_col(outfptr, TDOUBLE, colnum, 1, 1, nrows, &nulldouble, &(offs_subs[0]), &initflag, &status);

  offs_sub_step = offs_subs[0] * 2;

  for (i = 1; i < nrows; ++i) {
    offs_subs[i] = offs_subs[i-1] + offs_sub_step;
  }

  fits_write_col(outfptr, TDOUBLE, colnum, 1, 1, nrows, &(offs_subs[0]),&status);
  
  // TO DO:
  // Put the last HDUs back

  
  printf("Output file is %s\n",outname);
  
  if (buffer) free(buffer);
  
  printf("Closing the output file\n");
  // Now close the file
  fits_close_file(outfptr,  &status);
  if (status) fits_report_error(stderr, status);
  
}
