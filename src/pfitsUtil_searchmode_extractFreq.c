// Code to extract one or more frequency bands from a specified PSRFITS search mode file
// gcc -lm -o pfitsUtil_searchmode_extractFreq pfitsUtil_searchmode_extractFreq.c -lcfitsio



// Note:
// This code is horrible.  The reason is that if I simply copy the original file and then
// try and shrink the file to remove the channels that I don't want then cfitsio still allocates
// enough memory for the original file (which can be enormous).
//
// This code is me attempting to circumvent cfitsio memory allocation and there is possibly
// a much better way to do this.
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
  int i,j,k;
  int hdu=1;
  char inname[1024]; 
  char outname[1024];
  char usename[1024];
  char keyname[128],val[128],comment[128];
  int nchan,npol,nsblk,nsubint,nbit;
  int outChan1 = 10;
  int outChan2 = 2000;
  int newNchan;
  float newBW;
  char nullVal = 0;
  float nullVal_f = 0;
  int initflag=0;
  
  int colnum_data_in=0;
  int colnum_data_out=0;
  int colnum_out_datFreq=0;
  int colnum_in_datFreq=0;
  int colnum_out_datWts=0;
  int colnum_in_datWts=0;
  int colnum_out_datScl=0;
  int colnum_in_datScl=0;
  int colnum_out_datOffs=0;
  int colnum_in_datOffs=0;
  
  char *dataIn;
  char *dataOut;
  float *colArrayIn;
  float *colArrayOut;
  
  int writePos;
  int readPos;

  long naxes[4];
  int naxis=3;
  char tdim[16];

  float meanFreq;
  int itemp;
  int nhdus,nhdus_in;
  int hdutype;
  int found;

  char checkName1[1024];
  char checkName2[1024];
  char checkName3[1024];
  char checkName4[1024];
  char checkName5[1024];
  char checkName6[1024];
  char **hduNames;

  char card[1024];
  char newcard[1024];
  char newVal[1024];
  int keytype;
  int oldBytes;

  float freq1,freq2;
  int   setFreq=0;
  
  printf("Running pfitsUtil_searchmode_extractFreq\n");
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0) // input file
	strcpy(inname,argv[++i]);
      else if (strcmp(argv[i],"-o")==0) // output file
	strcpy(outname,argv[++i]);
      else if (strcmp(argv[i],"-c1")==0) // First channel number for output
	sscanf(argv[++i],"%d",&outChan1);
      else if (strcmp(argv[i],"-c2")==0) // Last channel number for output
	sscanf(argv[++i],"%d",&outChan2);
      else if (strcmp(argv[i],"-f1")==0) 
	{
	  sscanf(argv[++i],"%f",&freq1);
	  setFreq=1;
	}
      else if (strcmp(argv[i],"-f2")==0) 
	sscanf(argv[++i],"%f",&freq2);
    }

  if (setFreq==0)
    newNchan= outChan2-outChan1+1; // SHOULD CHECK FOR BAND INVERSION
    
  sprintf(usename,"!%s",outname,inname);
  if ( !fits_open_file(&infptr, inname, READONLY, &status) )
    {
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
      fits_read_key(infptr,TINT,"NCHAN",&nchan,NULL,&status);
      fits_read_key(infptr,TINT,"NPOL",&npol,NULL,&status);
      fits_read_key(infptr,TINT,"NSBLK",&nsblk,NULL,&status);
      fits_read_key(infptr,TINT,"NAXIS2",&nsubint,NULL,&status);
      fits_read_key(infptr,TINT,"NAXIS1",&oldBytes,NULL,&status);
      fits_read_key(infptr,TINT,"NBITS",&nbit,NULL,&status);          


      if (npol != 1 || nbit != 8)
	{
	  printf("Sorry - I can only work with 1 pol and 8 bit data at the moment\n");
	  exit(1);
	}

      fits_get_colnum(infptr, CASEINSEN, "DATA", &colnum_data_in, &status);
      fits_get_colnum(infptr, CASEINSEN, "DAT_FREQ", &colnum_in_datFreq, &status);  
      fits_get_colnum(infptr, CASEINSEN, "DAT_WTS", &colnum_in_datWts, &status);  
      fits_get_colnum(infptr, CASEINSEN, "DAT_SCL", &colnum_in_datScl, &status);  
      fits_get_colnum(infptr, CASEINSEN, "DAT_OFFS", &colnum_in_datOffs, &status);  

      if (setFreq==1)
	{
	  float *freq;
	  freq = (float *)malloc(sizeof(float)*nchan);
	  fits_read_col(infptr,TFLOAT,colnum_in_datFreq,1,1,nchan,&nullVal_f,freq,&initflag,&status);
	  for (i=0;i<nchan-1;i++)
	    {
	      if (freq1 > freq[i] && freq1 <= freq[i+1])
		outChan1 = i;
	      if (freq2 > freq[i] && freq2 <= freq[i+1])
		outChan2 = i;
	    }
	  newNchan= outChan2-outChan1+1; // SHOULD CHECK FOR BAND INVERSION
	  printf("New channels from %g to %g are %d to %d and number of channels = %d\n",freq1,freq2,outChan1,outChan2,newNchan);
	  free(freq);
	}
    
      sprintf(checkName1,"TFORM%d",colnum_data_in);
      sprintf(checkName2,"TDIM%d",colnum_data_in);
      sprintf(checkName3,"TFORM%d",colnum_in_datFreq);
      sprintf(checkName4,"TFORM%d",colnum_in_datOffs);
      sprintf(checkName5,"TFORM%d",colnum_in_datWts);
      sprintf(checkName6,"TFORM%d",colnum_in_datScl);

      fits_movabs_hdu(infptr,1,NULL,&status);
      dataIn = (char *)malloc(sizeof(char)*nchan*nsblk*npol);
      dataOut = (char *)malloc(sizeof(char)*newNchan*nsblk*npol);
      colArrayIn = (float *)malloc(sizeof(float)*nchan*npol);
      colArrayOut = (float *)malloc(sizeof(float)*newNchan*npol);
      
      if ( !fits_create_file(&outfptr, usename, &status) )
	{
	  // Go through each HDU in turn and get the suitable data and meta-data
	  while( !fits_movabs_hdu(infptr, hdu, NULL, &status) )
	    {
	      printf("hdu = %d\n",hdu);
	      if (hdu==1) strcpy(keyname,"PRIMARY");
	      else {
		fits_read_key(infptr,TSTRING,"EXTNAME",keyname,comment,&status);
	      }
	      if (strcmp(keyname,"SUBINT")==0)
		{
		  // I'm doing this line by line because if we copy_header then
		  // it allocates all the memory required for the original file (and that can be huge)		  
		  int keysexist=-1;
		  int morekeys=-1;
		  int len;
		  char val[1024],comment[1024],key[1024];
		  
		  printf("Processing subint\n");
		  
		  fits_get_hdrspace(infptr,&keysexist,&morekeys,&status);
		  printf("keys = %d %d\n",keysexist,morekeys);
		  fits_create_hdu(outfptr,&status);
		  
		  for (i=1;i<=keysexist;i++)
		    {
		      fits_read_record(infptr,i,card,&status);
		      if (status) fits_report_error(stderr, status);

		      fits_get_keyname(card,key,&len,&status);
		      fits_parse_value(card,val,comment,&status);
		   		    		      
  		      if (strcmp(key,checkName1)==0)
			{
			  sprintf(newVal,"%dB",newNchan*npol*nsblk); // Update if not 8 bit

			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
			  printf("newcard: >%s<\n",newcard);
			  if (*comment)
			    {
			      strcat(newcard, " / ");
			      strcat(newcard, comment);
			    }
			  fits_parse_template(newcard,card,&keytype,&status);
			  strcpy(newcard,card);
			}
		      else if (strcmp(key,checkName2)==0)
			{
			  sprintf(newVal,"(%d,%d,%d)",newNchan,npol,nsblk);
			  
			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
			  printf("newcard: >%s<\n",newcard);
			  if (*comment)
			    {
			      strcat(newcard, " / ");
			      strcat(newcard, comment);
			    }
			  fits_parse_template(newcard,card,&keytype,&status);
			  strcpy(newcard,card);
			}
		      else if (strcmp(key,checkName3)==0) // DAT_FREQ
			{
			  sprintf(newVal,"%dD",newNchan);

			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
			  printf("newcard: >%s<\n",newcard);
			  if (*comment)
			    {
			      strcat(newcard, " / ");
			      strcat(newcard, comment);
			    }
			  fits_parse_template(newcard,card,&keytype,&status);
			  strcpy(newcard,card);
			}
		      else if (strcmp(key,checkName4)==0 || strcmp(key,checkName5)==0 || strcmp(key,checkName6)==0) // DAT_WTS, DAT_OFFS,DAT_SCL
			{
			  sprintf(newVal,"%dE",newNchan);

			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
			  printf("newcard: >%s<\n",newcard);
			  if (*comment)
			    {
			      strcat(newcard, " / ");
			      strcat(newcard, comment);
			    }
			  fits_parse_template(newcard,card,&keytype,&status);
			  strcpy(newcard,card);
			}
		      else if (strcmp(key,"NAXIS1")==0)
			{
			  long newBytes;

			  // Note some mix of floats and doubles and bytesq
			  newBytes = oldBytes - nchan*nsblk*npol - nchan*3*4 -nchan*1*8 + newNchan*npol*nsblk + newNchan*3*4 + newNchan*1*8;
			  printf("oldBytes = %d (%d) [%d] %d %d\n",oldBytes,nchan,newBytes,sizeof(float),sizeof(double));
			  //			  strcpy(newVal,"525608");
			  sprintf(newVal,"%d",newBytes);
			  
			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
			  printf("newcard: >%s<\n",newcard);
			  if (*comment)
			    {
			      strcat(newcard, " / ");
			      strcat(newcard, comment);
			    }
			  fits_parse_template(newcard,card,&keytype,&status);
			  strcpy(newcard,card);
			}
		      fits_write_record(outfptr,card,&status);		     
		      if (status) {fits_report_error(stderr, status); exit(1);}

		      if (strcmp(key,checkName2)==0)
			fits_modify_vector_len(outfptr,colnum_data_in,(long)((long)(newNchan)*(long)npol*(long)nsblk),&status); // Should * bytespersample for non 8 bit
		      
		    }
		      
		}
	      else
		{
		  printf("Copying\n");
		  fits_copy_hdu(infptr, outfptr, 0, &status);		  	      	      
		}
	      hdu++;
	    }
	}
      if (status == END_OF_FILE) status = 0;	  
    }

  fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
  fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);

  fits_get_colnum(outfptr, CASEINSEN, "DATA", &colnum_data_out, &status);  
  fits_get_colnum(outfptr, CASEINSEN, "DAT_FREQ", &colnum_out_datFreq, &status);
  fits_get_colnum(outfptr, CASEINSEN, "DAT_WTS", &colnum_out_datWts, &status);
  fits_get_colnum(outfptr, CASEINSEN, "DAT_SCL", &colnum_out_datScl, &status);
  fits_get_colnum(outfptr, CASEINSEN, "DAT_OFFS", &colnum_out_datOffs, &status);
  if (status) {fits_report_error(stderr, status); exit(1);}

  for (i=0;i<nsubint;i++)
    {
      printf("Processing sub-integration %d of %d\n",i+1,nsubint);
      fits_read_col(infptr,TBYTE,colnum_data_in,i+1,1,nsblk*nchan*npol,&nullVal,dataIn,&initflag,&status);
      for (j=0;j<nsblk;j++)
	{
	  writePos = j*newNchan; // Note npol = 1 (assumed ** must fix)
	  readPos = j*nchan+outChan1;
	  memcpy(dataOut+writePos,dataIn+readPos,newNchan);  
	}
      fits_write_col(outfptr,TBYTE,colnum_data_out,i+1,1,newNchan*nsblk*npol,dataOut,&status);
      fits_read_col(infptr,TFLOAT,colnum_in_datFreq,i+1,1,nchan,&nullVal_f,colArrayIn,&initflag,&status);

      // Calculate central frequency
      // Seemingly we need to modify OBS_FREQ in the header to get
      // dspsr and pav to show the correct frequencies even though the
      // frequencies are correct in the SUBINT table
      meanFreq=0;
      for (j=0;j<newNchan;j++)
	meanFreq+=colArrayIn[outChan1+j];
      meanFreq/=(double)newNchan;
      newBW = fabs(colArrayIn[outChan1]-colArrayIn[outChan1+newNchan-1]);
      
      fits_write_col(outfptr,TFLOAT,colnum_out_datFreq,i+1,1,newNchan,colArrayIn+outChan1,&status);

      fits_read_col(infptr,TFLOAT,colnum_in_datWts,i+1,1,nchan,&nullVal_f,colArrayIn,&initflag,&status);
      fits_write_col(outfptr,TFLOAT,colnum_out_datWts,i+1,1,newNchan,colArrayIn+outChan1,&status);

      fits_read_col(infptr,TFLOAT,colnum_in_datScl,i+1,1,nchan*npol,&nullVal_f,colArrayIn,&initflag,&status);
      for (j=0;j<npol;j++)
	{
	  for (k=0;k<newNchan;k++)
	    colArrayOut[j*newNchan+k] = colArrayIn[j*nchan+k];	    
	}
      fits_write_col(outfptr,TFLOAT,colnum_out_datScl,i+1,1,newNchan*npol,colArrayOut,&status);
      fits_read_col(infptr,TFLOAT,colnum_in_datOffs,i+1,1,nchan*npol,&nullVal_f,colArrayIn,&initflag,&status);
      for (j=0;j<npol;j++)
	{
	  for (k=0;k<newNchan;k++)
	    colArrayOut[j*newNchan+k] = colArrayIn[j*nchan+k];	    
	}
      fits_write_col(outfptr,TFLOAT,colnum_out_datOffs,i+1,1,newNchan*npol,colArrayOut,&status);
    } 

  fits_update_key(outfptr, TINT, (char *)"NCHAN", &newNchan, NULL, &status );    
  fits_update_key(outfptr, TFLOAT, (char *)"REFFREQ", &meanFreq, NULL, &status );

  // Back to primary header
  fits_movabs_hdu(outfptr,1,NULL,&status);
  fits_update_key(outfptr, TFLOAT, (char *)"OBSFREQ", &meanFreq, NULL, &status );
  fits_update_key(outfptr, TFLOAT, (char *)"OBSBW", &newBW, NULL, &status );
  fits_update_key(outfptr, TINT, (char *)"OBSNCHAN", &newNchan, NULL, &status );    

  
  
  free(dataIn);
  free(dataOut);
  free(colArrayIn);
  free(colArrayOut);
  

  fits_close_file(infptr, &status);
  fits_close_file(outfptr, &status);
  
  printf("Completed\n");
}
