// gcc -lm -o pfitsUtil_searchmode_filter pfitsUtil_searchmode_filter.c -lcfitsio

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
  int i,j,k,k0,ii,jj;
  int hdu=1;
  char inname[1024]; 
  char outname[1024];
  char usename[1024];
  char keyname[128],val[128],comment[128];
  int nchan,npol,nsblk,nsubint,nbit;
  unsigned char nullVal = 0;
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
  
  unsigned char *dataIn;
  unsigned char *dataOut;
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

  int   samplesperbyte=1;

  int opol=-1;
  int outNpol;
  int outNbits;
  
  int setFilterHigh=0;
  int filterHighBound;
  int ival;

  int outSamplesPerByte=1;
  int obit = -1;
  
  printf("Running pfitsUtil_searchmode_extractFreq\n");
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)       // input file
	strcpy(inname,argv[++i]);
      else if (strcmp(argv[i],"-o")==0)  // output file
	strcpy(outname,argv[++i]);
      else if (strcmp(argv[i],"-op")==0) // output number of polarisation
	sscanf(argv[++i],"%d",&opol);
      else if (strcmp(argv[i],"-fh")==0) // filter high
	{
	  setFilterHigh = 1;
	  sscanf(argv[++i],"%d",&filterHighBound);
	}
      else if (strcmp(argv[i],"-ob")==0) // Output number of bits
	sscanf(argv[++i],"%d",&obit);
    }

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

      samplesperbyte = 8/nbit;
      printf("samplesperbyte = %d\n",samplesperbyte);
      printf("npol = %d\n",npol);
      if ((npol != 1 || nbit != 8) && (npol != 1 || nbit != 2) && (npol != 4 || nbit != 8))
	{
	  printf("Sorry - I can only work with (1 pol and 8 bit) or (1 pol and 2 bit) data or (4 pol and 8 bit) data at the moment\n");
	  exit(1);
	}
      if (opol==-1)
	  outNpol = npol;
      else
	outNpol = opol;

      if (obit==-1)
	outSamplesPerByte = samplesperbyte;
      else
	outSamplesPerByte = 8/obit;

      outNbits = 8/outSamplesPerByte;
      
      fits_get_colnum(infptr, CASEINSEN, "DATA", &colnum_data_in, &status);
      fits_get_colnum(infptr, CASEINSEN, "DAT_FREQ", &colnum_in_datFreq, &status);  
      fits_get_colnum(infptr, CASEINSEN, "DAT_WTS", &colnum_in_datWts, &status);  
      fits_get_colnum(infptr, CASEINSEN, "DAT_SCL", &colnum_in_datScl, &status);  
      fits_get_colnum(infptr, CASEINSEN, "DAT_OFFS", &colnum_in_datOffs, &status);  

    
      sprintf(checkName1,"TFORM%d",colnum_data_in);
      sprintf(checkName2,"TDIM%d",colnum_data_in);
      sprintf(checkName3,"TFORM%d",colnum_in_datFreq);
      sprintf(checkName4,"TFORM%d",colnum_in_datOffs);
      sprintf(checkName5,"TFORM%d",colnum_in_datWts);
      sprintf(checkName6,"TFORM%d",colnum_in_datScl);

      fits_movabs_hdu(infptr,1,NULL,&status);
      dataIn = (unsigned char *)malloc(sizeof(unsigned char)*nchan*nsblk*npol);
      dataOut = (unsigned char *)malloc(sizeof(unsigned char)*nchan*nsblk*outNpol); 
      colArrayIn = (float *)malloc(sizeof(float)*nchan*npol);
      colArrayOut = (float *)malloc(sizeof(float)*nchan*outNpol);

      printf("Status before creating new file = %d\n",status);
      
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
		  
		  
		  fits_get_hdrspace(infptr,&keysexist,&morekeys,&status);
		  //		  printf("keys = %d %d %d\n",keysexist,morekeys,status);
		  fits_create_hdu(outfptr,&status);
		  for (i=1;i<=keysexist;i++)
		    {
		      fits_read_record(infptr,i,card,&status);
		      if (status) fits_report_error(stderr, status);

		      fits_get_keyname(card,key,&len,&status);
		      fits_parse_value(card,val,comment,&status);
		   		    		      
  		      if (strcmp(key,checkName1)==0)
			{
			  sprintf(newVal,"%dB",nchan*outNpol*nsblk/outSamplesPerByte); 

			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
			  if (*comment)
			    {
			      strcat(newcard, " / ");
			      strcat(newcard, comment);
			    }
			  fits_parse_template(newcard,card,&keytype,&status);
			  strcpy(newcard,card);
			}
		      else if (strcmp(key,checkName2)==0) // TDIM
			{
			  sprintf(newVal,"(%d,%d,%d)",nchan/outSamplesPerByte,outNpol,nsblk);
			  
			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
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
			  sprintf(newVal,"%dD",nchan*outNpol); // THIS SHOULD NOT BE NPOL -- A MISTAKE IN UWL PSRFITS?

			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
			  if (*comment)
			    {
			      strcat(newcard, " / ");
			      strcat(newcard, comment);
			    }
			  fits_parse_template(newcard,card,&keytype,&status);
			  strcpy(newcard,card);
			}
		      else if (strcmp(key,checkName4)==0) // DAT_WTS
			{
			  sprintf(newVal,"%dE",nchan*outNpol); // THIS SHOULD NOT BE NPOL -- A MISTAKE IN UWL PSRFITS?

			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
			  if (*comment)
			    {
			      strcat(newcard, " / ");
			      strcat(newcard, comment);
			    }
			  fits_parse_template(newcard,card,&keytype,&status);
			  strcpy(newcard,card);
			}
		      else if  (strcmp(key,checkName5)==0 || strcmp(key,checkName6)==0) // DAT_OFFS,DAT_SCL
			{
			  sprintf(newVal,"%dE",nchan*outNpol);

			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
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

			  // Note some mix of floats and doubles and bytes
			  // DAT_FREQ = double (nchan)
			  // DAT_WTS = double 

			  //                                                   
			  newBytes = oldBytes
			    - nchan*nsblk*npol/samplesperbyte // DATA
			    - nchan*4  * npol // DAT_WTS (FLOAT) (the npol is wrong, but seems to be needed in current PSRFITS from UWL)
			    - nchan*npol*2*4 // DAT_OFFS/DAT_SCLS (FLOAT)
			    - nchan*1*8 * npol  // DAT_FREQ (DOUBLE) (the npol is wrong, but seems to be needed in current PSRFITS from UWL)
			    + nchan*outNpol*nsblk/outSamplesPerByte
			    + nchan*4 * outNpol // (the npol is wrong, but seems to be needed in current PSRFITS from UWL)
			    + nchan*outNpol*2*4
			    + nchan*1*8 * outNpol; // (the npol is wrong, but seems to be needed in current PSRFITS from UWL) 

			  printf("oldBytes = %d (%d) [%d] %d %d\n",oldBytes,nchan,newBytes,sizeof(float),sizeof(double));
			  //			  printf("sizeof = %d %d\n",sizeof(float),sizeof(double));
			  sprintf(newVal,"%d",newBytes);
			  
			  strcpy(newcard,key);
			  strcat(newcard," = ");
			  strcat(newcard,newVal);
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
		      if (strcmp(key,checkName2)==0) // TDIM
			{
			  printf("Doing the change %d (%d)\n",status,nchan*npol*nsblk/outSamplesPerByte);
			  fits_modify_vector_len(outfptr,colnum_data_in,(long)((long)(nchan)*(long)outNpol*(long)nsblk/outSamplesPerByte),&status); 
			  printf("Status = %d\n",status);
			}
		    }
		      
		}
	      else
		{
		  fits_copy_hdu(infptr, outfptr, 0, &status);		  	      	      
		}
	      hdu++;
	    }
	}
      if (status == END_OF_FILE) status = 0;	  
    }

  fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
  printf("Status = %d\n",status);
  fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);

  printf("Starting %d\n",status);
  fits_get_colnum(outfptr, CASEINSEN, "DATA", &colnum_data_out, &status);  
  printf("a %d\n",status);
  fits_get_colnum(outfptr, CASEINSEN, "DAT_FREQ", &colnum_out_datFreq, &status);
  fits_get_colnum(outfptr, CASEINSEN, "DAT_WTS", &colnum_out_datWts, &status);
  fits_get_colnum(outfptr, CASEINSEN, "DAT_SCL", &colnum_out_datScl, &status);
  fits_get_colnum(outfptr, CASEINSEN, "DAT_OFFS", &colnum_out_datOffs, &status);
  if (status) {fits_report_error(stderr, status); exit(1);}

  for (i=0;i<nsubint;i++)
    {
      printf("Processing sub-integration %d of %d (status = %d)\n",i+1,nsubint,status);
      fits_read_col(infptr,TBYTE,colnum_data_in,i+1,1,nsblk*nchan*npol/samplesperbyte,&nullVal,dataIn,&initflag,&status);
      for (k=0;k<outNpol;k++)
	{
	  for (j=0;j<nsblk;j++)
	    {
	      writePos = (k*nchan+j*outNpol*nchan)/outSamplesPerByte; 
	      readPos = (j*npol*nchan+k*nchan)/samplesperbyte;

	      // Should not do this is outSamplesPerByte != samplesperbyte
	      if (outSamplesPerByte == samplesperbyte)
		{
		  memcpy(dataOut+writePos,dataIn+readPos,nchan/samplesperbyte);  
		  if (setFilterHigh==1 && samplesperbyte==1)
		    {
		      for (ii=0;ii<nchan;ii++)
			{
			  ival = (int)dataOut[writePos+ii];
			  
			  if (dataOut[writePos+ii] > filterHighBound)
			    dataOut[writePos+ii] =  filterHighBound;
			}
		    }
		}
	      else
		{
		  unsigned char tc;
		  long n = 0;
		  int jj;
		  double bit_level = 127.5; // Should set carefully
		  
		  for (ii=0;ii<nchan/outSamplesPerByte;ii++)
		    {
		      tc=0;
		      for (jj=0;jj<8;jj++)
			{
			  if (dataIn[n+readPos] > bit_level)
			    tc = tc | (1 << (7-jj));
			  n++;
			}
		      dataOut[writePos+ii] = tc;
		    }
		}		  
	    }
	}
      fits_write_col(outfptr,TBYTE,colnum_data_out,i+1,1,nchan*nsblk*outNpol/outSamplesPerByte,dataOut,&status);
      fits_read_col(infptr,TFLOAT,colnum_in_datFreq,i+1,1,nchan,&nullVal_f,colArrayIn,&initflag,&status);
      fits_write_col(outfptr,TFLOAT,colnum_out_datFreq,i+1,1,nchan,colArrayIn,&status);
      fits_read_col(infptr,TFLOAT,colnum_in_datWts,i+1,1,nchan,&nullVal_f,colArrayIn,&initflag,&status);
      fits_write_col(outfptr,TFLOAT,colnum_out_datWts,i+1,1,nchan,colArrayIn,&status);
      fits_read_col(infptr,TFLOAT,colnum_in_datScl,i+1,1,nchan*outNpol,&nullVal_f,colArrayIn,&initflag,&status);
      for (j=0;j<outNpol;j++)
	{
	  for (k=0;k<nchan;k++)
	    colArrayOut[j*nchan+k] = colArrayIn[j*nchan+k];	    
	}
      fits_write_col(outfptr,TFLOAT,colnum_out_datScl,i+1,1,nchan*outNpol,colArrayOut,&status);
      fits_read_col(infptr,TFLOAT,colnum_in_datOffs,i+1,1,nchan*outNpol,&nullVal_f,colArrayIn,&initflag,&status);
      for (j=0;j<outNpol;j++)
	{
	  for (k=0;k<nchan;k++)
	    colArrayOut[j*nchan+k] = colArrayIn[j*nchan+k];	    
	}
      fits_write_col(outfptr,TFLOAT,colnum_out_datOffs,i+1,1,nchan*outNpol,colArrayOut,&status);
    }
  fits_update_key(outfptr, TINT, (char *)"NPOL", &outNpol, NULL, &status );
  fits_update_key(outfptr, TINT, (char *)"NBITS", &outNbits, NULL, &status );    
  if (status) {fits_report_error(stderr, status); exit(1);}
  if (status) {fits_report_error(stderr, status); exit(1);}
  // Back to primary header
  fits_movabs_hdu(outfptr,1,NULL,&status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  
  free(dataIn);
  free(dataOut);
  free(colArrayIn);
  free(colArrayOut);
  

  fits_close_file(infptr, &status);
  fits_close_file(outfptr, &status);
  
  printf("Completed\n");
}
