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
// Newton:
// gcc -lm -o pfits_crossCoupling pfits_crossCoupling.c pfits_loader.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include "fftw3.h"
#include "TKnonlinearFit.h"

int checkFreq(double freq,double f0,double bw);
void set1MHZ_flagging(int *flag,int nchan);
void scaleInput(double calAA,double calBB,double calReAB,double calImAb,double *aa,double *bb,double *reAB,double *imAB);
double nonlinearFunc(double x,const double *par,int obsNum);
void sdhdf_mult4x4(double src1[4][4], double src2[4][4], double dest[4][4]);
int sdhdf_inv4x4(double m[4][4],double invOut[4][4]);
double sdhdf_dotproduct(double *v1,double *v2);
void displayMatrix_4x4(double matrix[4][4]);
void sdhdf_setIdentity_4x4(double mat[4][4]);
void sdhdf_display_vec4(double *vec);
void sdhdf_multMat_vec_replace(double mat[4][4],double *vec);
void sdhdf_multMat_vec(double mat[4][4],double *vecIn,double *vecOut);
void sdhdf_copy_mat4(double in[4][4],double out[4][4]);
void sdhdf_mult4x4_replace(double src1[4][4], double src2[4][4]);
void sdhdf_copy_vec4(double *in,double *out);
void sdhdf_setGainPhase(double ma[4][4],double diffGain,double diffPhase);
void sdhdf_setGain2Phase(double ma[4][4],double diffGain,double diffPhase);
void sdhdf_setImperfection(double mcc[4][4],double eps,double phi);
void sdhdf_setImperfection_full(double mcc[4][4],double eps1,double phi1,double eps2,double phi2);
void sdhdf_setParallacticAngle(double msky[4][4],double pa);
void sdhdf_setAstroAngle(double mastro[4][4],double pa);
void sdhdf_setFeed(double mf[4][4],double gamma);
void sdhdf_setFeedAlphaChi(double mf[4][4],double alpha,double chi);
void sdhdf_convertStokes(double p1,double p2,double p3,double p4,double *stokesI,double *stokesQ,double *stokesU,double *stokesV);
void sdhdf_convertStokesVec(double p1,double p2,double p3,double p4,double *vec);
double calcRMS(double *in,int n);
void readPol(dSetStruct *dSet,int polNum,float *loadData);

#define MAX_BANDS 1024

// Global parameters
double gain[1024],phase[1024];
int nAct_global;
int model_global = 3;
int nBinFit = 15;

double fortran_mod( double a, double p);

#define MAX_FILES 64

int main(int argc,char *argv[])
{
  int i,j,k,ii,jj,kk,bval;
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int colnum_data,colnum_datOffs,colnum_datScl,colnum_datWts,colnum_para;
  long nchan,nbin,nsub,npol;   
  float *fChan;
  float n_fval=0;
  short int n_val =0;
  int initflag =0;
  float val;
  int *flag;
  char fname[MAX_FILES][128];
  char fnameCal[MAX_FILES][128];
  float calPhaseOn1[MAX_FILES],calPhaseOn2[MAX_FILES];
  float calPhaseOff1[MAX_FILES],calPhaseOff2[MAX_FILES];
  int nFiles=0;
  FILE *fout,*fin,*fout2;
  char oname[1024];
  
  
  double *para;
  float *loadAll_pol1,*loadAll_pol2,*loadAll_pol3,*loadAll_pol4;
  int polNum;
  int debugChan=-1;
  int debugBand=-1;
  int debugFile=-1;
  int binSelect=0;
  int binOff=-1;
  int offset;
  double paraVal;
  
  double fitX[MAX_FILES*3*16],fitY[MAX_FILES*3*16]; // 16 should be MAX_BINS
  int nFit; // =5 ; // 8;  // Set nFit up here as need to allocate the memory
  int nPts;
  double *fitParams;
  double *Im,*Qm,*Um,*Vm;

  lm_status_struct fit_status;
  lm_control_struct fit_control = lm_control_double;

  double stokesM[4];
  double stokesPreCalib[4];
  double **stokesA;
  double stokesOut[4];
  
  double mueller[4][4],pa,msky[4][4],minv[4][4];
  double mhand[4][4],mastro[4][4];
  double ma[4][4],mf[4][4],mcc[4][4];
  double eps,phi,deltaG2,alpha,psi,chi;

  double eps1,eps2,phi1,phi2;
  double aa,bb,reAB,imAB,reAB_prime,imAB_prime;
  
  double bandF1[MAX_BANDS],bandF2[MAX_BANDS];
  int bandFlag[MAX_BANDS];
  int nBands=0;

  double Ioff,Qoff,Uoff,Voff;
  double Ion,Qon,Uon,Von;
  double f0;
  double outQI,outUI,outVI;

  double calAA_on,calBB_on,calReAB_on,calImAB_on;
  double calAA_off,calBB_off,calReAB_off,calImAB_off;  
  int nOn,nOff;
  int nAdd;
  char fileListing[1024];


  double *calDiffGain;
  double *calDiffPhase;
  double *calAA,*calBB,*calReAB,*calImAB;

  
  // Hard code a nice clean band

  nBands=0;
  for (f0 = 704;f0 < 4032;f0+=8)
  //  for (f0 = 1344;f0 < 1497;f0+=16)
   //  for (f0 = 1344;f0 < 2368;f0+=32)
    {
      bandF1[nBands] = f0;
      //      bandF2[nBands] = f0+16;
      bandF2[nBands] = f0+8;
      bandFlag[nBands] = 0;
      nBands++;
    }
  printf("Nbands = %d\n",nBands);
  //  bandF1[0] = 1344+16;
  //  bandF2[0] = 1344+16*2;
  //  bandF1[0] = 1360;
  //  bandF2[0] = 1460;
  //  nBands = 1;
  
  //  bandF1[0] = 1241;

  // TRY HIGHER FREQUENCY BAND

  /*
    bandF1[0] = 1400;
    bandF2[0] = 1440;
    nBands=1;
  */
  //  bandF1[0] = 704;
  //  bandF2[0] = 4032;
  //  nBands=1;

  
  /*
  bandF1[0] = 3160;
  bandF2[0] = 3240;

  //  bandF2[0] = 1497;
  nBands=1;
  */
  
  // Read inputs
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-debugChan")==0)
	sscanf(argv[++i],"%d",&debugChan);
      else if (strcmp(argv[i],"-binOff")==0)
	sscanf(argv[++i],"%d",&binOff);
      else if (strcmp(argv[i],"-debugBand")==0)
	sscanf(argv[++i],"%d",&debugBand);
      else if (strcmp(argv[i],"-debugFile")==0)
	sscanf(argv[++i],"%d",&debugFile);
      else if (strcmp(argv[i],"-bin")==0)
	sscanf(argv[++i],"%d",&binSelect);
      else if (strcmp(argv[i],"-model")==0)
	sscanf(argv[++i],"%d",&model_global);
      else if (strcmp(argv[i],"-f")==0)
	strcpy(fileListing,argv[++i]);
      //	strcpy(fname[nFiles++],argv[i]);
    }
  printf("Starting off\n");
  // Obtain the filenames for the pulsar and the cal
  fin = fopen(fileListing,"r");
  nFiles=0;
  while (!feof(fin))
    {
      if (fscanf(fin,"%s %s %f %f %f %f",fname[nFiles],fnameCal[nFiles],&calPhaseOn1[nFiles],&calPhaseOn2[nFiles],
		 &calPhaseOff1[nFiles],&calPhaseOff2[nFiles]) == 6)
	nFiles++;								   
    }
  fclose(fin);

  // Process the cals
  for (i=0;i<nFiles;i++)
    {
      initialise(&dSet,debug);

      setFilename(fnameCal[i],dSet,debug);
      pfitsOpenFile(dSet,debug);
      pfitsLoadHeader(dSet,debug);
      nchan = dSet->head->nchan;
      nbin  = dSet->head->nbin;
      nsub  = dSet->head->nsub;
      npol  = dSet->head->npol;
      if (nsub != 1)
	{
	  printf("ERROR: Current pfits_crossCoupling requires time scrunched data\n");
	  exit(1);
	}

      if (i==0) // Asssuming all files have the same number of channels
	{
	  flag = (int *)malloc(sizeof(int)*nchan);
	  set1MHZ_flagging(flag,nchan);
	  
	  calDiffGain  = (double *)malloc(sizeof(double)*nFiles*nchan);
	  calDiffPhase = (double *)malloc(sizeof(double)*nFiles*nchan);
	  calAA = (double *)malloc(sizeof(double)*nFiles*nchan);
	  calBB = (double *)malloc(sizeof(double)*nFiles*nchan);
	  calReAB = (double *)malloc(sizeof(double)*nFiles*nchan);
	  calImAB = (double *)malloc(sizeof(double)*nFiles*nchan);
	}
      printf("Loaded %s %d %d %d %d\n",fname[i],nchan,nbin,nsub,npol);

      // Obtain Im, Qm, Um, Vm and parallactic angle

      loadAll_pol1 = (float *)malloc(sizeof(float)*nchan*nbin);
      loadAll_pol2 = (float *)malloc(sizeof(float)*nchan*nbin);
      loadAll_pol3 = (float *)malloc(sizeof(float)*nchan*nbin);
      loadAll_pol4 = (float *)malloc(sizeof(float)*nchan*nbin);
      
      if (i==0) // Assuming frequency channels are the same for all files
	{
	  int flagit;
	  fChan    = (float *)malloc(sizeof(float)*nchan);
	  for (j=0;j<nchan;j++)
	    fChan[j] = dSet->head->chanFreq[j];

	  // Which bands will we flag
	  for (j=0;j<nBands;j++)
	    {
	      flagit=0;
	      bandFlag[j]=0;
	      for (k=0;k<nchan;k++)
		{
		  if (fChan[k] >= bandF1[j] && fChan[k] <= bandF2[j] && flag[k]==1)
		    {
		      flagit=1;
		      break;
		    }
		}
	      printf("BAND %d %g %g %d\n",j,bandF1[j],bandF2[j],flagit);
	      if (flagit==1) bandFlag[j]=1;
	    }


	}
    
      // Read the data
      fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_data,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);
      
      // Polarisation 1
      readPol(dSet,0,loadAll_pol1);
      readPol(dSet,1,loadAll_pol2);
      readPol(dSet,2,loadAll_pol3);
      readPol(dSet,3,loadAll_pol4);

      // Get CAL ON and CAL OFF and then DIFF GAIN and DIFF PHASE for each channel
      sprintf(oname,"%s.process",fnameCal[i]);
      fout = fopen(oname,"w");
      for (j=0;j<nchan;j++)
	{
	  calAA_on=calBB_on=calReAB_on=calImAB_on = 0.0;
	  calAA_off=calBB_off=calReAB_off=calImAB_off = 0.0;
	  nOn = nOff =0.0;
	  for (k=0;k<nbin;k++)
	    {
	      if (k > calPhaseOn1[i]*nbin && k < calPhaseOn2[i]*nbin)
		{
		  calAA_on   += loadAll_pol1[j*nbin+k];
		  calBB_on   += loadAll_pol2[j*nbin+k];
		  calReAB_on += loadAll_pol3[j*nbin+k];
		  calImAB_on += loadAll_pol4[j*nbin+k];
		  nOn++;
		}
	      if (k > calPhaseOff1[i]*nbin && k < calPhaseOff2[i]*nbin)
		{
		  calAA_off   += loadAll_pol1[j*nbin+k];
		  calBB_off   += loadAll_pol2[j*nbin+k];
		  calReAB_off += loadAll_pol3[j*nbin+k];
		  calImAB_off += loadAll_pol4[j*nbin+k];
		  nOff++;
		}
	    }
	  calAA[i*nchan+j]   = calAA_on/nOn   - calAA_off/nOff;
	  calBB[i*nchan+j]   = calBB_on/nOn   - calBB_off/nOff;
	  calReAB[i*nchan+j] = calReAB_on/nOn - calReAB_off/nOff;
	  calImAB[i*nchan+j] = calImAB_on/nOn - calImAB_off/nOff;
	
	  calDiffGain[i*nchan+j]  = calBB[i*nchan+j]/calAA[i*nchan+j];
	  calDiffPhase[i*nchan+j] = atan2(calImAB[i*nchan+j],calReAB[i*nchan+j]);
	  if (calDiffPhase[i*nchan+j] > 2)
	    calDiffPhase[i*nchan+j] = calDiffPhase[i*nchan+j]-2*M_PI;
	  fprintf(fout,"CAL: %d %g %d %g %g %g %g %g %g %g %g %g %g\n",j,fChan[j],flag[j],calDiffGain[i*nchan+j],calDiffPhase[i*nchan+j],calAA[i*nchan+j],calBB[i*nchan+j],calReAB[i*nchan+j],calImAB[i*nchan+j],calAA_on/nOn,calBB_on/nOn,calReAB_on/nOn,calImAB_on/nOn);
	  if (j==debugChan)
	    printf("CAL_ON_OFF %d %g %g %g %g %g\n",i,fChan[j],calAA_on/nOn,calBB_on/nOn,calAA_off/nOff,calBB_off/nOff); 
	}
      fclose(fout);
      free(loadAll_pol1);
      free(loadAll_pol2);
      free(loadAll_pol3);
      free(loadAll_pol4);      
  
      deallocateMemory(&dSet,debug);
      fits_close_file(dSet->fp,&status);     
      fits_report_error(stderr,status);      
    }
  
  if (model_global==1)
    nFit = 8;
  else if (model_global==2)
    nFit = 5;
  else if (model_global==3)
    nFit = 3+2*nFiles;  // eps, phi, alpha and then for each file deltaG and psi
  else if (model_global==4)
    nFit = 3+3+2*nFiles;  // qa,ua,va,eps, phi, alpha and then for each file deltaG and psi
  else if (model_global==5)
    nFit = 6;  
  else if (model_global==6)
    nFit = 11;  

  printf("Loading %d files\n",nFiles);

  for (i=0;i<nFiles;i++)
    {
      // Initialise everything
      initialise(&dSet,debug);

      setFilename(fname[i],dSet,debug);
      pfitsOpenFile(dSet,debug);
      pfitsLoadHeader(dSet,debug);
      nchan = dSet->head->nchan;
      nbin  = dSet->head->nbin;
      nsub  = dSet->head->nsub;
      npol  = dSet->head->npol;
      if (nsub != 1)
	{
	  printf("ERROR: Current pfits_crossCoupling requires time scrunched data\n");
	  printf("(and dedispersed data)\n");
	}
      if (i==0) // Assuming all files will have the same number of channels
	{
	  Im   = (double *)malloc(sizeof(double)*nchan*nFiles*nBinFit);
	  Qm   = (double *)malloc(sizeof(double)*nchan*nFiles*nBinFit);
	  Um   = (double *)malloc(sizeof(double)*nchan*nFiles*nBinFit);
	  Vm   = (double *)malloc(sizeof(double)*nchan*nFiles*nBinFit);
	  para = (double *)malloc(sizeof(double)*nchan*nFiles);
	  fitParams = (double *)malloc(sizeof(double)*nchan*nFit);
	}

      printf("Loaded %s %d %d %d %d\n",fname[i],nchan,nbin,nsub,npol);

      // Obtain Im, Qm, Um, Vm and parallactic angle

      loadAll_pol1 = (float *)malloc(sizeof(float)*nchan*nbin);
      loadAll_pol2 = (float *)malloc(sizeof(float)*nchan*nbin);
      loadAll_pol3 = (float *)malloc(sizeof(float)*nchan*nbin);
      loadAll_pol4 = (float *)malloc(sizeof(float)*nchan*nbin);
      
      if (i==0) // Assuming frequency channels are the same for all files
	{
	  fChan    = (float *)malloc(sizeof(float)*nchan);
	  for (j=0;j<nchan;j++)
	    fChan[j] = dSet->head->chanFreq[j];
	  if (binOff==-1) // Note this should be set for correct usage
	    binOff = (int)(nbin/2.);
	}
    
      // Read the data
      fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_data,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);
      fits_get_colnum(dSet->fp,CASEINSEN,"PAR_ANG",&colnum_para,&status);
      
      // Read the parallactic angle
      fits_read_col(dSet->fp,TFLOAT,colnum_para,1,1,1,&n_fval,&val,&initflag,&status);
      para[i] = val;
          
      // Polarisation 1
      readPol(dSet,0,loadAll_pol1);
      readPol(dSet,1,loadAll_pol2);
      readPol(dSet,2,loadAll_pol3);
      readPol(dSet,3,loadAll_pol4);


      if (debugChan >= 0)
	{
	  fout = fopen("crossCoupling_profile.dat","w");
	  for (j=0;j<nbin;j++)
	    fprintf(fout,"PROFILE: %d %g %g %g %g\n",j,loadAll_pol1[debugChan*nbin+j],loadAll_pol2[debugChan*nbin+j],loadAll_pol3[debugChan*nbin+j],loadAll_pol4[debugChan*nbin+j]);
	  fclose(fout);
	}
    
      // Now set Im,Qm,Um,Vm
      //      k = debugBand;

      sprintf(oname,"crossCoupling_offPulse.%s.dat",fname[i]);
      fout = fopen(oname,"w");
      for (j=0;j<nchan;j++)
	{
	  aa   = loadAll_pol1[j*nbin + binOff];
	  bb   = loadAll_pol2[j*nbin + binOff];
	  reAB = loadAll_pol3[j*nbin + binOff];
	  imAB = loadAll_pol4[j*nbin + binOff];		  
	  fprintf(fout,"%g %g %g %g %g ",fChan[j],aa,bb,reAB,imAB);
	  scaleInput(calAA[i*nchan+j],calBB[i*nchan+j],calReAB[i*nchan+j],calImAB[i*nchan+j],&aa,&bb,&reAB,&imAB);
	  fprintf(fout,"%g %g %g %g\n",aa,bb,reAB,imAB);
	}
      fclose(fout);

      sprintf(oname,"crossCoupling_measuredStokes.%s.dat",fname[i]);
      fout = fopen(oname,"w");
      for (bval=0;bval<nBinFit;bval++)
	{
	  for (k=0;k<nBands;k++)
	    {
	      Im[bval*nBands*nFiles + i*nBands+k] = Qm[bval*nBands*nFiles + i*nBands+k] = Um[bval*nBands*nFiles + i*nBands+k] = Vm[bval*nBands*nFiles + i*nBands+k] = 0.0;
	      nAdd=0;
	      for (j=0;j<nchan;j++)
		{
		  //	      printf("FCHAN = %g %g %g\n",fChan[j],bandF1[k],bandF2[k]);
		  if (fChan[j] >= bandF1[k] && fChan[j] <= bandF2[k] && flag[j] != 1)
		    {
		      aa   = loadAll_pol1[j*nbin + binOff];
		      bb   = loadAll_pol2[j*nbin + binOff];
		      reAB = loadAll_pol3[j*nbin + binOff];
		      imAB = loadAll_pol4[j*nbin + binOff];		  
		      
		      // Do a complex division to remove the cal
		      scaleInput(calAA[i*nchan+j],calBB[i*nchan+j],calReAB[i*nchan+j],calImAB[i*nchan+j],&aa,&bb,&reAB,&imAB);
		      
		      Ioff =  aa+bb;
		      Qoff =  aa-bb;
		      Uoff =  2*reAB;
		      Voff =  2*imAB;

		      binSelect = bval; // ** CHECK THIS TO CHOOSE CORRECT BINS
		      
		      aa   = loadAll_pol1[j*nbin + binSelect];
		      bb   = loadAll_pol2[j*nbin + binSelect];
		      reAB = loadAll_pol3[j*nbin + binSelect];
		      imAB = loadAll_pol4[j*nbin + binSelect];
		      
		      // Do a complex division to remove the cal
		      scaleInput(calAA[i*nchan+j],calBB[i*nchan+j],calReAB[i*nchan+j],calImAB[i*nchan+j],&aa,&bb,&reAB,&imAB);		  
		      
		      Ion =  aa+bb;
		      Qon =  aa-bb;
		      Uon =  2*reAB;
		      Von =  2*imAB;
		      
		      Im[bval*nBands*nFiles + i*nBands+k] += (Ion - Ioff);
		      Qm[bval*nBands*nFiles + i*nBands+k] += (Qon - Qoff);
		      Um[bval*nBands*nFiles + i*nBands+k] += (Uon - Uoff);
		      Vm[bval*nBands*nFiles + i*nBands+k] += (Von - Voff);
		      nAdd++;
		    }
		}
	      Im[bval*nBands*nFiles + i*nBands+k] /= (double)nAdd;
	      Qm[bval*nBands*nFiles + i*nBands+k] /= (double)nAdd;
	      Um[bval*nBands*nFiles + i*nBands+k] /= (double)nAdd;
	      Vm[bval*nBands*nFiles + i*nBands+k] /= (double)nAdd;
	      
	      fprintf(fout,"STOKES_MEASURED: %d %d %g %g %g %g %g\n",k,bval,para[i],Im[bval*nBands*nFiles + i*nBands+k],Qm[bval*nBands*nFiles + i*nBands+k],Um[bval*nBands*nFiles + i*nBands+k],Vm[bval*nBands*nFiles + i*nBands+k]);
	    }
	}
      fclose(fout);
      free(loadAll_pol1);
      free(loadAll_pol2);
      free(loadAll_pol3);
      free(loadAll_pol4);      
  
      deallocateMemory(&dSet,debug);
      fits_close_file(dSet->fp,&status);     
      fits_report_error(stderr,status);
    }

  if (debugBand >= 0)
    {
      sprintf(oname,"crossCoupling_inputStokes.dat");
      fout = fopen(oname,"w");
      for (bval=0;bval<nBinFit;bval++)
	{
	  for (i=0;i<nFiles;i++)
	    fprintf(fout,"INPUT_STOKES: %d %g %g %g %g %g\n",i,para[i],Im[bval*nBands*nFiles + i*nBands+debugBand],Qm[bval*nBands*nFiles + i*nBands+debugBand],Um[bval*nBands*nFiles + i*nBands+debugBand],Vm[bval*nBands*nFiles + i*nBands+debugBand]);
	}
      fclose(fout);
    }

  
  // Now determine the cross coupling etc.
  //  j = debugChan;
  for (j=0;j<nBands;j++)
  {
    printf("BAND: %d %g %g\n",j,bandF1[j],bandF2[j]);


    // SHOULD ADD IN MORE BINS HERE
    for (bval=0;bval<nBinFit;bval++)
      {
	for (i=0;i<nFiles;i++)
	  {
	    // Set up the fitting values
	    fitX[bval*nFiles*3 + i]          = para[i];
	    fitY[bval*nFiles*3 + i]          = Qm[bval*nFiles*nBands + i*nBands+j]/Im[bval*nFiles*nBands + i*nBands+j];
	    fitX[bval*nFiles*3 + i+nFiles]   = para[i];
	    fitY[bval*nFiles*3 + i+nFiles]   = Um[bval*nFiles*nBands + i*nBands+j]/Im[bval*nFiles*nBands + i*nBands+j];
	    fitX[bval*nFiles*3 + i+2*nFiles] = para[i];
	    fitY[bval*nFiles*3 + i+2*nFiles] = Vm[bval*nFiles*nBands + i*nBands+j]/Im[bval*nFiles*nBands + i*nBands+j];	
	  }
      }
    
    nPts=nFiles*3*nBinFit;
    nAct_global = nFiles;
    
    // fitParams[0] = 0.1; // I_actual

    if (model_global == 1)
      {
	/*
	fitParams[0+j*nFit] = 0.1;  // Q_actual
	fitParams[1+j*nFit] = 0.1;  // U_actual
	fitParams[2+j*nFit] = 0.1;  // V_actual
	fitParams[3+j*nFit] = 0.01; // eps
	fitParams[4+j*nFit] = 0.0;  // phi
	fitParams[5+j*nFit] = 0.0;  // alpha
	fitParams[6+j*nFit] = 0.1;  // Delta G2
	fitParams[7+j*nFit] = 0.2;  // diff phase
	*/
	fitParams[0+j*nFit] = 0.0;  // Q_actual
	fitParams[1+j*nFit] = 0.0;  // U_actual
	fitParams[2+j*nFit] = 0.0;  // V_actual
	fitParams[3+j*nFit] = 0.0;  // eps
	fitParams[4+j*nFit] = 0.0; //-M_PI;  // M_PI/2.;  // phi
	//    fitParams[5+j*nFit] = -90.0*M_PI/180; // alpha
    	fitParams[5+j*nFit] = 0.0; // M_PI/4.;  // alpha
	fitParams[6+j*nFit] = 0.0;  // Delta G2
	//	fitParams[7+j*nFit] = 0.0;  // diff phase
	fitParams[7+j*nFit] = 0.0; //M_PI;  // diff phase
      }
    else if (model_global == 2)
      {
	// Result is dependent on the starting values
	fitParams[0+j*nFit] = 0.0; // eps
	fitParams[1+j*nFit] = 0.0;  // phi
	fitParams[2+j*nFit] = 0.0;  // alpha
	fitParams[3+j*nFit] = 0.0;  // Delta G2
	fitParams[4+j*nFit] = 0.0;  // diff phase
      }
    else if (model_global == 3)
      {
	// Result is dependent on the starting values
	fitParams[0+j*nFit] = 0.0; // eps
	fitParams[1+j*nFit] = 0.0;  // phi
	fitParams[2+j*nFit] = 0.0;  // alpha
      
	for (k=0;k<nFiles;k++)
	  {
	    fitParams[3+2*k+j*nFit] = 0.0;  // Delta G2
	    fitParams[4+2*k+j*nFit] = 0.0;  // diff phase
	  }
      }
    else if (model_global == 4)
      {
	// Result is dependent on the starting values
	fitParams[0+j*nFit] = 0.0; // qa
	fitParams[1+j*nFit] = 0.0; // ua
	fitParams[2+j*nFit] = 0.0; // va
	fitParams[3+j*nFit] = 0.0; // eps
	fitParams[4+j*nFit] = 0.0;  // phi
	fitParams[5+j*nFit] = 0.0; // 45.0*M_PI/180.0; //0.0;  // alpha
      
	for (k=0;k<nFiles;k++)
	  {
	    fitParams[6+2*k+j*nFit] = 0.0;  // Delta G2
	    fitParams[7+2*k+j*nFit] = 0.0;  // diff phase
	  }
      }
    else if (model_global == 5)
      {
	fitParams[0+j*nFit] = 0.0;  // Q_actual
	fitParams[1+j*nFit] = 0.0;  // U_actual
	fitParams[2+j*nFit] = 0.0;  // V_actual
	fitParams[3+j*nFit] = 0.0;  // eps
	fitParams[4+j*nFit] = 0.0;  // M_PI/2.;  // phi
	fitParams[5+j*nFit] = 0.0;  // alpha
      }
    else if (model_global == 6)
      {
	fitParams[0+j*nFit] = 0.0;  // Q_actual
	fitParams[1+j*nFit] = 0.0;  // U_actual
	fitParams[2+j*nFit] = 0.0;  // V_actual
	fitParams[3+j*nFit] = 0.0;  // eps1
	fitParams[4+j*nFit] = 0.0;  // phi1
	fitParams[5+j*nFit] = 0.0;  // eps2
	fitParams[6+j*nFit] = 0.0;  // phi2
    	fitParams[7+j*nFit] = 0.0; // 45*M_PI/180; // M_PI/4.;  // alpha
	fitParams[8+j*nFit] = 0.0; // M_PI/4.;  // chi
	fitParams[9+j*nFit] = 0.0;  // Delta G2
	fitParams[10+j*nFit] = 0.0; //M_PI;  // diff phase
      }
    
    lmcurve_fit(nFit,fitParams+j*nFit,nPts,fitX,fitY,nonlinearFunc,&fit_control,&fit_status);
    // Get angles +/- pi
    //    printf("I AM IN HERE %g\n",fitParams[4+j*nFit]);


    /*
    while (fitParams[4+j*nFit] < -M_PI || fitParams[4+j*nFit] > M_PI)
      {
	//	printf("I AM IN HERE\n");
	if (fitParams[4+j*nFit] < -M_PI) fitParams[4+j*nFit] += 2*M_PI;
	if (fitParams[4+j*nFit] >  M_PI) fitParams[4+j*nFit] -= 2*M_PI;
      }
    */
    
    if (model_global==1)
      {
	if (bandFlag[j]==0)
	  printf("Output_FIT: %d %g %g %g %g %g %g %g %g %g\n",j,0.5*(bandF1[j]+bandF2[j]),fitParams[0+j*nFit],fitParams[1+j*nFit],
		 fitParams[2+j*nFit],fitParams[3+j*nFit],fitParams[4+j*nFit],fitParams[5+j*nFit],fitParams[6+j*nFit],fitParams[7+j*nFit]);
      }
    else if (model_global==2)
      {
	printf("Output_FIT: %d %g %g %g %g %g %g\n",j,0.5*(bandF1[j]+bandF2[j]),fitParams[0+j*nFit],fitParams[1+j*nFit],
	       fitParams[2+j*nFit],fitParams[3+j*nFit],fitParams[4+j*nFit]);
      }
    else if (model_global==3)
      {
	printf("Output_FIT: %d %g %g %g %g ",j,0.5*(bandF1[j]+bandF2[j]),fitParams[0+j*nFit],fitParams[1+j*nFit],
	       fitParams[2+j*nFit]);
	for (k=0;k<nFiles;k++)
	  printf("%g %g ",fitParams[3+2*k+j*nFit],fitParams[4+2*k+j*nFit]);
	printf("\n");
      }
    else if (model_global==4)
      {
	printf("Output_FIT: %d %g %g %g %g %g %g %g ",j,0.5*(bandF1[j]+bandF2[j]),
	       fitParams[0+j*nFit],fitParams[1+j*nFit],fitParams[2+j*nFit],
	       fitParams[3+j*nFit],fitParams[4+j*nFit],fitParams[5+j*nFit]);
	for (k=0;k<nFiles;k++)
	  printf("%g %g ",fitParams[6+2*k+j*nFit],fitParams[7+2*k+j*nFit]);
	printf("\n");
      }
    else if (model_global==5)
      {
	printf("Output_FIT: %d %g %g %g %g %g %g %g\n",j,0.5*(bandF1[j]+bandF2[j]),fitParams[0+j*nFit],fitParams[1+j*nFit],
	       fitParams[2+j*nFit],fitParams[3+j*nFit],fitParams[4+j*nFit],fitParams[5+j*nFit]);
      }
    else if (model_global==6)
      {
	printf("Output_FIT: %d %g %g %g %g %g %g %g %g %g %g %g %g\n",j,0.5*(bandF1[j]+bandF2[j]),fitParams[0+j*nFit],fitParams[1+j*nFit],
	       fitParams[2+j*nFit],fitParams[3+j*nFit],fitParams[4+j*nFit],fitParams[5+j*nFit],fitParams[6+j*nFit],fitParams[7+j*nFit],fitParams[8+j*nFit],fitParams[9+j*nFit],fitParams[10+j*nFit]);
      }

    // Now output Q/I, U/I and V/I model
    if (debugBand >= 0)
      {
	fout = fopen("crossCoupling_model.dat","w");
	for (k=0;k<nFiles;k++)
	  {
	    outQI = nonlinearFunc(para[k],fitParams+debugBand*nFit,k);
	    outUI = nonlinearFunc(para[k],fitParams+debugBand*nFit,k+nFiles);
	    outVI = nonlinearFunc(para[k],fitParams+debugBand*nFit,k+2*nFiles);
	    fprintf(fout,"%d %g %g %g %g\n",k,para[k],outQI,outUI,outVI);
	  }
	fclose(fout);
      }


  }


  // Now re-process the first pulsar observation
  if (debugBand >= 0)
    {
      //      i=0;
      for (i=0;i<nFiles;i++)
      {
	// Initialise everything
	initialise(&dSet,debug);
	
	setFilename(fname[i],dSet,debug);
	pfitsOpenFile(dSet,debug);
	pfitsLoadHeader(dSet,debug);
		
	loadAll_pol1 = (float *)malloc(sizeof(float)*nchan*nbin);
	loadAll_pol2 = (float *)malloc(sizeof(float)*nchan*nbin);
	loadAll_pol3 = (float *)malloc(sizeof(float)*nchan*nbin);
	loadAll_pol4 = (float *)malloc(sizeof(float)*nchan*nbin);
	
	// Read the data
	fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
		
	readPol(dSet,0,loadAll_pol1);
	readPol(dSet,1,loadAll_pol2);
	readPol(dSet,2,loadAll_pol3);
	readPol(dSet,3,loadAll_pol4);
	
	// Now set Im,Qm,Um,Vm
	if (model_global==1)
	  {
	    eps     = fitParams[3+debugBand*nFit];
	    phi     = fitParams[4+debugBand*nFit];
	    alpha   = fitParams[5+debugBand*nFit];
	    deltaG2 = fitParams[6+debugBand*nFit];
	    psi     = fitParams[7+debugBand*nFit];
	  }
	else if (model_global==2)
	  {
	    eps     = fitParams[0+debugBand*nFit];
	    phi     = fitParams[1+debugBand*nFit];
	    alpha   = fitParams[2+debugBand*nFit];
	    deltaG2 = fitParams[3+debugBand*nFit];
	    psi     = fitParams[4+debugBand*nFit];
	  }
	else if (model_global==3)
	  {
	    eps     = fitParams[0+debugBand*nFit];
	    phi     = fitParams[1+debugBand*nFit];
	    alpha   = fitParams[2+debugBand*nFit];
	    deltaG2 = fitParams[3+2*i+debugBand*nFit];
	    psi     = fitParams[4+2*i+debugBand*nFit];
	  }
	else if (model_global==4)
	  {
	    eps     = fitParams[3+debugBand*nFit];
	    phi     = fitParams[4+debugBand*nFit];
	    alpha   = fitParams[5+debugBand*nFit];
	    deltaG2 = fitParams[6+2*i+debugBand*nFit];
	    psi     = fitParams[7+2*i+debugBand*nFit];
	  }
	else if (model_global==5)
	  {
	    eps     = fitParams[3+debugBand*nFit];
	    phi     = fitParams[4+debugBand*nFit];
	    alpha   = fitParams[5+debugBand*nFit];
	    deltaG2 = 0;
	    psi     = 0;
	  }
	else if (model_global==6)
	  {
	    eps1     = fitParams[3+debugBand*nFit];
	    phi1     = fitParams[4+debugBand*nFit];
	    eps2     = fitParams[5+debugBand*nFit];
	    phi2     = fitParams[6+debugBand*nFit];
	    alpha   = fitParams[7+debugBand*nFit];
	    chi   = fitParams[8+debugBand*nFit];
	    deltaG2 = fitParams[9+debugBand*nFit];
	    psi     = fitParams[10+debugBand*nFit];
	  }
      
	stokesA = (double **)malloc(sizeof(double *)*nbin*2);
	for (j=0;j<nbin*2;j++)
	  {
	    stokesA[j] = (double *)malloc(sizeof(double)*4);
	    stokesA[j][0] = stokesA[j][1] =stokesA[j][2] =stokesA[j][3] = 0;
	  }

	
	sdhdf_setIdentity_4x4(mueller);
	sdhdf_setIdentity_4x4(mhand);
	sdhdf_setIdentity_4x4(mcc);
	sdhdf_setGain2Phase(ma,deltaG2,psi);
	if (model_global == 6)
	  {
	    sdhdf_setImperfection_full(mcc,eps1,phi1,eps2,phi2);
	    sdhdf_setFeedAlphaChi(mf,alpha,chi);
	  }
	else
	  {
	    sdhdf_setImperfection(mcc,eps,phi);
	    sdhdf_setFeed(mf,alpha);
	  }
	pa      = para[i]*M_PI/180.0;        
	sdhdf_setParallacticAngle(msky,pa);
	sdhdf_setAstroAngle(mastro,M_PI/4.); // NOTE HARDCODE
	//			sdhdf_setAstroAngle(mastro,M_PI/2.); // NOTE HARDCODE			
	
	sdhdf_mult4x4_replace(mueller,ma);
	sdhdf_mult4x4_replace(mueller,mcc);
	sdhdf_mult4x4_replace(mueller,mf);		
	sdhdf_mult4x4_replace(mueller,msky);
	sdhdf_mult4x4_replace(mueller,mastro); 
	if (i==0)
	  {
	    printf("m_lna matrix:\n");
	    displayMatrix_4x4(ma);
	    printf("---------------------------------\n");
	    printf("m_coupling matrix:\n");
	    displayMatrix_4x4(mcc);
	    printf("---------------------------------\n");
	    printf("m_feed matrix:\n");
	    displayMatrix_4x4(mf);
	    printf("---------------------------------\n");
	    printf("m_sky matrix:\n");
	    displayMatrix_4x4(msky);
	    printf("---------------------------------\n");
	    printf("m_astro matrix:\n");
	    displayMatrix_4x4(mastro);
	    printf("---------------------------------\n");
	    printf("Mueller matrix:\n");
	    displayMatrix_4x4(mueller);
	    printf("---------------------------------\n");
	  }
	sdhdf_inv4x4(mueller,minv);
      
	for (offset=0;offset<2;offset++)
	  {
	    for (j=0;j<nbin;j++)
	      {
		stokesM[0] = stokesM[1] = stokesM[2] = stokesM[3] = 0.0;
		nAdd=0;
		for (jj=0;jj<nchan;jj++)
		  {
		    if (fChan[jj] >= bandF1[debugBand] && fChan[jj] <= bandF2[debugBand] && flag[jj] == 0)
		      {
			aa   = loadAll_pol1[jj*nbin + j];
			bb   = loadAll_pol2[jj*nbin + j];
			reAB = loadAll_pol3[jj*nbin + j];
			imAB = loadAll_pol4[jj*nbin + j];

			if (i==debugFile && j==binSelect)
			  {
			    printf("astroPreproc: %d %d %g %g %g %g ",j,jj,aa,bb,reAB,imAB);
			  }

			// Do a complex division to remove the cal
			scaleInput(calAA[i*nchan+jj],calBB[i*nchan+jj],calReAB[i*nchan+jj],calImAB[i*nchan+jj],&aa,&bb,&reAB,&imAB);
			if (i==debugFile && j==binSelect)
			  {
			    printf(" %g %g %g %g\n",aa,bb,reAB,imAB);
			  }

			Ion =  aa+bb;
			Qon =  aa-bb;
			Uon =  2*reAB;
			Von =  2*imAB;
			

			stokesM[0] = Ion;
			stokesM[1] = Qon;
			stokesM[2] = Uon;
			stokesM[3] = Von;
			
			//			sdhdf_multMat_vec(minv,stokesM,stokesA[j+offset*nbin]);
			sdhdf_multMat_vec(minv,stokesM,stokesOut);
			stokesA[j+offset*nbin][0] += stokesOut[0];
			stokesA[j+offset*nbin][1] += stokesOut[1];
			stokesA[j+offset*nbin][2] += stokesOut[2];
			stokesA[j+offset*nbin][3] += stokesOut[3];
			
			/*
			  stokesM[0] += Ion;
			  stokesM[1] += Qon;
			  stokesM[2] += Uon;
			  stokesM[3] += Von;
			*/
			nAdd++;
		      }
		  }
		stokesA[j+offset*nbin][0]/=(double)nAdd;
		stokesA[j+offset*nbin][1]/=(double)nAdd;
		stokesA[j+offset*nbin][2]/=(double)nAdd;
		stokesA[j+offset*nbin][3]/=(double)nAdd;

		//		stokesM[0] /= (double)nAdd;
		//		stokesM[1] /= (double)nAdd;
		//		stokesM[2] /= (double)nAdd;
		//		stokesM[3] /= (double)nAdd;
		

		// SHOULD AVERAGE AFTER DETERMINING ACTUAL STOKES
		
		// Model FD_HAND = -1
		//		mhand[1][1]=-1; // Q
		//		mhand[3][3]=-1; // V
		
		//		sdhdf_mult4x4_replace(mueller,mhand); // Should this be before or after msky? or before after other parts
		
	      }
	  }
	sprintf(oname,"crossCoupling_calibratedProfile.%s.dat",fname[i]);
	fout = fopen(oname,"w");

	for (j=0;j<nbin*2;j++)
	  {
	    /*
	    fprintf(fout,"CALIBRATED_PROFILE: %d %d %g %g %g %g %g %g %g %g %g\n",j,debugBand,(bandF1[debugBand]+bandF2[debugBand])/2.,
		   stokesA[j][0]-stokesA[binOff][0],stokesA[j][1]-stokesA[binOff][1],stokesA[j][2]-stokesA[binOff][2],
		   stokesA[j][3]-stokesA[binOff][3],stokesM[0],stokesM[1],stokesM[2],stokesM[3]);
	    */
	    fprintf(fout,"CALIBRATED_PROFILE: %d %d %g %g %g %g %g %g %g %g %g\n",j,debugBand,(bandF1[debugBand]+bandF2[debugBand])/2.,
		    (stokesA[j][0]-stokesA[binOff][0])/(stokesA[0][0]-stokesA[binOff][0]),(stokesA[j][1]-stokesA[binOff][1])/(stokesA[0][0]-stokesA[binOff][0]),(stokesA[j][2]-stokesA[binOff][2])/(stokesA[0][0]-stokesA[binOff][0]),
		    (stokesA[j][3]-stokesA[binOff][3])/(stokesA[0][0]-stokesA[binOff][0]),stokesM[0],stokesM[1],stokesM[2],stokesM[3]);
	  }
	fclose(fout);
	for (j=0;j<nbin*2;j++)
	  free(stokesA[j]);
	free(stokesA);

      }



    }  
  free(flag);
  free(calDiffGain);
  free(calDiffPhase);
  free(calAA); free(calBB); free(calReAB); free(calImAB);
  free(fChan);
  free(fitParams);
  free(Im);  free(Qm); free(Um); free(Vm); free(para);
}


double calcRMS(double *in,int n)
{
  int i;
  double x=0,x2=0;
  double sdev;
  for (i=0;i<n;i++)
    {
      x+=in[i];
      x2+=pow(in[i],2);
    }

  sdev = sqrt(1/(double)n*x2 - pow(1.0/(double)n*x2,2));
  return sdev;
}

void sdhdf_copy_mat4(double in[4][4],double out[4][4])
{
  int i,j;
  for (i=0;i<4;i++)
    {
      for (j=0;j<4;j++)
	out[i][j] = in[i][j];
    }
}


void sdhdf_copy_vec4(double *in,double *out)
{
  int i;
  for (i=0;i<4;i++)
    out[i]=in[i];
}

void sdhdf_display_vec4(double *vec)
{
  printf("%g\n",vec[0]);
  printf("%g\n",vec[1]);
  printf("%g\n",vec[2]);
  printf("%g\n",vec[3]);
}

void sdhdf_multMat_vec_replace(double mat[4][4],double *vec)
{
  double newVec[4];
  int i;
  
  newVec[0] = mat[0][0]*vec[0] + mat[1][0]*vec[1] + mat[2][0]*vec[2] + mat[3][0]*vec[3];
  newVec[1] = mat[0][1]*vec[0] + mat[1][1]*vec[1] + mat[2][1]*vec[2] + mat[3][1]*vec[3];
  newVec[2] = mat[0][2]*vec[0] + mat[1][2]*vec[1] + mat[2][2]*vec[2] + mat[3][2]*vec[3];
  newVec[3] = mat[0][3]*vec[0] + mat[1][3]*vec[1] + mat[2][3]*vec[2] + mat[3][3]*vec[3];
  for (i=0;i<4;i++)
    vec[i] = newVec[i];
}

void sdhdf_multMat_vec(double mat[4][4],double *vec,double *vecOut)
{
  vecOut[0] = mat[0][0]*vec[0] + mat[1][0]*vec[1] + mat[2][0]*vec[2] + mat[3][0]*vec[3];
  vecOut[1] = mat[0][1]*vec[0] + mat[1][1]*vec[1] + mat[2][1]*vec[2] + mat[3][1]*vec[3];
  vecOut[2] = mat[0][2]*vec[0] + mat[1][2]*vec[1] + mat[2][2]*vec[2] + mat[3][2]*vec[3];
  vecOut[3] = mat[0][3]*vec[0] + mat[1][3]*vec[1] + mat[2][3]*vec[2] + mat[3][3]*vec[3];
}

void sdhdf_setIdentity_4x4(double mat[4][4])
{
  int i,j;
  for (i=0;i<4;i++)
    {
      for (j=0;j<4;j++)
	mat[i][j]=0;
    }
  mat[0][0]=mat[1][1]=mat[2][2]=mat[3][3]=1;    
}

void sdhdf_mult4x4_replace(double src1[4][4], double src2[4][4])
{
  double temp[4][4];
  int i,j;
  sdhdf_mult4x4(src1,src2,temp);
  for (i=0;i<4;i++)
    {
      for (j=0;j<4;j++)
	src1[i][j] = temp[i][j];
    }
}

// Multiply 4x4 matrix
void sdhdf_mult4x4(double src1[4][4], double src2[4][4], double dest[4][4])
{

  dest[0][0] = src1[0][0]*src2[0][0] + src1[1][0]*src2[0][1] + src1[2][0]*src2[0][2] + src1[3][0]*src2[0][3];
  dest[0][1] = src1[0][1]*src2[0][0] + src1[1][1]*src2[0][1] + src1[2][1]*src2[0][2] + src1[3][1]*src2[0][3];
  dest[0][2] = src1[0][2]*src2[0][0] + src1[1][2]*src2[0][1] + src1[2][2]*src2[0][2] + src1[3][2]*src2[0][3];
  dest[0][3] = src1[0][3]*src2[0][0] + src1[1][3]*src2[0][1] + src1[2][3]*src2[0][2] + src1[3][3]*src2[0][3];



  dest[1][0] = src1[0][0]*src2[1][0] + src1[1][0]*src2[1][1] + src1[2][0]*src2[1][2] + src1[3][0]*src2[1][3];
  dest[1][1] = src1[0][1]*src2[1][0] + src1[1][1]*src2[1][1] + src1[2][1]*src2[1][2] + src1[3][1]*src2[1][3];
  dest[1][2] = src1[0][2]*src2[1][0] + src1[1][2]*src2[1][1] + src1[2][2]*src2[1][2] + src1[3][2]*src2[1][3];
  dest[1][3] = src1[0][3]*src2[1][0] + src1[1][3]*src2[1][1] + src1[2][3]*src2[1][2] + src1[3][3]*src2[1][3];

  dest[2][0] = src1[0][0]*src2[2][0] + src1[1][0]*src2[2][1] + src1[2][0]*src2[2][2] + src1[3][0]*src2[2][3];
  dest[2][1] = src1[0][1]*src2[2][0] + src1[1][1]*src2[2][1] + src1[2][1]*src2[2][2] + src1[3][1]*src2[2][3];
  dest[2][2] = src1[0][2]*src2[2][0] + src1[1][2]*src2[2][1] + src1[2][2]*src2[2][2] + src1[3][2]*src2[2][3];
  dest[2][3] = src1[0][3]*src2[2][0] + src1[1][3]*src2[2][1] + src1[2][3]*src2[2][2] + src1[3][3]*src2[2][3];

  dest[3][0] = src1[0][0]*src2[3][0] + src1[1][0]*src2[3][1] + src1[2][0]*src2[3][2] + src1[3][0]*src2[3][3];
  dest[3][1] = src1[0][1]*src2[3][0] + src1[1][1]*src2[3][1] + src1[2][1]*src2[3][2] + src1[3][1]*src2[3][3];
  dest[3][2] = src1[0][2]*src2[3][0] + src1[1][2]*src2[3][1] + src1[2][2]*src2[3][2] + src1[3][2]*src2[3][3];
  dest[3][3] = src1[0][3]*src2[3][0] + src1[1][3]*src2[3][1] + src1[2][3]*src2[3][2] + src1[3][3]*src2[3][3];

  
  



  
  /*
  dest[0][0] = src1[0][0] * src2[0][0] + src1[0][1] * src2[1][0] + src1[0][2] * src2[2][0] + src1[0][3] * src2[3][0];
  dest[0][1] = src1[0][0] * src2[0][1] + src1[0][1] * src2[1][1] + src1[0][2] * src2[2][1] + src1[0][3] * src2[3][1];
  dest[0][2] = src1[0][0] * src2[0][2] + src1[0][1] * src2[1][2] + src1[0][2] * src2[2][2] + src1[0][3] * src2[3][2];
  dest[0][3] = src1[0][0] * src2[0][3] + src1[0][1] * src2[1][3] + src1[0][2] * src2[2][3] + src1[0][3] * src2[3][3];

  dest[1][0] = src1[1][0] * src2[0][0] + src1[1][1] * src2[1][0] + src1[1][2] * src2[2][0] + src1[1][3] * src2[3][0];
  dest[1][1] = src1[1][0] * src2[0][1] + src1[1][1] * src2[1][1] + src1[1][2] * src2[2][1] + src1[1][3] * src2[3][1];
  dest[1][2] = src1[1][0] * src2[0][2] + src1[1][1] * src2[1][2] + src1[1][2] * src2[2][2] + src1[1][3] * src2[3][2];
  dest[1][3] = src1[1][0] * src2[0][3] + src1[1][1] * src2[1][3] + src1[1][2] * src2[2][3] + src1[1][3] * src2[3][3];

  dest[2][0] = src1[2][0] * src2[0][0] + src1[2][1] * src2[1][0] + src1[2][2] * src2[2][0] + src1[2][3] * src2[3][0];
  dest[2][1] = src1[2][0] * src2[0][1] + src1[2][1] * src2[1][1] + src1[2][2] * src2[2][1] + src1[2][3] * src2[3][1];
  dest[2][2] = src1[2][0] * src2[0][2] + src1[2][1] * src2[1][2] + src1[2][2] * src2[2][2] + src1[2][3] * src2[3][2];
  dest[2][3] = src1[2][0] * src2[0][3] + src1[2][1] * src2[1][3] + src1[2][2] * src2[2][3] + src1[2][3] * src2[3][3];

  dest[3][0] = src1[3][0] * src2[0][0] + src1[3][1] * src2[1][0] + src1[3][2] * src2[2][0] + src1[3][3] * src2[3][0];
  dest[3][1] = src1[3][0] * src2[0][1] + src1[3][1] * src2[1][1] + src1[3][2] * src2[2][1] + src1[3][3] * src2[3][1];
  dest[3][2] = src1[3][0] * src2[0][2] + src1[3][1] * src2[1][2] + src1[3][2] * src2[2][2] + src1[3][3] * src2[3][2];
  dest[3][3] = src1[3][0] * src2[0][3] + src1[3][1] * src2[1][3] + src1[3][2] * src2[2][3] + src1[3][3] * src2[3][3];
  */
}


// Invert 4x4 matrix
int sdhdf_inv4x4(double m[4][4],double invOut[4][4])
{
  double det;
  double inv[4][4];
  int i,j;

  inv[0][0] = m[1][1]  * m[2][2] * m[3][3] -
    m[1][1]  * m[2][3] * m[3][2] -
    m[2][1]  * m[1][2]  * m[3][3] +
    m[2][1]  * m[1][3]  * m[3][2] +
    m[3][1]  * m[1][2]  * m[2][3] -
    m[3][1]  * m[1][3]  * m[2][2];
  
  inv[1][0] = -m[1][0]  * m[2][2] * m[3][3] +
    m[1][0]  * m[2][3] * m[3][2] +
    m[2][0]  * m[1][2]  * m[3][3] -
    m[2][0]  * m[1][3]  * m[3][2] -
    m[3][0] * m[1][2]  * m[2][3] +
    m[3][0] * m[1][3]  * m[2][2];
  
  inv[2][0] = m[1][0]  * m[2][1] * m[3][3] -
    m[1][0]  * m[2][3] * m[3][1] -
    m[2][0]  * m[1][1] * m[3][3] +
    m[2][0]  * m[1][3] * m[3][1] +
    m[3][0] * m[1][1] * m[2][3] -
    m[3][0] * m[1][3] * m[2][1];
  
  inv[3][0] = -m[1][0]  * m[2][1] * m[3][2] +
    m[1][0]  * m[2][2] * m[3][1] +
    m[2][0]  * m[1][1] * m[3][2] -
    m[2][0]  * m[1][2] * m[3][1] -
    m[3][0] * m[1][1] * m[2][2] +
    m[3][0] * m[1][2] * m[2][1];
  
  inv[0][1] = -m[0][1]  * m[2][2] * m[3][3] +
    m[0][1]  * m[2][3] * m[3][2] +
    m[2][1]  * m[0][2] * m[3][3] -
    m[2][1]  * m[0][3] * m[3][2] -
    m[3][1] * m[0][2] * m[2][3] +
    m[3][1] * m[0][3] * m[2][2];
  
  inv[1][1] = m[0][0]  * m[2][2] * m[3][3] -
    m[0][0]  * m[2][3] * m[3][2] -
    m[2][0]  * m[0][2] * m[3][3] +
    m[2][0]  * m[0][3] * m[3][2] +
    m[3][0] * m[0][2] * m[2][3] -
    m[3][0] * m[0][3] * m[2][2];
  
  inv[2][1] = -m[0][0]  * m[2][1] * m[3][3] +
    m[0][0]  * m[2][3] * m[3][1] +
    m[2][0]  * m[0][1] * m[3][3] -
    m[2][0]  * m[0][3] * m[3][1] -
    m[3][0] * m[0][1] * m[2][3] +
    m[3][0] * m[0][3] * m[2][1];
  
  inv[3][1] = m[0][0]  * m[2][1] * m[3][2] -
    m[0][0]  * m[2][2] * m[3][1] -
    m[2][0]  * m[0][1] * m[3][2] +
    m[2][0]  * m[0][2] * m[3][1] +
    m[3][0] * m[0][1] * m[2][2] -
    m[3][0] * m[0][2] * m[2][1];
  
  inv[0][2] = m[0][1]  * m[1][2] * m[3][3] -
    m[0][1]  * m[1][3] * m[3][2] -
    m[1][1]  * m[0][2] * m[3][3] +
    m[1][1]  * m[0][3] * m[3][2] +
    m[3][1] * m[0][2] * m[1][3] -
    m[3][1] * m[0][3] * m[1][2];
  
  inv[1][2] = -m[0][0]  * m[1][2] * m[3][3] +
    m[0][0]  * m[1][3] * m[3][2] +
    m[1][0]  * m[0][2] * m[3][3] -
    m[1][0]  * m[0][3] * m[3][2] -
    m[3][0] * m[0][2] * m[1][3] +
    m[3][0] * m[0][3] * m[1][2];
  
  inv[2][2] = m[0][0]  * m[1][1] * m[3][3] -
    m[0][0]  * m[1][3] * m[3][1] -
    m[1][0]  * m[0][1] * m[3][3] +
    m[1][0]  * m[0][3] * m[3][1] +
    m[3][0] * m[0][1] * m[1][3] -
    m[3][0] * m[0][3] * m[1][1];
  
  inv[3][2] = -m[0][0]  * m[1][1] * m[3][2] +
    m[0][0]  * m[1][2] * m[3][1] +
    m[1][0]  * m[0][1] * m[3][2] -
    m[1][0]  * m[0][2] * m[3][1] -
    m[3][0] * m[0][1] * m[1][2] +
    m[3][0] * m[0][2] * m[1][1];
  
  inv[0][3] = -m[0][1] * m[1][2] * m[2][3] +
    m[0][1] * m[1][3] * m[2][2] +
    m[1][1] * m[0][2] * m[2][3] -
    m[1][1] * m[0][3] * m[2][2] -
    m[2][1] * m[0][2] * m[1][3] +
    m[2][1] * m[0][3] * m[1][2];
  
  inv[1][3] = m[0][0] * m[1][2] * m[2][3] -
    m[0][0] * m[1][3] * m[2][2] -
    m[1][0] * m[0][2] * m[2][3] +
    m[1][0] * m[0][3] * m[2][2] +
    m[2][0] * m[0][2] * m[1][3] -
    m[2][0] * m[0][3] * m[1][2];
  
  inv[2][3] = -m[0][0] * m[1][1] * m[2][3] +
    m[0][0] * m[1][3] * m[2][1] +
    m[1][0] * m[0][1] * m[2][3] -
    m[1][0] * m[0][3] * m[2][1] -
    m[2][0] * m[0][1] * m[1][3] +
    m[2][0] * m[0][3] * m[1][1];
  
  inv[3][3] = m[0][0] * m[1][1] * m[2][2] -
    m[0][0] * m[1][2] * m[2][1] -
    m[1][0] * m[0][1] * m[2][2] +
    m[1][0] * m[0][2] * m[2][1] +
    m[2][0] * m[0][1] * m[1][2] -
    m[2][0] * m[0][2] * m[1][1];
  
  det = m[0][0] * inv[0][0] + m[0][1] * inv[1][0] + m[0][2] * inv[2][0] + m[0][3] * inv[3][0];
  //  printf("det = %g\n",det);
  if (det == 0)
    return -1;

  det = 1.0 / det;
  
  for (i = 0; i < 4; i++)
    {
      for (j=0;j<4;j++)
	{
	  invOut[i][j] = inv[i][j] * det;
	  //	  printf("output: %g\n",invOut[i][j]);
	}
    }
  return 0;
}

/* Computes the dot product of two vectors */
double sdhdf_dotproduct(double *v1,double *v2)
{
  double dot;
  dot = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
  return dot;
}

void displayMatrix_4x4(double matrix[4][4])
{
  printf("%g %g %g %g\n",matrix[0][0],matrix[1][0],matrix[2][0],matrix[3][0]);
  printf("%g %g %g %g\n",matrix[0][1],matrix[1][1],matrix[2][1],matrix[3][1]);
  printf("%g %g %g %g\n",matrix[0][2],matrix[1][2],matrix[2][2],matrix[3][2]);
  printf("%g %g %g %g\n",matrix[0][3],matrix[1][3],matrix[2][3],matrix[3][3]);
}

// Delta G and Delta Phi (radians)
void sdhdf_setGainPhase(double ma[4][4],double diffGain,double diffPhase)
{
  sdhdf_setIdentity_4x4(ma);
  
  ma[1][0] = ma[0][1] = diffGain/2;
  ma[2][2] = ma[3][3] = cos(diffPhase);
  ma[3][2] = -sin(diffPhase);
  ma[2][3] = sin(diffPhase);
  
}

// Delta G/2 and Delta Phi (radians)
void sdhdf_setGain2Phase(double ma[4][4],double diffGain,double diffPhase)
{
  sdhdf_setIdentity_4x4(ma);
  
  ma[1][0] = ma[0][1] = diffGain;
  ma[2][2] = ma[3][3] = cos(diffPhase);
  ma[3][2] = -sin(diffPhase);
  ma[2][3] = sin(diffPhase);
  
}

// 
void sdhdf_setImperfection(double mcc[4][4],double eps,double phi)
{
  sdhdf_setIdentity_4x4(mcc);
  
  mcc[2][0] =2*eps*cos(phi);
  mcc[3][0] =2*eps*cos(phi);
  mcc[0][2] = 2*eps*cos(phi);
  mcc[0][3] = 2*eps*sin(phi);
}

void sdhdf_setImperfection_full(double mcc[4][4],double eps1,double phi1,double eps2,double phi2)
{
  double a,b,c,d;

  a = eps1*cos(phi1)+eps2*cos(phi2);
  b = eps1*sin(phi1)+eps2*sin(phi2);
  c = eps1*cos(phi1)-eps2*cos(phi2);
  d = eps1*sin(phi1)-eps2*sin(phi2);
    
  sdhdf_setIdentity_4x4(mcc);

  mcc[2][0] = a;
  mcc[3][0] = b;
  mcc[2][1] = c;
  mcc[3][1] = d;

  mcc[0][2] = a;
  mcc[1][2] = -c;
  mcc[0][3] = b;
  mcc[1][3] = -d;

  
  
}

// theta in radians **
void sdhdf_setAstroAngle(double mastro[4][4],double theta)
{
  sdhdf_setIdentity_4x4(mastro);
  mastro[1][1] = mastro[2][2] = cos(2*theta);

  mastro[2][1] = sin(2*theta);
  mastro[1][2] = -sin(2*theta);       
}

// pa in radians **
void sdhdf_setParallacticAngle(double msky[4][4],double pa)
{
  sdhdf_setIdentity_4x4(msky);
  msky[1][1] = msky[2][2] = cos(2*pa);

  msky[2][1] = sin(2*pa);
  msky[1][2] = -sin(2*pa);       
}


void sdhdf_setFeedAlphaChi(double mf[4][4],double alpha,double chi)
{
  sdhdf_setIdentity_4x4(mf);
  mf[1][1] = cos(2*alpha);
  mf[2][1] = sin(2*alpha)*cos(chi);
  mf[3][1] = sin(2*alpha)*sin(chi);

  mf[1][2] = -sin(2*alpha)*cos(chi);
  mf[2][2] = pow(cos(alpha),2) - pow(sin(alpha),2)*cos(2*chi);
  mf[3][2] = -pow(sin(alpha),2)*sin(2*chi);

  mf[1][3] = -sin(2*alpha)*sin(chi);
  mf[2][3] = -pow(sin(alpha),2)*sin(2*chi);
  mf[3][3] = pow(cos(alpha),2) + pow(sin(alpha),2)*cos(2*chi);    
}

// gamma in radians **
void sdhdf_setFeed(double mf[4][4],double gamma)
{
  sdhdf_setIdentity_4x4(mf);

  // CHAGING BACK TO HEILES EQ 15  // This is from Equation 25 of Robishaw, which is different from the Handbook of pulsar astronomy
   mf[1][1] = mf[2][2] = cos(2*gamma);
   mf[2][1] = sin(2*gamma);
   mf[1][2] = -sin(2*gamma);

  
  //  mf[1][1] = mf[3][3] = cos(2*gamma);
  //  mf[3][1] = sin(2*gamma);
  //  mf[1][3] = -sin(2*gamma);
  
}

void sdhdf_convertStokes(double p1,double p2,double p3,double p4,double *stokesI,double *stokesQ,double *stokesU,double *stokesV)
{
  *stokesI = p1 + p2;
  *stokesQ = p1 - p2;
  *stokesU = 2*p3;
  *stokesV = 2*p4;
}

void sdhdf_convertStokesVec(double p1,double p2,double p3,double p4,double *vec)
{
  vec[0] = p1 + p2;
  vec[1] = p1 - p2;
  vec[2] = 2*p3;
  vec[3] = 2*p4;
}

double nonlinearFunc(double x,const double *par,int obsNum)
{
  int k;
  double result=0.0;
  double msky[4][4];
  double mastro[4][4];
  double mhand[4][4];
  double mf[4][4];
  double mcc[4][4];
  double deltaG2,phi,pa,gamma;
  double cc_a,cc_b,cc_c,cc_d;
  double psi;
  double alpha,eps,chi;
  double eps1,eps2,phi1,phi2;
  double stokesMeasured[4];
  double stokesActual[4];
  double mueller[4][4],mueller_inv[4][4];
  double ma[4][4];

  int binNum;
  
  int n,nSeg;

  //  stokesActual[0] = par[0];

  if (model_global==1)
    {
      stokesActual[0] = 1;
      stokesActual[1] = par[0];
      stokesActual[2] = par[1];
      stokesActual[3] = par[2];
      eps             = par[3];
      phi             = par[4];
      alpha           = par[5];
      
      deltaG2         = par[6];
      psi             = par[7];
    }
  else if (model_global==2)
    {
      stokesActual[0] = 1;
      stokesActual[1] = 0.002272; // 0.00187; //0.0018062; /// HARDCODED!!!!!
      stokesActual[2] = -0.17725; // -0.2397; //-0.261829;
      stokesActual[3] = -0.08849; // -0.1526; // -0.18081;
      eps             = par[0];
      phi             = par[1];
      alpha           = par[2];
      
      deltaG2         = par[3];
      psi             = par[4];
    }
  else if (model_global==3)
    {
      stokesActual[0] = 1;
      stokesActual[1] = 0.002272; // 0.00187; //0.0018062; /// HARDCODED!!!!!
      stokesActual[2] = -0.17725; // -0.2397; //-0.261829;
      stokesActual[3] = -0.08849; // -0.1526; // -0.18081;
      eps             = par[0];
      phi             = par[1];
      alpha           = par[2];
      
      deltaG2         = par[3+2*obsNum];
      psi             = par[4+2*obsNum];
    }
  else if (model_global==4)
    {
      stokesActual[0] = 1;
      stokesActual[1] = par[0];
      stokesActual[2] = par[1];
      stokesActual[3] = par[2];
      eps             = par[3];
      phi             = par[4];
      alpha           = par[5];
      
      deltaG2         = par[6+2*obsNum];
      psi             = par[7+2*obsNum];
    }
  else if (model_global==5)
    {
      stokesActual[0] = 1;
      stokesActual[1] = par[0];
      stokesActual[2] = par[1];
      stokesActual[3] = par[2];
      eps             = par[3];
      phi             = par[4];
      alpha           = par[5];
      
      deltaG2         = 0;
      psi             = 0;
    }
  else if (model_global==6)
    {
      stokesActual[0] = 1;
      stokesActual[1] = par[0];
      stokesActual[2] = par[1];
      stokesActual[3] = par[2];

      eps1             = par[3];
      phi1             = par[4];
      eps2             = par[5];
      phi2             = par[6];

      alpha            = par[7];
      chi              = par[8];
      
      deltaG2          = par[9];
      psi              = par[10];
    }

  // Which bin number
  binNum = (int)((double)obsNum/(double)(3*nAct_global));

  
  if (obsNum - binNum*nAct_global*3 < nAct_global)        {n = obsNum - 3*binNum*nAct_global; nSeg = 1;}
  else if (obsNum - binNum*nAct_global*3 < 2*nAct_global) {n = obsNum - 3*binNum*nAct_global - nAct_global; nSeg = 2;}
  else if (obsNum - binNum*nAct_global*3 < 3*nAct_global) {n = obsNum - 3*binNum*nAct_global - 2*nAct_global; nSeg = 3;}

  //  printf("FIT HERE obsNum = %d binNum = %d n = %d nSeg = %d nAct_global = %d nBinFit = %d\n",obsNum,binNum,n,nSeg,nAct_global,nBinFit);

  
  //  printf("Entering with obsNum %d gain = %g %g x = %g, n = %d\n",obsNum,gain[n],phase[n],x,n);
  sdhdf_setIdentity_4x4(mueller);
  sdhdf_setIdentity_4x4(mhand);
  sdhdf_setIdentity_4x4(mcc);

  sdhdf_setGain2Phase(ma,deltaG2,psi);

  if (model_global==6)
    {
      sdhdf_setImperfection_full(mcc,eps1,phi1,eps2,phi2);
      sdhdf_setFeedAlphaChi(mf,alpha,chi);
    }
  else
    {
        sdhdf_setImperfection(mcc,eps,phi);
	sdhdf_setFeed(mf,alpha);
    }
  //  sdhdf_setFeedAlphaChi(mf,alpha,chi);

  // Model FD_HAND = -1
  //  mhand[1][1]=-1; // Q
  //  mhand[3][3]=-1; // V

  sdhdf_mult4x4_replace(mueller,mhand); // Should this be before or after msky? or before after other parts  
  sdhdf_mult4x4_replace(mueller,ma);
  sdhdf_mult4x4_replace(mueller,mcc);
  sdhdf_mult4x4_replace(mueller,mf);
  /*
  mueller[0][0] = 1;
  mueller[1][0] = -2*eps*sin(phi)*sin(2*alpha) + deltaG2*cos(2*alpha);   /// WARNING THIS WAS WRONG IN A PREVIOUS VERSION *******
  mueller[2][0] = 2*eps*cos(phi);
  mueller[3][0] = 2*eps*sin(phi)*cos(2*alpha) + deltaG2*sin(2*alpha);

  mueller[0][1] = deltaG2;
  mueller[1][1] = cos(2*alpha);
  mueller[2][1] = 0;
  mueller[3][1] = sin(2*alpha);

  mueller[0][2] = 2*eps*cos(phi+psi);
  mueller[1][2] = sin(2*alpha)*sin(psi);
  mueller[2][2] = cos(psi);
  mueller[3][2] = -cos(2*alpha)*sin(psi);

  mueller[0][3] = 2*eps*sin(phi+psi);
  mueller[1][3] = -sin(2*alpha)*cos(psi);
  mueller[2][3] = sin(psi);
  mueller[3][3] = cos(2*alpha)*cos(psi);
  */
  
  //  gamma = 45*M_PI/180.0;
  //  sdhdf_setFeed(mf,gamma);
		  
  
  pa      = x*M_PI/180.0;             
  sdhdf_setParallacticAngle(msky,pa);

    sdhdf_setAstroAngle(mastro,M_PI/4.); // NOTE HARDCODE
  //  sdhdf_setAstroAngle(mastro,M_PI/2.); // NOTE HARDCODE


  
  sdhdf_mult4x4_replace(mueller,msky);
  sdhdf_mult4x4_replace(mueller,mastro); 

  
  sdhdf_multMat_vec(mueller,stokesActual,stokesMeasured);

  if (nSeg==1)
    result = stokesMeasured[1]/stokesMeasured[0];
  else if (nSeg==2)
    result = stokesMeasured[2]/stokesMeasured[0];
  else if (nSeg==3)
    result = stokesMeasured[3]/stokesMeasured[0];
    
  return result;
}

void readPol(dSetStruct *dSet,int polNum,float *loadData)
{
  int nchan =  dSet->head->nchan;
  int nbin  =  dSet->head->nbin;
  int nsub  =  dSet->head->nsub;  
  int colnum_data,colnum_datOffs,colnum_datScl,colnum_datWts;
  int status=0;
  float n_fval=0;
  short int n_val =0;
  int initflag =0;
  short int *sval;
  int kk,jj;
  float *datScl;
  float *datOffs;
  float *datWts;
  
  sval         = (short int *)calloc(nchan*nbin,sizeof(short int));
  datScl       = (float *)malloc(sizeof(float)*nchan*nsub);
  datOffs      = (float *)malloc(sizeof(float)*nchan*nsub);
  datWts       = (float *)malloc(sizeof(float)*nchan*nsub);
  
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_data,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);

  fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,1,1+polNum*nchan,nchan,&n_fval,datOffs,&initflag,&status); // Would need updating if more than one subint
  fits_read_col(dSet->fp,TFLOAT,colnum_datScl,1,1+polNum*nchan,nchan,&n_fval,datScl,&initflag,&status);
  fits_read_col(dSet->fp,TSHORT,colnum_data,1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);
  for (kk=0;kk<nchan;kk++)
    {
      for (jj=0;jj<nbin;jj++)
	loadData[kk*nbin+jj] = (sval[kk*nbin+jj]*datScl[kk]+datOffs[kk]);
    }
  free(sval);
  free(datScl);
  free(datOffs);
  free(datWts);
}

void scaleInput(double calAA,double calBB,double calReAB,double calImAB,double *aa,double *bb,double *reAB,double *imAB)
{
  double phase = atan2(calImAB,calReAB);
  double in_aa = *aa;
  double in_bb = *bb;
  double in_reab = *reAB;
  double in_imab = *imAB;


  //  double len = sqrt(calAA*calBB);
  double len = sqrt(pow(calReAB,2)+pow(calImAB,2));

  //  printf("In scale with %g %g %g %g %g %g\n",in_aa,in_bb,in_reab,in_imab,len,phase);
  
  *aa = in_aa/calAA;
  *bb = in_bb/calBB;

  *reAB = (in_reab*cos(phase)+in_imab*sin(phase))/len;
  *imAB = (in_imab*cos(phase)-in_reab*sin(phase))/len;
  //  printf("In scale output = %g %g %g %g\n",*aa,*bb,*reAB,*imAB);
}

void set1MHZ_flagging(int *flag,int nchan)
{
  double freq;
  int i,j;
    
  for (i=0;i<nchan;i++)
    {
      freq = 704+i;
      flag[i]=0;
      if (checkFreq(freq,1857.5,25)==1) flag[i]=1; // OPTUS
      else if (checkFreq(freq,1815,20)==1) flag[i]=1; // TELSTRA
      else if (checkFreq(freq,1728.1,8)==1) flag[i]=1; //Aliased signal
      else if (checkFreq(freq,1982.35,8.6)==1) flag[i]=1; //Aliased signal
      else if (checkFreq(freq,2115,10)==1) flag[i]=1; // Vodafone
      else if (checkFreq(freq,1720,20)==1) flag[i]=1; // Mobile handset
      else if (checkFreq(freq,2142.5,5)==1) flag[i]=1; // Optus
      else if (checkFreq(freq,2147.5,5)==1) flag[i]=1; // Optus
      else if (checkFreq(freq,2167.5,5)==1) flag[i]=1; // Optus
      else if (checkFreq(freq,2312,20)==1) flag[i]=1; // NBN
      else if (checkFreq(freq,2332,20)==1) flag[i]=1; // NBN
      else if (checkFreq(freq,2352,20)==1) flag[i]=1; // NBN
      else if (checkFreq(freq,1890,20)==1) flag[i]=1; // Unknown (Vodafone 4G?)
      else if (checkFreq(freq,1925,10)==1) flag[i]=1; // Unknown (Vodafone 4G?)


      // Satellites
      if (checkFreq(freq,1622,10)==1) flag[i]=1; // Iridium

      // UNEXPLAINED -- JUST DONE BY HAND
      if (checkFreq(freq,2401,256)==1) flag[i]=1;  // BAD SUBBANDS
      if (checkFreq(freq,2681,20)==1) flag[i]=1; 
      if (checkFreq(freq,3453,20)==1) flag[i]=1; 
      
      // Sub-band boundaries
      for (j=0;j<26;j++)
	{
	  
	  if (checkFreq(freq,704+j*128,20)==1) flag[i]=1; 
	}
    }
  
}

int checkFreq(double freq,double f0,double bw)
{
  if (freq >= f0-bw/2. && freq <= f0+bw/2.) return 1;
  return 0;
}
