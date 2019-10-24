// Code to zap pulse profiles
//
//gcc -lm -o pfits_zapProfile pfits_zapProfile.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -O3

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include <stdint.h>
#include "fitsio.h"

#define MAX_ZAP_LIST 512

void calcMinMax(float *timeAve,int nbin,int nchan,float *minVal,float *maxVal,int *imax,int showSub,int selectSB,float f0,float chanbw,int type);
void saveFile(float *datWts, int nsub, int nchan,dSetStruct *dSet,int autoRun);
void get_dm_tbin(dSetStruct *dSet,double *dm,double *tbin);
void formFrequencyTimeAverages(int nsub,int nchan,int nbin,double dm,double tbin,float fref,float *fChan,float *loadAll_pol1,float *loadAll_pol2,float *timeAve,float *freqAve,float *profile,float *datWts,int doDM,int subNum,int selectSB,float f0,float chanbw);
void help();
void displayZapCommand(int *zapFreq,int nchan,int *zapSub,int nsub,float *fChan,char *fname);
int modOp(int a,int b);

//===================== Method 1: =============================================
//Algorithm from N. Wirth’s book Algorithms + data structures = programs of 1976

typedef float elem_type;

#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }

elem_type kth_smallest(elem_type a[], uint16_t n, uint16_t k)
{
  uint64_t i,j,l,m ;
  elem_type x ;
  l=0 ; m=n-1 ;
  while (l<m) {
    x=a[k] ;
    i=l ;
    j=m ;
    do {
      while (a[i]<x) i++ ;
      while (x<a[j]) j-- ;
      if (i<=j) {
	ELEM_SWAP(a[i],a[j]) ;
	i++ ; j-- ;
      }
    } while (i<=j) ;
    if (j<k) l=i ;
    if (k<i) m=j ;
  }
  return a[k] ;
}

#define wirth_median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

//===================== Method 2: =============================================
//This is the faster median determination method.
//Algorithm from Numerical recipes in C of 1992

elem_type quick_select_median(elem_type arr[], uint16_t n)
{
  uint16_t low, high ;
  uint16_t median;
  uint16_t middle, ll, hh;
  low = 0 ; high = n-1 ; median = (low + high) / 2;
  for (;;) {
    if (high <= low) /* One element only */
      return arr[median] ;
    if (high == low + 1) { /* Two elements only */
      if (arr[low] > arr[high])
	ELEM_SWAP(arr[low], arr[high]) ;
      return arr[median] ;
    }
    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])
      ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])
      ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])
      ELEM_SWAP(arr[middle], arr[low]) ;
    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;
    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do ll++; while (arr[low] > arr[ll]) ;
      do hh--; while (arr[hh] > arr[low]) ;
      if (hh < ll)
	break;
      ELEM_SWAP(arr[ll], arr[hh]) ;
    }
    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;
    /* Re-set active partition */
    if (hh <= median)
      low = ll;
    if (hh >= median)
      high = hh - 1;
  }
  return arr[median] ;
}


void help()
{
  printf("pfits_zapProfile -- program to remove RFI from pulsar fold-mode data sets\n");
  printf("Command line arguments when starting software\n");
  printf("-f <filename>\n");
  printf("-tbin <tbin (sec)> - bin time in seconds - should be read from PSRFITS header, but sometimes missing\n");
  printf("\n");
  printf("Key presses:\n");
  printf("1: plot frequency-phase\n");
  printf("2: plot sub-integration-phase\n");
  printf("3: plot profile\n");
  printf("\n");
  printf("b: select specific subband\n");
  printf("+: increase selected subband\n");
  printf("-: decrease selected subband\n");
  printf("h: this help\n");
  printf("m: set the maximum value for the colour-scale in the plots\n");
  printf("p: display paz zap command\n");
  printf("q: quit - note that the file is not automatically saved - must press 's' to save first\n");
  printf("r: recalculate scalings\n");
  printf("s: save zapped file\n");
  printf("u: un-zoom\n");
  printf("X: (right mouse button) delete individual channel or sub-integration\n");
  printf("z: zoom into a region on the plot\n");
  printf("Z: zap a range of channels or subintegrations\n");  
}


int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,j;

  float *loadAll_pol1;
  float *loadAll_pol2;
  
  float *freqAve;
  float *timeAve;
  float *profile;
  float *profileX;
  
  short int *sval;
  short int n_val =0;
  float n_fval=0;
  long nchan,nsub,nbin,npol;
  long ii,jj,kk;
  int colnum_data;
  int colnum_datOffs;
  int colnum_datScl;
  int colnum_datWts;
  int initflag = 0;
  int polNum = 0;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float tr[6];
  int   imaxVal_freqAve,imaxVal_timeAve,imaxVal_profile;
  float minVal_freqAve,maxVal_freqAve;
  float minVal_timeAve,maxVal_timeAve;
  float minVal_profile,maxVal_profile;
  
  float *datScl;
  float *datOffs;
  float *datWts;
  float minx,maxx,miny,maxy;
  float *fChan,fref;
  
  double dm=0;
  double tbin=0;
  
  float mx,my;
  char key;
  int didDelete=0;
  int doDM=1;
  int selectSB=-1;
  
  float val;
  
  double mean1,mean2;
  int plotType=1;
  int dedispBin;
  double tdiff;
  int binOff;
  int autoScale=1;
  int helpFindBright=1;
  float fx[2],fy[2];
  float maxValChan;
  float f0,f1,chanbw;
  
  int *zapFreq;
  int *zapSub;

  float threshold=-1;
  int autoRun=0;
  int autoSave=0;

  int show1Sub=-1;
  int showSub=-1;

  char zapListFile[128];
  int useZapListFile=0;
  float fs[MAX_ZAP_LIST],fe[MAX_ZAP_LIST];
  int  nZapList=0;
  FILE *fin;
  
  // Initialise everything
  initialise(&dSet,debug);
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-tbin")==0)
	sscanf(argv[++i],"%lf",&tbin);
      else if (strcmp(argv[i],"-autoSave")==0)
	autoSave=1;
      else if (strcmp(argv[i],"-nodm")==0)
	doDM=0;
      else if (strcmp(argv[i],"-zapList")==0)
	{
	  strcpy(zapListFile,argv[++i]);
	  useZapListFile=1;
	}
      else if (strcmp(argv[i],"-autoRun")==0)
	{
	  sscanf(argv[++i],"%f",&threshold);
	  autoRun=1;
	}
    }
       
  if (useZapListFile==1)
    {
      printf("Got here\n");
      if (!(fin = fopen(zapListFile,"r")))
	{
	  printf("Unable to open file %s\n",zapListFile);
	  exit(1);
	}
      while (!feof(fin))
	{
	  if (fscanf(fin,"%f %f",&fs[nZapList],&fe[nZapList])==2)
	    nZapList++;
	}
      fclose(fin);
      printf("Loaded %d zap commands\n",nZapList);
    }

  
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  // Now load the data to plot
  // Average in polarisation
  // Produce data files also time averaged and frequency averaged

  nchan = dSet->head->nchan;
  nbin  = dSet->head->nbin;
  nsub  = dSet->head->nsub;
  npol  = dSet->head->npol;

  get_dm_tbin(dSet,&dm,&tbin);
  
  printf("Loaded header\n");
  printf("Number of channels = %d\n",nchan);
  printf("Number of bins = %d\n",nbin);
  printf("Number of sub-integrations = %d\n",nsub);
  printf("Number of polarisations = %d\n",npol);
  printf("Dispersion measure = %g\n",dm);

  loadAll_pol1 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
  loadAll_pol2 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);

  timeAve  = (float *)calloc(nchan*nbin,sizeof(float));
  freqAve  = (float *)calloc(nsub*nbin,sizeof(float));
  profile  = (float *)calloc(nbin,sizeof(float));
  profileX = (float *)calloc(nbin,sizeof(float));
  fChan    = (float *)malloc(sizeof(float)*nchan);
  sval     = (short int *)calloc(nchan*nbin,sizeof(short int));

  zapFreq = (int *)calloc(nchan,sizeof(int));
  zapSub = (int *)calloc(nsub,sizeof(int));
    
  for (ii=0;ii<nbin;ii++)
    profileX[ii] = ii;
  
  datScl = (float *)malloc(sizeof(float)*nchan*nsub);
  datOffs = (float *)malloc(sizeof(float)*nchan*nsub);
  datWts = (float *)malloc(sizeof(float)*nchan*nsub);
  
  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_data,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);

  printf("Loading all data\n");
  for (ii=0;ii<nchan;ii++)
    fChan[ii] = dSet->head->chanFreq[ii];
  fref = dSet->head->freq;
    
  f0 = fChan[0];
  f1 = fChan[nchan-1];
  chanbw = fChan[1]-fChan[0];
  
  printf("Reference freq = %g\n",fref);
  
  for (ii=0;ii<nsub;ii++)
    {
      // Polarisation 1
      polNum=0;
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datWts,ii+1,1,nchan,&n_fval,datWts+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);
	   
      for (kk=0;kk<nchan;kk++)
	{
	  mean1=0;
	  for (jj=0;jj<nbin;jj++)
	    {
 	      val = (sval[kk*nbin+jj]*datScl[ii*nchan+kk]+datOffs[ii*nchan+kk]);
	      loadAll_pol1[ii*nbin*nchan+kk*nbin+jj] = val;
	      mean1 += val;
	    }
	  for (jj=0;jj<nbin;jj++)
	    loadAll_pol1[ii*nbin*nchan+kk*nbin+jj] -= mean1/(double)nbin;
	}
      
      // Polarisation 2
      polNum=1;
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);
      // Don't need to read the weights as not polarisation dependent
      fits_read_col(dSet->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);

      for (kk=0;kk<nchan;kk++)
	{
	  mean2=0;
	  for (jj=0;jj<nbin;jj++)
	    {
 	      val = (sval[kk*nbin+jj]*datScl[ii*nchan+kk]+datOffs[ii*nchan+kk]);
	      loadAll_pol2[ii*nbin*nchan+kk*nbin+jj] = val;
	      mean2+=val;
	    }	      
	  for (jj=0;jj<nbin;jj++)
	    loadAll_pol2[ii*nbin*nchan+kk*nbin+jj] -= mean2/(double)nbin;
	}
    }

  formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
  printf("Completed loading all data\n");
  
  // Get minimum and maximum values
  calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
  calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
  calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);
  printf("Frequency averaging: Minimum value = %g, maximum value = %g (subint for max = %d)\n",minVal_freqAve,maxVal_freqAve,imaxVal_freqAve);
  printf("Time averaging: Minimum value = %g, maximum value = %g (freq. channel for max = %d)\n",minVal_timeAve,maxVal_timeAve,imaxVal_timeAve);

  if (autoRun==1)
    {
      int i0,i1,ii,jj;

      i0 = 0;
      i1 = nchan;
      
      for (i=i0;i<i1;i++)
	{
	  for (j=0;j<nbin;j++)
	    {
	      if (timeAve[i*nbin+j] > threshold)
		{
		  zapFreq[i] = 1;
		  for (jj=0;jj<nsub;jj++)
		    datWts[jj*nchan+(int)i] = 0;
		  
		  for (jj=0;jj<nbin;jj++)
		    timeAve[(int)i*nbin+jj] = 0;
		}
	    }
	}
      saveFile(datWts,nsub,nchan,dSet,autoRun);      
    }
  else
    {
      // Do the plot
      
      tr[0] = 0;  tr[1] = 1;  tr[2] = 0;
      tr[3] = 0;  tr[4] = 0;  tr[5] = 1;
      
      //  maxVal = 2000;
      
      cpgbeg(0,"/xs",1,1);
      cpgask(0);
      
      minx = 0;  maxx = nbin;  miny = 0;  maxy = nchan;
      
      help();
      
      
      do {
	cpgeras();
	//	cpgenv(minx,maxx,miny,maxy,0,1);
	//	cpgsvp(0.2,0.9,0.2,0.9);
	if (plotType==1)
	  {
	    cpgswin(minx,maxx,f0+miny*chanbw,f0+maxy*chanbw);
	    cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
	    	    
	    cpgswin(minx,maxx,miny,maxy);
	    cpgbox("ABCTSN",0,0,"ABCTSM",0,0);
	  }
	else
	  {
	    cpgswin(minx,maxx,miny,maxy);
	    cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
	  }
	    //
	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
	if (plotType==1)
	  {
	    int i0,i1;
		
	    if (autoScale==1)
	      {
		int tt=1;

		
		if (miny < maxy) {i0 = (int)miny; i1 = (int)maxy;}
		else {i0 = (int)maxy; i1 = (int)miny;}
		if (i0<0)       i0=0;
		if (i1 > nchan) i1=nchan;

		for (i=i0;i<i1;i++)
		  {
		    for (j=0;j<nbin;j++)
		      {
			if (tt==1)
			  {minVal_timeAve = maxVal_timeAve = timeAve[(i)*nbin+j]; tt=0;}
			else if (minVal_timeAve > timeAve[(i)*nbin+j]) minVal_timeAve = timeAve[(i)*nbin+j];
			else if (maxVal_timeAve < timeAve[(i)*nbin+j]) {maxVal_timeAve = timeAve[(i)*nbin+j]; maxValChan = (i);}
		      }
		  }
	      }
	    //	    maxVal_timeAve = 4000;
	    printf("At this plotting point with %g %g\n",minVal_timeAve,maxVal_timeAve);

	    tr[0] = 0;  tr[1] = 1;  tr[2] = 0;
	    tr[3] = 0;  tr[4] = 0;  tr[5] = 1;
      

	    cpgimag(timeAve,nbin,nchan,1,nbin,1,nchan,minVal_timeAve,maxVal_timeAve,tr);
	    if (helpFindBright==1)
	      {
		fx[0] = minx; fx[1] = maxx;
		fy[0] = fy[1] = maxValChan+1;
		printf("Drawing line at %f %f %f %f %f\n",maxValChan,fx[0],fy[0],fx[1],fy[1]);
		cpgsls(1);cpgsci(3);
		cpgline(2,fx,fy);
		cpgsls(1); cpgsci(1);
	      }
	    if (selectSB > -1)
	      {
		char sb[128];

		sprintf(sb,"File: %s Sub-band %d",dSet->fileName,selectSB);
		if (show1Sub == 1)
		  {
		    sprintf(sb,"File: %s Sub-band %d, sub-int = %d",dSet->fileName,selectSB,showSub);
		  }

		cpglab("Bin number","Frequency (MHz)",sb);
	      }
	    else
	      {
		char sb[128];
		sprintf(sb,"File: %s",dSet->fileName);

		cpglab("Bin number","Frequency (MHz)",sb);
	      }
	    cpgmtxt("R",2.5,0.5,0.5,"Frequency channel");
	  }
	else if (plotType==2)
	  {
	    char sb[128];
	    
	    cpgimag(freqAve,nbin,nsub,1,nbin,1,nsub,minVal_freqAve,maxVal_freqAve,tr);
	    if (selectSB > -1)
	      sprintf(sb,"File: %s Sub-band %d",dSet->fileName,selectSB);
	    else
	      sprintf(sb,"File: %s",dSet->fileName);
	    cpglab("Bin number","Subint number",sb);
	  }
	else if (plotType==3)
	  {
	    cpgline(nbin,profileX,profile);
	    cpglab("Bin number","","");
	  }

	// Overplot areas to zap
	if (useZapListFile==1)
	  {
	    cpgsfs(3);
	    for (i=0;i<nZapList;i++)
	      cpgrect(minx,maxx,(fs[i]-f0)/chanbw+1,(fe[i]-f0)/chanbw+1);
	    cpgsfs(1);
	  }




	cpgcurs(&mx,&my,&key);
	if (key=='z')
	  {
	    float mx2,my2;
	    cpgband(3,0,mx,my,&mx2,&my2,&key);
	    if (my > my2) {miny = my2; maxy = my;}
	    if (my2 > my) {miny = my; maxy = my2;}
	  }
	else if (key=='k' && useZapListFile==1) // Kill using specified zap list in region
	  {
	    int j;
	    cpgsfs(3);
	    printf("In here with %d zaps (miny = %g, maxy = %g)\n",nZapList,miny,maxy);
	  
	    for (i=0;i<nZapList;i++)
	      {
		for (j=(int)((fs[i]-f0)/chanbw+1);j<(int)((fe[i]-f0)/chanbw+1);j++)
		  {
		    if (j>=miny && j<=maxy)
		      {
			zapFreq[j] = 1;
			
			for (jj=0;jj<nsub;jj++)
			  datWts[jj*nchan+j] = 0;
			for (jj=0;jj<nbin;jj++)
			  timeAve[j*nbin+jj] = 0;
		      }
		  }
	      }
	    calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
	       
	  }
	else if (key=='a') // Auto zap using median zapping
	  {
	    int ii,jj,kk,nk;
	    float *b;
	    float *median,med,s1,s2;
	    float *min,*max,*mean,*sdev;
	    float threshold1,threshold2,maxVal,minVal;	   
	    double stat;
	    int b0,b1;

	    b = (float *)malloc(sizeof(float)*nbin);
	    median = (float *)malloc(sizeof(float)*nsub*nchan);
	    min = (float *)malloc(sizeof(float)*nsub*nchan);
	    max = (float *)malloc(sizeof(float)*nsub*nchan);
	    mean = (float *)malloc(sizeof(float)*nsub*nchan);
	    sdev = (float *)malloc(sizeof(float)*nsub*nchan);
	    
	    printf("Select baseline bin range (0 %d) ",nbin);
	    scanf("%d %d",&b0,&b1);
	    
	    printf("Select threshold1 (~5000 is a good number) ");
	    scanf("%f",&threshold1);
	    printf("Select threshold2 (~5 is a good number) ");
	    scanf("%f",&threshold2);
	    printf("Automatic zapping\n");
	    nk=0;
	    for (jj=0;jj<nsub;jj++)
	      {
		for (ii=0;ii<nchan;ii++)
		  {
		    s1=s2=0;
		    nk=0;
		    for (kk=b0;kk<b1;kk++)
		      {
			b[nk] = loadAll_pol1[jj*nbin*nchan+ii*nbin+kk]+loadAll_pol2[jj*nbin*nchan+ii*nbin+kk];
			if (nk==0)
			  maxVal=minVal=b[nk];
			s1+=b[nk];
			s2+=pow(b[nk],2);
			if (maxVal < b[nk]) maxVal = b[nk];
			if (minVal > b[nk]) minVal = b[nk];
			nk++;
		      }
		    mean[jj*nchan+ii]   = s1/(double)nk;
		    sdev[jj*nchan+ii]   = sqrt(s2/(double)nk-(pow(s1/(double)nk,2)));
		    median[jj*nchan+ii] = quick_select_median(b,nk);
		    max[jj*nchan+ii]    = maxVal;
		    min[jj*nchan+ii]    = minVal;
		  }
	      }
	    med = quick_select_median(median,nchan*nk);
	    printf("Median = %g for all data\n",med);

	    for (jj=0;jj<nsub;jj++)
	      {
		for (ii=0;ii<nchan;ii++)
		  {
		    stat = (max[jj*nchan+ii]-median[jj*nchan+ii])/sdev[jj*nchan+ii];
		    // if (fabs(median[jj*nchan+ii]) > fabs(threshold*med))
		    if (fabs(max[jj*nchan+ii]) > fabs(threshold1*med))
		      datWts[jj*nchan+ii] = 0;
		    if (fabs(min[jj*nchan+ii]) > fabs(threshold1*med))
		      datWts[jj*nchan+ii] = 0;
		       //		    printf("stat = %g %g %g %g\n",stat,max[jj*nchan+ii],median[jj*nchan+ii],sdev[jj*nchan+ii]);
		    if (fabs(stat) > threshold2)
		      datWts[jj*nchan+ii] = 0;
		  }
		  
	      }
	    formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
	    
	    // Get minimum and maximum values
	    calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
	    calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
	    calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);

	    free(b);
	    free(min);
	    free(max);
	    free(mean);
	    free(sdev);
	    free(median);

	  }
	else if (key=='A') // Left mouse click
	  {
	    int ii,kk,jj;
	    int iclose = (int)(my+0.5);
	    float maxVal,minVal;
	    float b[nbin];
	    
	    if (plotType==1)
	      {
		printf("Closest channel number is %d\n",iclose);

		if (iclose >0 && iclose < nchan)
		  {
		    for (jj=0;jj<nsub;jj++)
		      {
			
			for (kk=0;kk<nbin;kk++)
			  {
			    b[kk] = loadAll_pol1[jj*nbin*nchan+iclose*nbin+kk]+loadAll_pol2[jj*nbin*nchan+iclose*nbin+kk];
			    if (kk==0)
			      {
				maxVal = b[kk];
				minVal = b[kk];
			      }
			    if (maxVal < b[kk]) maxVal = b[kk];
			    if (minVal > b[kk]) minVal = b[kk];
			  }
			printf("Subint: %d, Minimum value = %g Maximum value = %g, median = %g weighting = %g\n",jj,minVal,maxVal,quick_select_median(b,nbin),datWts[jj*nchan+iclose]);
		      }
		  }
		
	      }
	  }
	else if (key=='t') // Text output
	  {
	    int ii,jj;
	    float mean=0;
	    long nc=0;
	    for (jj=(int)minx;jj<=(int)maxx;jj++)
	      {
		for (ii=(int)miny;ii<=(int)maxy;ii++)
		  {
		    if (plotType==1)
		      {
			printf("bin: %d chan: %d %g\n",ii,jj,timeAve[ii*nbin+jj]);
			mean+=timeAve[ii*nbin+jj];
			nc++;
		      }
		  }
	      }
	    printf("Mean = %g\n",mean/(double)nc);
	  }
	  
	else if (key=='l')  // For Lawrence
	  {
	    zapFreq[(int)maxValChan] = 1;
	    for (jj=0;jj<nsub;jj++)
	      datWts[jj*nchan+(int)maxValChan] = 0;
	    
	    for (jj=0;jj<nbin;jj++)
	      timeAve[(int)maxValChan*nbin+jj] = 0;

	    calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);

	  }
	else if (key=='L')
	  {
	    int tt=1;
	    int i0,i1;
	    
	    printf("Maximum value = %g (minimum value = %g)\n",maxVal_timeAve,minVal_timeAve);
	    printf("Enter threshold (for zap if maxVal > threshold*median in each channel and subintegration) ");
	    scanf("%f",&threshold);
	    
	    if (miny < maxy) {i0 = (int)miny; i1 = (int)maxy;}
	    else {i0 = (int)maxy; i1 = (int)miny;}
	    if (i0<0) i0=0;
	    if (i1 > nchan) i1=nchan;
	    
	    for (i=i0;i<i1;i++)
	      {
		for (j=0;j<nbin;j++)
		  {
		    if (timeAve[i*nbin+j] > threshold)
		      {
			zapFreq[i] = 1;
			for (jj=0;jj<nsub;jj++)
			  datWts[jj*nchan+(int)i] = 0;
			
			for (jj=0;jj<nbin;jj++)
			  timeAve[(int)i*nbin+jj] = 0;
		      }
		  }
	      }
	    
	  }
	else if (key=='h')
	  help();
	else if (key=='r')
	  {
	    formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
	    calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
	    calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
	    calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);
	    
	  }
	else if (key=='m')
	  {
	    if (plotType==1)printf("Current maximum value is %g\n",maxVal_timeAve);
	    else if (plotType==2)printf("Current maximum value is %g\n",maxVal_freqAve);
	    printf("Please enter new maximum value ");
	    if (plotType==1)scanf("%f",&maxVal_timeAve);
	    else if (plotType==2)scanf("%f",&maxVal_freqAve);	
	  }
	else if (key=='1') {
	  plotType=1;   minx = 0;  maxx = nbin;  miny = 0;  maxy = nchan;
	  if (selectSB > -1)
	    {
	      int i0,i1;
	      
	      i0 = (704+selectSB*128-f0)/chanbw;
	      i1 = (704+(selectSB+1)*128-f0)/chanbw;
	      
	      miny = i0; // Channel number
	      maxy = i1;
	    }

	  if (didDelete==1)
	    formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
	  didDelete=0;
	}
	else if (key=='2') {
	  plotType=2;   minx = 0;  maxx = nbin;  miny = 0;  maxy = nsub+1;
	  //	  if (didDelete==1)
	  formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
	  calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
	  calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
	  calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);		
	  
	    //	  didDelete=0;
	}    
	else if (key=='3') {
	  plotType=3;   minx = 0;  maxx = nbin;  miny = minVal_profile;  maxy = maxVal_profile;
	  if (didDelete==1)
	    formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
	  didDelete=0;
	}
	else if (key=='p')
	  {
	    displayZapCommand(zapFreq,nchan,zapSub,nsub,fChan,"");
	  }
	else if (key=='u')
	  {
	    printf("Unzoom\n");
	    if (selectSB > -1 && plotType==1)
	      {
		int i0,i1;
		
		i0 = (704+selectSB*128-f0)/chanbw;
		i1 = (704+(selectSB+1)*128-f0)/chanbw;

		miny = i0;
		maxy = i1;
	      }
	    else
	      {
		if (plotType==1)
		  { miny = 0; maxy = nchan; minx = 0; maxx = nbin;}
		else if (plotType==2)
		  { miny = 0; maxy = nsub; minx = 0; maxx = nbin;}
		else if (plotType==3)
		  { miny = minVal_profile; maxy = maxVal_profile; minx = 0; maxx = nbin;}
	      }
	  }
	else if (key=='X' && plotType==1)
	  {
	    int zapChannel;
	    didDelete=1;
	    zapChannel = (int)(my-0.5); // Check -0.5
	    printf("Zapping channel %d\n",zapChannel);
	    zapFreq[zapChannel] = 1;
	    for (jj=0;jj<nsub;jj++)
	      datWts[jj*nchan+zapChannel] = 0;
	    
	    for (jj=0;jj<nbin;jj++)
	      timeAve[zapChannel*nbin+jj] = 0;
	    //	    formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);	    
	    calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
	    printf("Minimum value = %g, maximum value = %g, channel number for maximum value = %d plotType=1\n",minVal_timeAve,maxVal_timeAve,imaxVal_timeAve);
	  }
	else if (key=='X' && plotType==2)
	  {
	    int zap_sub;
	    didDelete=1;
	    zap_sub = (int)(my-0.5); // Check -0.5
	    zapSub[zap_sub] = 1;
	    printf("Zapping sub-int %d\n",zap_sub);
	    for (jj=0;jj<nchan;jj++)
	      datWts[zap_sub*nchan+jj] = 0;
	    
	    for (jj=0;jj<nbin;jj++)
	      freqAve[zap_sub*nbin+jj] = 0;
	    calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
	    printf("Minimum value = %g, maximum value = %g, subint number for maximum value = %d\n",minVal_freqAve,maxVal_freqAve,imaxVal_freqAve);
	  }
	else if (key=='X' && plotType==3)
	  printf("Sorry -- cannot delete bins in a profile\n");
	else if (key=='b') // Select specific subband
	  {
	    int i0,i1;
	    
	    if (selectSB==-1)
	      {
		printf("Selecting sub-band 0\n");
		selectSB=0;
		i0 = (704-f0)/chanbw;
		i1 = (704+128-f0)/chanbw;
		if (plotType==1)
		  {
		    miny = i0; // Channel number
		    maxy = i1;
		  }
		else
		  formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
		calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
		calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
		calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);		
	      }
	    else
	      {
		printf("De-selecting sub-band\n");
		selectSB=-1;
		if (plotType==1)
		  {
		    miny = 0;
		    maxy = nchan;
		  }
		formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
		calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
		calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
		calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);		

	      }
	  }
	else if (key=='+') // Increasing subband
	  {
	    int i0,i1;
	    selectSB++;
	    if (selectSB > 25)
	      selectSB = 25;
	    printf("Selecting sub-band %d\n",selectSB);
	    i0 = (704+selectSB*128-f0)/chanbw;
	    i1 = (704+(selectSB+1)*128-f0)/chanbw;
	    if (plotType==1)
	      {
		miny = i0; // Channel number
		maxy = i1;
	      }
	    //	    else
	    //	    formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
	    //	    calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
	    //	    calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
	    //	    calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);
	    
	  }
	else if (key=='i') // Show only a single sub-integration
	  {
	    show1Sub*=-1;
	    if (show1Sub == -1)
	      showSub=-1;
	    else
	      showSub=0;
	    
	    formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
		// Get minimum and maximum values
	    calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
	    calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
	    calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);

	  }
	else if (key==']' && show1Sub==1) // Increase subint
	  {
	    showSub++;
	    if (showSub == nsub) showSub=nsub-1;
	    formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
		// Get minimum and maximum values
	    calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
	    calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
	    calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);

	  }
	else if (key=='[' && show1Sub==1) // Decrease subint
	  {
	    showSub--;
	    if (showSub < 0) showSub=0;
	    formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
		// Get minimum and maximum values
	    calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
	    calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
	    calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);

	  }
	else if (key=='-') // Decreasing subband
	  {
	    int i0,i1;
	    selectSB--;
	    if (selectSB < 0)
	      selectSB = 0;
	    printf("Selecting sub-band %d\n",selectSB);
	    i0 = (704+selectSB*128-f0)/chanbw;
	    i1 = (704+(selectSB+1)*128-f0)/chanbw;
	    if (plotType==1)
	      {
		miny = i0; // Channel number
		maxy = i1;
	      }
	    //	    else
	    
	  }
	else if (key=='s') // Save new file
	  {
	    saveFile(datWts,nsub,nchan,dSet,autoRun);	    
	  }
	else if (key=='Z' && plotType==1)
	  {
	    float mx2,my2;
	    int zapChannel;
	    float min,max;
	    int kk0,kk1;
	    
	    didDelete=1;
	    cpgband(3,0,mx,my,&mx2,&my2,&key);
	    if (my > my2) {min = my2; max = my;}
	    if (my2 > my) {min = my; max = my2;}
	    
	    if (selectSB > -1)
	      {
		printf("Selecting specific subband\n");
		kk0 = (704+selectSB*128-f0)/chanbw;
		kk1 = (704+(selectSB+1)*128-f0)/chanbw;
	      }
	    else
	      {		    
		kk0 = 0;
		kk1 = nchan;
	      }
	    
	    
	    for (zapChannel = (int)(min-0.5);zapChannel <= (int)(max-0.5);zapChannel++)
	      {
		if (zapChannel > 0 && zapChannel < nchan)
		  {
		    zapFreq[zapChannel] = 1;
		    printf("Zapping channel %d\n",zapChannel);
		    for (jj=0;jj<nsub;jj++)
		      datWts[jj*nchan+(zapChannel)] = 0;
		    
		    for (jj=0;jj<nbin;jj++)
		      timeAve[(zapChannel)*nbin+jj] = 0;
		  }
	      }
	    
	    //	    formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts,doDM,showSub,selectSB,f0,chanbw);
	    calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
	    printf("Minimum value = %g, maximum value = %g, channel number for maximum value = %d plotType=1\n",minVal_timeAve,maxVal_timeAve,imaxVal_timeAve);
	  }
	else if (key=='Z' && plotType==2)
	  {
	    float mx2,my2;
	    int zap_sub;
	    float min,max;
	    didDelete=1;
	    cpgband(3,0,mx,my,&mx2,&my2,&key);
	    if (my > my2) {min = my2; max = my;}
	    if (my2 > my) {min = my; max = my2;}
	    
	    for (zap_sub = (int)(min-0.5);zap_sub <= (int)(max-0.5);zap_sub++)
	      {
		if (zap_sub > 0 && zap_sub < nsub)
		  {
		    zapSub[zap_sub] = 1;
		    printf("Zapping subint %d\n",zap_sub);
		    for (jj=0;jj<nchan;jj++)
		      datWts[zap_sub*nchan+jj] = 0;
		    
		    for (jj=0;jj<nbin;jj++)
		      freqAve[zap_sub*nbin+jj] = 0;
		  }
	      }
	    calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
	    printf("Minimum value = %g, maximum value = %g, subint number for maximum value = %d\n",minVal_freqAve,maxVal_freqAve,imaxVal_freqAve);
	    
	  }
	else if (key!='q')
	  printf("Unknown key press %c\n",key);
      } while (key != 'q');
      displayZapCommand(zapFreq,nchan,zapSub,nsub,fChan,"");
      if (autoSave==1)
	{
	  char outname[1024];
	  
	  saveFile(datWts,nsub,nchan,dSet,1);
	  sprintf(outname,"%s.zapCmd",dSet->fileName);
	  displayZapCommand(zapFreq,nchan,zapSub,nsub,fChan,outname);
	}
      cpgend();  
      if (status)
	{
	  fits_report_error(stderr,status);
	  exit(1);
	}
    }
  // De-allocate the memory
  //  pfitsCloseFile(dSet,debug);
  deallocateMemory(&dSet,debug);


  free(zapFreq);
  free(zapSub);
  free(loadAll_pol1);
  free(loadAll_pol2);
  free(timeAve);
  free(freqAve);
  free(profile);
  free(profileX);
  free(sval);
  free(fChan);
  free(datScl);
  free(datOffs);
  free(datWts);

}

//  calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve,showSub,selectSB,f0,chanbw,1);
//  calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve,showSub,selectSB,f0,chanbw,2);
//  calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile,showSub,selectSB,f0,chanbw,3);

void calcMinMax(float *timeAve,int nbin,int nchan,float *minVal,float *maxVal,int *imax,int showSub,int selectSB,float f0,float chanbw,int type)
{
  int ii,jj;
  int ii0,ii1;

  if (type==2)
    {
      if (selectSB == -1)
	{
	  ii0 = 0; ii1 = nchan;
	}
      else
	{
	  ii0 = (704+selectSB*128-f0)/chanbw;
	  ii1 = (704+(selectSB+1)*128-f0)/chanbw;
	}
    }
  else
    {
      ii0=0;
      ii1 = nchan;
    }
  
  for (ii=ii0;ii<ii1;ii++)
    {
      for (jj=0;jj<nbin;jj++)
	{
	  if (jj==0 && ii==ii0)
	    {
	      *imax = 0;
	      *minVal = *maxVal = timeAve[(ii-ii0)*nbin+jj];
	    }
	  else if (*minVal > timeAve[(ii-ii0)*nbin+jj]) *minVal = timeAve[(ii-ii0)*nbin+jj];
	  else if (*maxVal < timeAve[(ii-ii0)*nbin+jj]) {*maxVal = timeAve[(ii-ii0)*nbin+jj]; *imax = (ii-ii0);}
	}
    }
  printf("minimum value = %g, maximum value = %g\n",*minVal,*maxVal);
}

void saveFile(float *datWts, int nsub, int nchan, dSetStruct *dSet,int autoRun)
{
  char outname[1024];
  char temp[1024];
  fitsfile *outfptr;
  int status=0;
  int hdu=1;
  int colnum_datWts;
  long ii;

  if (autoRun==0)
    {
      printf("Enter output filename ");
      scanf("%s",temp);
      sprintf(outname,"!%s",temp);
      printf("Output file name is %s\n",outname);
    }
  else
    {
      sprintf(outname,"!%s.autoZap",dSet->fileName);
    }
  fits_create_file(&outfptr, outname, &status);
  if (status) {fits_report_error(stderr, status); printf("NOT WRITING FILE\n"); status=0; return;}
  
  /* Copy every HDU until we get an error */
  while( !fits_movabs_hdu(dSet->fp, hdu++, NULL, &status) )
    {
      fits_copy_hdu(dSet->fp, outfptr, 0, &status);
    }
  /* Reset status after normal error */
  if (status == END_OF_FILE) status = 0;

  // Update the DAT_WTS column
  fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);
  for (ii=0;ii<nsub;ii++)
    fits_write_col(outfptr,TFLOAT,colnum_datWts,ii+1,1,nchan,datWts+ii*nchan,&status);
  
  fits_close_file(outfptr, &status);
  if (status) {fits_report_error(stderr, status); printf("FILE MAY BE CORRUPT\n"); status=0; return;}
  printf("Data file written (%s)\n",outname);
  
}

void get_dm_tbin(dSetStruct *dSet,double *dm,double *tbin)
{
  int status=0;
  float fdm,ftbin;

  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_read_key(dSet->fp,TFLOAT,"DM",&fdm,NULL,&status);
  *dm = fdm;
  if (*tbin == 0)
    {
      fits_read_key(dSet->fp,TFLOAT,"TBIN",&ftbin,NULL,&status);
      *tbin = ftbin;
    }
  printf("Loaded DM = %g\n",*dm);
  printf("Loaded TBIN = %g\n",*tbin);
}

void formFrequencyTimeAverages(int nsub,int nchan,int nbin,double dm,double tbin,float fref,float *fChan,float *loadAll_pol1,float *loadAll_pol2,float *timeAve,float *freqAve,float *profile,float *datWts,int doDM,int showSub,int selectSB,float f0,float chanbw)
{
  long int ii,jj,kk;
  double tdiff;
  int binOff;
  int dedispBin;
  int t=0;
  int ii0,ii1;
  int kk0,kk1;
  
  if (showSub==-1)
    {
      ii0 = 0;
      ii1 = nsub;
    }
  else
    {
      ii0 = showSub;
      ii1 = ii0+1;
    }

  if (selectSB > -1)
    {
      kk0 = (704+selectSB*128-f0)/chanbw;
      kk1 = (704+(selectSB+1)*128-f0)/chanbw;
      printf("Selecting specific subband have %d and %d\n",kk0,kk1);

    }
  else
    {

      kk0 = 0;
      kk1 = nchan;
    }
  
  printf("Forming freqTime\n");
  printf("Resetting arrays to zero\n");
  // Must be quicker way to reset all to 0
  for (ii=0;ii<nsub*nbin;ii++)
      freqAve[ii]=0;
  for (ii=0;ii<nchan*nbin;ii++)
      timeAve[ii]=0;
  for (ii=0;ii<nbin;ii++)
    profile[ii] = 0;

  printf("Forming frequency and time averaged data sets\n");
  // Form frequency and time averaged data sets
  for (ii=0;ii<nsub;ii++)
    {
      if (showSub==-1 || ii==showSub)
	{
	  for (kk=0;kk<nchan;kk++)
	    {
	      if (doDM==1)
		{
		  tdiff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(fChan[kk]/1000.0,-2));
		  binOff = (int)(tdiff/tbin+0.5);
		}
	      else
		binOff=0;
	      
	      
	      if (datWts[ii*nchan+kk] > 0)
		{
		  for (jj=0;jj<nbin;jj++)
		    {
		      //
		      //		  dedispBin = jj+binOff;
		      //		  if (dedispBin >= nbin)
		      dedispBin = modOp(jj+binOff,nbin);
		      
		      //		  else
		      //		    while (dedispBin < 0) dedispBin+=nbin;
		      //		  dedispBin = jj+binOff;
		      //
		      //		  while (dedispBin >= nbin) dedispBin-=nbin; // DO THIS BETTER
		      timeAve[(kk)*nbin+dedispBin] += (loadAll_pol1[(ii)*nbin*nchan+kk*nbin+jj]+loadAll_pol2[(ii)*nbin*nchan+kk*nbin+jj]);
		      freqAve[(ii)*nbin+dedispBin] += (loadAll_pol1[(ii)*nbin*nchan+kk*nbin+jj]+loadAll_pol2[(ii)*nbin*nchan+kk*nbin+jj]);
		      //		  if (selectSB == 1)
		      //		    printf("Setting: %d %d %d %g\n",ii,jj,kk,timeAve[(kk-kk0)*nbin+dedispBin]);
		      profile[dedispBin] += (loadAll_pol1[(ii)*nbin*nchan+kk*nbin+jj]+loadAll_pol2[(ii)*nbin*nchan+kk*nbin+jj]);	      
		    }
		}
	    }
	}
    }
  printf("Complete forming data sets\n");
}

void displayZapCommand(int *zapFreq,int nchan,int *zapSub,int nsub,float *fChan,char *fname)
{
  int i;
  int haveZapChannel=0;
  int haveZapSub=0;
  int t=0;
  double chanbw = fabs(fChan[1]-fChan[0]);
  FILE *fout;

  if (strlen(fname)>0)
    fout = fopen(fname,"w");
  
  printf("\n");
  printf("paz command with individual subintegrations and channel numbers:\n\n");
  
  printf("paz");
  for (i=0;i<nchan;i++)
    {
      if (zapFreq[i] == 1)
	haveZapChannel=1;
    }
  for (i=0;i<nsub;i++)
    {
      if (zapSub[i] == 1)
	haveZapSub=1;
    }
  
  if (haveZapChannel==1)
    {
      
      printf(" -z \"");
      t=0;
      for (i=0;i<nchan;i++)
	{
	  if (zapFreq[i]==1)
	    {
	      if (strlen(fname)>0)
		fprintf(fout,"ch: %d\n",i);
	      
	      if (t==0){printf("%d",i); t=1;}
	      else {printf(" %d",i);}
	    }
	}
      printf("\"");
    }

    if (haveZapSub==1)
    {
      printf(" -w \"");
      t=0;
      for (i=0;i<nsub;i++)
	{
	  if (zapSub[i]==1)
	    {
	      if (strlen(fname)>0)
		fprintf(fout,"sb: %d\n",i);
	      if (t==0){printf("%d",i); t=1;}
	      else {printf(" %d",i);}
	    }
	}
      printf("\"");
    }

  printf("\n");

  printf("------------\n");
  printf("paz commands with ranges and frequencies instead of channels\n\n");
  printf("paz ");
  if (haveZapChannel==1)
    {
      int s1,s2;
      int pos=0;

      for (i=0;i<nchan;i++)
	{
	  if (zapFreq[i]==1 && pos==0)
	    {
	      if (pos==0){s1 = i; pos = 1;}		
	    }
	 else if (zapFreq[i]!=1 && pos==1)
	    {
	      s2=i-1;
	      if (fChan[s1] < fChan[s2])
		{printf(" -F \"%g %g\"",fChan[s1]-chanbw/2.,fChan[s2]+chanbw/2.); pos=0;}
	      else
		{printf(" -F \"%g %g\"",fChan[s2]-chanbw/2.,fChan[s1]+chanbw/2.); pos=0;}
	    }
	   }
    }
  
  if (haveZapSub==1)
    {
      int s1,s2;
      int pos=0;
      
      for (i=0;i<nsub;i++)
	{
	  if (zapSub[i]==1 && pos==0)
	    {
	      if (pos==0){s1 = i; pos = 1;}		
	    }
	  else if (zapSub[i]!=1 && pos==1)
	    {
	      s2=i-1;
	      if (s2==s1)
		printf(" -w %d",s1);
	      else
		printf(" -W \"%d %d\"",s1,s2);
	      pos=0;
	    }
	}
    }
  printf("\n");

  if (strlen(fname)>0)
    fclose(fout);
  
}

int modOp(int a,int b)
{
  int r=a%b;
  return r<0 ? r+b:r;
}


