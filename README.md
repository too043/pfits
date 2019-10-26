# pfits

Author:    George.Hobbs@csiro.au  
Copyright: CSIRO 2019

**pfits** is a software package to read, manipulate and process PSRFITS format search- and fold-mode pulsar astronomy data files

## Compile
./bootstrap  
./configure  
make  
make install

## Routines

### pfits_describe

+ summarises the header information in a PSRFITS file
+ The output can be written to the screen or a file

### pfits_fv

+ code that reproduces some of the "fv" (fits viewer) functionality
+ allows the user to obtain detailed information about the file.

### pfits_plot

+ determines whether the data is search mode or fold mode
+ for fold mode:

if 1 channel then it plots the profile (with polarisation information if present)
if multiple channels then it plots of colour scale image for I, Q, U and V (unless just I is requested)

+ for search mode:

the user can select a time range, or nsamp range, or subint range. If > N points (where N is typically 4096) then the data gets averaged and the minimum, mean and maximum value in the average is recorded. The following options are available:

- sum in frequency
- dedisperse in frequency (prior to summing if requested)
- form stokes I
- plot 4-pol data underneath each other
- plot frequency-time plots

### pfits_dedisperse

+ read in a search mode file, dedisperses, frequency-sums (if requested)
+ a new PSRFITS file is output
+ the option to output multiple dispersed data files if requested

### pfits_fold

+ folds a search mode file with a given pulse period or using a tempo2 predictor
+ outputs a new PSRFITS fold-mode file

### pfits_singlepulse

+ extracts single pulses from a search mode file
+ can output PSRFITS fold-mode files for each pulse or calculate various parameters (e.g., flux density) and output those instead

### pfits_makeSpectrum
### pfitsUtil_copyToNew
### pfits_addFRB_1bit   
### pfits_makeTimeMovie
### pfitsUtil_createFoldTemplate
### pfits_addFRB
### pfits_merge
### pfitsUtil_foldmode_combineFreq
### pfits_output
### pfitsUtil_foldmode_combineTime
### pfits_cal
### pfitsUtil_foldmode_extractFreq
### pfits_change
### pfits_outValues
### pfitsUtil_foldmode_extractTime
### pfits_convertText    
### pfitsUtil_foldmode_useTemplate
### pfitsUtil_headerInformation
### pfits_plotfold  
### pfitsUtil_searchmode_combineFreq
### pfits_extract1bit
### pfits_plotsearch           
### pfitsUtil_searchmode_combineTime
### pfits_extractImage
### pfitsUtil_searchmode_extractFreq
### pfits_fftSearch
### pfits_pulseExtraction      
### pfitsUtil_searchMode_extractPart
### pfits_frb
### pfits_read_1bitExtraction
### pfitsUtil_updateMetaData
### pfits_zapProfile
### pfits_getZeroDM                 
### pfits_loader
### pfits_sim
### pfits_zapUWL
### pfits_makeImages
### pfits_simRealData
### pfits_makeMovie
### pfits_statistics

