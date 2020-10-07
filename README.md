# pfits
---
Author:    George.Hobbs@csiro.au  
Copyright: CSIRO 2020 

**pfits** is a software package to read, manipulate and process PSRFITS format search- and fold-mode pulsar astronomy data files 

## Dependencies
cfitsio
pgplot

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

### fitsProcess
### pfits_addFRB
### pfits_addFRB_1bit
### pfits_cal
### pfits_change
### pfits_convertText
### pfits_extract1bit
### pfits_extractImage
### pfits_fftSearch
### pfits_fluxcal
### pfits_foldStatistics
### pfits_formBandpass
### pfits_frb
### pfits_getZeroDM
### pfits_imageFold
### pfits_makeImages
### pfits_makeMovie
### pfits_makeSpectrum
### pfits_makeTimeMovie
### pfits_merge
### pfits_output
### pfits_output1
### pfits_outValues
### pfits_plotfold
### pfits_plotsearch
### pfits_read_1bitExtraction
### pfits_searchMode_offScls
### pfits_sim
### pfits_sim1
### pfits_simRealData
### pfits_statistics
### pfits_studyRFI
### pfitsUtil_copyToNew
### pfitsUtil_createFoldTemplate
### pfitsUtil_foldmode_combineFreq
### pfitsUtil_foldmode_combineTime
### pfitsUtil_foldmode_extractFreq
### pfitsUtil_foldmode_extractTime
### pfitsUtil_foldmode_useTemplate
### pfitsUtil_headerInformation
### pfitsUtil_searchmode_combineFreq
### pfitsUtil_searchmode_combineTime
### pfitsUtil_searchmode_extractFreq
### pfitsUtil_searchMode_extractPart
### pfitsUtil_searchmode_filter
### pfitsUtil_updateMetaData
### pfits_zapProfile
### pfits_zapUWL
### pfits_zapUWL_auto
