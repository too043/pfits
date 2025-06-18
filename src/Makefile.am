LIBS = -lerfa -lcalceph

noinst_LTLIBRARIES = libpfits.la

libpfits_la_headers = pfits.h

libpfits_la_SOURCES = pfits.c

bin_PROGRAMS = fitsProcess \
  pfits_fv \
  pfitsUtil_copyToNew \
  pfitsUtil_headerInformation \
  pfitsUtil_foldmode_combineFreq \
  pfitsUtil_foldmode_combineTime \
  pfitsUtil_foldmode_extractFreq \
  pfitsUtil_foldmode_extractTime \
  pfitsUtil_foldmode_useTemplate \
  pfitsUtil_searchmode_combineFreq \
  pfitsUtil_searchmode_combineTime \
  pfitsUtil_searchMode_extractPart \
  pfitsUtil_searchmode_extractFreq \
  pfitsUtil_updateMetaData \
  pfitsUtil_searchmode_filter \
  pfitsUtil_createFoldTemplate

LDADD = $(LIBS) libpfits.la @CFITSIO_LIBS@ @FFTW_LIBS@

if HAVE_PGPLOT

  bin_PROGRAMS += pfits_plot \
    pfits_plot \
    pfits_plotFold \
    pfits_zapProfile \
    pfits_zapUWL \
    pfits_zapUWL_auto \
    pfits_studyRFI \
    pfits_statistics \
    pfits_simRealData \
    pfits_sim \
    pfits_sim1 \
    pfits_searchMode_offScls \
    pfits_read_1bitExtraction \
    pfits_cal \
    pfits_change \
    pfits_dedisperse \
  	pfits_describe \
    pfits_extract1bit \
    pfits_fftSearch \
    pfits_fluxcal \
    pfits_foldStatistics \
    pfits_formBandpass \
    pfits_frb \
    pfits_getZeroDM \
    pfits_imageFold \
    pfits_makeImages \
    pfits_makeMovie \
    pfits_makeSpectrum \
    pfits_makeTimeMovie \
    pfits_merge \
    pfits_output \
    pfits_output1 \
    pfits_outValues \
    pfits_checkRandom

  pfits_plot_SOURCES = pfits_plot.c pfits_setup.c pfits_loader.c
  pfits_plotFold_SOURCES = pfits_plotFold.c pfits_setup.c pfits_loader.c
  pfits_zapProfile_SOURCES = pfits_zapProfile.c pfits_setup.c 
  pfits_zapUWL_SOURCES = pfits_zapUWL.c pfits_setup.c 
  pfits_zapUWL_auto_SOURCES = pfits_zapUWL_auto.c pfits_setup.c 
  pfits_studyRFI_SOURCES = pfits_studyRFI.c pfits_setup.c pfits_loader.c
  pfits_statistics_SOURCES = pfits_statistics.c pfits_setup.c pfits_loader.c
  pfits_simRealData_SOURCES = pfits_simRealData.c pfits_setup.c pfits_loader.c T2toolkit.c
  pfits_sim_SOURCES = pfits_sim.c pfits_setup.c pfits_loader.c T2toolkit.c
  pfits_sim1_SOURCES = pfits_sim1.c pfits_setup.c pfits_loader.c T2toolkit.c
  pfits_searchMode_offScls_SOURCES = pfits_searchMode_offScls.c pfits_setup.c 
  pfits_read_1bitExtraction_SOURCES = pfits_read_1bitExtraction.c pfits_setup.c pfits_loader.c
  pfits_cal_SOURCES = pfits_cal.c pfits_setup.c pfits_loader.c
  pfits_change_SOURCES = pfits_change.c pfits_setup.c
  pfits_fftSearch_SOURCES = pfits_fftSearch.c pfits_setup.c
  pfits_makeImages_SOURCES = pfits_makeImages.c pfits_setup.c pfits_loader.c
  pfits_makeMovie_SOURCES = pfits_makeMovie.c pfits_setup.c pfits_loader.c
  pfits_makeTimeMovie_SOURCES = pfits_makeTimeMovie.c pfits_setup.c pfits_loader.c
  pfits_foldStatistics_SOURCES = pfits_foldStatistics.c pfits_setup.c
  pfits_formBandpass_SOURCES = pfits_formBandpass.c pfits_setup.c
  pfits_dedisperse_SOURCES = pfits_dedisperse.c pfits_setup.c pfits_loader.c
  pfits_describe_SOURCES = pfits_describe.c pfits_setup.c pfits_loader.c
  pfits_extract1bit_SOURCES = pfits_extract1bit.c pfits_setup.c pfits_loader.c
  pfits_fluxcal_SOURCES = pfits_fluxcal.c pfits_setup.c pfits_loader.c
  pfits_frb_SOURCES = pfits_frb.c pfits_setup.c pfits_loader.c
  pfits_getZeroDM_SOURCES = pfits_getZeroDM.c pfits_setup.c pfits_loader.c
  pfits_imageFold_SOURCES = pfits_imageFold.c pfits_setup.c pfits_loader.c
  pfits_makeSpectrum_SOURCES = pfits_makeSpectrum.c pfits_setup.c pfits_loader.c
  pfits_merge_SOURCES = pfits_merge.c pfits_setup.c pfits_loader.c
  pfits_output_SOURCES = pfits_output.c pfits_setup.c pfits_loader.c
  pfits_output1_SOURCES = pfits_output1.c pfits_setup.c pfits_loader.c
  pfits_outValues_SOURCES = pfits_outValues.c pfits_setup.c pfits_loader.c
  pfits_checkRandom_SOURCES = pfits_checkRandom.c pfits_setup.c pfits_loader.c
  
  AM_CFLAGS = @PGPLOT_CFLAGS@
  LDADD += @PGPLOT_LIBS@

endif
