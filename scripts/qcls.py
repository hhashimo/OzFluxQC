"""
    OzFlux QC v2.1.2 26 October 2012
"""
#    Version History:
#    <<v1.0: 21 July 2011, code diversion reconciliation>>
#    <<v1.4 30 Sep 2011, final version arrising from OzFlux Black Mountain data workshop>>
#    <<v1.9.9a 8 June 2012, version arrising from conclusion of OzFlux UTS data workshop>>
#    <<v2.0    July 2012, version arrising from merge of PI and JC code after Methven workshop>>
#    <<v2.1.1  September 2012, various PI modifications (see ChangeLog.txt)>>
#    <<v2.1.2 (current) October 2012, changes to metadata, fn2nc for Tumba (see ChangeLog.txt)>>

import sys
import logging
import ast
import constants as c
import copy
import numpy
import qcck
import qcgf
import qcio
import qcts
import qcutils
import time
import xlrd
import meteorologicalfunctions as mf

log = logging.getLogger('qc.ls')

def l2qc(cf,ds1):
    """
        Perform initial QA/QC on flux data
        Generates L2 from L1 data
        * check parameters specified in control file
        
        Functions performed:
            qcck.do_rangecheck*
            qcck.do_CSATcheck
            qcck.do_7500check
            qcck.do_diurnalcheck*
            qcck.do_excludedates*
            qcck.do_excludehours*
            qcts.albedo
        """
    # make a copy of the L1 data
    ds2 = copy.deepcopy(ds1)
    ds2.globalattributes['nc_level'] = 'L2'
    ds2.globalattributes['EPDversion'] = sys.version
    ds2.globalattributes['Functions'] = ''
    # put the control file name into the global attributes
    ds2.globalattributes['controlfile_name'] = cf['controlfile_name']
    # apply the quality control checks (range, diurnal, exclude dates and exclude hours
    qcck.do_qcchecks(cf,ds2)
    # do the CSAT diagnostic check
    qcck.do_CSATcheck(cf,ds2)
    # do the LI-7500 diagnostic check
    qcck.do_7500check(cf,ds2)
    # constrain albedo estimates to full sun angles
    qcts.albedo(cf,ds2)
    log.info(' Finished the albedo constraints')    # apply linear corrections to the data
    log.info(' Applying linear corrections ...')
    qcck.do_linear(cf,ds2)
    # write series statistics to file
    qcio.get_seriesstats(cf,ds2)
    # write the percentage of good data as a variable attribute
    qcutils.get_coverage_individual(ds2)
    
    return ds2

def l3qc(cf,ds2):
    """
        Corrections
        Generates L3 from L2 data
        
        Functions performed:
            qcts.AddMetVars (optional)
            qcts.CorrectSWC (optional*)
            qcck.do_linear (all sites)
            qcutils.GetMergeList + qcts.MergeSeries Ah_EC (optional)x
            qcts.TaFromTv (optional)
            qcutils.GetMergeList + qcts.MergeSeries Ta_EC (optional)x
            qcts.CoordRotation2D (all sites)
            qcts.MassmanApprox (optional*)y
            qcts.Massman (optional*)y
            qcts.CalculateFluxes (used if Massman not optioned)x
            qcts.CalculateFluxesRM (used if Massman optioned)y
            qcts.FhvtoFh (all sites)
            qcts.Fe_WPL (WPL computed on fluxes, as with Campbell algorithm)+x
            qcts.Fc_WPL (WPL computed on fluxes, as with Campbell algorithm)+x
            qcts.Fe_WPLcov (WPL computed on kinematic fluxes (ie, covariances), as with WPL80)+y
            qcts.Fc_WPLcov (WPL computed on kinematic fluxes (ie, covariances), as with WPL80)+y
            qcts.CalculateNetRadiation (optional)
            qcutils.GetMergeList + qcts.MergeSeries Fsd (optional)
            qcutils.GetMergeList + qcts.MergeSeries Fn (optional*)
            qcts.InterpolateOverMissing (optional)
            AverageSeriesByElements (optional)
            qcts.CorrectFgForStorage (all sites)
            qcts.Average3SeriesByElements (optional)
            qcts.CalculateAvailableEnergy (optional)
            qcck.do_qcchecks (all sites)
            qcck.gaps (optional)
            
            *:  requires ancillary measurements for paratmerisation
            +:  each site requires one pair, either Fe_WPL & Fc_WPL (default) or Fe_WPLCov & FcWPLCov
            x:  required together in option set
            y:  required together in option set
        """
    # make a copy of the L2 data
    ds3 = copy.deepcopy(ds2)
    ds3.globalattributes['nc_level'] = 'L3'
    ds3.globalattributes['EPDversion'] = sys.version
    # initialise the global attribute to document the functions used
    ds3.globalattributes['Functions'] = ''
    # put the control file name into the global attributes
    ds3.globalattributes['controlfile_name'] = cf['controlfile_name']
    # check to see if the user wants to gap fill anything
    qcgf.GapFill_L2(cf,ds2,ds3)
    # correct measured soil water content using empirical relationship to collected samples
    #qcts.CorrectSWC(cf,ds3)
    # apply linear corrections to the data
    qcck.do_linear(cf,ds3)
    # merge the HMP and corrected 7500 data
    qcts.MergeSeries(cf,ds3,'Ah',[0,10])
    # get the air temperature from the CSAT virtual temperature
    qcts.TaFromTv(cf,ds3)
    # merge the HMP and corrected CSAT data
    qcts.MergeSeries(cf,ds3,'Ta',[0,10])
    # merge the 7500 CO2 concentration
    qcts.MergeSeries(cf,ds3,'Cc',[0,10])
    # add relevant meteorological values to L3 data
    qcts.CalculateMeteorologicalVariables(ds3)
    # do the 2D coordinate rotation
    qcts.CoordRotation2D(cf,ds3)
    # do the Massman frequency attenuation correction
    qcts.MassmanStandard(cf,ds3)
    # calculate the fluxes
    qcts.CalculateFluxes(cf,ds3)
    # approximate wT from virtual wT using wA (ref: Campbell OPECSystem manual)
    qcts.FhvtoFh(cf,ds3)
    # correct the H2O & CO2 flux due to effects of flux on density measurements
    qcts.Fe_WPL(cf,ds3)
    qcts.Fc_WPL(cf,ds3)
    # convert CO2 units if required
    qcutils.ConvertCO2Units(cf,ds3,Cc='Cc')
    # calculate Fc storage term - single height only at present
    qcts.CalculateFcStorage(cf,ds3)
    # convert Fc and Fc_storage units if required
    qcutils.ConvertFcUnits(cf,ds3,Fc='Fc',Fc_storage='Fc_storage')
    # correct Fc for storage term - only recommended if storage calculated from profile available
    qcts.CorrectFcForStorage(cf,ds3)
    # create a series of filtered ustar
    qcts.FilterUstar(cf,ds3)
    # merge the incoming shortwave radiation
    qcts.MergeSeries(cf,ds3,'Fsd',[0,10])
    # calculate the net radiation from the Kipp and Zonen CNR1
    qcts.CalculateNetRadiation(cf,ds3,'Fn_KZ','Fsd','Fsu','Fld','Flu')
    qcts.MergeSeries(cf,ds3,'Fn',[0,10])
    # combine wind speed from the Wind Sentry and  the CSAT
    qcts.MergeSeries(cf,ds3,'Ws',[0,10])
    # combine wind direction from the Wind Sentry and  the CSAT
    qcts.MergeSeries(cf,ds3,'Wd',[0,10])
    # correct soil heat flux for storage
    #    ... either average the raw ground heat flux, soil temperature and moisture
    #        and then do the correction (OzFlux "standard")
    if qcutils.cfoptionskey(cf,Key='CorrectIndividualFg'):
        #    ... or correct the individual ground heat flux measurements (James' method)
            qcts.CorrectIndividualFgForStorage(cf,ds3)
    else:
        qcts.AverageSeriesByElements(cf,ds3,'Fg')
        qcts.AverageSeriesByElements(cf,ds3,'Ts')
        qcts.AverageSeriesByElements(cf,ds3,'Sws')
        qcts.CorrectFgForStorage(cf,ds3,Fg_out='Fg',Fg_in='Fg',Ts_in='Ts',SWC_in='Sws')
    # calculate the available energy
    qcts.CalculateAvailableEnergy(ds3,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg')
    # create new series using MergeSeries or AverageSeries
    qcck.CreateNewSeries(cf,ds3)
    # create a series of daily averaged soil moisture interpolated back to the time step
    #qcts.DailyAverageSws_Interpolated(cf,ds3,Sws_out='Sws_daily',Sws_in='Sws')
    # re-apply the quality control checks (range, diurnal and rules)
    qcck.do_qcchecks(cf,ds3)
    # coordinate gaps in the three main fluxes
    qcck.CoordinateFluxGaps(cf,ds3)
    # coordinate gaps in Ah_7500_Av with Fc
    qcck.CoordinateAh7500AndFcGaps(cf,ds3)
    # get the statistics for the QC flags and write these to an Excel spreadsheet
    qcio.get_seriesstats(cf,ds3)
    # write the percentage of good data as a variable attribute
    qcutils.get_coverage_individual(ds3)
    # write the percentage of good data for groups
    qcutils.get_coverage_groups(ds3)

    return ds3

def l4qc(cf,ds3):

    # !!! code here to use existing L4 file
    # logic
    # if the L4 doesn't exist
    #  - create ds4 by using copy.deepcopy(ds3)
    # if the L4 does exist and the "UseExistingL4File" option is False
    #  - create ds4 by using copy.deepcopy(ds3)
    # if the L4 does exist and the "UseExistingL4File" option is True
    #  - read the contents of the L4 netCDF file
    #  - check the start and end dates of the L3 and L4 data
    #     - if these are the same then tell the user there is nothing to do
    #  - copy the L3 data to the L4 data structure
    #  - replace the L3 data with the L4 data
    #ds4 = copy.deepcopy(ds3)
    ds4 = qcio.copy_datastructure(cf,ds3)
    # ds4 will be empty (logical false) if an error occurs in copy_datastructure
    # return from this routine if this is the case
    if not ds4: return ds4
    # set some attributes for this level    
    ds4.globalattributes['nc_level'] = 'L4'
    ds4.globalattributes['EPDversion'] = sys.version
    # put the control file name into the global attributes
    ds4.globalattributes['controlfile_name'] = cf['controlfile_name']
    # now do the meteorological driver gap filling
    for ThisOne in cf['Drivers'].keys():
        # interpolate over any gaps up to 1 hour in length
        qcts.InterpolateOverMissing(ds4,series=ThisOne,maxlen=2)
        # gap fill using data from an alternative site
        qcgf.GapFillFromAlternate(cf,ds4,series=ThisOne)
        # gap fill using climatology
        qcgf.GapFillFromClimatology(cf,ds4,series=ThisOne)
        # re-apply the quality control checks (range, diurnal and rules)
        qcck.do_qcchecks_oneseries(cf,ds4,series=ThisOne)
        # interpolate over any remaining gaps up to 3 hours in length
        qcts.InterpolateOverMissing(ds4,series=ThisOne,maxlen=6)
        # fill any remaining gaps climatology
        qcgf.GapFillFromClimatology(cf,ds4,series=ThisOne)
    # re-calculate the meteorological variables
    qcts.CalculateMeteorologicalVariables(ds4)
    # now do the flux gap filling methods
    for ThisOne in cf['Targets'].keys():
        # interpolate over any gaps up to 1 hour in length
        qcts.InterpolateOverMissing(ds4,series=ThisOne,maxlen=2)
        qcgf.GapFillFluxUsingMDS(cf,ds4,series=ThisOne)
        # _namecollector will put a list of series names in ds.soloserieslist for later processing
        qcgf.GapFillUsingSOLO_namecollector(cf,ds4,series=ThisOne)
        qcgf.GapFillFluxFromDayRatio(cf,ds4,series=ThisOne)
        qcgf.GapFillFromClimatology(cf,ds4,series=ThisOne)
    # do the gap filling using SOLO on all series identified by _namecollector
    qcgf.GapFillUsingSOLO(ds3,ds4)
    # now re-apply the quality control checks
    for ThisOne in cf['Drivers'].keys():
        # re-apply the quality control checks (range, diurnal and rules)
        qcck.do_qcchecks_oneseries(cf,ds4,series=ThisOne)
        # interpolate over any remaining gaps up to 3 hours in length
        qcts.InterpolateOverMissing(ds4,series=ThisOne,maxlen=3)
        # fill any remaining gaps climatology
        qcgf.GapFillFromClimatology(cf,ds4,series=ThisOne)
    for ThisOne in cf['Targets'].keys():
        # gap fill the observations using gap fill data in the order specified in cf
        qcts.MergeSeries(cf,ds4,ThisOne,[0,10,20,30,40,50])
        # re-apply the quality control checks (range, diurnal and rules)
        qcck.do_qcchecks_oneseries(cf,ds4,series=ThisOne)
        # interpolate over any remaining gaps up to 3 hours in length
        qcts.InterpolateOverMissing(ds4,series=ThisOne,maxlen=3)
        # fill any remaining gaps climatology
        qcgf.GapFillFromClimatology(cf,ds4,series=ThisOne)
    ## estimate Reco using the methods specified in the control file
    #qcgf.EstimateReco(cf,ds4)
    ## get a single series of NEE
    #qcgf.CalculateNEE(cf,ds4)
    ## ... and partition NEE into GPP and Reco
    #qcgf.PartitionNEE(cf,ds4)
    # write the percentage of good data as a variable attribute
    qcutils.get_coverage_individual(ds4)
    # write the percentage of good data for groups
    qcutils.get_coverage_groups(ds4)

    return ds4