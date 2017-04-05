import sys
import logging
import ast
import constants as c
import copy
import numpy
import os
import qcck
import qcgf
import qcio
import qcrp
import qcts
import qcutils
import time
import xlrd
import meteorologicalfunctions as mf

log = logging.getLogger('qc.ls')

def l1qc(cf):
    # get the data series from the Excel file
    in_filename = qcio.get_infilenamefromcf(cf)
    if not qcutils.file_exists(in_filename,mode="quiet"):
        msg = " Input file "+in_filename+" not found ..."
        log.error(msg)
        ds1 = qcio.DataStructure()
        ds1.returncodes = {"value":1,"message":msg}
        return ds1
    file_name,file_extension = os.path.splitext(in_filename)
    if "csv" in file_extension.lower():
        ds1 = qcio.csv_read_series(cf)
        if ds1.returncodes["value"] != 0:
            return ds1
        # get a series of Excel datetime from the Python datetime objects
        qcutils.get_xldatefromdatetime(ds1)
    else:
        ds1 = qcio.xl_read_series(cf)
        if ds1.returncodes["value"] != 0:
            return ds1
        # get a series of Python datetime objects from the Excel datetime
        qcutils.get_datetimefromxldate(ds1)
    # get the netCDF attributes from the control file
    qcts.do_attributes(cf,ds1)
    # round the Python datetime to the nearest second
    qcutils.round_datetime(ds1,mode="nearest_second")
    #check for gaps in the Python datetime series and fix if present
    fixtimestepmethod = qcutils.get_keyvaluefromcf(cf,["options"],"FixTimeStepMethod",default="round")
    if qcutils.CheckTimeStep(ds1): qcutils.FixTimeStep(ds1,fixtimestepmethod=fixtimestepmethod)
    # recalculate the Excel datetime
    qcutils.get_xldatefromdatetime(ds1)
    # get the Year, Month, Day etc from the Python datetime
    qcutils.get_ymdhmsfromdatetime(ds1)
    # write the processing level to a global attribute
    ds1.globalattributes['nc_level'] = str("L1")
    # get the start and end date from the datetime series unless they were
    # given in the control file
    if 'start_date' not in ds1.globalattributes.keys():
        ds1.globalattributes['start_date'] = str(ds1.series['DateTime']['Data'][0])
    if 'end_date' not in ds1.globalattributes.keys():
        ds1.globalattributes['end_date'] = str(ds1.series['DateTime']['Data'][-1])
    # calculate variances from standard deviations and vice versa
    qcts.CalculateStandardDeviations(cf,ds1)
    # create new variables using user defined functions
    qcts.DoFunctions(cf,ds1)
    # create a series of synthetic downwelling shortwave radiation
    qcts.get_synthetic_fsd(ds1)

    return ds1

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
    # set some attributes for this level    
    qcutils.UpdateGlobalAttributes(cf,ds2,"L2")
    ds2.globalattributes['Functions'] = ''
    # put the control file name into the global attributes
    ds2.globalattributes['controlfile_name'] = cf['controlfile_name']
    # apply the quality control checks (range, diurnal, exclude dates and exclude hours
    qcck.do_qcchecks(cf,ds2)
    # do the CSAT diagnostic check
    qcck.do_CSATcheck(cf,ds2)
    # do the IRGA diagnostic check
    qcck.do_IRGAcheck(cf,ds2)
    # constrain albedo estimates to full sun angles
    #qcts.albedo(cf,ds2)
    #log.info(' Finished the albedo constraints')    # apply linear corrections to the data
    #log.info(' Applying linear corrections ...')
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
    # set some attributes for this level    
    qcutils.UpdateGlobalAttributes(cf,ds3,"L3")
    # initialise the global attribute to document the functions used
    ds3.globalattributes['Functions'] = ''
    # put the control file name into the global attributes
    ds3.globalattributes['controlfile_name'] = cf['controlfile_name']
    # check to see if we have any imports
    qcgf.ImportSeries(cf,ds3)
    # correct measured soil water content using empirical relationship to collected samples
    qcts.CorrectSWC(cf,ds3)
    # apply linear corrections to the data
    qcck.do_linear(cf,ds3)
    # merge whatever humidities are available
    qcts.MergeHumidities(cf,ds3,convert_units=True)
    # get the air temperature from the CSAT virtual temperature
    qcts.TaFromTv(cf,ds3)
    # merge the HMP and corrected CSAT data
    qcts.MergeSeries(cf,ds3,'Ta',[0,10],convert_units=True)
    qcutils.CheckUnits(ds3,"Ta","C",convert_units=True)
    # calculate humidities (absolute, specific and relative) from whatever is available
    qcts.CalculateHumidities(ds3)
    # merge the 7500 CO2 concentration
    qcts.MergeSeries(cf,ds3,'Cc',[0,10],convert_units=True)
    # PRI - disable CO2 units conversion from whatever to mg/m3
    #     - this step is, as far as I can see, redundant, see qcts.Fc_WPL()
    #qcutils.CheckUnits(ds3,"Cc","mg/m3",convert_units=True)
    # add relevant meteorological values to L3 data
    qcts.CalculateMeteorologicalVariables(ds3)
    # check to see if the user wants to use the fluxes in the L2 file
    if not qcutils.cfoptionskeylogical(cf,Key="UseL2Fluxes",default=False):
        # check the covariancve units and change if necessary
        qcts.CheckCovarianceUnits(ds3)
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
    # merge the incoming shortwave radiation
    qcts.MergeSeries(cf,ds3,'Fsd',[0,10])
    # calculate the net radiation from the Kipp and Zonen CNR1
    qcts.CalculateNetRadiation(cf,ds3,Fn_out='Fn_KZ',Fsd_in='Fsd',Fsu_in='Fsu',Fld_in='Fld',Flu_in='Flu')
    qcts.MergeSeries(cf,ds3,'Fn',[0,10])
    # combine wind speed from the Wind Sentry and  the CSAT
    qcts.MergeSeries(cf,ds3,'Ws',[0,10])
    # combine wind direction from the Wind Sentry and  the CSAT
    qcts.MergeSeries(cf,ds3,'Wd',[0,10])
    # correct soil heat flux for storage
    #    ... either average the raw ground heat flux, soil temperature and moisture
    #        and then do the correction (OzFlux "standard")
    qcts.AverageSeriesByElements(cf,ds3,'Ts')
    qcts.AverageSeriesByElements(cf,ds3,'Sws')
    if qcutils.cfoptionskeylogical(cf,Key='CorrectIndividualFg'):
        #    ... or correct the individual ground heat flux measurements (James' method)
            qcts.CorrectIndividualFgForStorage(cf,ds3)
            qcts.AverageSeriesByElements(cf,ds3,'Fg')
    else:
        qcts.AverageSeriesByElements(cf,ds3,'Fg')
        qcts.CorrectFgForStorage(cf,ds3,Fg_out='Fg',Fg_in='Fg',Ts_in='Ts',Sws_in='Sws')
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
    qcutils.UpdateGlobalAttributes(cf,ds4,"L4")
    ds4.cf = cf
    # calculate the available energy
    if "Fa" not in ds4.series.keys():
        qcts.CalculateAvailableEnergy(ds4,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg')
    # create a dictionary to hold the gap filling data
    ds_alt = {}
    # check to see if we have any imports
    qcgf.ImportSeries(cf,ds4)
    # re-apply the quality control checks (range, diurnal and rules)
    qcck.do_qcchecks(cf,ds4)
    # now do the meteorological driver gap filling
    for ThisOne in cf["Drivers"].keys():
        if ThisOne not in ds4.series.keys(): log.error("Series "+ThisOne+" not in data structure"); continue
        # parse the control file for information on how the user wants to do the gap filling
        qcgf.GapFillParseControlFile(cf,ds4,ThisOne,ds_alt)
    # *** start of the section that does the gap filling of the drivers ***
    # fill short gaps using interpolation
    qcgf.GapFillUsingInterpolation(cf,ds4)
    # gap fill using climatology
    qcgf.GapFillFromClimatology(ds4)
    # do the gap filling using the ACCESS output
    qcgf.GapFillFromAlternate(cf,ds4,ds_alt)
    if ds4.returncodes["alternate"]=="quit": return ds4
    # gap fill using SOLO
    qcgf.GapFillUsingSOLO(cf,ds3,ds4)
    if ds4.returncodes["solo"]=="quit": return ds4
    # merge the first group of gap filled drivers into a single series
    qcts.MergeSeriesUsingDict(ds4,merge_order="prerequisite")
    # re-calculate the ground heat flux but only if requested in control file
    opt = qcutils.get_keyvaluefromcf(cf,["Options"],"CorrectFgForStorage",default="No",mode="quiet")
    if opt.lower()!="no":
        qcts.CorrectFgForStorage(cf,ds4,Fg_out='Fg',Fg_in='Fg_Av',Ts_in='Ts',Sws_in='Sws')
    # re-calculate the net radiation
    qcts.CalculateNetRadiation(cf,ds4,Fn_out='Fn',Fsd_in='Fsd',Fsu_in='Fsu',Fld_in='Fld',Flu_in='Flu')
    # re-calculate the available energy
    qcts.CalculateAvailableEnergy(ds4,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg')
    # merge the second group of gap filled drivers into a single series
    qcts.MergeSeriesUsingDict(ds4,merge_order="standard")
    # re-calculate the water vapour concentrations
    qcts.CalculateHumiditiesAfterGapFill(ds4)
    # re-calculate the meteorological variables
    qcts.CalculateMeteorologicalVariables(ds4)
    # the Tumba rhumba
    qcts.CalculateComponentsFromWsWd(ds4)
    # check for any missing data
    qcutils.get_missingingapfilledseries(ds4)
    # write the percentage of good data as a variable attribute
    qcutils.get_coverage_individual(ds4)
    # write the percentage of good data for groups
    qcutils.get_coverage_groups(ds4)

    return ds4

def l5qc(cf,ds4):
    ds5 = qcio.copy_datastructure(cf,ds4)
    # ds4 will be empty (logical false) if an error occurs in copy_datastructure
    # return from this routine if this is the case
    if not ds5: return ds5
    # set some attributes for this level
    qcutils.UpdateGlobalAttributes(cf,ds5,"L5")
    ds5.cf = cf
    # create a dictionary to hold the gap filling data
    ds_alt = {}
    # check to see if we have any imports
    qcgf.ImportSeries(cf,ds5)
    # re-apply the quality control checks (range, diurnal and rules)
    qcck.do_qcchecks(cf,ds5)
    # now do the flux gap filling methods
    label_list = qcutils.get_label_list_from_cf(cf)
    for ThisOne in label_list:
        # parse the control file for information on how the user wants to do the gap filling
        qcgf.GapFillParseControlFile(cf,ds5,ThisOne,ds_alt)
    # *** start of the section that does the gap filling of the fluxes ***
    # apply the turbulence filter (if requested)
    qcck.ApplyTurbulenceFilter(cf,ds5)
    # fill short gaps using interpolation
    #qcgf.GapFillUsingInterpolation(cf,ds5)
    # do the gap filling using SOLO
    qcgf.GapFillUsingSOLO(cf,ds4,ds5)
    if ds5.returncodes["solo"]=="quit": return ds5
    ## gap fill using marginal distribution sampling
    #qcgf.GapFillFluxUsingMDS(cf,ds5)
    ## gap fill using ratios
    #qcgf.GapFillFluxFromDayRatio(cf,ds5)
    # gap fill using climatology
    qcgf.GapFillFromClimatology(ds5)
    # merge the gap filled drivers into a single series
    qcts.MergeSeriesUsingDict(ds5,merge_order="standard")
    # write the percentage of good data as a variable attribute
    qcutils.get_coverage_individual(ds5)
    # write the percentage of good data for groups
    qcutils.get_coverage_groups(ds5)

    return ds5

def l6qc(cf,ds5):
    ds6 = qcio.copy_datastructure(cf,ds5)
    # ds6 will be empty (logical false) if an error occurs in copy_datastructure
    # return from this routine if this is the case
    if not ds6: return ds6
    # set some attributes for this level
    qcutils.UpdateGlobalAttributes(cf,ds6,"L6")
    # parse the control file
    qcrp.ParseL6ControlFile(cf,ds6)
    # check to see if we have any imports
    qcgf.ImportSeries(cf,ds6)
    # check units
    qcutils.CheckUnits(ds6,"Fc","umol/m2/s",convert_units=True)
    ## filter Fc for night time and ustar threshold, write to ds as "ER"
    #result = qcrp.GetERFromFc(cf,ds6)
    #if result==0: return
    # apply the turbulence filter (if requested)
    qcck.ApplyTurbulenceFilter(cf,ds6)
    qcrp.GetERFromFc2(cf,ds6)
    # estimate ER using SOLO
    qcrp.ERUsingSOLO(cf,ds6)
    # estimate ER using FFNET
    qcrp.ERUsingFFNET(cf,ds6)
    # estimate ER using Lloyd-Taylor
    qcrp.ERUsingLloydTaylor(cf,ds6)
    # estimate ER using Lasslop et al
    qcrp.ERUsingLasslop(cf,ds6)
    # merge the estimates of ER with the observations
    qcts.MergeSeriesUsingDict(ds6,merge_order="standard")
    # calculate NEE from Fc and ER
    qcrp.CalculateNEE(cf,ds6)
    # calculate NEP from NEE
    qcrp.CalculateNEP(cf,ds6)
    # calculate ET from Fe
    qcrp.CalculateET(ds6)
    # partition NEE into GPP and ER
    qcrp.PartitionNEE(cf,ds6)
    # write the percentage of good data as a variable attribute
    qcutils.get_coverage_individual(ds6)
    # write the percentage of good data for groups
    qcutils.get_coverage_groups(ds6)
    # do the L6 summary
    qcrp.L6_summary(cf,ds6)

    return ds6
