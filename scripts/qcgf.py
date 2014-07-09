import ast
from calendar import isleap
from configobj import ConfigObj
import constants as c
import csv
import datetime
import logging
import numpy
import matplotlib.dates as mdt
import matplotlib.pyplot as plt
import os
import qcio
import qcts
import qcutils
from scipy import interpolate
import shutil
import statsmodels.api as sm
import subprocess
import sys
import Tkinter
import xlrd

log = logging.getLogger('qc.gf')

def EstimateReco(cf,ds):
    # check to see if there is an Reco section in the control file
    section = qcutils.get_cfsection(cf,series='Reco',mode='quiet')
    if len(section)==0: return
    # step through the methods given in the control file
    method_list = cf[section]['Reco'].keys()
    for ThisMethod in method_list:
        if ThisMethod == 'LloydTaylor': EstimateRecoUsingLloydTaylor(cf,ds)
        if ThisMethod == 'SOLO': EstimateRecoUsingSOLO(cf,ds)

def EstimateRecoUsingLloydTaylor(cf,ds):
    # need to check that all information required is in the control file
    section = qcutils.get_cfsection(cf,series='Reco',mode='quiet')
    if len(section)==0:
        log.error('Reco section not found in control file')
        return
    if 'LloydTaylor' not in cf[section]['Reco']:
        log.error('LloydTaylor section not found in control file')
        return
    # get the driver
    if 'drivers' not in cf[section]['Reco']['LloydTaylor']:
        log.error('drivers key not found in LloydTaylor section of control file')
        return
    else:
        driver_list = eval(cf[section]['Reco']['LloydTaylor']['drivers'])
        driver_label = driver_list[0]
    # get the monthly values for the activation energy, E0
    if 'E0' not in cf[section]['Reco']['LloydTaylor']:
        log.error('E0 key not found in LloydTaylor section of control file')
        return
    else:
        E0_monthly = numpy.array(eval(cf[section]['Reco']['LloydTaylor']['E0']))
    # get the monthly values for the base respiration, rb
    if 'rb' not in cf[section]['Reco']['LloydTaylor']:
        log.error('rb key not found in LloydTaylor section of control file')
        return
    else:
        rb_monthly = numpy.array(eval(cf[section]['Reco']['LloydTaylor']['rb']))
    # get the output label
    if 'output' not in cf[section]['Reco']['LloydTaylor']:
        log.error('output key not found in LloydTaylor section of control file')
        return
    else:
        out_label = cf[section]['Reco']['LloydTaylor']['output']
    # ... and make an array of values for each month
    nRecs = int(ds.globalattributes['nc_nrecs'])
    E0 = numpy.ma.ones(nRecs)
    rb = numpy.ma.zeros(nRecs)
    month = ds.series['Month']['Data']
    lwr = numpy.max(numpy.array([numpy.min(month),1]))
    upr = numpy.min(numpy.array([numpy.max(month),12]))
    for m in range(lwr,upr+1):
        index = numpy.where(month==m)[0]
        E0[index] = E0_monthly[m-1]
        rb[index] = rb_monthly[m-1]
    # get the driver data
    Ts,flag = qcutils.GetSeriesasMA(ds,driver_label)
    # estimate Reco using the Lloyd-Taylor expression
    t1 = 1/(c.Tref-c.T0)
    t2 = 1/(Ts-c.T0)
    Reco = rb*numpy.exp(E0*(t1-t2))
    # put the estimated Reco into the data structure
    units=qcutils.GetUnitsFromds(ds, 'Fc')
    attr = qcutils.MakeAttributeDictionary(long_name='Reco estimated using Lloyd-Taylor',units=units)
    qcutils.CreateSeries(ds,out_label,Reco,Flag=flag,Attr=attr)
    
def EstimateRecoUsingSOLO(cf,ds):
    log.info('Estimating Reco using SOLO is not implemented yet')
    pass

def CalculateNEE(cf,ds):
    # now get a single series of NEE using the following rules:
    #  - when Fsd > 10 W/m2 (day time)
    #    - use Fc gap filled by chosen method
    #  - when Fsd < 10 W/m2 (night time)
    #    - use observed Fc (Fc_flag=0) when u* > ustar_threshold
    #    - use modelled Reco when Fc_flag!=0 or u* < ustar_threshold
    Fsd,Fsd_flag = qcutils.GetSeriesasMA(ds,'Fsd')
    ustar,ustar_flag = qcutils.GetSeriesasMA(ds,'ustar')
    Fc_label = str(cf['Derived']['NEE']['Fc'])
    Fc,Fc_flag = qcutils.GetSeriesasMA(ds,Fc_label)
    Reco_label = str(cf['Derived']['NEE']['Reco'])
    Reco,Reco_flag = qcutils.GetSeriesasMA(ds,Reco_label)
    Fsd_threshold = float(cf['Params']['Fsd_threshold'])
    ustar_threshold = float(cf['Params']['ustar_threshold'])
    # create a series for NEE and the NEE QC flag
    nRecs = int(ds.globalattributes['nc_nrecs'])
    NEE = numpy.ma.array([-9999]*nRecs,dtype=numpy.float64)
    NEE_flag = numpy.ma.array([0]*nRecs,dtype=numpy.int32)
    # fill the day time NEE
    index = numpy.ma.where(Fsd>Fsd_threshold)[0]
    NEE[index] = Fc[index]
    # fill the night time NEE when ustar is above the threshold
    index = numpy.ma.where((Fsd<=Fsd_threshold)&(ustar>ustar_threshold))[0]
    NEE[index] = Fc[index]
    # fill the night time NEE when ustar is below the threshold
    index = numpy.ma.where((Fsd<=Fsd_threshold)&(ustar<=ustar_threshold))[0]
    NEE[index] = Reco[index]
    # create the NEE series in the data structure
    units=qcutils.GetUnitsFromds(ds,Fc_label)
    attr = qcutils.MakeAttributeDictionary(long_name='Net ecosystem exchange',units=units)
    qcutils.CreateSeries(ds,'NEE',NEE,Flag=NEE_flag,Attr=attr)

def PartitionNEE(cf,ds):
    # check that there is a GPP section in the control file
    section = qcutils.get_cfsection(cf,series='GPP',mode='quiet')
    if len(section)==0:
        log.info('PartitionNEE: GPP section not found in control file')
        return
    # now check that the GPP section in the control file contains the
    # information needed.
    if 'NEE' not in cf[section]['GPP']:
        log.info('PartitionNEE: NEE key not in GPP section of control file')
        return
    else:
        NEE_label = str(cf[section]['GPP']['NEE'])
    if 'Reco' not in cf[section]['GPP']:
        log.info('PartitionNEE: Reco key not in GPP section of control file')
        return
    else:
        Reco_label = str(cf[section]['GPP']['Reco'])
    # now check that the requested series are in the data structure
    if NEE_label not in ds.series.keys():
        log.info('PartitionNEE: requested series '+NEE_label+' is not in the data structure')
        return
    else:
        NEE,f = qcutils.GetSeriesasMA(ds,NEE_label)
    if Reco_label not in ds.series.keys():
        log.info('PartitionNEE: requested series '+Reco_label+' is not in the data structure')
        return
    else:
        Reco,f = qcutils.GetSeriesasMA(ds,Reco_label)
    # check the units
    NEE_units=qcutils.GetUnitsFromds(ds,NEE_label)
    Reco_units=qcutils.GetUnitsFromds(ds,Reco_label)
    if NEE_units!=Reco_units:
        log.error('Units of NEE ('+NEE_label+') and Reco ('+Reco_label+') are different')
        return
    # ... at last, we can do the partitioning
    GEP = NEE - Reco
    # create the GEP QC flag
    nRecs = int(ds.globalattributes['nc_nrecs'])
    GEP_flag = numpy.ma.array([0]*nRecs,dtype=numpy.int32)
    # now check for non-zero values of GPP at night
    attr = qcutils.MakeAttributeDictionary(long_name='Gross ecosystem productivity',units=NEE_units)
    qcutils.CreateSeries(ds,'GEP',GEP,Flag=GEP_flag,Attr=attr)

def GapFillFromAlternate(cf,ds,series=''):
    """
        Gap fill using data from alternate sites specified in the control file
        """
    section = qcutils.get_cfsection(cf,series=series,mode='quiet')
    if len(section)==0: return
    if 'GapFillFromAlternate' not in cf[section][series].keys(): return
    log.info(' Gap filling '+series+' using data from alternate sites')
    ts = ds.globalattributes['time_step']
    # Gap fill using data from alternate sites specified in the control file
    ds_alt = {}               # create a dictionary for the data from alternate sites
    open_ncfiles = []         # create an empty list of open netCDF files
    # loop over the entries in the GapFillFromAlternate section
    for Alt in cf[section][series]['GapFillFromAlternate'].keys():
        log.info(' Gap filling '+ThisOne+' by replacing with alternate site data')
        # get the file name for the alternate site
        alt_filename = cf[section][series]['GapFillFromAlternate'][Alt]['FileName']
        # get the variable name for the alternate site data if specified, otherwise use the same name
        if 'AltVarName' in cf[section][series]['GapFillFromAlternate'][Alt].keys():
            alt_varname = cf[section][series]['GapFillFromAlternate'][Alt]['AltVarName']
        else:
            alt_varname = series
        # check to see if the alternate site file is already open
        if alt_filename not in open_ncfiles:
            # open and read the file if it is not already open
            n = len(open_ncfiles)
            open_ncfiles.append(alt_filename)
            ds_alt[n] = qcio.nc_read_series_file(alt_filename)
        else:
            # get the file index number if it is already open
            n = open_ncfiles.index(alt_filename)
        # check to see if alternate site data needs transform
        if 'Transform' in cf[section][series]['GapFillFromAlternate'][Alt].keys():
            # get the datetime series for the alternate site
            AltDateTime = ds_alt[n].series['DateTime']['Data']
            # get the data for the alternate site
            AltSeriesData = ds_alt[n].series[alt_varname]['Data']
            # evaluate the list of start dates, end dates and transform coefficients
            TList = ast.literal_eval(cf[section][section]['GapFillFromAlternate'][Alt]['Transform'])
            # loop over the datetime ranges for the transform
            for TListEntry in TList:
                qcts.TransformAlternate(TListEntry,AltDateTime,AltSeriesData,ts=ts)
        qcts.ReplaceWhereMissing(ds.series[series],ds.series[series],ds_alt[n].series[alt_varname],FlagOffset=100)
    if 'GapFillFromAlternate' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', GapFillFromAlternate'

def GapFill_L2(cf,ds2,ds3):
    for series in ds3.series.keys():
        section = qcutils.get_cfsection(cf,series=series,mode='quiet')
        if len(section)==0: continue
        for gftype in cf[section][series].keys():
            if gftype=="GapFillFromAlternate":
                log.error("GapFill_L2: Gap filling from alternate data source not implemented yet")
                pass
            if gftype=="GapFillFromClimatology":
                gfClimatology_oneseries(cf,ds3,series)
            if gftype=="GapFillUsingSOLO":
                GapFillUsingSOLO_namecollector(cf,ds3,series=series)
    if len(ds3.soloserieslist)!=0:
        GapFillUsingSOLO(ds2,ds3)
        for series in ds3.soloserieslist:
            qcts.MergeSeries(cf,ds3,series,[0,10,20,30,40,50])

def GapFillFromClimatology(ds):
    '''
    Gap fill missing data using data from the climatology spreadsheet produced by
    the climatology.py script.
    '''
    if "climatology" not in dir(ds): return
    # get the list of climatology data files and take out the duplicates
    climatology_file_list = list(set(ds.climatology["file_name"]))
    # more than 1 climatology data file not supported at present
    if len(climatology_file_list)>1: raise NotImplementedError
    # get the name of the climatology file
    cli_filename = climatology_file_list[0]
    # return to calling routine if the file doesn't exist
    if not os.path.exists(cli_filename):
        log.error(" GapFillFromClimatology: Climatology file "+cli_filename+" doesn't exist")
        return
    # open the climatology file
    cli_xlbook = xlrd.open_workbook(cli_filename)
    # loop over the series to gap filled using climatology
    for label,output,method in zip(ds.climatology["series"],ds.climatology["output"],ds.climatology["method"]):
        # check to see if there are any gaps in "series"
        index = numpy.where(abs(ds.series[label]['Data']-float(-9999))<c.eps)[0]
        if len(index)==0: continue                      # no gaps found in "series"
        # do the gap filling
        log.info(" Gap filling "+label+" using climatology")
        # choose the gap filling method
        if method=="monthly":
            gfClimatology_monthly(ds,label,output,cli_xlbook)
        elif method=="interpolated daily":
            gfClimatology_interpolateddaily(ds,label,output,cli_xlbook)
        else:
            log.error(" GapFillFromClimatology: unrecognised method option for "+series)
            continue
    if 'GapFillFromClimatology' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', GapFillFromClimatology'
    # remove the "climatology" attribute from ds
    del ds.climatology

def gfClimatology_monthly(ds,series,output,xlbook):
    thissheet = xlbook.sheet_by_name(series)
    val1d = numpy.zeros_like(ds.series[series]['Data'])
    values = numpy.zeros([48,12])
    for month in range(1,13):
        xlCol = (month-1)*5 + 2
        values[:,month-1] = thissheet.col_values(xlCol)[2:50]
    for i in range(len(ds.series[series]['Data'])):
        h = numpy.int(2*ds.series['Hdh']['Data'][i])
        m = numpy.int(ds.series['Month']['Data'][i])
        val1d[i] = values[h,m-1]
    ds.series[output]['Data'][index] = val1d[index]
    ds.series[output]['Flag'][index] = numpy.int32(40)
    
def gfClimatology_interpolateddaily(ds,series,output,xlbook):
    # gap fill from interpolated 30 minute data
    ts = ds.globalattributes["time_step"]
    ldt = ds.series["DateTime"]["Data"]
    thissheet = xlbook.sheet_by_name(series+'i(day)')
    datemode = xlbook.datemode
    basedate = datetime.datetime(1899, 12, 30)
    nts = thissheet.ncols - 1
    ndays = thissheet.nrows - 2
    # read the time stamp values from the climatology worksheet
    tsteps = thissheet.row_values(1,start_colx=1,end_colx=nts+1)
    # read the data from the climatology workbook
    val1d = numpy.ma.zeros(ndays*nts,dtype=numpy.float64)
    # initialise an array for the datetime of the climatological values
    cdt = [None]*nts*ndays
    # loop over the rows (days) of data
    for xlRow in range(ndays):
        # get the Excel datetime value
        xldatenumber = thissheet.cell_value(xlRow+2,0)
        # convert this to a Python Datetime
        xldatetime = basedate+datetime.timedelta(days=xldatenumber+1462*datemode)
        # fill the climatology datetime array
        cdt[xlRow*nts:(xlRow+1)*nts] = [xldatetime+datetime.timedelta(hours=hh) for hh in tsteps]
        # fill the climatological value array
        val1d[xlRow*nts:(xlRow+1)*nts] = thissheet.row_values(xlRow+2,start_colx=1,end_colx=nts+1)
    # get the data to be filled with climatological values
    data,flag = qcutils.GetSeriesasMA(ds,series)
    # get an index of missing values
    idx = numpy.ma.where(data.mask==True)[0]
    # there must be a better way to do this ...
    # simply using the index (idx) to set a slice of the data array to the gap filled values in val1d
    # does not seem to work (mask stays true on replaced values in data), the work around is to
    # step through the indices, find the time of the missing value in data, find the same time in the
    # gap filled values val1d and set the missing element of data to this element of val1d
    li = 0
    for ii in idx:
        jj = cdt.index(ldt[ii],li)
        li = jj
        data[ii] = val1d[jj]
        flag[ii] = numpy.int32(40)
    # get the attribute dictionary of the series being gap filled
    attr = qcutils.GetAttributeDictionary(ds,series)
    # put the gap filled data back into the data structure
    qcutils.CreateSeries(ds,output,data,Flag=flag,Attr=attr)

def GapFillFluxFromDayRatio(cf,ds,series=''):
    section = qcutils.get_cfsection(cf,series=series,mode='quiet')
    if len(section)==0: return
    if 'GapFillFluxFromDayRatio' not in cf[section][series].keys(): return
    ndays = 365
    if isleap(ds.series['Year']['Data'][0]): ndays = 366
    dt = int(ds.globalattributes['time_step'])    # time step in minutes
    nts = int(24/(float(dt)/float(60)))           # number of time steps in a day
    nmn = int(12)                                 # number of months in year
    # if no series list passed in then create one
    log.info(' Gap filling '+series+' using daily ratio (day) and climatology (night)')
    # get the details from the control file
    alt_filename = cf[section][series]['GapFillFluxFromDayRatio']['file_name']
    if not os.path.exists(alt_filename):
        log.error(" GapFillFromDayRatio: Climatology file "+alt_filename+" doesn't exist")
        return
    xlbook = xlrd.open_workbook(alt_filename)
    ratio_label = cf[section][series]['GapFillFluxFromDayRatio']['ratio_xlSheet']
    ratio_xlsheet = xlbook.sheet_by_name(ratio_label)
    driver_label = cf[section][series]['GapFillFluxFromDayRatio']['drivers']
    flux_xlsheet = xlbook.sheet_by_name(cf[section][series]['GapFillFluxFromDayRatio']['flux_xlSheet'])
    out_label = cf[section][series]['GapFillFluxFromDayRatio']['output']
    # get the data series as masked arrays from the data structure
    Fsd, Fsd_flag = qcutils.GetSeriesasMA(ds,'Fsd')
    us, us_flag = qcutils.GetSeriesasMA(ds,'ustar')
    driver, driver_flag = qcutils.GetSeriesasMA(ds, driver_label)
    nRecs = len(driver)
    # get the flux series to be gap filled
    flux, flux_flag = qcutils.GetSeriesasMA(ds,series)
    # initialise arrays for the ratios and the fluxes
    ratio_ts = numpy.zeros(ndays*nts,dtype=numpy.float64)
    flux_monthly = numpy.zeros((nts,nmn),dtype=numpy.float64)
    # read the ratio from the climatology workbook
    for xlRow in range(ndays):
        ratio_ts[xlRow*nts:(xlRow+1)*nts] = ratio_xlsheet.row_values(xlRow+2,start_colx=1,end_colx=nts+1)
    # add an extra constraint (WUE must be <0 during the day) if we are doing Fc
    # PRI - relocated from original position in code (just before calls to CreateSeries
    #       and eventually deprecated as it may be dangerous
    #if 'Fc' in ThisOne: index = numpy.ma.where((Fsd>50)&(ratio_ts<0))
    # read the monthly flux data from the climatology workbook
    for xlCol in range(nmn):
        flux_monthly[:,xlCol] = flux_xlsheet.col_values(xlCol+1)[2:nts+2]
    # now we interpolate the climatological fluxes from monthly to a daily time step
    flux_monthly_tiled = numpy.tile(flux_monthly,(3,3))                              # tile the climatological flux into a 3 by 3 matrix
    nx = numpy.shape(flux_monthly_tiled)[1]; ny = numpy.shape(flux_monthly_tiled)[0] # get the dimensions of the original flux array
    nxi = ndays*3; nyi = numpy.shape(flux_monthly_tiled)[0]                          # get the dimensions of the daily array (still tiled)
    x = numpy.linspace(1,nx,nx); y = numpy.linspace(1,ny,ny)                         # 1d array of x and y dimension indices, tiled original
    xi = numpy.linspace(1,nx,nxi); yi = numpy.linspace(1,ny,nyi)                     # 1d array of x and y dimension indices, tiled interpolated
    # interpolate from monthly to daily time step
    nk = interpolate.RectBivariateSpline(y,x,flux_monthly_tiled,kx=2,ky=2)           # quadratic spline
    flux_daily_tiled = nk(yi,xi)                                                     # do the interpolation
    flux_daily = flux_daily_tiled[nyi/3:2*nyi/3,nxi/3:2*nxi/3]                       # pick out the central tile (still a 2D array)
    flux_ic = numpy.ravel(flux_daily,order='F')                                      # convert 2D array to 1D array, same length as data series, correct order
    # calculate the flux from the ratio (EF, BR or WUE) and the driver (Fa, Fe or Fe)
    flux_rd = ratio_ts * driver
    # create the gap filling series and initialise it to the interpolated climatological fluxes
    flux_gf = flux_ic.copy()
    # now figure out when its daytime and night time
    index = numpy.ma.where(Fsd>50)[0]
    # replace the interpolated climatological fluxes with the daytime values
    # calculated from the ratio (EF, BR or WUE) and the driver (Fa, Fe or Fe)
    flux_gf[index] = flux_rd[index]
    # generate a QC flag to indicate data has been estimated from ratios (day time)
    # and interpolated climatology (night time).
    flag = numpy.array([60]*nRecs,dtype=numpy.int32)
    # put the gap filled data into the data structure
    units=qcutils.GetUnitsFromds(ds, series)
    attr = qcutils.MakeAttributeDictionary(long_name='gap filled using ratio (daily)',units=units)
    qcutils.CreateSeries(ds,out_label,flux_gf,Flag=flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='interpolated '+ratio_label,units='None')
    qcutils.CreateSeries(ds,ratio_label,ratio_ts,Flag=flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='ratio times driver',units=units)
    qcutils.CreateSeries(ds,series+'_rd',flux_rd,Flag=flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='interpolated climatological fluxes',units=units)
    qcutils.CreateSeries(ds,series+'_ic',flux_ic,Flag=flag,Attr=attr)
    if 'GapFillFluxFromDayRatio' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', GapFillFluxFromDayRatio'

def GapFillFluxUsingMDS(cf,ds,series=""):
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    if len(section)==0: return
    if "GapFillFluxUsingMDS" in cf[section][series].keys():
        log.info(" GapFillFluxUsingMDS: not implemented yet")
        return

def GapFillParseControlFile(cf,ds,series=""):
    # find the section containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # return empty handed if the series is not in a section
    if len(section)==0: return
    if series not in ds.series.keys(): return
    if "GapFillFromACCESS" in cf[section][series].keys():
        # create the ACCESS directory in the data structure
        if "access" not in dir(ds): ds.access = {}
        # create a dictionary for this series
        ds.access[series] = {}
        ds.access[series]["file_name"] = cf[section][series]["GapFillFromACCESS"]["file_name"]
        ds.access[series]["output"] = cf[section][series]["GapFillFromACCESS"]["output"]
        ds.access[series]["var@maxr"] = ""
        ds.access[series]["results"] = {"startdate":[],"enddate":[],"n":[],"r_max":[],"bias":[],
                                        "rmse":[],"var_tow":[],"var_acc":[],"max_lag":[],
                                        "m_rlm":[],"b_rlm":[],"m_ols":[],"b_ols":[]}
        data,flag = qcutils.MakeEmptySeries(ds,ds.access[series]["output"])
        attr = qcutils.GetAttributeDictionary(ds,ds.access[series]["output"])
        qcutils.CreateSeries(ds,ds.access[series]["output"],data,Flag=flag,Attr=attr)
    if "GapFillUsingSOLO" in cf[section][series].keys():
        # create the SOLO directory in the data structure
        if "solo" not in dir(ds): ds.solo = {}
        # create a dictionary for this series
        ds.solo[series] = {}
        ds.solo[series]["drivers"] = ast.literal_eval(cf[section][series]["GapFillUsingSOLO"]["drivers"])
        ds.solo[series]["output"] = cf[section][series]["GapFillUsingSOLO"]["output"]
        ds.solo[series]["results"] = {"startdate":[],"enddate":[],"n":[],"r_max":[],"bias":[],
                                      "rmse":[],"var_obs":[],"var_mod":[],"m_ols":[],"b_ols":[]}
        data,flag = qcutils.MakeEmptySeries(ds,ds.solo[series]["output"])
        attr = qcutils.GetAttributeDictionary(ds,ds.solo[series]["output"])
        qcutils.CreateSeries(ds,ds.solo[series]["output"],data,Flag=flag,Attr=attr)
    if "GapFillFromClimatology" in cf[section][series].keys():
        if "climatology" not in dir(ds):
            ds.climatology = {}
            ds.climatology["series"] = []
            ds.climatology["output"] = []
            ds.climatology["file_name"] = []
            ds.climatology["method"] = []
        ds.climatology["series"].append(series)
        file_name = cf[section][series]["GapFillFromClimatology"]["file_name"]
        ds.climatology["file_name"].append(file_name)
        output = cf[section][series]["GapFillFromClimatology"]["output"]
        ds.climatology["output"].append(output)
        data,flag = qcutils.MakeEmptySeries(ds,output)
        attr = qcutils.GetAttributeDictionary(ds,output)
        qcutils.CreateSeries(ds,output,data,Flag=flag,Attr=attr)
        if "method" not in cf[section][series]["GapFillFromClimatology"].keys():
            # default if "method" missing is "interpolated_daily"
            ds.climatology["method"].append("interpolated_daily")
        else:
            ds.climatology["method"].append(cf[section][series]["GapFillFromClimatology"]["method"])
    if "MergeSeries" in cf[section][series].keys():
        if "merge" not in dir(ds):
            ds.merge = {}
            ds.merge["series"] = []
            ds.merge["source"] = []
        ds.merge["series"].append(series)
        source = ast.literal_eval(cf[section][series]["MergeSeries"]["Source"])
        ds.merge["source"].append(source)

def GapFillFromACCESS(ds4):
    '''
    This is the gap fill from ACCESS GUI.
    The ACCESS GUI is displayed separately from the main OzFluxQC GUI.
    It consists of text to display the start and end datetime of the file,
    two entry boxes for the start and end datetimes of the ACCESS gap fill and
    a button to insert the gap fill data ("Use ACCESS") and a button to exit
    the ACCESS GUI when we are done.  On exit, the OzFluxQC main GUI continues
    and eventually writes the gap filled data to file.
    '''
    if "access" not in dir(ds4): return
    # get the list of ACCESS data files and take out the duplicates
    access_file_list = list(set([ds4.access[x]["file_name"] for x in ds4.access.keys()]))
    #access_file_list = list(set(ds4.access["file_name"]))
    # more than 1 ACCESS data file not supported at present
    if len(access_file_list)>1: raise NotImplementedError
    ldt_tower = ds4.series["DateTime"]["Data"]
    # we need the start and end time of the period of overlap between the ACCESS and tower data
    # read the ACCESS data file
    ds_access = qcio.nc_read_series(access_file_list[0])
    ldt_access = ds_access.series["DateTime"]["Data"]
    # check that the ACCESS and tower data overlap
    if ldt_access[0]>ldt_tower[-1]:
        log.error("ACCESS data starts after tower data finishes")
        return
    if ldt_access[-1]<ldt_tower[0]:
        log.error("ACCESS data end before tower data starts")
        return
    # store the start and end date of the overlap period
    access_dates = {"overlap_startdate":max([ldt_tower[0],ldt_access[0]]),
                    "overlap_enddate":min([ldt_tower[-1],ldt_access[-1]])}
    # make the GUI
    access_gui = Tkinter.Toplevel()
    access_gui.wm_title("ACCESS GUI : "+str(ds4.access.keys()))
    access_gui.grid()
    # top row
    nrow = 0
    access_gui.filestartLabel = Tkinter.Label(access_gui,text="Overlap start date")
    access_gui.filestartLabel.grid(row=nrow,column=0,columnspan=3)
    access_gui.fileendLabel = Tkinter.Label(access_gui,text="Overlap end date")
    access_gui.fileendLabel.grid(row=nrow,column=3,columnspan=3)
    # second row
    nrow = nrow + 1
    access_gui.filestartValue = Tkinter.Label(access_gui,text=str(access_dates["overlap_startdate"]))
    access_gui.filestartValue.grid(row=nrow,column=0,columnspan=3)
    access_gui.fileendValue = Tkinter.Label(access_gui,text=str(access_dates["overlap_enddate"]))
    access_gui.fileendValue.grid(row=nrow,column=3,columnspan=3)
    # third row
    nrow = nrow + 1
    access_gui.startLabel = Tkinter.Label(access_gui, text="Start date (YYYY-MM-DD)")
    access_gui.startLabel.grid(row=nrow,column=0,columnspan=3)
    access_gui.startEntry = Tkinter.Entry(access_gui)
    access_gui.startEntry.grid(row=nrow,column=3,columnspan=3)
    # fourth row
    nrow = nrow + 1
    access_gui.endLabel = Tkinter.Label(access_gui, text="End date   (YYYY-MM-DD)")
    access_gui.endLabel.grid(row=nrow,column=0,columnspan=3)
    access_gui.endEntry = Tkinter.Entry(access_gui)
    access_gui.endEntry.grid(row=nrow,column=3,columnspan=3)
    # bottom row
    nrow = nrow + 1
    access_gui.runButton = Tkinter.Button(access_gui,text="Run",command=lambda:gfACCESS_main(ds4,ds_access,access_gui,access_dates))
    access_gui.runButton.grid(row=nrow,column=0,columnspan=2)
    access_gui.acceptButton = Tkinter.Button(access_gui,text="Done",command=lambda:gfACCESS_finish(ds4,access_gui))
    access_gui.acceptButton.grid(row=nrow,column=2,columnspan=2)
    access_gui.doneButton = Tkinter.Button(access_gui,text="Close plots",command=lambda:gfACCESS_closeplots())
    access_gui.doneButton.grid(row=nrow,column=4,columnspan=2)
    access_gui.wait_window(access_gui)

def gfACCESS_finish(ds,access_gui):
    # close any open plot windows
    for i in plt.get_fignums(): plt.close(i)
    # destroy the ACCESS GUI
    access_gui.destroy()
    # write Excel spreadsheet with fit statistics
    qcio.xl_write_ACCESSStats(ds)
    # remove the ACCESS dictionary from the data structure
    del ds.access

def gfACCESS_closeplots():
    for i in plt.get_fignums(): plt.close(i)

def gfACCESS_main(ds_tower,ds_access,access_gui,access_dates):
    '''
    This is the main routine for using ACCESS data to gap fill drivers.
    '''
    log.info(" Gap filling "+str(ds_tower.access.keys())+" using ACCESS data")
    # read the control file again, this allows the contents of the control file to
    # be changed with the ACCESS GUI still displayed
    cfname = ds_tower.globalattributes["controlfile_name"]
    cf = qcio.get_controlfilecontents(cfname)
    # get the site name
    site_name = ds_tower.globalattributes["site_name"]
    # get the time step and a local pointer to the datetime series
    ts = int(ds_tower.globalattributes["time_step"])
    # check the time step of the tower and ACCESS data are the same
    if ts!=int(ds_access.globalattributes["time_step"]):
        log.error("gfACCESS_main: ACCESS and tower time steps are different")
        return
    # get a local datetime series
    ldt_tower = ds_tower.series["DateTime"]["Data"]
    ldt_access = ds_access.series["DateTime"]["Data"]
    xldt_tower = ds_tower.series["xlDateTime"]["Data"]
    # get the start and end datetimes entered in the ACCESS GUI
    startdate = access_gui.startEntry.get()
    if len(startdate)==0: startdate = access_dates["overlap_startdate"]
    enddate = access_gui.endEntry.get()
    if len(enddate)==0: enddate = access_dates["overlap_enddate"]
    # get the indices of the start and end datetimes in the tower and the ACCESS data.
    si_tower = qcutils.GetDateIndex(ldt_tower,str(startdate),ts=ts,match="startnextday",default=0)
    ei_tower = qcutils.GetDateIndex(ldt_tower,str(enddate),ts=ts,match="endpreviousday",default=-1)
    si_access = qcutils.GetDateIndex(ldt_access,str(startdate),ts=ts,match="startnextday",default=0)
    ei_access = qcutils.GetDateIndex(ldt_access,str(enddate),ts=ts,match="endpreviousday",default=-1)
    # now make sure the start and end datetime match
    startdate = max([ldt_tower[si_tower],ldt_access[si_access]])
    enddate = min([ldt_tower[ei_tower],ldt_access[ei_access]])
    # and get the start and end indices again in case they didn't
    si_tower = qcutils.GetDateIndex(ldt_tower,str(startdate),ts=ts,match="startnextday")
    ei_tower = qcutils.GetDateIndex(ldt_tower,str(enddate),ts=ts,match="endpreviousday")
    si_access = qcutils.GetDateIndex(ldt_access,str(startdate),ts=ts,match="startnextday")
    ei_access = qcutils.GetDateIndex(ldt_access,str(enddate),ts=ts,match="endpreviousday")
    # get the datetime series for the overlap period
    odt_tower = ldt_tower[si_tower:ei_tower+1]
    odt_access = ldt_access[si_access:ei_access+1]
    # get the number of samples per day and the number of days in the period
    nrecs = ei_tower - si_tower + 1
    nperhr = int(float(60)/ts+0.5)
    nperday = int(float(24)*nperhr+0.5)
    ndays = nrecs/nperday
    # now loop over the variables to be gap filled using the ACCESS data
    fig_num = 0
    for label in ds_tower.access.keys():
        ds_tower.access[label]["results"]["startdate"].append(xldt_tower[si_tower])
        ds_tower.access[label]["results"]["enddate"].append(xldt_tower[ei_tower])
        output = ds_tower.access[label]["output"]
        fig_num = fig_num + 1
        title = site_name+' : Comparison of tower and ACCESS data for '+label
        # get the tower data
        data_tower,flag_tower = qcutils.GetSeriesasMA(ds_tower,label,si=si_tower,ei=ei_tower)
        # get the units
        attr_tower = qcutils.GetAttributeDictionary(ds_tower,label)
        units_tower = attr_tower["units"]
        # get a list of ACCESS variables for this tower variable
        access_var_list = [item for item in ds_access.series.keys() if label in item]
        # check the series in the ACCESS data
        if len(access_var_list)==0:
            log.error("gfACCESS_main: series "+label+" not in ACCESS data file")
            continue
        # get the ACCESS series that has the highest correlation with the tower data
        r = numpy.zeros(len(access_var_list))
        for idx,var in enumerate(access_var_list):
            data_access,flag_access = qcutils.GetSeriesasMA(ds_access,var,si=si_access,ei=ei_access)
            attr_access = qcutils.GetAttributeDictionary(ds_access,var)
            units_access = attr_access["units"]
            r[idx] = numpy.ma.corrcoef(data_tower,data_access)[0,1]
        maxidx = numpy.argmax(r)
        # save the name of the ACCESS variable that has the highest correlation with the tower data
        var_maxr = access_var_list[maxidx]
        ds_tower.access[label]["var_maxr"] = var_maxr
        data_access,flag_access = qcutils.GetSeriesasMA(ds_access,var_maxr,si=si_access,ei=ei_access)
        # plot the data for this period
        pd = gfACCESS_initplot(site_name=site_name,label=label,fig_num=fig_num,title=title,
                               ts=ts,ndays=ndays,nperday=nperday,
                               units_tower=units_tower,units_access=units_access)
        gfACCESS_plot(pd,ds_tower,odt_tower,data_tower,odt_access,data_access,r)
        # mask both series when either one is missing
        data_access.mask = numpy.ma.mask_or(data_access.mask,data_tower.mask)
        data_tower.mask = numpy.ma.mask_or(data_access.mask,data_tower.mask)
        # get non-masked versions of the data, these are used with the robust statistics module
        data_access_nm = numpy.ma.compressed(data_access)
        data_tower_nm = numpy.ma.compressed(data_tower)
        # best fit of ACCESS data to tower data by robust least squares
        resrlm = sm.RLM(data_tower_nm,sm.add_constant(data_access_nm,prepend=False),M=sm.robust.norms.TukeyBiweight()).fit()
        m_rlm = resrlm.params[0]; b_rlm = resrlm.params[1]
        data_access_rlm = m_rlm*data_access+b_rlm
        # put the best-fit adjusted ACCESS data into the output series
        ds_tower.series[output]["Data"][si_tower:ei_tower+1] = data_access_rlm
        ds_tower.series[output]["Flag"][si_tower:ei_tower+1] = numpy.int32(20)
    # make sure this processing step gets written to the global attribute "Functions"
    if "GapFillFromACCESS" not in ds_tower.globalattributes["Functions"]:
        ds_tower.globalattributes["Functions"] = ds_tower.globalattributes["Functions"]+", GapFillFromACCESS"

def gfACCESS_initplot(**kwargs):
    pd = {"margin_bottom":0.05,"margin_top":0.05,"margin_left":0.075,"margin_right":0.05,
          "xy_height":0.25,"xy_width":0.20,"xyts_space":0.05,"xyxy_space":0.05,
          "ts_width":0.9}
   # calculate bottom of the first time series and the height of the time series plots
    pd["ts_bottom"] = pd["margin_bottom"]+pd["xy_height"]+pd["xyxy_space"]+pd["xy_height"]+pd["xyts_space"]
    pd["ts_height"] = (1.0 - pd["margin_top"] - pd["ts_bottom"])
    for key, value in kwargs.iteritems():
        pd[key] = value
    return pd

def gfACCESS_plot(pd,ds_tower,odt_tower,data_tower,odt_access,data_access,r):
    # mask both series when either one is missing
    data_access.mask = numpy.ma.mask_or(data_access.mask,data_tower.mask)
    data_tower.mask = numpy.ma.mask_or(data_access.mask,data_tower.mask)
    # get non-masked versions of the data, these are used with the robust statistics module
    data_access_nm = numpy.ma.compressed(data_access)
    data_tower_nm = numpy.ma.compressed(data_tower)
    # turn on interactive plotting
    plt.ion()
    # create the figure canvas
    fig=plt.figure(pd["fig_num"],figsize=(13,9))
    fig.clf()
    plt.figtext(0.5,0.96,pd["title"],ha='center',size=16)
    # top row of XY plots
    # correlation coefficients
    rect1 = [0.10,pd["margin_bottom"]+pd["margin_bottom"]+pd["xy_height"],
             pd["xy_width"],pd["xy_height"]]
    ax1 = plt.axes(rect1)
    ind=numpy.arange(len(r))
    ax1.bar(ind,r,0.35)
    ax1.set_ylabel('r')
    ax1.set_xlabel('Grid')
    plt.draw()
    # lagged correlation
    rect2 = [0.40,pd["margin_bottom"]+pd["margin_bottom"]+pd["xy_height"],
             pd["xy_width"],pd["xy_height"]]
    ax2 = plt.axes(rect2)
    l=ax2.xcorr(data_tower_nm,data_access_nm,maxlags=pd["nperday"])
    ax2.set_ylabel('r')
    ax2.set_xlabel('Lags')
    plt.draw()
    max_lag = pd["ts"]*(numpy.argmax(l[1])-pd["nperday"])
    ds_tower.access[pd["label"]]["results"]["max_lag"].append(max_lag)
    # bottom row of XY plots
    # scatter plot of 30 minute data
    rect3 = [0.10,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax3 = plt.axes(rect3)
    ax3.plot(data_access,data_tower,'b.')
    ax3.set_ylabel('Tower ('+pd["units_tower"]+')')
    ax3.set_xlabel('ACCESS ('+pd["units_access"]+')')
    resrlm = sm.RLM(data_tower_nm,sm.add_constant(data_access_nm,prepend=False),M=sm.robust.norms.TukeyBiweight()).fit()
    pd["m_rlm"] = resrlm.params[0]
    pd["b_rlm"] = resrlm.params[1]
    ds_tower.access[pd["label"]]["results"]["m_rlm"].append(resrlm.params[0])
    ds_tower.access[pd["label"]]["results"]["b_rlm"].append(resrlm.params[1])
    pd["rlm"] = pd["m_rlm"]*data_access+pd["b_rlm"]
    eqnstr = 'y = %.3fx + %.3f'%(pd["m_rlm"],pd["b_rlm"])
    ax3.plot(data_access_nm,resrlm.fittedvalues,'r--',linewidth=3)
    ax3.text(0.5,0.915,eqnstr,fontsize=8,horizontalalignment='center',transform=ax3.transAxes,color='red')
    resols = sm.OLS(data_tower_nm,sm.add_constant(data_access_nm,prepend=False)).fit()
    pd["m_ols"] = resols.params[0]
    pd["b_ols"] = resols.params[1]
    ds_tower.access[pd["label"]]["results"]["m_ols"].append(resols.params[0])
    ds_tower.access[pd["label"]]["results"]["b_ols"].append(resols.params[1])
    pd["ols"] = pd["m_ols"]*data_access+pd["b_ols"]
    eqnstr = 'y = %.3fx + %.3f'%(pd["m_ols"],pd["b_ols"])
    ax3.plot(data_access_nm,resrlm.fittedvalues,'g--',linewidth=3)
    ax3.text(0.5,0.85,eqnstr,fontsize=8,horizontalalignment='center',transform=ax3.transAxes,color='green')
    ax3.text(0.6,0.075,'30 minutes',fontsize=10,horizontalalignment='left',transform=ax3.transAxes)
    plt.draw()
    # scatter plot of daily averages
    rect4 = [0.40,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax4 = plt.axes(rect4)
    data_tower_2d = numpy.ma.reshape(data_tower,[pd["ndays"],pd["nperday"]])
    data_tower_daily_avg = numpy.ma.average(data_tower_2d,axis=1)
    data_access_2d = numpy.ma.reshape(data_access,[pd["ndays"],pd["nperday"]])
    data_access_daily_avg = numpy.ma.average(data_access_2d,axis=1)
    ax4.plot(data_access_daily_avg,data_tower_daily_avg,'b.')
    ax4.set_ylabel('Tower ('+pd["units_tower"]+')')
    ax4.set_xlabel('ACCESS ('+pd["units_access"]+')')
    data_access_daily_avg.mask = numpy.ma.mask_or(data_access_daily_avg.mask,data_tower_daily_avg.mask)
    data_tower_daily_avg.mask = numpy.ma.mask_or(data_access_daily_avg.mask,data_tower_daily_avg.mask)
    data_access_daily_avg_nm = numpy.ma.compressed(data_access_daily_avg)
    data_tower_daily_avg_nm = numpy.ma.compressed(data_tower_daily_avg)
    resrlm = sm.RLM(data_tower_daily_avg_nm,sm.add_constant(data_access_daily_avg_nm,prepend=False),M=sm.robust.norms.TukeyBiweight()).fit()
    eqnstr = 'y = %.3fx + %.3f'%(resrlm.params[0],resrlm.params[1])
    ax4.plot(data_access_daily_avg_nm,resrlm.fittedvalues,'r--',linewidth=3)
    ax4.text(0.5,0.915,eqnstr,fontsize=8,horizontalalignment='center',transform=ax4.transAxes,color='red')
    resrlm = sm.OLS(data_tower_daily_avg_nm,sm.add_constant(data_access_daily_avg_nm,prepend=False)).fit()
    eqnstr = 'y = %.3fx + %.3f'%(resrlm.params[0],resrlm.params[1])
    ax4.plot(data_access_daily_avg_nm,resrlm.fittedvalues,'g--',linewidth=3)
    ax4.text(0.5,0.85,eqnstr,fontsize=8,horizontalalignment='center',transform=ax4.transAxes,color='green')
    ax4.text(0.6,0.075,'Daily average',fontsize=10,horizontalalignment='left',transform=ax4.transAxes)
    plt.draw()
    # diurnal average plot
    rect5 = [0.70,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax5 = plt.axes(rect5)
    data_tower_hourly_avg = numpy.ma.average(data_tower_2d,axis=0)
    data_access_hourly_avg = numpy.ma.average(data_access_2d,axis=0)
    data_access_hourly_rlm = numpy.ma.average(pd["m_rlm"]*data_access_2d+pd["b_rlm"],axis=0)
    data_access_hourly_ols = numpy.ma.average(pd["m_ols"]*data_access_2d+pd["b_ols"],axis=0)
    ind = numpy.arange(len(data_tower_hourly_avg))*float(pd["ts"])/float(60)
    ax5.plot(ind,data_tower_hourly_avg,'ro',label='Tower')
    ax5.plot(ind,data_access_hourly_avg,'b-',label='ACCESS-A')
    ax5.plot(ind,data_access_hourly_rlm,'r-',label='ACCESS-A (RLM)')
    ax5.plot(ind,data_access_hourly_ols,'g-',label='ACCESS-A (OLS)')
    ax5.set_ylabel(pd["label"]+' ('+pd["units_tower"]+')')
    ax5.set_xlim(0,24)
    ax5.xaxis.set_ticks([0,6,12,18,24])
    ax5.set_xlabel('Hour')
    ax5.legend(loc='upper right',frameon=False,prop={'size':8})
    plt.draw()
    # time series
    rect_ts = [pd["margin_left"],pd["ts_bottom"],pd["ts_width"],pd["ts_height"]]
    axes_ts = plt.axes(rect_ts)
    axes_ts.plot(odt_tower,data_tower,'ro',label="Tower")
    axes_ts.plot(odt_access,data_access,'b-',label="ACCESS-A")
    axes_ts.plot(odt_access,pd["rlm"],'r-',label="ACCESS-A (RLM)")
    axes_ts.plot(odt_access,pd["ols"],'g-',label="ACCESS-A (OLS)")
    axes_ts.set_ylabel(pd["label"]+' ('+pd["units_tower"]+')')
    axes_ts.legend(loc='upper right',frameon=False,prop={'size':8})
    plt.draw()
    # now get some statistics
    numpoints = numpy.ma.count(data_tower)
    ds_tower.access[pd["label"]]["results"]["n"].append(numpoints)
    diff = data_tower - data_access
    bias = numpy.ma.average(diff)
    ds_tower.access[pd["label"]]["results"]["bias"].append(bias)
    rmse = numpy.ma.sqrt(numpy.ma.average(diff*diff))
    ds_tower.access[pd["label"]]["results"]["rmse"].append(rmse)
    var_tow = numpy.ma.var(data_tower)
    ds_tower.access[pd["label"]]["results"]["var_tow"].append(var_tow)
    var_acc = numpy.ma.var(data_access)
    ds_tower.access[pd["label"]]["results"]["var_acc"].append(var_acc)
    r_max = numpy.ma.maximum(r)
    ds_tower.access[pd["label"]]["results"]["r_max"].append(r_max)
    text_left = 0.70
    num_left = 0.80
    row_bottom = 0.375
    row_space = 0.030
    i = 0
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'Lag (minutes)')
    plt.figtext(num_left,row_posn,'%.4g'%(max_lag))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'Var (ACCESS)')
    plt.figtext(num_left,row_posn,'%.4g'%(var_acc))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'Var (tower)')
    plt.figtext(num_left,row_posn,'%.4g'%(var_tow))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'RMSE')
    plt.figtext(num_left,row_posn,'%.4g'%(rmse))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'Bias')
    plt.figtext(num_left,row_posn,'%.4g'%(bias))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'r')
    plt.figtext(num_left,row_posn,'%.4g'%(r_max))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'No. points')
    plt.figtext(num_left,row_posn,str(numpoints))
    # save a hard copy of the plot
    sdt = odt_tower[0].strftime("%Y%m%d")
    edt = odt_tower[-1].strftime("%Y%m%d")
    figname = "plots/"+pd["site_name"].replace(" ","")+"_ACCESS_"+pd["label"]
    figname = figname+"_"+sdt+"_"+edt+'.png'
    fig.savefig(figname,format='png')
    # show the plot
    plt.draw()

def GapFillUsingSOLO(dsa,dsb):
    '''
    This is the "Run SOLO" GUI.
    The SOLO GUI is displayed separately from the main OzFluxQC GUI.
    It consists of text to display the start and end datetime of the file,
    two entry boxes for the start and end datetimes of the SOLO run and
    a button to run SOLO ("Run SOLO") and a button to exit the SOLO GUI
    when we are done.  On exit, the OzFluxQC main GUI continues and eventually
    writes the gap filled data to file.
    '''
    if "solo" not in dir(dsb): return
    # local pointer to the datetime series
    ldt = dsb.series["DateTime"]["Data"]
    # set up the GUI
    solo_gui = Tkinter.Toplevel()
    solo_gui.wm_title("SOLO GUI : "+str(dsb.solo.keys()))
    solo_gui.grid()
    # top row
    nrow = 0
    solo_gui.nodesLabel = Tkinter.Label(solo_gui,text="Nodes")
    solo_gui.nodesLabel.grid(row=nrow,column=0,columnspan=1,sticky="E")
    solo_gui.nodesEntry = Tkinter.Entry(solo_gui,width=6)
    solo_gui.nodesEntry.grid(row=nrow,column=1,columnspan=1,sticky="W")
    solo_gui.nodesEntry.insert(0,"Auto")
    #solo_gui.nodesEntry.insert(0,'5')
    solo_gui.trainingLabel = Tkinter.Label(solo_gui,text="Training")
    solo_gui.trainingLabel.grid(row=nrow,column=2,columnspan=1,sticky="E")
    solo_gui.trainingEntry = Tkinter.Entry(solo_gui,width=6)
    solo_gui.trainingEntry.grid(row=nrow,column=3,columnspan=1,sticky="W")
    solo_gui.trainingEntry.insert(0,"500")
    solo_gui.factorLabel = Tkinter.Label(solo_gui,text="Nda factor")
    solo_gui.factorLabel.grid(row=nrow,column=4,columnspan=1,sticky="E")
    solo_gui.factorEntry = Tkinter.Entry(solo_gui,width=6)
    solo_gui.factorEntry.grid(row=nrow,column=5,columnspan=1,sticky="W")
    solo_gui.factorEntry.insert(0,"5")
    # second row
    nrow = nrow + 1
    solo_gui.learningrateLabel = Tkinter.Label(solo_gui,text="Learning")
    solo_gui.learningrateLabel.grid(row=nrow,column=2,columnspan=1,sticky="E")
    solo_gui.learningrateEntry = Tkinter.Entry(solo_gui,width=6)
    solo_gui.learningrateEntry.grid(row=nrow,column=3,columnspan=1,sticky="W")
    solo_gui.learningrateEntry.insert(0,"0.01")
    solo_gui.iterationsLabel = Tkinter.Label(solo_gui,text="Iterations")
    solo_gui.iterationsLabel.grid(row=nrow,column=4,columnspan=1,sticky="E")
    solo_gui.iterationsEntry = Tkinter.Entry(solo_gui,width=6)
    solo_gui.iterationsEntry.grid(row=nrow,column=5,columnspan=1,sticky="W")
    solo_gui.iterationsEntry.insert(0,"500")
    # third row
    nrow = nrow + 1
    solo_gui.filestartLabel = Tkinter.Label(solo_gui,text="File start date")
    solo_gui.filestartLabel.grid(row=nrow,column=0,columnspan=3)
    solo_gui.fileendLabel = Tkinter.Label(solo_gui,text="File end date")
    solo_gui.fileendLabel.grid(row=nrow,column=3,columnspan=3)
    # fourth row
    nrow = nrow + 1
    solo_gui.filestartValue = Tkinter.Label(solo_gui,text=str(ldt[0]))
    solo_gui.filestartValue.grid(row=nrow,column=0,columnspan=3)
    solo_gui.fileendValue = Tkinter.Label(solo_gui,text=str(ldt[-1]))
    solo_gui.fileendValue.grid(row=nrow,column=3,columnspan=3)
    # fifth row
    nrow = nrow + 1
    solo_gui.startLabel = Tkinter.Label(solo_gui, text="Start date (YYYY-MM-DD)")
    solo_gui.startLabel.grid(row=nrow,column=0,columnspan=3)
    solo_gui.startEntry = Tkinter.Entry(solo_gui)
    solo_gui.startEntry.grid(row=nrow,column=3,columnspan=3)
    # sixth row
    nrow = nrow + 1
    solo_gui.endLabel = Tkinter.Label(solo_gui, text="End date   (YYYY-MM-DD)")
    solo_gui.endLabel.grid(row=nrow,column=0,columnspan=3)
    solo_gui.endEntry = Tkinter.Entry(solo_gui)
    solo_gui.endEntry.grid(row=nrow,column=3,columnspan=3)
    # bottom row
    nrow = nrow + 1
    solo_gui.runButton = Tkinter.Button (solo_gui, text="Run",command=lambda:gfSOLO_main(dsa,dsb,solo_gui))
    solo_gui.runButton.grid(row=nrow,column=0,columnspan=2)
    solo_gui.doneButton = Tkinter.Button (solo_gui, text="Done",command=lambda:gfSOLO_finish(dsb,solo_gui))
    solo_gui.doneButton.grid(row=nrow,column=2,columnspan=2)
    solo_gui.doneButton = Tkinter.Button (solo_gui, text="Close plots",command=lambda:gfSOLO_closeplots())
    solo_gui.doneButton.grid(row=nrow,column=4,columnspan=2)
    solo_gui.wait_window(solo_gui)

def gfSOLO_finish(ds,solo_gui):
    # close any open plot windows
    for i in plt.get_fignums(): plt.close(i)
    # destroy the SOLO GUI
    solo_gui.destroy()
    # write Excel spreadsheet with fit statistics
    qcio.xl_write_SOLOStats(ds)
    # remove the ACCESS dictionary from the data structure
    del ds.solo

def gfSOLO_closeplots():
    for i in plt.get_fignums(): plt.close(i)

def gfSOLO_main(dsa,dsb,solo_gui):
    '''
    This is the main routine for running SOLO, an artifical neural network for gap filling fluxes.
    '''
    log.info(" Gap filling "+str(dsb.solo.keys())+" using SOLO ANN")
    # read the control file again, this allows the contents of the control file to
    # be changed with the SOLO GUI still displayed
    cfname = dsb.globalattributes["controlfile_name"]
    cf = qcio.get_controlfilecontents(cfname)
    # get some useful things
    site_name = dsa.globalattributes["site_name"]
    # get the time step and a local pointer to the datetime series
    ts = dsb.globalattributes["time_step"]
    ldt = dsb.series["DateTime"]["Data"]
    xldt = dsb.series["xlDateTime"]["Data"]
    startdate = solo_gui.startEntry.get()
    enddate = solo_gui.endEntry.get()
    # get the start and end datetime indices
    si = qcutils.GetDateIndex(ldt,startdate,ts=ts,default=0,match="exact")
    ei = qcutils.GetDateIndex(ldt,enddate,ts=ts,default=-1,match="exact")
    # check the start and end indices
    if si >= ei:
        print " GapFillUsingSOLO: end datetime index smaller that start datetime index"
        return
    if si==0 and ei==-1:
        print " GapFillUsingSOLO: no start and end datetime specified, using all data"
        nRecs = int(dsb.globalattributes["nc_nrecs"])
    else:
        nRecs = ei - si + 1
    # loop over the series to be gap filled using solo
    fig_num = 0
    for series in dsb.solo.keys():
        dsb.solo[series]["results"]["startdate"].append(xldt[si])
        dsb.solo[series]["results"]["enddate"].append(xldt[ei])
        drivers = dsb.solo[series]["drivers"]
        output = dsb.solo[series]["output"]
        # set the number of nodes for the inf files
        nodesAuto = gfSOLO_setnodesEntry(solo_gui,drivers)
        # write the inf files for sofm, solo and seqsolo
        gfSOLO_writeinffiles(solo_gui)
        # run SOFM
        result = gfSOLO_runsofm(dsa,dsb,solo_gui,drivers,series,nRecs,si=si,ei=ei)
        if result!=1: return
        # run SOLO
        result = gfSOLO_runsolo(dsa,dsb,drivers,series,nRecs,si=si,ei=ei)
        if result!=1: return
        # run seqsolo and put the solo_modelled data into the ds series
        result = gfSOLO_runseqsolo(dsa,dsb,drivers,series,output,nRecs,si=si,ei=ei)
        if result!=1: return
        # plot the results
        fig_num = fig_num + 1
        title = site_name+' : Comparison of tower and SOLO data for '+series
        pd = gfSOLO_initplot(site_name=site_name,label=series,fig_num=fig_num,title=title,
                             nDrivers=len(drivers))
        gfSOLO_plot(pd,dsa,dsb,drivers,series,output,solo_gui,si=si,ei=ei)
        # reset the nodesEntry in the solo_gui
        if nodesAuto: gfSOLO_resetnodesEntry(solo_gui)
    if 'GapFillUsingSOLO' not in dsb.globalattributes['Functions']:
        dsb.globalattributes['Functions'] = dsb.globalattributes['Functions']+', GapFillUsingSOLO'

def gfSOLO_setnodesEntry(solo_gui,drivers):
    nodesAuto = False
    if str(solo_gui.nodesEntry.get()).lower()=="auto":
        nodesAuto = True
        solo_gui.nodesEntry.delete(0,Tkinter.END)
        solo_gui.nodesEntry.insert(0,str(len(drivers)+1))
    return nodesAuto

def gfSOLO_resetnodesEntry(solo_gui):
    solo_gui.nodesEntry.delete(0,Tkinter.END)
    solo_gui.nodesEntry.insert(0,"Auto")

def gfSOLO_writeinffiles(solo_gui):
    # sofm inf file
    f = open('solo/inf/sofm.inf','w')
    f.write(str(solo_gui.nodesEntry.get())+'\n')
    f.write(str(solo_gui.trainingEntry.get())+'\n')
    f.write(str(20)+'\n')
    f.write(str(0.01)+'\n')
    f.write(str(1234)+'\n')
    f.write('solo/input/sofm_input.csv'+'\n')
    f.write('solo/output/sofm_1.out'+'\n')
    f.write('solo/output/sofm_2.out'+'\n')
    f.write('solo/output/sofm_3.out'+'\n')
    f.write('solo/output/sofm_4.out'+'\n')
    f.write(str(50)+'\n')
    f.write('### Comment lines ###\n')
    f.write('Line 1: No. of nodes - default is the number of drivers plus 1 (changeable via GUI if used)\n')
    f.write('Line 2: No. of training iterations - default is 500 (changeable via GUI if used)\n')
    f.write('Line 3: No. of iterations per screen output - default is 20\n')
    f.write('Line 4: Spacing between initial weights - default is 0.01\n')
    f.write('Line 5: Seed for random number generator - default is 1234\n')
    f.write('Line 6: input data filename with path relative to current directory\n')
    f.write('Line 7: first output filename with path relative to current directory\n')
    f.write('Line 8: second output filename with path relative to current directory\n')
    f.write('Line 9: third output filename with path relative to current directory\n')
    f.write('Line 10: fourth output filename with path relative to current directory (used by SOLO)\n')
    f.write('Line 11: No. iterations per write of weights to screen - default is 50\n')
    f.close()
    # solo inf file
    f = open('solo/inf/solo.inf','w')
    f.write(str(solo_gui.nodesEntry.get())+'\n')
    f.write(str(solo_gui.factorEntry.get())+'\n')
    f.write('solo/output/sofm_4.out'+'\n')
    f.write('solo/input/solo_input.csv'+'\n')
    f.write('training'+'\n')
    f.write(str(5678)+'\n')
    f.write(str(0)+'\n')
    f.write('solo/output/eigenValue.out'+'\n')
    f.write('solo/output/eigenVector.out'+'\n')
    f.write('solo/output/accumErr.out'+'\n')
    f.write('solo/output/accumRR.out'+'\n')
    f.write('solo/output/trainProcess.out'+'\n')
    f.write('solo/output/freqTable.out'+'\n')
    f.write('solo/output/hidOutputWt.out'+'\n')
    f.write('solo/output/errorMap.out'+'\n')
    f.write('solo/output/finResult.out'+'\n')
    f.write('solo/output/trainWin.out'+'\n')
    f.write('solo/output/trainWout.out'+'\n')
    f.write('### Comment lines ###\n')
    f.write('Line 1: No. of nodes - default is the number of drivers plus 1 (changeable via GUI if used)\n')
    f.write('Line 2: multiplier for minimum number of points per node (NdaFactor) - default is 5 (ie 5*(no. of drivers+1) (changeable via GUI if used)\n')
    f.write('Line 3: fourth output file from SOFM, used as input to SOLO\n')
    f.write('Line 4: input data filename with path relative to current directory\n')
    f.write('Line 5: type of run ("training" or "simulation", always "training" for SOLO)\n')
    f.write('Line 6: seed for random number generator - default is 5678\n')
    f.write('Line 7: "calThreshold", not used by SOLO\n')
    f.write('Lines 8 to 18: output files from SOLO with path relative to current directory\n')
    f.close()
    # seqsolo inf file
    f = open('solo/inf/seqsolo.inf','w')
    f.write(str(solo_gui.nodesEntry.get())+'\n')
    f.write(str(0)+'\n')
    f.write(str(solo_gui.learningrateEntry.get())+'\n')
    f.write(str(solo_gui.iterationsEntry.get())+'\n')
    f.write('solo/output/sofm_4.out'+'\n')
    f.write('solo/input/seqsolo_input.csv'+'\n')
    f.write('simulation'+'\n')
    f.write(str(9100)+'\n')
    f.write(str(0)+'\n')
    f.write('solo/output/eigenValue.out'+'\n')
    f.write('solo/output/eigenVector.out'+'\n')
    f.write('solo/output/trainWout.out'+'\n')
    f.write('solo/output/freqTable.out'+'\n')
    f.write('solo/output/errorMap.out'+'\n')
    f.write('solo/output/finResult.out'+'\n')
    f.write('solo/output/trainingRMSE.out'+'\n')
    f.write('solo/output/seqOut0.out'+'\n')
    f.write('solo/output/seqOut1.out'+'\n')
    f.write('solo/output/seqOut2.out'+'\n')
    f.write('solo/output/seqHidOutW.out'+'\n')
    f.write('solo/output/seqFreqMap.out'+'\n')
    f.write(str(-9999.0)+'\n')
    f.write('### Comment lines ###\n')
    f.write('Line 1: No. of nodes - default is the number of drivers plus 1 (changeable via GUI if used)\n')
    f.write('Line 2: NdaFactor - not used by SEQSOLO, default value is 0\n')
    f.write('Line 3: learning rate - default value 0.01 (must be between 0.0 1nd 1.0, changeable via GUI if used)\n')
    f.write('Line 4: number of iterations for sequential training, default value is 500 (changeable via GUI if used)\n')
    f.write('Line 5: fourth output file from SOFM, used as input file by SEQSOLO\n')
    f.write('Line 6: input data filename with path relative to current directory\n')
    f.write('Line 7: type of run ("training" or "simulation", always "simulation" for SEQSOLO)\n')
    f.write('Line 8: seed for random number generator - default is 9100\n')
    f.write('Line 9: "calThreshold" - minimum number of data points for SOLO node to be used in simulation, default value is 0 (use all nodes)\n')
    f.write('Lines 10 to 21: output files from SEQSOLO with path relative to current directory\n')
    f.write('Line 22: missing data value, default value is -9999\n')
    f.close()

def gfSOLO_runsofm(dsa,dsb,solo_gui,driverlist,targetlabel,nRecs,si=0,ei=-1):
    '''
    Run sofm, the pre-processor for SOLO.
    '''
    # check to see if we need to run sofm again
    # construct the sofm output file name
    ldt = dsb.series['DateTime']['Data']
    sofmoutname=dsb.globalattributes['site_name'].replace(' ','')     # site name
    for x in driverlist: sofmoutname = sofmoutname + x                # drivers
    sofmoutname = sofmoutname + ldt[si].strftime('%Y%m%d%H%M')        # start datetime
    sofmoutname = sofmoutname + ldt[ei].strftime('%Y%m%d%H%M')        # end datetime
    sofmoutname = sofmoutname + str(solo_gui.nodesEntry.get())        # nodes
    sofmoutname = sofmoutname + str(solo_gui.trainingEntry.get())     # sofm training iterations
    sofmoutname = sofmoutname + str(solo_gui.factorEntry.get())       # Nda factor
    sofmoutname = sofmoutname + str(solo_gui.learningrateEntry.get()) # learning rate
    sofmoutname = sofmoutname + str(solo_gui.iterationsEntry.get())   # seqsolo training iterations
    sofmoutname = sofmoutname + '.out'
    sofminfname = sofmoutname + '.inf'
    ## check to see if the user wants to use any existing sofm output
    #if qcutils.cfoptionskey(cf,'UseExistingSOFMOutput'):
        #if gfSOLO_checkforprevioussofmrun(sofmoutname,sofminfname): return
    # get the number of drivers
    ndrivers = len(driverlist)
    # add an extra column for the target data
    sofminputdata = numpy.zeros((nRecs,ndrivers))
    # now fill the driver data array
    i = 0
    badlines = []
    for TheseOnes in driverlist:
        driver,flag = qcutils.GetSeries(dsb,TheseOnes,si=si,ei=ei)
        index = numpy.where(abs(driver-float(-9999)<c.eps))[0]
        if len(index)!=0:
            log.error(' GapFillUsingSOLO: -9999 found in driver '+TheseOnes+' at lines '+str(index))
            badlines = badlines+index.tolist()
        sofminputdata[:,i] = driver[:]
        i = i + 1
    if len(badlines)!=0:
        nBad = len(badlines)
        goodlines = [x for x in range(0,nRecs) if x not in badlines]
        sofminputdata = sofminputdata[goodlines,:]
        log.info(' GapFillUsingSOLO: removed '+str(nBad)+' lines from sofm input file')
        nRecs = len(goodlines)
    # now write the drivers to the SOFM input file
    sofmfile = open('solo/input/sofm_input.csv','wb')
    wr = csv.writer(sofmfile,delimiter=',')
    for i in range(sofminputdata.shape[0]):
        wr.writerow(sofminputdata[i,0:ndrivers])
    sofmfile.close()
    # if the output file from a previous run exists, delete it
    if os.path.exists('solo/output/sofm_4.out'): os.remove('solo/output/sofm_4.out')
    # now run SOFM
    #log.info(' GapFillUsingSOLO: running SOFM')
    sofmlogfile = open('solo/log/sofm.log','wb')
    subprocess.call(['./solo/bin/sofm','solo/inf/sofm.inf'],stdout=sofmlogfile)
    sofmlogfile.close()
    # check to see if the sofm output file exists, this is used to indicate that sofm ran correctly
    if os.path.exists('solo/output/sofm_4.out'):
        # write out the details of this run
        gfSOLO_savethissofmrun(sofmoutname,sofminfname)
        return 1
    else:
        log.error(' gfSOLO_runsofm: SOFM did not run correctly, check the SOLO GUI and the log files')
        return 0

def gfSOLO_savethissofmrun(sofmoutname,sofminfname):
    '''
    Writes the output file from sofm to the archive directory so that it is
    available for future runs if needed.
    The sofm output file "sofm_4.out" is saved to the "solo/archive" directory
    so that it can be re-used in a later run.  The file is saved using a new
    file name that is constructed as follows:
     <site_name><target><driverlist><start_datetime><end_datetime>
    Running sofm is slow.  Saving the putput files can speed up the process
    of gap filling using sofm/solo.
    '''
    try:
        os.makedirs('solo/archive')
    except OSError:
        if not os.path.isdir('solo/archive'): raise
    shutil.copy2('solo/output/sofm_4.out','solo/archive/'+sofmoutname)
    shutil.copy2('solo/inf/sofm.inf','solo/archive/'+sofminfname)

def gfSOLO_checkforprevioussofmrun(sofmoutname,sofminfname):
    '''
    Checks to see if there is a sofm output file from a previous run in
    the "solo/archive" directory.
    If a sofm output file from the current site, target, driver list, start
    and end date exists then we will copy the file from the archive to the
    "solo/output" directory and then skip running sofm.  This will save time
    because running sofm is the slowest part of the SOLO gap filling process.
    USEAGE: if gfSOLO_checkforpreviousrun(sofmoutname):
    INPUT: sofmoutname - the name of the sofm output file built from the
                         site name, the target, the driver list and the start
                         and end datetimes.
    RETURNS: logical true if a suitable file from a previous run exists
             logical false if it doesn't
    '''
    if os.path.exists('solo/archive/'+sofmoutname) and os.path.exists('solo/archive/'+sofminfname):
        log.info(' GapFillUsingSOLO: Using matching output file from previous sofm run')
        shutil.copy2('solo/archive/'+sofmoutname,'solo/output/sofm_4.out')
        shutil.copy2('solo/archive/'+sofminfname,'solo/inf/sofm.inf')
        return True
    else:
        return False

def gfSOLO_runsolo(dsa,dsb,driverlist,targetlabel,nRecs,si=0,ei=-1):
    '''
    Run SOLO.
    Note that although we pass in <targetlabel>, we will use <targetlabel>_L3 as the
    real target (<targetlabel>_L3 is the data at L3, prior to any gap filling).  This
    is to avoid using gap filled data (in <targetlabel> from a previous run) to train
    the neural network.
    '''
    ndrivers = len(driverlist)
    # add an extra column for the target data
    soloinputdata = numpy.zeros((nRecs,ndrivers+1))
    # now fill the driver data array, drivers come from the modified ds
    i = 0
    for TheseOnes in driverlist:
        driver,flag = qcutils.GetSeries(dsb,TheseOnes,si=si,ei=ei)
        soloinputdata[:,i] = driver[:]
        i = i + 1
    # a clean copy of the target is pulled from the unmodified ds each time
    target,flag = qcutils.GetSeries(dsa,targetlabel,si=si,ei=ei)
    # now load the target data into the data array
    soloinputdata[:,ndrivers] = target[:]
    # now strip out the bad data
    cind = numpy.zeros(nRecs)
    for i in range(ndrivers+1):
        index = numpy.where(soloinputdata[:,i]==-9999)[0]
        if len(index!=0): cind[index] = 1
    index = numpy.where(cind==0)[0]
    nRecs_good = len(index)
    gooddata = numpy.zeros((nRecs_good,ndrivers+1))
    for i in range(ndrivers+1):
        gooddata[:,i] = soloinputdata[:,i][index]
    # and then write the solo input file, the name is assumed by the solo.inf control file
    solofile = open('solo/input/solo_input.csv','wb')
    wr = csv.writer(solofile,delimiter=',')
    for i in range(gooddata.shape[0]):
        wr.writerow(gooddata[i,0:ndrivers+1])
    solofile.close()
    # if the output file from a previous run exists, delete it
    if os.path.exists('solo/output/eigenValue.out'): os.remove('solo/output/eigenValue.out')
    # now run SOLO
    #log.info(' GapFillUsingSOLO: running SOLO')
    solologfile = open('solo/log/solo.log','wb')
    if 'win' in sys.platform:
        subprocess.call(['./solo/bin/solo.exe','solo/inf/solo.inf'],stdout=solologfile)
    else:
        subprocess.call(['./solo/bin/solo','solo/inf/solo.inf'],stdout=solologfile)
    solologfile.close()
    # check to see if the solo output file exists, this is used to indicate that solo ran correctly
    if os.path.exists('solo/output/eigenValue.out'):
        return 1
    else:
        log.error(' gfSOLO_runsolo: SOLO did not run correctly, check the SOLO GUI and the log files')
        return 0

def gfSOLO_runseqsolo(dsa,dsb,driverlist,targetlabel,outputlabel,nRecs,si=0,ei=-1):
    '''
    Run SEQSOLO.
    '''
    # get the number of drivers    
    ndrivers = len(driverlist)
    # add an extra column for the target data
    seqsoloinputdata = numpy.zeros((nRecs,ndrivers+1))
    # now fill the driver data array
    i = 0
    for TheseOnes in driverlist:
        driver,flag = qcutils.GetSeries(dsb,TheseOnes,si=si,ei=ei)
        seqsoloinputdata[:,i] = driver[:]
        i = i + 1
    # a clean copy of the target is pulled from the unmodified ds each time
    target,flag = qcutils.GetSeries(dsa,targetlabel,si=si,ei=ei)
    # now load the target data into the data array
    seqsoloinputdata[:,ndrivers] = target[:]
    # and then write the seqsolo input file
    seqsolofile = open('solo/input/seqsolo_input.csv','wb')
    wr = csv.writer(seqsolofile,delimiter=',')
    for i in range(seqsoloinputdata.shape[0]):
        wr.writerow(seqsoloinputdata[i,0:ndrivers+1])
    seqsolofile.close()
    # if the output file from a previous run exists, delete it
    if os.path.exists('solo/output/seqOut2.out'): os.remove('solo/output/seqOut2.out')
    # now run SEQSOLO
    #log.info(' GapFillUsingSOLO: running SEQSOLO')
    seqsolologfile = open('solo/log/seqsolo.log','wb')
    if 'win' in sys.platform:
        subprocess.call(['./solo/bin/seqsolo.exe','solo/inf/seqsolo.inf'],stdout=seqsolologfile)
    else:
        subprocess.call(['./solo/bin/seqsolo','solo/inf/seqsolo.inf'],stdout=seqsolologfile)
    seqsolologfile.close()
    # check to see if the solo output file exists, this is used to indicate that solo ran correctly
    if os.path.exists('solo/output/seqOut2.out'):
        # now read in the seqsolo results, use the seqOut2 file so that the learning capability of
        # seqsolo can be used via the "learning rate" and "Iterations" GUI options
        seqdata = numpy.genfromtxt('solo/output/seqOut2.out')
        # if no output series exists then create one now
        if outputlabel not in dsb.series.keys():
            # qcutils.GetSeriesasMA will create a blank series if it doesn't already exist
            flux, flux_flag = qcutils.GetSeriesasMA(dsb,outputlabel)
            attr = qcutils.MakeAttributeDictionary(long_name=targetlabel+' modelled by SOLO',
                                                   units=dsa.series[targetlabel]['Attr']['units'])
            qcutils.CreateSeries(dsb,outputlabel,flux,Flag=flux_flag,Attr=attr)
        # put the SOLO modelled data back into the data series
        if ei==-1:
            dsb.series[outputlabel]['Data'][si:] = seqdata[:,1]
            dsb.series[outputlabel]['Flag'][si:] = numpy.int32(30)
        else:
            dsb.series[outputlabel]['Data'][si:ei+1] = seqdata[:,1]
            dsb.series[outputlabel]['Flag'][si:ei+1] = numpy.int32(30)
        return 1
    else:
        log.error(' gfSOLO_runseqsolo: SEQSOLO did not run correctly, check the SOLO GUI and the log files')
        return 0

def gfSOLO_initplot(**kwargs):
    # set the margins, heights, widths etc
    pd = {"margin_bottom":0.075,"margin_top":0.075,"margin_left":0.05,"margin_right":0.05,
          "xy_height":0.20,"xy_width":0.20,"xyts_space":0.05,"xyts_space":0.05,
          "ts_width":0.9}
    # set the keyword arguments
    for key, value in kwargs.iteritems():
        pd[key] = value
    # calculate bottom of the first time series and the height of the time series plots
    pd["ts_bottom"] = pd["margin_bottom"]+pd["xy_height"]+pd["xyts_space"]
    pd["ts_height"] = (1.0 - pd["margin_top"] - pd["ts_bottom"])/float(pd["nDrivers"]+1)
    return pd

def gfSOLO_plot(pd,dsa,dsb,driverlist,targetlabel,outputlabel,solo_gui,si=0,ei=-1):
    """ Plot the results of the SOLO run. """
    # get the time step
    ts = int(dsb.globalattributes['time_step'])
    # get a local copy of the datetime series
    xdt = numpy.array(dsb.series['DateTime']['Data'][si:ei+1])
    Hdh,f = qcutils.GetSeriesasMA(dsb,'Hdh',si=si,ei=ei)
    # get the observed and modelled values
    obs,f = qcutils.GetSeriesasMA(dsa,targetlabel,si=si,ei=ei)
    mod,f = qcutils.GetSeriesasMA(dsb,outputlabel,si=si,ei=ei)
    # make the figure
    plt.ion()
    fig = plt.figure(pd["fig_num"],figsize=(13,9))
    fig.clf()
    plt.figtext(0.5,0.95,pd["title"],ha='center',size=16)
    # XY plot of the diurnal variation
    rect1 = [0.10,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax1 = plt.axes(rect1)
    # get the diurnal stats of the observations
    Hr1,Av1,Sd1,Mx1,Mn1 = gf_getdiurnalstats(Hdh,obs,ts)
    ax1.plot(Hr1,Av1,'b-',label="Obs")
    # get the diurnal stats of all SOLO predictions
    Hr2,Av2,Sd2,Mx2,Mn2 = gf_getdiurnalstats(Hdh,mod,ts)
    ax1.plot(Hr2,Av2,'r-',label="SOLO(all)")
    if numpy.ma.count_masked(obs)!=0:
        index = numpy.ma.where(obs.mask==False)[0]
        # get the diurnal stats of SOLO predictions when observations are present
        Hr3,Av3,Sd3,Mx3,Mn3=gf_getdiurnalstats(Hdh[index],mod[index],ts)
        ax1.plot(Hr3,Av3,'g-',label="SOLO(obs)")
    plt.xlim(0,24)
    plt.xticks([0,6,12,18,24])
    ax1.set_ylabel(targetlabel)
    ax1.set_xlabel('Hour')
    ax1.legend(loc='upper right',frameon=False,prop={'size':8})
    # XY plot of the 30 minute data
    rect2 = [0.40,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax2 = plt.axes(rect2)
    ax2.plot(mod,obs,'b.')
    ax2.set_ylabel(targetlabel+'_obs')
    ax2.set_xlabel(targetlabel+'_SOLO')
    # plot the best fit line
    coefs = numpy.ma.polyfit(mod,obs,1)
    xfit = numpy.ma.array([numpy.ma.minimum(mod),numpy.ma.maximum(mod)])
    yfit = numpy.polyval(coefs,xfit)
    r = numpy.ma.corrcoef(mod,obs)
    ax2.plot(xfit,yfit,'r--',linewidth=3)
    eqnstr = 'y = %.3fx + %.3f, r = %.3f'%(coefs[0],coefs[1],r[0][1])
    ax2.text(0.5,0.875,eqnstr,fontsize=8,horizontalalignment='center',transform=ax2.transAxes)
    # write the fit statistics to the plot
    numpoints = numpy.ma.count(obs)
    numfilled = numpy.ma.count(mod)-numpy.ma.count(obs)
    diff = mod - obs
    bias = numpy.ma.average(diff)
    dsb.solo[targetlabel]["results"]["bias"].append(bias)
    rmse = numpy.ma.sqrt(numpy.ma.mean((obs-mod)*(obs-mod)))
    plt.figtext(0.65,0.225,'No. points')
    plt.figtext(0.75,0.225,str(numpoints))
    dsb.solo[targetlabel]["results"]["n"].append(numpoints)
    plt.figtext(0.65,0.200,'Nodes')
    plt.figtext(0.75,0.200,str(solo_gui.nodesEntry.get()))
    plt.figtext(0.65,0.175,'Training')
    plt.figtext(0.75,0.175,str(solo_gui.trainingEntry.get()))
    plt.figtext(0.65,0.150,'Nda factor')
    plt.figtext(0.75,0.150,str(solo_gui.factorEntry.get()))
    plt.figtext(0.65,0.125,'Learning rate')
    plt.figtext(0.75,0.125,str(solo_gui.learningrateEntry.get()))
    plt.figtext(0.65,0.100,'Iterations')
    plt.figtext(0.75,0.100,str(solo_gui.iterationsEntry.get()))
    plt.figtext(0.815,0.225,'No. filled')
    plt.figtext(0.915,0.225,str(numfilled))
    plt.figtext(0.815,0.200,'Slope')
    plt.figtext(0.915,0.200,str(qcutils.round2sig(coefs[0],sig=4)))
    dsb.solo[targetlabel]["results"]["m_ols"].append(coefs[0])
    plt.figtext(0.815,0.175,'Offset')
    plt.figtext(0.915,0.175,str(qcutils.round2sig(coefs[1],sig=4)))
    dsb.solo[targetlabel]["results"]["b_ols"].append(coefs[1])
    plt.figtext(0.815,0.150,'r')
    plt.figtext(0.915,0.150,str(qcutils.round2sig(r[0][1],sig=4)))
    dsb.solo[targetlabel]["results"]["r_max"].append(r[0][1])
    plt.figtext(0.815,0.125,'RMSE')
    plt.figtext(0.915,0.125,str(qcutils.round2sig(rmse,sig=4)))
    dsb.solo[targetlabel]["results"]["rmse"].append(rmse)
    var_obs = numpy.ma.var(obs)
    dsb.solo[targetlabel]["results"]["var_obs"].append(var_obs)
    var_mod = numpy.ma.var(mod)
    dsb.solo[targetlabel]["results"]["var_mod"].append(var_mod)
    plt.draw()
    # time series of drivers and target
    ts_axes = []
    rect = [pd["margin_left"],pd["ts_bottom"],pd["ts_width"],pd["ts_height"]]
    ts_axes.append(plt.axes(rect))
    ts_axes[0].plot(xdt,obs,'b.',xdt,mod,'r-')
    ts_axes[0].set_xlim(xdt[0],xdt[-1])
    TextStr = targetlabel+'_obs ('+dsa.series[targetlabel]['Attr']['units']+')'
    ts_axes[0].text(0.05,0.85,TextStr,color='b',horizontalalignment='left',transform=ts_axes[0].transAxes)
    TextStr = outputlabel+'('+dsb.series[outputlabel]['Attr']['units']+')'
    ts_axes[0].text(0.85,0.85,TextStr,color='r',horizontalalignment='right',transform=ts_axes[0].transAxes)
    plt.draw()
    for ThisOne,i in zip(driverlist,range(1,pd["nDrivers"]+1)):
        this_bottom = pd["ts_bottom"] + i*pd["ts_height"]
        rect = [pd["margin_left"],this_bottom,pd["ts_width"],pd["ts_height"]]
        ts_axes.append(plt.axes(rect,sharex=ts_axes[0]))
        data,flag = qcutils.GetSeriesasMA(dsb,ThisOne,si=si,ei=ei)
        ts_axes[i].plot(xdt,data)
        plt.setp(ts_axes[i].get_xticklabels(),visible=False)
        TextStr = ThisOne+'('+dsb.series[ThisOne]['Attr']['units']+')'
        ts_axes[i].text(0.05,0.85,TextStr,color='b',horizontalalignment='left',transform=ts_axes[i].transAxes)
        plt.draw()
    # save a hard copy of the plot
    sdt = xdt[0].strftime("%Y%m%d")
    edt = xdt[-1].strftime("%Y%m%d")
    figname = "plots/"+pd["site_name"].replace(" ","")+"_SOLO_"+pd["label"]
    figname = figname+"_"+sdt+"_"+edt+'.png'
    fig.savefig(figname,format='png')

def gf_getdiurnalstats(DecHour,Data,dt):
    nInts = 24*int((60/dt)+0.5)
    Hr = numpy.array([-9999]*nInts,dtype=numpy.float64)
    Av = numpy.array([-9999]*nInts,dtype=numpy.float64)
    Sd = numpy.array([-9999]*nInts,dtype=numpy.float64)
    Mx = numpy.array([-9999]*nInts,dtype=numpy.float64)
    Mn = numpy.array([-9999]*nInts,dtype=numpy.float64)
    for i in range(nInts):
        Hr[i] = float(i)*dt/60.
        li = numpy.where((abs(DecHour-Hr[i])<c.eps)&(abs(Data-float(-9999))>c.eps))
        if numpy.size(li)!=0:
            Av[i] = numpy.mean(Data[li])
            Sd[i] = numpy.std(Data[li])
            Mx[i] = numpy.max(Data[li])
            Mn[i] = numpy.min(Data[li])
    return Hr, Av, Sd, Mx, Mn

def gf_getdateticks(start, end):
    from datetime import timedelta as td
    delta = end - start
    if delta <= td(minutes=10):
        loc = mdt.MinuteLocator()
        fmt = mdt.DateFormatter('%H:%M')
    elif delta <= td(minutes=30):
        loc = mdt.MinuteLocator(byminute=range(0,60,5))
        fmt = mdt.DateFormatter('%H:%M')
    elif delta <= td(hours=1):
        loc = mdt.MinuteLocator(byminute=range(0,60,15))
        fmt = mdt.DateFormatter('%H:%M')
    elif delta <= td(hours=6):
        loc = mdt.HourLocator()
        fmt = mdt.DateFormatter('%H:%M')
    elif delta <= td(days=1):
        loc = mdt.HourLocator(byhour=range(0,24,3))
        fmt = mdt.DateFormatter('%H:%M')
    elif delta <= td(days=3):
        loc = mdt.HourLocator(byhour=range(0,24,12))
        fmt = mdt.DateFormatter('%d/%m %H')
    elif delta <= td(weeks=2):
        loc = mdt.DayLocator()
        fmt = mdt.DateFormatter('%d/%m')
    elif delta <= td(weeks=12):
        loc = mdt.WeekdayLocator()
        fmt = mdt.DateFormatter('%d/%m')
    elif delta <= td(weeks=104):
        loc = mdt.MonthLocator()
        fmt = mdt.DateFormatter('%d/%m')
    elif delta <= td(weeks=208):
        loc = mdt.MonthLocator(interval=3)
        fmt = mdt.DateFormatter('%d/%m/%y')
    else:
        loc = mdt.MonthLocator(interval=6)
        fmt = mdt.DateFormatter('%d/%m/%y')
    return loc,fmt