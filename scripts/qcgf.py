import ast
from calendar import isleap
from configobj import ConfigObj
import constants as c
import csv
import datetime
import dateutil
import logging
import numpy
import matplotlib as mpl
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
    Ts,flag,attr = qcutils.GetSeriesasMA(ds,driver_label)
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
    Fsd,Fsd_flag,Fsd_attr = qcutils.GetSeriesasMA(ds,'Fsd')
    ustar,ustar_flag,ustar_attr = qcutils.GetSeriesasMA(ds,'ustar')
    Fc_label = str(cf['Derived']['NEE']['Fc'])
    Fc,Fc_flag,Fc_attr = qcutils.GetSeriesasMA(ds,Fc_label)
    Reco_label = str(cf['Derived']['NEE']['Reco'])
    Reco,Reco_flag,Reco_attr = qcutils.GetSeriesasMA(ds,Reco_label)
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
        NEE,f,a = qcutils.GetSeriesasMA(ds,NEE_label)
    if Reco_label not in ds.series.keys():
        log.info('PartitionNEE: requested series '+Reco_label+' is not in the data structure')
        return
    else:
        Reco,f,a = qcutils.GetSeriesasMA(ds,Reco_label)
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
    climatology_file_list = list(set([ds.climatology[x]["file_name"] for x in ds.climatology.keys()]))
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
    # loop over the series to be gap filled using climatology
    for label in ds.climatology.keys():
        # check to see if there are any gaps in "series"
        index = numpy.where(abs(ds.series[label]['Data']-float(-9999))<c.eps)[0]
        if len(index)==0: continue                      # no gaps found in "series"
        # local pointers to output name and method
        output = ds.climatology[label]["output"]
        method = ds.climatology[label]["method"]
        # do the gap filling
        log.info(" Gap filling "+label+" using climatology")
        # choose the gap filling method
        if method=="monthly":
            gfClimatology_monthly(ds,label,output,cli_xlbook)
        elif method=="interpolated daily":
            gfClimatology_interpolateddaily(ds,label,cli_xlbook)
        else:
            log.error(" GapFillFromClimatology: unrecognised method option for "+label)
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
    
def gfClimatology_interpolateddaily(ds,series,xlbook):
    # gap fill from interpolated 30 minute data
    ts = ds.globalattributes["time_step"]
    ldt = ds.series["DateTime"]["Data"]
    output = ds.climatology[series]["output"]
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
    data,flag,attr = qcutils.GetSeriesasMA(ds,series)
    # get an index of missing values
    idx = numpy.ma.where(data.mask==True)[0]
    # there must be a better way to do this ...
    # simply using the index (idx) to set a slice of the data array to the gap filled values in val1d
    # does not seem to work (mask stays true on replaced values in data), the work around is to
    # step through the indices, find the time of the missing value in data, find the same time in the
    # gap filled values val1d and set the missing element of data to this element of val1d
    # actually ...
    # this may not be the fastest but it may be the most robust because it matches dates of missing data
    # to dates in the climatology file
    li = 0
    for ii in idx:
        try:
            jj = cdt.index(ldt[ii],li)
            li = jj
            data[ii] = val1d[jj]
            flag[ii] = numpy.int32(40)
        except ValueError:
            data[ii] = numpy.float64(-9999)
            flag[ii] = numpy.int32(41)
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
    Fsd,Fsd_flag,Fsd_attr = qcutils.GetSeriesasMA(ds,'Fsd')
    us,us_flag,us_attr = qcutils.GetSeriesasMA(ds,'ustar')
    driver,driver_flag,driver_attr = qcutils.GetSeriesasMA(ds, driver_label)
    nRecs = len(driver)
    # get the flux series to be gap filled
    flux,flux_flag,flux_attr = qcutils.GetSeriesasMA(ds,series)
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

def gfACCESS_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about the ACCESS data used
        to gap fill the tower data."""
    # get the section of the control file containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # return without doing anything if the series isn't in a control file section
    if len(section)==0:
        log.error("gfACCESS_createACCESSdict: Series "+series+" not found in control file, skipping ...")
        return
    # create the ACCESS directory in the data structure
    if "access" not in dir(ds): ds.access = {}
    # create the dictionary keys for this series
    ds.access[series] = {}
    # site name
    ds.access[series]["site_name"] = ds.globalattributes["site_name"]
    # ACCESS file name
    ds.access[series]["file_name"] = cf[section][series]["GapFillFromACCESS"]["file_name"]
    # ACCESS variable name if different from name used in control file
    if "access_name" in cf[section][series]["GapFillFromACCESS"]:
        ds.access[series]["access_name"] = cf[section][series]["GapFillFromACCESS"]["access_name"]
    else:
        ds.access[series]["access_name"] = series
    # name of ACCESS output series in ds
    ds.access[series]["output"] = cf[section][series]["GapFillFromACCESS"]["output"]
    # results of best fity for plotting later on
    ds.access[series]["results"] = {"startdate":[],"middate":[],"enddate":[],"n":[],"r_max":[],
                                    "bias":[],"rmse":[],"fb":[],"nmse":[],
                                    "avg_tow":[],"avg_acc":[],
                                    "var_tow":[],"var_acc":[],"var_ratio":[],
                                    "lag_maxr":[],"m_rlm":[],"b_rlm":[],"m_ols":[],"b_ols":[]}

def gfClimatology_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about the climatological data used
        to gap fill the tower data."""
    # get the section of the control file containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # return without doing anything if the series isn't in a control file section
    if len(section)==0:
        log.error("GapFillFromClimatology: Series "+series+" not found in control file, skipping ...")
        return
    # create the climatology directory in the data structure
    if "climatology" not in dir(ds): ds.climatology = {}
    # create the dictionary keys for this series
    ds.climatology[series] = {}
    # site name
    ds.climatology[series]["site_name"] = ds.globalattributes["site_name"]
    # Climatology file name
    ds.climatology[series]["file_name"] = cf[section][series]["GapFillFromClimatology"]["file_name"]
    # name of climatology output series in ds
    ds.climatology[series]["output"] = cf[section][series]["GapFillFromClimatology"]["output"]
    # climatology gap filling method
    if "method" not in cf[section][series]["GapFillFromClimatology"].keys():
        # default if "method" missing is "interpolated_daily"
        ds.climatology[series]["method"] = "interpolated_daily"
    else:
        ds.climatology[series]["method"] = cf[section][series]["GapFillFromClimatology"]["method"]

def GapFillParseControlFile(cf,ds,series=""):
    # find the section containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # return empty handed if the series is not in a section
    if len(section)==0: return
    if "GapFillFromACCESS" in cf[section][series].keys():
        # create the ACCESS dictionary in ds
        gfACCESS_createdict(cf,ds,series)
        # create an empty series in ds if the ACCESS output series doesn't exist yet
        if ds.access[series]["output"] not in ds.series.keys():
            data,flag,attr = qcutils.MakeEmptySeries(ds,ds.access[series]["output"])
            qcutils.CreateSeries(ds,ds.access[series]["output"],data,Flag=flag,Attr=attr)
    if "GapFillUsingSOLO" in cf[section][series].keys():
        # create the SOLO directory in the data structure
        if "solo" not in dir(ds): ds.solo = {}
        # create a dictionary for this series
        ds.solo[series] = {}
        ds.solo[series]["drivers"] = ast.literal_eval(cf[section][series]["GapFillUsingSOLO"]["drivers"])
        ds.solo[series]["output"] = cf[section][series]["GapFillUsingSOLO"]["output"]
        ds.solo[series]["results"] = {"startdate":[],"middate":[],"enddate":[],"n":[],"r_max":[],
                                      "bias":[],"rmse":[],"var_obs":[],"var_mod":[],"m_ols":[],"b_ols":[]}
        data,flag,attr = qcutils.MakeEmptySeries(ds,ds.solo[series]["output"])
        qcutils.CreateSeries(ds,ds.solo[series]["output"],data,Flag=flag,Attr=attr)
    if "GapFillFromClimatology" in cf[section][series].keys():
        # create the climatology dictionary in the data structure
        gfClimatology_createdict(cf,ds,series)
        # create an empty series in ds if the ACCESS output series doesn't exist yet
        if ds.climatology[series]["output"] not in ds.series.keys():
            data,flag,attr = qcutils.MakeEmptySeries(ds,ds.climatology[series]["output"])
            qcutils.CreateSeries(ds,ds.climatology[series]["output"],data,Flag=flag,Attr=attr)
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
    # check the time step of the tower and ACCESS data are the same
    if int(ds4.globalattributes["time_step"])!=int(ds_access.globalattributes["time_step"]):
        log.error(" ACCESS and tower time steps are different")
        return
    ldt_access = ds_access.series["DateTime"]["Data"]
    # check that the ACCESS and tower data overlap
    if ldt_access[0]>ldt_tower[-1]:
        log.error(" ACCESS data starts after tower data finishes")
        return
    if ldt_access[-1]<ldt_tower[0]:
        log.error(" ACCESS data end before tower data starts")
        return
    # store the start and end date of the overlap period
    startdate = max([ldt_tower[0],ldt_access[0]])
    enddate = min([ldt_tower[-1],ldt_access[-1]])
    access_info = {"overlap_startdate":startdate.strftime("%Y-%m-%d %H:%M"),
                    "overlap_enddate":enddate.strftime("%Y-%m-%d %H:%M")}
    # make the GUI
    access_gui = Tkinter.Toplevel()
    access_gui.wm_title("ACCESS GUI")
    access_gui.grid()
    # top row
    nrow = 0
    access_gui.filestartLabel = Tkinter.Label(access_gui,text="Overlap start date")
    access_gui.filestartLabel.grid(row=nrow,column=0,columnspan=3)
    access_gui.fileendLabel = Tkinter.Label(access_gui,text="Overlap end date")
    access_gui.fileendLabel.grid(row=nrow,column=3,columnspan=3)
    # second row
    nrow = nrow + 1
    access_gui.filestartValue = Tkinter.Label(access_gui,text=access_info["overlap_startdate"])
    access_gui.filestartValue.grid(row=nrow,column=0,columnspan=3)
    access_gui.fileendValue = Tkinter.Label(access_gui,text=access_info["overlap_enddate"])
    access_gui.fileendValue.grid(row=nrow,column=3,columnspan=3)
    # third row
    nrow = nrow + 1
    access_gui.startLabel = Tkinter.Label(access_gui, text="Start date (YYYY-MM-DD)")
    access_gui.startLabel.grid(row=nrow,column=0,columnspan=3)
    access_gui.startEntry = Tkinter.Entry(access_gui,width=15)
    access_gui.startEntry.grid(row=nrow,column=3,columnspan=3)
    # fourth row
    nrow = nrow + 1
    access_gui.endLabel = Tkinter.Label(access_gui, text="End date   (YYYY-MM-DD)")
    access_gui.endLabel.grid(row=nrow,column=0,columnspan=3)
    access_gui.endEntry = Tkinter.Entry(access_gui,width=15)
    access_gui.endEntry.grid(row=nrow,column=3,columnspan=3)
    # fifth row
    nrow = nrow + 1
    access_gui.peropt = Tkinter.IntVar()
    access_gui.peropt.set(1)
    access_gui.manualperiod = Tkinter.Radiobutton(access_gui,text="Manual",variable=access_gui.peropt,value=1)
    access_gui.manualperiod.grid(row=nrow,column=0,columnspan=3,sticky="W")
    access_gui.doallplots = Tkinter.Radiobutton(access_gui,text="Monthly",variable=access_gui.peropt,value=2)
    access_gui.doallplots.grid(row=nrow,column=3,columnspan=3,sticky="W")
    # sixth row
    nrow = nrow + 1
    access_gui.daysperiod = Tkinter.Radiobutton(access_gui,text="No. days",variable=access_gui.peropt,value=3)
    access_gui.daysperiod.grid(row=nrow,column=0,sticky="W")
    access_gui.daysentry = Tkinter.Entry(access_gui,width=5)
    access_gui.daysentry.grid(row=nrow,column=1,columnspan=1,sticky="W")
    access_gui.pointsperiod = Tkinter.Radiobutton(access_gui,text="No. pts",variable=access_gui.peropt,value=4)
    access_gui.pointsperiod.grid(row=nrow,column=3,sticky="W")
    access_gui.pointsentry = Tkinter.Entry(access_gui,width=5)
    access_gui.pointsentry.grid(row=nrow,column=4,columnspan=1,sticky="W")
    # seventh row
    nrow = nrow + 1
    access_gui.minptsLabel = Tkinter.Label(access_gui,text="Min points")
    access_gui.minptsLabel.grid(row=nrow,column=0,columnspan=1,sticky="E")
    access_gui.minpts = Tkinter.Entry(access_gui,width=5)
    access_gui.minpts.grid(row=nrow,column=1,columnspan=1,sticky="W")
    access_gui.minpts.insert(0,"100")
    access_gui.owopt = Tkinter.IntVar()
    access_gui.owopt.set(1)
    access_gui.overwrite = Tkinter.Checkbutton(access_gui, text="Overwrite", variable=access_gui.owopt)
    access_gui.overwrite.grid(row=nrow,column=3,columnspan=2,sticky="w")
    # eighth row
    nrow = nrow + 1
    access_gui.runButton = Tkinter.Button(access_gui,text="Run",command=lambda:gfACCESS_run(ds4,ds_access,access_gui,access_info))
    access_gui.runButton.grid(row=nrow,column=0,columnspan=3)
    access_gui.acceptButton = Tkinter.Button(access_gui,text="Done",command=lambda:gfACCESS_done(ds4,access_gui))
    access_gui.acceptButton.grid(row=nrow,column=3,columnspan=3)
    # ninth row
    nrow = nrow + 1
    access_gui.progress_row = nrow
    access_gui.progress = Tkinter.Label(access_gui, text='Waiting for input ...')
    access_gui.progress.grid(row=nrow,column=0,columnspan=6,sticky="W")

    access_gui.wait_window(access_gui)

def gfACCESS_progress(access_gui,text):
    """
        Update progress message in ACCESS GUI
        """
    access_gui.progress.destroy()
    access_gui.progress = Tkinter.Label(access_gui, text=text)
    access_gui.progress.grid(row=9,column=0,columnspan=2,sticky="W")
    access_gui.update()

def gfACCESS_done(ds,access_gui):
    # plot the summary statistics if required
    if access_gui.peropt.get()==1: gfACCESS_plotsummary(ds)
    # destroy the ACCESS GUI
    access_gui.destroy()
    # write Excel spreadsheet with fit statistics
    qcio.xl_write_ACCESSStats(ds)
    # remove the ACCESS dictionary from the data structure
    del ds.access

def gfACCESS_plotsummary(ds):
    """ Plot single pages of summary results for groups of variables. """
    # get a list of variables for which ACCESS data was available
    label_list = ds.access.keys()
    if len(ds.access[label_list[0]]["results"]["startdate"])==0:
        log.info("gfACCESS: no summary data to plot")
        return
    # get the Excel datemode, needed to convert the Excel datetime to Python datetimes
    datemode = int(ds.globalattributes['xl_datemode'])
    # site name for titles
    site_name = ds.globalattributes["site_name"]
    # datetimes are stored in ds.access as Excel datetimes, here we convert to Python datetimes
    # for ease of handling and plotting.
    # start datetimes of the periods compared first
    basedate = datetime.datetime(1899, 12, 30)
    dt_start = []
    for xldt in ds.access[label_list[0]]["results"]["startdate"]:
        dt_start.append(basedate+datetime.timedelta(days=xldt+1462*datemode))
    startdate = min(dt_start)
    # and then the end datetimes
    dt_end = []
    for xldt in ds.access[label_list[0]]["results"]["enddate"]:
        dt_end.append(basedate+datetime.timedelta(days=xldt+1462*datemode))
    enddate = max(dt_end)
    # get the major tick locator and label format
    MTLoc = mdt.AutoDateLocator(minticks=3,maxticks=5)
    MTFmt = mdt.DateFormatter('%b')
    # group lists of the resuts to be plotted
    result_list = ["r_max","bias","rmse","var_ratio","lag_maxr","m_ols","b_ols"]
    ylabel_list = ["r","Bias","RMSE","Var ratio","Lag","Slope","Offset"]
    # turn on interactive plotting
    plt.ion()
    # now loop over the group lists
    for nFig in ds.cf["ACCESS_Summary"].keys():
        plot_title = ds.cf["ACCESS_Summary"][str(nFig)]["Title"]
        var_list = ast.literal_eval(ds.cf["ACCESS_Summary"][str(nFig)]["Variables"])
        # set up the subplots on the page
        fig,axs = plt.subplots(len(result_list),len(var_list),figsize=(13,9))
        fig.canvas.set_window_title("ACCESS summary: "+plot_title)
        # make a title string for the plot and render it
        title_str = "ACCESS: "+plot_title+"; "+site_name+" "+datetime.datetime.strftime(startdate,"%Y-%m-%d")
        title_str = title_str+" to "+datetime.datetime.strftime(enddate,"%Y-%m-%d")
        fig.suptitle(title_str, fontsize=14, fontweight='bold')
        # initialise a string to take the concatenated variable names, used in the name of the hard-copy of the plot
        figlab = ""
        # now loop over the variables in the group list
        for col,label in enumerate(var_list):
            # append the variable name to the variable name string
            figlab = figlab+label
            # and loop over rows in plot
            for row,rlabel,ylabel in zip(range(len(result_list)),result_list,ylabel_list):
                # get the results to be plotted
                result = numpy.ma.masked_equal(ds.access[label]["results"][rlabel],float(-9999))
                # put the data into the right order to be plotted
                dt,data = gfACCESS_plotsummary_getdata(dt_start,dt_end,result)
                # plot the results
                axs[row,col].plot(dt,data)
                # put in the major ticks
                axs[row,col].xaxis.set_major_locator(MTLoc)
                # if this is the left-most column, add the Y axis labels
                if col==0: axs[row,col].set_ylabel(ylabel,visible=True)
                # if this is not the last row, hide the tick mark labels
                if row<len(result_list)-1: plt.setp(axs[row,col].get_xticklabels(),visible=False)
                # if this is the first row, add the column title
                if row==0: axs[row,col].set_title(label)
                # if this is the last row, add the major tick mark and axis labels
                if row==len(result_list)-1:
                    axs[row,col].xaxis.set_major_formatter(MTFmt)
                    axs[row,col].set_xlabel('Month',visible=True)
        # draw the plot
        plt.draw()
        # make the hard-copy file name and save the plot as a PNG file
        sdt = startdate.strftime("%Y%m%d")
        edt = enddate.strftime("%Y%m%d")
        figname = "plots/"+site_name.replace(" ","")+"_ACCESS_FitStatistics_"+figlab
        figname = figname+"_"+sdt+"_"+edt+".png"
        fig.savefig(figname,format="png")

def gfACCESS_plotsummary_getdata(dt_start,dt_end,result):
    dt = []
    data = []
    for s,e,r in zip(dt_start,dt_end,result):
        dt.append(s)
        data.append(r)
        dt.append(e)
        data.append(r)
    return dt,data

def gfACCESS_run(ds_tower,ds_access,access_gui,access_info):
    # get some settings from the ACCESS GUI
    access_info["peropt"] = access_gui.peropt.get()
    access_info["overwrite"] = True
    if access_gui.owopt.get()==1: access_info["overwrite"] = True
    access_info["min_points"] = int(access_gui.minpts.get())
    log.info(" Gap filling "+str(ds_tower.access.keys())+" using ACCESS data")
    if access_gui.peropt.get()==1:
        gfACCESS_progress(access_gui,"Starting manual run ...")
        # get the start and end datetimes entered in the ACCESS GUI
        access_info["startdate"] = access_gui.startEntry.get()
        if len(access_info["startdate"])==0: access_info["startdate"] = access_info["overlap_startdate"]
        access_info["enddate"] = access_gui.endEntry.get()
        if len(access_info["enddate"])==0: access_info["enddate"] = access_info["overlap_enddate"]
        gfACCESS_main(ds_tower,ds_access,access_info)
        gfACCESS_progress(access_gui,"Finished manual run ...")
    elif access_gui.peropt.get()==2:
        gfACCESS_progress(access_gui,"Starting auto (monthly) run ...")
        # get the start datetime entered in the ACCESS GUI
        access_info["startdate"] = access_gui.startEntry.get()
        if len(access_info["startdate"])==0: access_info["startdate"] = access_info["overlap_startdate"]
        startdate = dateutil.parser.parse(access_info["startdate"])
        overlap_startdate = dateutil.parser.parse(access_info["overlap_startdate"])
        overlap_enddate = dateutil.parser.parse(access_info["overlap_enddate"])
        enddate = startdate+dateutil.relativedelta.relativedelta(months=1)
        enddate = min([overlap_enddate,enddate])
        access_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<overlap_enddate:
            gfACCESS_main(ds_tower,ds_access,access_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(months=1)
            access_info["startdate"] = startdate.strftime("%Y-%m-%d")
            access_info["enddate"] = enddate.strftime("%Y-%m-%d")
        # plot the summary statistics
        gfACCESS_plotsummary(ds_tower)
        gfACCESS_progress(access_gui,"Finished auto (monthly) run ...")
    elif access_gui.peropt.get()==3:
        pass
    elif access_gui.peropt.get()==4:
        pass

def gfACCESS_main(ds_tower,ds_access,access_info):
    '''
    This is the main routine for using ACCESS data to gap fill drivers.
    '''
    # read the control file again, this allows the contents of the control file to
    # be changed with the ACCESS GUI still displayed
    cfname = ds_tower.globalattributes["controlfile_name"]
    cf = qcio.get_controlfilecontents(cfname,mode="quiet")
    # get the site name
    site_name = ds_tower.globalattributes["site_name"]
    # get the time step and a local pointer to the datetime series
    ts = int(ds_tower.globalattributes["time_step"])
    ldt_tower = ds_tower.series["DateTime"]["Data"]
    ldt_access = ds_access.series["DateTime"]["Data"]
    xldt_tower = ds_tower.series["xlDateTime"]["Data"]
    # get the start and end datetimes
    startdate = access_info["startdate"]
    enddate = access_info["enddate"]
    if len(enddate)==0: enddate = access_info["overlap_enddate"]
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
    msg = " Gap filling with ACCESS: "+odt_tower[0].strftime("%Y-%m-%d")+" to "+odt_tower[-1].strftime("%Y-%m-%d")
    log.info(msg)
    # get the number of samples per day and the number of days in the period
    nrecs = ei_tower - si_tower + 1
    nperhr = int(float(60)/ts+0.5)
    nperday = int(float(24)*nperhr+0.5)
    ndays = nrecs/nperday
    # now loop over the variables to be gap filled using the ACCESS data close any open plot windows
    #for i in plt.get_fignums(): plt.close(i)
    if len(plt.get_fignums())!=0:
        for i in plt.get_fignums(): plt.close(i)
    fig_num = 0
    for label in ds_tower.access.keys():
        acc_name = ds_tower.access[label]["access_name"]
        # save the start, mid and end datetimes for later output
        ds_tower.access[label]["results"]["startdate"].append(xldt_tower[si_tower])
        ds_tower.access[label]["results"]["enddate"].append(xldt_tower[ei_tower])
        ds_tower.access[label]["results"]["middate"].append(numpy.average([xldt_tower[si_tower],
                                                                           xldt_tower[ei_tower]]))
        # get the tower data
        if label not in ds_tower.series.keys():
            log.error("gfACCESS_main: series "+label+" not in tower data file")
            continue
        data_tower,flag_tower,attr_tower = qcutils.GetSeriesasMA(ds_tower,label,si=si_tower,ei=ei_tower)
        # get the units
        units_tower = attr_tower["units"]
        # get the output series label
        output = ds_tower.access[label]["output"]
        # increment the figure number
        fig_num = fig_num + 1
        title = site_name+' : Comparison of tower and ACCESS data for '+label
        # get a list of ACCESS variables for this tower variable
        access_var_list = [item for item in ds_access.series.keys() if acc_name in item]
        if label=="Fn":
            access_var_list = [item for item in access_var_list if "lw" not in item]
            access_var_list = [item for item in access_var_list if "sw" not in item]
        # check the series in the ACCESS data
        if len(access_var_list)==0:
            log.error("gfACCESS_main: series "+label+" not in ACCESS data file")
            continue
        # get the ACCESS series that has the highest correlation with the tower data
        r = numpy.ones(len(access_var_list))*float(-9999)
        if numpy.ma.count(data_tower)>access_info["min_points"]:
            for idx,var in enumerate(access_var_list):
                data_access,flag_access,attr_access = qcutils.GetSeriesasMA(ds_access,var,si=si_access,ei=ei_access)
                units_access = attr_access["units"]
                r[idx] = numpy.ma.corrcoef(data_tower,data_access)[0,1]
        maxidx = numpy.ma.argmax(r)
        # save the name of the ACCESS variable that has the highest correlation with the tower data
        var_maxr = access_var_list[maxidx]
        data_access,flag_access,attr_access = qcutils.GetSeriesasMA(ds_access,var_maxr,si=si_access,ei=ei_access)
        # plot the data for this period
        pd = gfACCESS_initplot(site_name=site_name,label=label,var_maxr=var_maxr,
                               fig_num=fig_num,title=title,
                               si_access=si_access,ei_access=ei_access,
                               si_tower=si_tower,ei_tower=ei_tower,
                               startdate=startdate,enddate=enddate,
                               r=r,ts=ts,nrecs=nrecs,ndays=ndays,nperday=nperday)
        gfACCESS_plotdetailed(pd,ds_tower,ds_access,access_info)
        # put the ordinary least-squares adjusted ACCESS data into the output series
        if access_info["overwrite"]==False:
            ind = numpy.where(abs(ds_tower.series[output]["Data"][pd["si_tower"]:pd["ei_tower"]+1]-float(-9999))<c.eps)[0]
            ds_tower.series[output]["Data"][pd["si_tower"]:pd["ei_tower"]+1][ind] = pd["ols"][ind]
        else:
            ds_tower.series[output]["Data"][pd["si_tower"]:pd["ei_tower"]+1] = pd["ols"]
        # make a QC flag for the gap filled data
        # default value is 20 for tower data replaced by ACCESS data
        flag = numpy.ones(len(pd["ols"]))*numpy.int32(20)
        # check for missing ACCESS data (eg no tower data for this period so no OLS statistics)
        ind = numpy.where(pd["ols"]==numpy.float64(-9999))[0]
        # set the QC flag for these times to 21
        flag[ind] = numpy.int32(21)
        # set the flag
        ds_tower.series[output]["Flag"][pd["si_tower"]:pd["ei_tower"]+1] = flag
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

def gfACCESS_plotdetailed(pd,ds_tower,ds_access,access_info):
    """ Needs an effing good refactor!"""
    # get local pointer to the datetime series
    ldt_tower = ds_tower.series["DateTime"]["Data"]
    ldt_access = ds_access.series["DateTime"]["Data"]
    xldt_tower = ds_tower.series["xlDateTime"]["Data"]
    # get the datetime series for the overlap period
    odt_tower = ldt_tower[pd["si_tower"]:pd["ei_tower"]+1]
    odt_access = ldt_access[pd["si_access"]:pd["ei_access"]+1]
    # get the tower and ACCESS data
    data_tower,flag,attr_tower = qcutils.GetSeriesasMA(ds_tower,pd["label"],si=pd["si_tower"],ei=pd["ei_tower"])
    units_tower = attr_tower["units"]
    data_access,flag,attr_access = qcutils.GetSeriesasMA(ds_access,pd["var_maxr"],si=pd["si_access"],ei=pd["ei_access"])
    units_access = attr_access["units"]
    access_uncorr = numpy.ma.array(data_access)
    access_uncorr_2d = numpy.ma.reshape(access_uncorr,[pd["ndays"],pd["nperday"]])
    # get copies of the data for modifying with masked "or"
    tower = numpy.ma.array(data_tower)
    access = numpy.ma.array(data_access)
    # mask both series when either one is missing
    access.mask = numpy.ma.mask_or(access.mask,tower.mask)
    tower.mask = numpy.ma.mask_or(access.mask,tower.mask)
    # get non-masked versions of the data, these are used with the robust statistics module
    access_nm = numpy.ma.compressed(access)
    tower_nm = numpy.ma.compressed(tower)
    # turn on interactive plotting
    plt.ion()
    # create the figure canvas
    fig = plt.figure(pd["fig_num"],figsize=(13,9))
    fig.canvas.set_window_title(pd["label"])
    plt.figtext(0.5,0.96,pd["title"],ha='center',size=16)
    # top row of XY plots
    # correlation coefficients
    rect1 = [0.10,pd["margin_bottom"]+pd["margin_bottom"]+pd["xy_height"],
             pd["xy_width"],pd["xy_height"]]
    ax1 = plt.axes(rect1)
    ind=numpy.arange(len(pd["r"]))
    ax1.bar(ind,pd["r"],0.35)
    ax1.set_ylabel("r")
    ax1.set_xlabel('Grid')
    # lagged correlation
    rect2 = [0.40,pd["margin_bottom"]+pd["margin_bottom"]+pd["xy_height"],
             pd["xy_width"],pd["xy_height"]]
    ax2 = plt.axes(rect2)
    ax2.set_ylabel('r')
    ax2.set_xlabel('Lags')
    if numpy.ma.count(data_tower)>=access_info["min_points"]:
        lag_uncorr = ax2.xcorr(tower_nm,access_nm,maxlags=int(pd["nperday"]/2),color="b")
        # get the lag at maximum correlation as the number of time steps
        pd["lag_uncorr_maxr_nots"] = numpy.argmax(lag_uncorr[1])-int(pd["nperday"]/2)
        # get the lag in minutes
        pd["lag_uncorr_maxr_minutes"] = pd["ts"]*(numpy.argmax(lag_uncorr[1])-int(pd["nperday"]/2))
        ds_tower.access[pd["label"]]["results"]["lag_maxr"].append(pd["lag_uncorr_maxr_minutes"])
        # correct for any lag between the ACCESS data and the tower data
        if pd["lag_uncorr_maxr_minutes"]!=0:
            data_tower,data_access = gfACCESS_lagaccess(pd,ds_tower,ds_access)
            # get copies for use with "or"'d mask
            tower = numpy.ma.array(data_tower)
            access = numpy.ma.array(data_access)
            # mask both series when either one is missing
            access.mask = numpy.ma.mask_or(access.mask,tower.mask)
            tower.mask = numpy.ma.mask_or(access.mask,tower.mask)
            # get non-masked versions of the data, these are used with the robust statistics module
            access_nm = numpy.ma.compressed(access)
            tower_nm = numpy.ma.compressed(tower)
        # plot the lagged correlatioin for the corrected data
        lag_corr = ax2.xcorr(tower_nm,access_nm,maxlags=pd["nperday"]/12,color="r")
        # get the lag in minutes
        pd["lag_corr_maxr_minutes"] = pd["ts"]*(numpy.argmax(lag_corr[1])-pd["nperday"]/12)
    else:
        log.error("Less than 100 points available for series "+pd["label"]+" ...")
        pd["lag_uncorr_maxr_nots"] = float(-9999)
        pd["lag_uncorr_maxr_minutes"] = float(-9999)
        pd["lag_corr_maxr_minutes"] = float(-9999)
    # bottom row of XY plots
    # scatter plot of 30 minute data
    rect3 = [0.10,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax3 = plt.axes(rect3)
    ax3.set_ylabel('Tower ('+units_tower+')')
    ax3.set_xlabel('ACCESS ('+units_access+')')
    ax3.text(0.6,0.075,'30 minutes',fontsize=10,horizontalalignment='left',transform=ax3.transAxes)
    if numpy.ma.count(data_tower)>=access_info["min_points"]:
        ax3.plot(access,tower,'b.')
        resrlm = sm.RLM(tower_nm,sm.add_constant(access_nm,prepend=False),M=sm.robust.norms.TukeyBiweight()).fit()
        pd["m_rlm"] = resrlm.params[0]
        pd["b_rlm"] = resrlm.params[1]
        eqnstr = 'y = %.3fx + %.3f'%(pd["m_rlm"],pd["b_rlm"])
        ax3.plot(access_nm,resrlm.fittedvalues,'r--',linewidth=3)
        ax3.text(0.5,0.915,eqnstr,fontsize=8,horizontalalignment='center',transform=ax3.transAxes,color='red')
        resols = sm.OLS(tower_nm,sm.add_constant(access_nm,prepend=False)).fit()
        pd["m_ols"] = resols.params[0]
        pd["b_ols"] = resols.params[1]
        eqnstr = 'y = %.3fx + %.3f'%(pd["m_ols"],pd["b_ols"])
        ax3.plot(access_nm,resols.fittedvalues,'g--',linewidth=3)
        ax3.text(0.5,0.85,eqnstr,fontsize=8,horizontalalignment='center',transform=ax3.transAxes,color='green')
        pd["ols"] = pd["m_ols"]*data_access+pd["b_ols"]
        pd["rlm"] = pd["m_rlm"]*data_access+pd["b_rlm"]
    else:
        pd["m_rlm"] = float(-9999); pd["b_rlm"] = float(-9999)
        pd["m_ols"] = float(-9999); pd["b_ols"] = float(-9999)
        pd["ols"] = numpy.ma.ones(len(data_access))*numpy.float64(-9999)
        pd["rlm"] = numpy.ma.ones(len(data_access))*numpy.float64(-9999)
    ds_tower.access[pd["label"]]["results"]["m_rlm"].append(pd["m_rlm"])
    ds_tower.access[pd["label"]]["results"]["b_rlm"].append(pd["b_rlm"])
    ds_tower.access[pd["label"]]["results"]["m_ols"].append(pd["m_ols"])
    ds_tower.access[pd["label"]]["results"]["b_ols"].append(pd["b_ols"])
    # scatter plot of daily averages
    rect4 = [0.40,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax4 = plt.axes(rect4)
    ax4.set_ylabel('Tower ('+units_tower+')')
    ax4.set_xlabel('ACCESS ('+units_access+')')
    ax4.text(0.6,0.075,'Daily average',fontsize=10,horizontalalignment='left',transform=ax4.transAxes)
    if numpy.ma.count(data_tower)>=access_info["min_points"]:
        tower_2d = numpy.ma.reshape(tower,[pd["ndays"],pd["nperday"]])
        tower_daily_avg = numpy.ma.average(tower_2d,axis=1)
        access_2d = numpy.ma.reshape(access,[pd["ndays"],pd["nperday"]])
        access_daily_avg = numpy.ma.average(access_2d,axis=1)
        ax4.plot(access_daily_avg,tower_daily_avg,'b.')
        access_daily_avg.mask = numpy.ma.mask_or(access_daily_avg.mask,tower_daily_avg.mask)
        tower_daily_avg.mask = numpy.ma.mask_or(access_daily_avg.mask,tower_daily_avg.mask)
        access_daily_avg_nm = numpy.ma.compressed(access_daily_avg)
        tower_daily_avg_nm = numpy.ma.compressed(tower_daily_avg)
        resrlm = sm.RLM(tower_daily_avg_nm,sm.add_constant(access_daily_avg_nm,prepend=False),M=sm.robust.norms.TukeyBiweight()).fit()
        eqnstr = 'y = %.3fx + %.3f'%(resrlm.params[0],resrlm.params[1])
        ax4.plot(access_daily_avg_nm,resrlm.fittedvalues,'r--',linewidth=3)
        ax4.text(0.5,0.915,eqnstr,fontsize=8,horizontalalignment='center',transform=ax4.transAxes,color='red')
        resrlm = sm.OLS(tower_daily_avg_nm,sm.add_constant(access_daily_avg_nm,prepend=False)).fit()
        eqnstr = 'y = %.3fx + %.3f'%(resrlm.params[0],resrlm.params[1])
        ax4.plot(access_daily_avg_nm,resrlm.fittedvalues,'g--',linewidth=3)
        ax4.text(0.5,0.85,eqnstr,fontsize=8,horizontalalignment='center',transform=ax4.transAxes,color='green')
    # diurnal average plot
    rect5 = [0.70,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax5 = plt.axes(rect5)
    access_hourly_avg = numpy.ma.average(access_uncorr_2d,axis=0)
    ind = numpy.arange(len(access_hourly_avg))*float(pd["ts"])/float(60)
    ax5.plot(ind,access_hourly_avg,'b-',label='ACCESS-A')
    ax5.set_ylabel(pd["label"]+' ('+units_tower+')')
    ax5.set_xlim(0,24)
    ax5.xaxis.set_ticks([0,6,12,18,24])
    ax5.set_xlabel('Hour')
    if numpy.ma.count(data_tower)>=access_info["min_points"]:
        tower_hourly_avg = numpy.ma.average(tower_2d,axis=0)
        access_hourly_rlm = numpy.ma.average(pd["m_rlm"]*access_2d+pd["b_rlm"],axis=0)
        access_hourly_ols = numpy.ma.average(pd["m_ols"]*access_2d+pd["b_ols"],axis=0)
        ind = numpy.arange(len(tower_hourly_avg))*float(pd["ts"])/float(60)
        ax5.plot(ind,tower_hourly_avg,'ro',label='Tower')
        ax5.plot(ind,access_hourly_rlm,'r-',label='ACCESS-A (RLM)')
        ax5.plot(ind,access_hourly_ols,'g-',label='ACCESS-A (OLS)')
    ax5.legend(loc='upper right',frameon=False,prop={'size':8})
    # time series
    rect_ts = [pd["margin_left"],pd["ts_bottom"],pd["ts_width"],pd["ts_height"]]
    axes_ts = plt.axes(rect_ts)
    axes_ts.plot(odt_access,access_uncorr,'b-',label="ACCESS-A")
    axes_ts.set_ylabel(pd["label"]+' ('+units_tower+')')
    if numpy.ma.count(data_tower)>=access_info["min_points"]:
        axes_ts.plot(odt_tower,tower,'ro',label="Tower")
        axes_ts.plot(odt_access,pd["rlm"],'r-',label="ACCESS-A (RLM)")
        axes_ts.plot(odt_access,pd["ols"],'g-',label="ACCESS-A (OLS)")
    axes_ts.legend(loc='upper right',frameon=False,prop={'size':8})
    # now get some statistics
    numpoints = numpy.ma.count(tower)
    ds_tower.access[pd["label"]]["results"]["n"].append(numpoints)
    if numpoints>=access_info["min_points"]:
        diff = tower - access
        bias = numpy.ma.average(diff)
        fb = numpy.ma.average((tower-access)/(0.5*(tower+access)))
        rmse = numpy.ma.sqrt(numpy.ma.average(diff*diff))
        nmse = rmse/(numpy.ma.maximum(tower)-numpy.ma.minimum(tower))
        var_tow = numpy.ma.var(tower)
        var_acc = numpy.ma.var(access)
        var_ratio = var_tow/var_acc
        avg_tow = numpy.ma.average(tower)
        avg_acc = numpy.ma.average(access)
    else:
        bias = float(-9999)
        fb = float(-9999)
        rmse = float(-9999)
        nmse = float(-9999)
        var_tow = float(-9999)
        var_acc = float(-9999)
        var_ratio = float(-9999)
        avg_tow = float(-9999)
        avg_acc = float(-9999)
    ds_tower.access[pd["label"]]["results"]["bias"].append(bias)
    ds_tower.access[pd["label"]]["results"]["fb"].append(fb)
    ds_tower.access[pd["label"]]["results"]["rmse"].append(rmse)
    ds_tower.access[pd["label"]]["results"]["nmse"].append(nmse)
    ds_tower.access[pd["label"]]["results"]["var_tow"].append(var_tow)
    ds_tower.access[pd["label"]]["results"]["var_acc"].append(var_acc)
    ds_tower.access[pd["label"]]["results"]["var_ratio"].append(var_ratio)
    ds_tower.access[pd["label"]]["results"]["avg_tow"].append(avg_tow)
    ds_tower.access[pd["label"]]["results"]["avg_acc"].append(avg_acc)
    r_max = numpy.ma.maximum(pd["r"])
    ds_tower.access[pd["label"]]["results"]["r_max"].append(r_max)
    text_left = 0.675
    num_left = 0.825
    row_bottom = 0.35
    row_space = 0.030
    i = 0
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'Lag (corrected)')
    plt.figtext(num_left,row_posn,'%.4g'%(pd["lag_corr_maxr_minutes"]))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'Lag (uncorrected)')
    plt.figtext(num_left,row_posn,'%.4g'%(pd["lag_uncorr_maxr_minutes"]))
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
    plt.draw()
    # save a hard copy of the plot
    sdt = odt_tower[0].strftime("%Y%m%d")
    edt = odt_tower[-1].strftime("%Y%m%d")
    figname = "plots/"+pd["site_name"].replace(" ","")+"_ACCESS_"+pd["label"]
    figname = figname+"_"+sdt+"_"+edt+'.png'
    fig.savefig(figname,format='png')

def gfACCESS_lagaccess(pd,ds_tower,ds_access):
    # local pointers to the datetime series
    ldt_tower = ds_tower.series["DateTime"]["Data"]
    ldt_access = ds_access.series["DateTime"]["Data"]
    xldt_tower = ds_tower.series["xlDateTime"]["Data"]
    # lag the ACCESS datetime series
    ldt_access_lagcorr = [x+datetime.timedelta(minutes=int(pd["lag_uncorr_maxr_minutes"])) for x in ldt_access]
    # get the indices of the start and end datetimes in the tower and the ACCESS data.
    pd["si_tower"] = qcutils.GetDateIndex(ldt_tower,str(pd["startdate"]),ts=pd["ts"],match="startnextday",default=0)
    pd["ei_tower"] = qcutils.GetDateIndex(ldt_tower,str(pd["enddate"]),ts=pd["ts"],match="endpreviousday",default=-1)
    pd["si_access"] = qcutils.GetDateIndex(ldt_access_lagcorr,str(pd["startdate"]),ts=pd["ts"],match="startnextday",default=0)
    pd["ei_access"] = qcutils.GetDateIndex(ldt_access_lagcorr,str(pd["enddate"]),ts=pd["ts"],match="endpreviousday",default=-1)
    # now make sure the start and end datetime match
    pd["startdate"] = max([ldt_tower[pd["si_tower"]],ldt_access_lagcorr[pd["si_access"]]])
    pd["enddate"] = min([ldt_tower[pd["ei_tower"]],ldt_access_lagcorr[pd["ei_access"]]])
    # and get the start and end indices again in case they didn't
    pd["si_tower"] = qcutils.GetDateIndex(ldt_tower,str(pd["startdate"]),ts=pd["ts"],match="startnextday")
    pd["ei_tower"] = qcutils.GetDateIndex(ldt_tower,str(pd["enddate"]),ts=pd["ts"],match="endpreviousday")
    pd["si_access"] = qcutils.GetDateIndex(ldt_access_lagcorr,str(pd["startdate"]),ts=pd["ts"],match="startnextday")
    pd["ei_access"] = qcutils.GetDateIndex(ldt_access_lagcorr,str(pd["enddate"]),ts=pd["ts"],match="endpreviousday")
    # get the datetime series for the overlap period
    odt_tower = ldt_tower[pd["si_tower"]:pd["ei_tower"]+1]
    odt_access = ldt_access_lagcorr[pd["si_access"]:pd["ei_access"]+1]
    # get the number of samples per day and the number of days in the period
    pd["nrecs"] = pd["ei_tower"] - pd["si_tower"] + 1
    nperhr = int(float(60)/pd["ts"]+0.5)
    pd["nperday"] = int(float(24)*nperhr+0.5)
    pd["ndays"] = pd["nrecs"]/pd["nperday"]
    # save the start, mid and end datetimes for later output
    ds_tower.access[pd["label"]]["results"]["startdate"][-1] = xldt_tower[pd["si_tower"]]
    ds_tower.access[pd["label"]]["results"]["enddate"][-1] = xldt_tower[pd["ei_tower"]]
    ds_tower.access[pd["label"]]["results"]["middate"][-1] = numpy.average([xldt_tower[pd["si_tower"]],
                                                                       xldt_tower[pd["ei_tower"]]])
    # get the data again, this time with the lag corrected
    data_tower,flag,attr_tower = qcutils.GetSeriesasMA(ds_tower,pd["label"],si=pd["si_tower"],ei=pd["ei_tower"])
    data_access,flag,attr_access = qcutils.GetSeriesasMA(ds_access,pd["var_maxr"],si=pd["si_access"],ei=pd["ei_access"])
    return data_tower,data_access

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
    solo_gui.wm_title("SOLO GUI")
    solo_gui.grid()
    # top row
    nrow = 0
    solo_gui.nodesLabel = Tkinter.Label(solo_gui,text="Nodes")
    solo_gui.nodesLabel.grid(row=nrow,column=0,columnspan=1,sticky="E")
    solo_gui.nodesEntry = Tkinter.Entry(solo_gui,width=6)
    solo_gui.nodesEntry.grid(row=nrow,column=1,columnspan=1,sticky="W")
    solo_gui.nodesEntry.insert(0,"Auto")
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
    solo_gui.runButton.grid(row=nrow,column=0,columnspan=3)
    solo_gui.doneButton = Tkinter.Button (solo_gui, text="Done",command=lambda:gfSOLO_finish(dsb,solo_gui))
    solo_gui.doneButton.grid(row=nrow,column=3,columnspan=3)
    solo_gui.wait_window(solo_gui)

def gfSOLO_finish(ds,solo_gui):
    # destroy the SOLO GUI
    solo_gui.destroy()
    # write Excel spreadsheet with fit statistics
    qcio.xl_write_SOLOStats(ds)
    # remove the ACCESS dictionary from the data structure
    del ds.solo

def gfSOLO_getserieslist(cf):
    series_list = []
    for series in cf["Drivers"].keys():
        if "GapFillUsingSOLO" in cf["Drivers"][series]: series_list.append(series)
    for series in cf["Fluxes"].keys():
        if "GapFillUsingSOLO" in cf["Fluxes"][series]: series_list.append(series)
    return series_list

def gfSOLO_main(dsa,dsb,solo_gui):
    '''
    This is the main routine for running SOLO, an artifical neural network for gap filling fluxes.
    '''
    log.info(" Gap filling "+str(dsb.solo.keys())+" using SOLO ANN")
    # read the control file again, this allows the contents of the control file to
    # be changed with the SOLO GUI still displayed
    cfname = dsb.globalattributes["controlfile_name"]
    cf = qcio.get_controlfilecontents(cfname)
    solo_series = gfSOLO_getserieslist(cf)
    #for series in dsb.solo.keys():
    for series in solo_series:
        section = qcutils.get_cfsection(cf,series=series,mode="quiet")
        if len(section)==0: continue
        if series not in dsb.series.keys(): continue
        dsb.solo[series]["drivers"] = ast.literal_eval(cf[section][series]["GapFillUsingSOLO"]["drivers"])
        dsb.solo[series]["output"] = cf[section][series]["GapFillUsingSOLO"]["output"]
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
        print " GapFillUsingSOLO: end datetime index ("+str(ei)+") smaller that start ("+str(si)+")"
        return
    if si==0 and ei==-1:
        print " GapFillUsingSOLO: no start and end datetime specified, using all data"
        nRecs = int(dsb.globalattributes["nc_nrecs"])
    else:
        nRecs = ei - si + 1
    # loop over the series to be gap filled using solo
    fig_num = 0
    for series in solo_series:
        dsb.solo[series]["results"]["startdate"].append(xldt[si])
        dsb.solo[series]["results"]["enddate"].append(xldt[ei])
        dsb.solo[series]["results"]["middate"].append(numpy.average([xldt[si],xldt[ei]]))
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
    # get the number of drivers
    ndrivers = len(driverlist)
    # add an extra column for the target data
    sofminputdata = numpy.zeros((nRecs,ndrivers))
    # now fill the driver data array
    i = 0
    badlines = []
    for TheseOnes in driverlist:
        driver,flag,attr = qcutils.GetSeries(dsb,TheseOnes,si=si,ei=ei)
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
        driver,flag,attr = qcutils.GetSeries(dsb,TheseOnes,si=si,ei=ei)
        soloinputdata[:,i] = driver[:]
        i = i + 1
    # a clean copy of the target is pulled from the unmodified ds each time
    target,flag,attr = qcutils.GetSeries(dsa,targetlabel,si=si,ei=ei)
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
        driver,flag,attr = qcutils.GetSeries(dsb,TheseOnes,si=si,ei=ei)
        seqsoloinputdata[:,i] = driver[:]
        i = i + 1
    # a clean copy of the target is pulled from the unmodified ds each time
    target,flag,attr = qcutils.GetSeries(dsa,targetlabel,si=si,ei=ei)
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
            flux,flag,attr = qcutils.GetSeriesasMA(dsb,outputlabel)
            attr = qcutils.MakeAttributeDictionary(long_name=targetlabel+' modelled by SOLO',
                                                   units=dsa.series[targetlabel]['Attr']['units'])
            qcutils.CreateSeries(dsb,outputlabel,flux,Flag=flag,Attr=attr)
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
    Hdh,f,a = qcutils.GetSeriesasMA(dsb,'Hdh',si=si,ei=ei)
    # get the observed and modelled values
    obs,f,a = qcutils.GetSeriesasMA(dsa,targetlabel,si=si,ei=ei)
    mod,f,a = qcutils.GetSeriesasMA(dsb,outputlabel,si=si,ei=ei)
    # make the figure
    #plt.ion()
    fig = plt.figure(pd["fig_num"],figsize=(13,9))
    fig.clf()
    fig.canvas.set_window_title(targetlabel)
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
    #plt.draw()
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
    #plt.draw()
    for ThisOne,i in zip(driverlist,range(1,pd["nDrivers"]+1)):
        this_bottom = pd["ts_bottom"] + i*pd["ts_height"]
        rect = [pd["margin_left"],this_bottom,pd["ts_width"],pd["ts_height"]]
        ts_axes.append(plt.axes(rect,sharex=ts_axes[0]))
        data,flag,attr = qcutils.GetSeriesasMA(dsb,ThisOne,si=si,ei=ei)
        ts_axes[i].plot(xdt,data)
        plt.setp(ts_axes[i].get_xticklabels(),visible=False)
        TextStr = ThisOne+'('+dsb.series[ThisOne]['Attr']['units']+')'
        ts_axes[i].text(0.05,0.85,TextStr,color='b',horizontalalignment='left',transform=ts_axes[i].transAxes)
        #plt.draw()
    plt.show(block=False)
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