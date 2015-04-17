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
import platform
import qcck
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

def GapFillFromClimatology(ds):
    '''
    Gap fill missing data using data from the climatology spreadsheet produced by
    the climatology.py script.
    '''
    if "climatology" not in dir(ds): return
    # loop over the series to be gap filled using climatology
    cli_xlbooks = {}
    for output in ds.climatology.keys():
        # check to see if there are any gaps in "series"
        #index = numpy.where(abs(ds.series[label]['Data']-float(c.missing_value))<c.eps)[0]
        #if len(index)==0: continue                      # no gaps found in "series"
        cli_filename = ds.climatology[output]["file_name"]
        if not os.path.exists(cli_filename):
            log.error(" GapFillFromClimatology: Climatology file "+cli_filename+" doesn't exist")
            continue
        if cli_filename not in cli_xlbooks: cli_xlbooks[cli_filename] = xlrd.open_workbook(cli_filename)
        # local pointers to the series name and climatology method
        label = ds.climatology[output]["label_tower"]
        method = ds.climatology[output]["method"]
        # do the gap filling
        log.info(" Gap filling "+label+" using climatology")
        # choose the gap filling method
        if method=="monthly":
            gfClimatology_monthly(ds,label,output,cli_xlbooks)
        elif method=="interpolated daily":
            gfClimatology_interpolateddaily(ds,label,output,cli_xlbooks)
        else:
            log.error(" GapFillFromClimatology: unrecognised method option for "+label)
            continue
    if 'GapFillFromClimatology' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', GapFillFromClimatology'
    # remove the "climatology" attribute from ds
    #del ds.climatology

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
    
def gfClimatology_interpolateddaily(ds,series,output,xlbooks):
    # gap fill from interpolated 30 minute data
    xlfilename = ds.climatology[output]["file_name"]
    ts = ds.globalattributes["time_step"]
    ldt = ds.series["DateTime"]["Data"]
    thissheet = xlbooks[xlfilename].sheet_by_name(series+'i(day)')
    datemode = xlbooks[xlfilename].datemode
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
        xldatenumber = int(thissheet.cell_value(xlRow+2,0))
        # convert this to a Python Datetime
        xldatetime = basedate+datetime.timedelta(days=xldatenumber+1462*datemode)
        # fill the climatology datetime array
        cdt[xlRow*nts:(xlRow+1)*nts] = [xldatetime+datetime.timedelta(hours=hh) for hh in tsteps]
        # fill the climatological value array
        val1d[xlRow*nts:(xlRow+1)*nts] = thissheet.row_values(xlRow+2,start_colx=1,end_colx=nts+1)
    # get the data to be filled with climatological values
    data,flag,attr = qcutils.GetSeriesasMA(ds,series)
    # get an index of missing values
    idx = numpy.ma.where(numpy.ma.getmaskarray(data)==True)[0]
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
            data[ii] = numpy.float64(c.missing_value)
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

def gfalternate_createdict(cf,ds,series,ds_alt):
    """
    Purpose:
     Creates a dictionary in ds to hold information about the alternate data used to gap fill the tower data.
    Usage:
    Side effects:
    Author: PRI
    Date: August 2014
    """
    # get the section of the control file containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # return without doing anything if the series isn't in a control file section
    if len(section)==0:
        log.error("GapFillFromAlternate: Series "+series+" not found in control file, skipping ...")
        return
    # create the alternate directory in the data structure
    if "alternate" not in dir(ds): ds.alternate = {}
    # name of alternate output series in ds
    output_list = cf[section][series]["GapFillFromAlternate"].keys()
    # loop over the outputs listed in the control file
    for output in output_list:
        # create the dictionary keys for this output
        ds.alternate[output] = {}
        ds.alternate[output]["label_tower"] = series
        # source name
        ds.alternate[output]["source"] = cf[section][series]["GapFillFromAlternate"][output]["source"]
        # site name
        ds.alternate[output]["site_name"] = ds.globalattributes["site_name"]
        # alternate data file name
        ds.alternate[output]["file_name"] = cf[section][series]["GapFillFromAlternate"][output]["file_name"]
        # if the file has not already been read, do it now
        if ds.alternate[output]["file_name"] not in ds_alt:
            ds_alt[ds.alternate[output]["file_name"]] = qcio.nc_read_series(ds.alternate[output]["file_name"])
        # get the type of fit
        ds.alternate[output]["fit_type"] = "OLS"
        if "fit" in cf[section][series]["GapFillFromAlternate"][output]:
            if cf[section][series]["GapFillFromAlternate"][output]["fit"].lower() in ["ols","mrev","replace"]:
                ds.alternate[output]["fit_type"] = cf[section][series]["GapFillFromAlternate"][output]["fit"]
            else:
                log.info("gfAlternate: unrecognised fit option for series "+output)
        # force the fit through the origin
        ds.alternate[output]["thru0"] = "no"
        if "thru0" in cf[section][series]["GapFillFromAlternate"][output]:
            if cf[section][series]["GapFillFromAlternate"][output]["thru0"].lower() in ["yes","true"]:
                ds.alternate[output]["thru0"] = "yes"
            else:
                log.info("gfAlternate: unrecognised thru0 option for series "+output)
        # correct for lag?
        ds.alternate[output]["lag"] = "yes"
        if "lag" in cf[section][series]["GapFillFromAlternate"][output]:
            if cf[section][series]["GapFillFromAlternate"][output]["lag"].lower() in ["no","false"]:
                ds.alternate[output]["lag"] = "no"
            else:
                log.info("gfAlternate: unrecognised lag option for series "+output)
        # alternate data variable name if different from name used in control file
        if "alternate_name" in cf[section][series]["GapFillFromAlternate"][output]:
            ds.alternate[output]["alternate_name"] = cf[section][series]["GapFillFromAlternate"][output]["alternate_name"]
        else:
            ds.alternate[output]["alternate_name"] = series
        # results of best fit for plotting later on
        ds.alternate[output]["results"] = {"startdate":[],"enddate":[],"No. points":[],"r":[],
                                           "Bias":[],"RMSE":[],"Frac Bias":[],"NMSE":[],
                                           "Avg (tower)":[],"Avg (alternate)":[],
                                           "Var (tower)":[],"Var (alternate)":[],"Var ratio":[],
                                           "Lag (uncorrected)":[],"Lag (corrected)":[],
                                           "Slope":[],"Offset":[]}
        # create an empty series in ds if the alternate output series doesn't exist yet
        if output not in ds.series.keys():
            data,flag,attr = qcutils.MakeEmptySeries(ds,output)
            qcutils.CreateSeries(ds,output,data,Flag=flag,Attr=attr)

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
    # name of alternate output series in ds
    output_list = cf[section][series]["GapFillFromClimatology"].keys()
    # loop over the outputs listed in the control file
    for output in output_list:
        # create the dictionary keys for this output
        ds.climatology[output] = {}
        ds.climatology[output]["label_tower"] = series
        # site name
        ds.climatology[output]["site_name"] = ds.globalattributes["site_name"]
        # Climatology file name
        ds.climatology[output]["file_name"] = cf[section][series]["GapFillFromClimatology"][output]["file_name"]
        # climatology variable name if different from name used in control file
        if "climatology_name" in cf[section][series]["GapFillFromClimatology"][output]:
            ds.climatology[output]["climatology_name"] = cf[section][series]["GapFillFromClimatology"][output]["climatology_name"]
        else:
            ds.climatology[output]["climatology_name"] = series
        # climatology gap filling method
        if "method" not in cf[section][series]["GapFillFromClimatology"][output].keys():
            # default if "method" missing is "interpolated_daily"
            ds.climatology[output]["method"] = "interpolated_daily"
        else:
            ds.climatology[output]["method"] = cf[section][series]["GapFillFromClimatology"][output]["method"]
        # create an empty series in ds if the climatology output series doesn't exist yet
        if output not in ds.series.keys():
            data,flag,attr = qcutils.MakeEmptySeries(ds,output)
            qcutils.CreateSeries(ds,output,data,Flag=flag,Attr=attr)

def gfMergeSeries_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about the merging of gap filled
        and tower data."""
    merge_prereq_list = ["Fsd","Fsu","Fld","Flu","Ts","Sws"]
    # get the section of the control file containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # create the merge directory in the data structure
    if "merge" not in dir(ds): ds.merge = {}
    # check to see if this series is in the "merge first" list
    # series in the "merge first" list get merged first so they can be used with existing tower
    # data to re-calculate Fg, Fn and Fa
    merge_order = "standard"
    if series in merge_prereq_list: merge_order = "prerequisite"
    if merge_order not in ds.merge.keys(): ds.merge[merge_order] = {}
    # create the dictionary keys for this series
    ds.merge[merge_order][series] = {}
    # output series name
    ds.merge[merge_order][series]["output"] = series
    # site name
    ds.merge[merge_order][series]["source"] = ast.literal_eval(cf[section][series]["MergeSeries"]["Source"])
    # create an empty series in ds if the output series doesn't exist yet
    if ds.merge[merge_order][series]["output"] not in ds.series.keys():
        data,flag,attr = qcutils.MakeEmptySeries(ds,ds.merge[merge_order][series]["output"])
        qcutils.CreateSeries(ds,ds.merge[merge_order][series]["output"],data,Flag=flag,Attr=attr)

def gfSOLO_createdict(cf,ds,series):
    """ Creates a dictionary in ds to hold information about the SOLO data used
        to gap fill the tower data."""
    # get the section of the control file containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # return without doing anything if the series isn't in a control file section
    if len(section)==0:
        log.error("GapFillUsingSOLO: Series "+series+" not found in control file, skipping ...")
        return
    # create the solo directory in the data structure
    if "solo" not in dir(ds): ds.solo = {}
    # create the dictionary keys for this series
    ds.solo[series] = {}
    # site name
    ds.solo[series]["site_name"] = ds.globalattributes["site_name"]
    # list of drivers
    ds.solo[series]["drivers"] = ast.literal_eval(cf[section][series]["GapFillUsingSOLO"]["drivers"])
    # name of SOLO output series in ds
    ds.solo[series]["output"] = cf[section][series]["GapFillUsingSOLO"]["output"]
    # results of best fit for plotting later on
    ds.solo[series]["results"] = {"startdate":[],"enddate":[],"No. points":[],"r":[],
                                  "Bias":[],"RMSE":[],"Frac Bias":[],"NMSE":[],
                                  "Avg (obs)":[],"Avg (SOLO)":[],
                                  "Var (obs)":[],"Var (SOLO)":[],"Var ratio":[],
                                  "m_ols":[],"b_ols":[]}

def GapFillParseControlFile(cf,ds,series,ds_alt):
    # find the section containing the series
    section = qcutils.get_cfsection(cf,series=series,mode="quiet")
    # return empty handed if the series is not in a section
    if len(section)==0: return
    if "GapFillFromAlternate" in cf[section][series].keys():
        # create the alternate dictionary in ds
        gfalternate_createdict(cf,ds,series,ds_alt)
    if "GapFillUsingSOLO" in cf[section][series].keys():
        # create the SOLO dictionary in ds
        gfSOLO_createdict(cf,ds,series)
        # create an empty series in ds if the SOLO output series doesn't exist yet
        if ds.solo[series]["output"] not in ds.series.keys():
            data,flag,attr = qcutils.MakeEmptySeries(ds,ds.solo[series]["output"])
            qcutils.CreateSeries(ds,ds.solo[series]["output"],data,Flag=flag,Attr=attr)
    if "GapFillFromClimatology" in cf[section][series].keys():
        # create the climatology dictionary in the data structure
        gfClimatology_createdict(cf,ds,series)
    if "MergeSeries" in cf[section][series].keys():
        # create the merge series dictionary in the data structure
        gfMergeSeries_createdict(cf,ds,series)

def GapFillFromAlternate(ds4,ds_alt):
    '''
    This is the gap fill from alternate data GUI.
    The alternate data gap fill GUI is displayed separately from the main OzFluxQC GUI.
    It consists of text to display the start and end datetime of the file,
    two entry boxes for the start and end datetimes of the alternate data gap fill and
    a button to insert the gap fill data ("Run") and a button to exit ("Done")
    the GUI when we are done.  On exit, the OzFluxQC main GUI continues
    and eventually writes the gap filled data to file.
    '''
    # set the default return code
    ds4.returncodes["alternate"] = "normal"
    if "alternate" not in dir(ds4): return
    # get a local pointer to the tower datetime series
    ldt_tower = ds4.series["DateTime"]["Data"]
    # get the start and end datetime of the tower data
    startdate = ldt_tower[0]
    enddate = ldt_tower[-1]
    alternate_info = {"overlap_startdate":startdate.strftime("%Y-%m-%d %H:%M"),
                    "overlap_enddate":enddate.strftime("%Y-%m-%d %H:%M")}
    # make the GUI
    alt_gui = Tkinter.Toplevel()
    alt_gui.wm_title("Gap fill using alternate data")
    alt_gui.grid()
    # top row
    nrow = 0
    alt_gui.filestartLabel = Tkinter.Label(alt_gui,text="Overlap start date")
    alt_gui.filestartLabel.grid(row=nrow,column=0,columnspan=3)
    alt_gui.fileendLabel = Tkinter.Label(alt_gui,text="Overlap end date")
    alt_gui.fileendLabel.grid(row=nrow,column=3,columnspan=3)
    # second row
    nrow = nrow + 1
    alt_gui.filestartValue = Tkinter.Label(alt_gui,text=alternate_info["overlap_startdate"])
    alt_gui.filestartValue.grid(row=nrow,column=0,columnspan=3)
    alt_gui.fileendValue = Tkinter.Label(alt_gui,text=alternate_info["overlap_enddate"])
    alt_gui.fileendValue.grid(row=nrow,column=3,columnspan=3)
    # third row
    nrow = nrow + 1
    alt_gui.startLabel = Tkinter.Label(alt_gui, text="Start date (YYYY-MM-DD)")
    alt_gui.startLabel.grid(row=nrow,column=0,columnspan=3)
    alt_gui.startEntry = Tkinter.Entry(alt_gui,width=15)
    alt_gui.startEntry.grid(row=nrow,column=3,columnspan=3)
    # fourth row
    nrow = nrow + 1
    alt_gui.endLabel = Tkinter.Label(alt_gui, text="End date   (YYYY-MM-DD)")
    alt_gui.endLabel.grid(row=nrow,column=0,columnspan=3)
    alt_gui.endEntry = Tkinter.Entry(alt_gui,width=15)
    alt_gui.endEntry.grid(row=nrow,column=3,columnspan=3)
    # fifth row
    nrow = nrow + 1
    alt_gui.peropt = Tkinter.IntVar()
    alt_gui.peropt.set(1)
    alt_gui.manualperiod = Tkinter.Radiobutton(alt_gui,text="Manual",variable=alt_gui.peropt,value=1)
    alt_gui.manualperiod.grid(row=nrow,column=0,columnspan=1,sticky="W")
    alt_gui.minptsLabel = Tkinter.Label(alt_gui,text="Min. pts (%)")
    alt_gui.minptsLabel.grid(row=nrow,column=1,columnspan=2,sticky="E")
    alt_gui.minptsEntry = Tkinter.Entry(alt_gui,width=5)
    alt_gui.minptsEntry.grid(row=nrow,column=3,columnspan=2,sticky="W")
    alt_gui.minptsEntry.insert(0,"25")
    # sixth row
    nrow = nrow + 1
    alt_gui.automonthly = Tkinter.Radiobutton(alt_gui,text="Monthly",variable=alt_gui.peropt,value=2)
    alt_gui.automonthly.grid(row=nrow,column=0,columnspan=1,sticky="W")
    alt_gui.daysLabel = Tkinter.Radiobutton(alt_gui,text="No. days",variable=alt_gui.peropt,value=3)
    alt_gui.daysLabel.grid(row=nrow,column=1,columnspan=2,sticky="W")
    alt_gui.daysEntry = Tkinter.Entry(alt_gui,width=5)
    alt_gui.daysEntry.grid(row=nrow,column=3,columnspan=2,sticky="W")
    alt_gui.daysEntry.insert(0,"120")
    # seventh row
    nrow = nrow + 1
    alt_gui.pltopt = Tkinter.IntVar()
    alt_gui.pltopt.set(1)
    alt_gui.showplots = Tkinter.Checkbutton(alt_gui, text="Show plots", variable=alt_gui.pltopt)
    alt_gui.showplots.grid(row=nrow,column=0,columnspan=3,sticky="w")
    alt_gui.owopt = Tkinter.IntVar()
    alt_gui.owopt.set(1)
    alt_gui.overwrite = Tkinter.Checkbutton(alt_gui, text="Overwrite", variable=alt_gui.owopt)
    alt_gui.overwrite.grid(row=nrow,column=3,columnspan=3,sticky="w")
    # eighth row
    nrow = nrow + 1
    alt_gui.doneButton = Tkinter.Button(alt_gui,text="Done",command=lambda:gfalternate_done(ds4,alt_gui))
    alt_gui.doneButton.grid(row=nrow,column=0,columnspan=2)
    alt_gui.runButton = Tkinter.Button(alt_gui,text="Run",command=lambda:gfalternate_run(ds4,ds_alt,alt_gui,alternate_info))
    alt_gui.runButton.grid(row=nrow,column=2,columnspan=2)
    alt_gui.quitButton = Tkinter.Button(alt_gui,text="Quit",command=lambda:gfalternate_quit(ds4,alt_gui))
    alt_gui.quitButton.grid(row=nrow,column=4,columnspan=2)
    # ninth row
    nrow = nrow + 1
    alt_gui.progress_row = nrow
    alt_gui.progress = Tkinter.Label(alt_gui, text='Waiting for input ...')
    alt_gui.progress.grid(row=nrow,column=0,columnspan=6,sticky="W")

    alt_gui.wait_window(alt_gui)

def gfalternate_progress(alt_gui,text):
    """
        Update progress message in alternate GUI
        """
    alt_gui.progress.destroy()
    alt_gui.progress = Tkinter.Label(alt_gui, text=text)
    alt_gui.progress.grid(row=9,column=0,columnspan=6,sticky="W")
    alt_gui.update()

def gfalternate_done(ds,alt_gui):
    # plot the summary statistics
    gfalternate_plotsummary(ds)
    # destroy the alternate GUI
    alt_gui.destroy()
    # write Excel spreadsheet with fit statistics
    qcio.xl_write_AlternateStats(ds)
    # put the return code into ds.alternate
    ds.returncodes["alternate"] = "normal"

def gfalternate_quit(ds,alt_gui):
    # destroy the alternate GUI
    alt_gui.destroy()
    # put the return code into ds.alternate
    ds.returncodes["alternate"] = "quit"

def gfalternate_plotsummary(ds):
    """ Plot single pages of summary results for groups of variables. """
    # get a list of variables for which alternate data is available
    output_list = ds.alternate.keys()
    series_list = [ds.alternate[item]["label_tower"] for item in ds.alternate.keys()]
    if len(ds.alternate[output_list[0]]["results"]["startdate"])==0:
        log.info("gfalternate: no summary data to plot")
        return
    # get the Excel datemode, needed to convert the Excel datetime to Python datetimes
    datemode = int(ds.globalattributes['xl_datemode'])
    # site name for titles
    site_name = ds.globalattributes["site_name"]
    # datetimes are stored in ds.alternate as Excel datetimes, here we convert to Python datetimes
    # for ease of handling and plotting.
    # start datetimes of the periods compared first
    basedate = datetime.datetime(1899, 12, 30)
    dt_start = []
    for xldt in ds.alternate[output_list[0]]["results"]["startdate"]:
        dt_start.append(basedate+datetime.timedelta(days=xldt+1462*datemode))
    startdate = min(dt_start)
    # and then the end datetimes
    dt_end = []
    for xldt in ds.alternate[output_list[0]]["results"]["enddate"]:
        dt_end.append(basedate+datetime.timedelta(days=xldt+1462*datemode))
    enddate = max(dt_end)
    # get the major tick locator and label format
    MTLoc = mdt.AutoDateLocator(minticks=3,maxticks=5)
    MTFmt = mdt.DateFormatter('%b')
    # group lists of the resuts to be plotted
    result_list = ["r","Bias","RMSE","Var ratio","Lag (uncorrected)","Slope","Offset"]
    ylabel_list = ["r","Bias","RMSE","Var ratio","Lag","Slope","Offset"]
    # turn on interactive plotting
    plt.ion()
    # now loop over the group lists
    for nFig in ds.cf["Alternate_Summary"].keys():
        plot_title = ds.cf["Alternate_Summary"][str(nFig)]["Title"]
        var_list = ast.literal_eval(ds.cf["Alternate_Summary"][str(nFig)]["Variables"])
        # set up the subplots on the page
        fig,axs = plt.subplots(len(result_list),len(var_list),figsize=(13,9))
        fig.canvas.set_window_title("Alternate summary: "+plot_title)
        # make a title string for the plot and render it
        title_str = "Alternate: "+plot_title+"; "+site_name+" "+datetime.datetime.strftime(startdate,"%Y-%m-%d")
        title_str = title_str+" to "+datetime.datetime.strftime(enddate,"%Y-%m-%d")
        fig.suptitle(title_str, fontsize=14, fontweight='bold')
        # initialise a string to take the concatenated variable names, used in the name of the hard-copy of the plot
        figlab = ""
        # now loop over the variables in the group list
        for col,output in enumerate(var_list):
            if output not in output_list:
                log.error("Series "+output+" requested for summary plot is not available")
                continue
            source = ds.alternate[output]["source"]
            # append the variable name to the variable name string
            figlab = figlab+output
            # and loop over rows in plot
            for row,rlabel,ylabel in zip(range(len(result_list)),result_list,ylabel_list):
                # if this is the first row, add the column title
                #if row==0: axs[row,col].set_title(output+" ("+source+")")
                if row==0: axs[row,col].set_title(output)
                # if this is the left-most column, add the Y axis labels
                if col==0: axs[row,col].set_ylabel(ylabel,visible=True)
                # get the results to be plotted
                result = numpy.ma.masked_equal(ds.alternate[output]["results"][rlabel],float(c.missing_value))
                if numpy.ma.count(result)==0: result = numpy.ma.ones(len(dt_start),dtype=numpy.float32)*float(c.large_value)
                # put the data into the right order to be plotted
                dt,data = gfalternate_plotsummary_getdata(dt_start,dt_end,result)
                # plot the results
                axs[row,col].plot(dt,data)
                # put in the major ticks
                axs[row,col].xaxis.set_major_locator(MTLoc)
                # if this is not the last row, hide the tick mark labels
                if row<len(result_list)-1: plt.setp(axs[row,col].get_xticklabels(),visible=False)
                # if this is the last row, add the major tick mark and axis labels
                if row==len(result_list)-1:
                    axs[row,col].xaxis.set_major_formatter(MTFmt)
                    axs[row,col].set_xlabel('Month',visible=True)
        # draw the plot
        plt.draw()
        # make the hard-copy file name and save the plot as a PNG file
        sdt = startdate.strftime("%Y%m%d")
        edt = enddate.strftime("%Y%m%d")
        figname = "plots/"+site_name.replace(" ","")+"_Alternate_FitStatistics_"+figlab
        figname = figname+"_"+sdt+"_"+edt+".png"
        fig.savefig(figname,format="png")

def gfalternate_plotsummary_getdata(dt_start,dt_end,result):
    dt = []
    data = []
    for s,e,r in zip(dt_start,dt_end,result):
        dt.append(s)
        data.append(r)
        dt.append(e)
        data.append(r)
    return dt,data

def gfalternate_run(ds_tower,ds_alt,alt_gui,alternate_info):
    # populate the alternate_info dictionary with things that will be useful
    alternate_info["peropt"] = alt_gui.peropt.get()
    alternate_info["overwrite"] = True
    if alt_gui.owopt.get()==0: alternate_info["overwrite"] = False
    alternate_info["show_plots"] = True
    if alt_gui.pltopt.get()==0: alternate_info["show_plots"] = False
    alternate_info["min_percent"] = int(alt_gui.minptsEntry.get())
    alternate_info["site_name"] = ds_tower.globalattributes["site_name"]
    alternate_info["time_step"] = int(ds_tower.globalattributes["time_step"])
    alternate_info["nperhr"] = int(float(60)/alternate_info["time_step"]+0.5)
    alternate_info["nperday"] = int(float(24)*alternate_info["nperhr"]+0.5)
    alternate_info["maxlags"] = int(float(12)*alternate_info["nperhr"]+0.5)
    alternate_info["tower"] = {}
    alternate_info["alternate"] = {}
    series_list = [ds_tower.alternate[item]["label_tower"] for item in ds_tower.alternate.keys()]
    alternate_info["series_list"] = series_list
    log.info(" Gap filling "+str(list(set(series_list)))+" using alternate data")
    if alt_gui.peropt.get()==1:
        gfalternate_progress(alt_gui,"Starting manual run ...")
        # get the start and end datetimes entered in the alternate GUI
        alternate_info["startdate"] = alt_gui.startEntry.get()
        if len(alternate_info["startdate"])==0: alternate_info["startdate"] = alternate_info["overlap_startdate"]
        alternate_info["enddate"] = alt_gui.endEntry.get()
        if len(alternate_info["enddate"])==0: alternate_info["enddate"] = alternate_info["overlap_enddate"]
        gfalternate_main(ds_tower,ds_alt,alternate_info)
        gfalternate_progress(alt_gui,"Finished manual run ...")
    elif alt_gui.peropt.get()==2:
        gfalternate_progress(alt_gui,"Starting auto (monthly) run ...")
        # get the start datetime entered in the alternate GUI
        alternate_info["startdate"] = alt_gui.startEntry.get()
        if len(alternate_info["startdate"])==0: alternate_info["startdate"] = alternate_info["overlap_startdate"]
        startdate = dateutil.parser.parse(alternate_info["startdate"])
        overlap_startdate = dateutil.parser.parse(alternate_info["overlap_startdate"])
        overlap_enddate = dateutil.parser.parse(alternate_info["overlap_enddate"])
        enddate = startdate+dateutil.relativedelta.relativedelta(months=1)
        enddate = min([overlap_enddate,enddate])
        alternate_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<overlap_enddate:
            gfalternate_main(ds_tower,ds_alt,alternate_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(months=1)
            alternate_info["startdate"] = startdate.strftime("%Y-%m-%d")
            alternate_info["enddate"] = enddate.strftime("%Y-%m-%d")
        gfalternate_progress(alt_gui,"Finished auto (monthly) run ...")
    elif alt_gui.peropt.get()==3:
        gfalternate_progress(alt_gui,"Starting auto (days) run ...")
        # get the start datetime entered in the alternate GUI
        alternate_info["startdate"] = alt_gui.startEntry.get()
        if len(alternate_info["startdate"])==0: alternate_info["startdate"] = alternate_info["overlap_startdate"]
        startdate = dateutil.parser.parse(alternate_info["startdate"])
        overlap_startdate = dateutil.parser.parse(alternate_info["overlap_startdate"])
        overlap_enddate = dateutil.parser.parse(alternate_info["overlap_enddate"])
        nDays = int(alt_gui.daysEntry.get())
        enddate = startdate+dateutil.relativedelta.relativedelta(days=nDays)
        enddate = min([overlap_enddate,enddate])
        alternate_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<overlap_enddate:
            gfalternate_main(ds_tower,ds_alt,alternate_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(days=nDays)
            alternate_info["startdate"] = startdate.strftime("%Y-%m-%d")
            alternate_info["enddate"] = enddate.strftime("%Y-%m-%d")
        gfalternate_progress(alt_gui,"Finished auto (days) run ...")
    else:
        log.error("GapFillFromAlternate: unrecognised period option")

def gfalternate_getdateindices(ldt_tower,ldt_alternate,alternate_info,si_match,ei_match):
   #gfalternate_getdateindices(ldt_tower,ldt_alternate,alternate_info,"exact","exact")
    startdate = alternate_info["startdate"]
    enddate = alternate_info["enddate"]
    ts = alternate_info["time_step"]
    # get the indices of the start and end datetimes in the tower and the alternate data.
    si_tower = qcutils.GetDateIndex(ldt_tower,startdate,ts=ts,match=si_match,default=0)
    ei_tower = qcutils.GetDateIndex(ldt_tower,enddate,ts=ts,match=ei_match,default=-1)
    si_alternate = qcutils.GetDateIndex(ldt_alternate,startdate,ts=ts,match=si_match,default=0)
    ei_alternate = qcutils.GetDateIndex(ldt_alternate,enddate,ts=ts,match=ei_match,default=-1)
    # now make sure the start and end datetime match
    sdt = max([ldt_tower[si_tower],ldt_alternate[si_alternate]])
    edt = min([ldt_tower[ei_tower],ldt_alternate[ei_alternate]])
    # and get the start and end indices again in case they didn't
    si_tower = qcutils.GetDateIndex(ldt_tower,str(sdt),ts=ts,match=si_match,default=0)
    ei_tower = qcutils.GetDateIndex(ldt_tower,str(edt),ts=ts,match=ei_match,default=-1)
    si_alternate = qcutils.GetDateIndex(ldt_alternate,str(sdt),ts=ts,match=si_match,default=0)
    ei_alternate = qcutils.GetDateIndex(ldt_alternate,str(edt),ts=ts,match=ei_match,default=-1)
    if ei_tower==-1: ei_tower = len(ldt_tower)-1
    if ei_alternate==-1: ei_alternate = len(ldt_alternate)-1
    return {"si":si_tower,"ei":ei_tower},{"si":si_alternate,"ei":ei_alternate}

def gfalternate_getalternatevaratmaxr(alternate_var_list,data_tower,ds_alternate,alternate_info):
    """ Get the name of the alternate variable that has the highest correlation with the tower data."""
    # local pointers to the start and end indices
    si = alternate_info["alternate"]["exact"]["si"]
    ei = alternate_info["alternate"]["exact"]["ei"]
    # create an array for the correlations and a list for the alternate variables in order of decreasing correlation
    r = numpy.zeros(len(alternate_var_list))
    altvarlist = [None]*len(alternate_var_list)
    # check that the number of tower data points is greater than the minimum
    if numpy.ma.count(data_tower)>alternate_info["min_points"]:
        # loop over the variables in the alternate file
        for idx,var in enumerate(alternate_var_list):
            # get the alternate data
            data_alternate,flag,attr = qcutils.GetSeriesasMA(ds_alternate,var,si=si,ei=ei)
            if numpy.ma.count(data_alternate)>alternate_info["min_points"]:
                # check the lengths of the tower and alternate data are the same
                if len(data_alternate)!=len(data_tower):
                    msg = "gfalternate_getalternatevaratmaxr: alternate data length is "+str(len(data_alternate))
                    log.info(msg)
                    msg = "gfalternate_getalternatevaratmaxr: tower data length is "+str(len(data_tower))
                    log.info(msg)
                    raise ValueError('gfalternate_getalternatevaratmaxr: data_tower and data_alternate lengths differ')
                # put the correlation into the r array
                rval = numpy.ma.corrcoef(data_tower,data_alternate)[0,1]
                if rval=="nan": rval = float(0)
            else:
                #msg = "gfalternate_getalternatevaratmaxr: less than "+str(alternate_info["min_points"])+" for "+var
                #log.info(msg)
                rval = float(0)
            r[idx] = rval
            altvarlist[idx] = var
    ## get the index of the maximum r value
    #maxidx = numpy.ma.argmax(numpy.abs(r))
    # save the correlation array for later plotting
    alternate_info["r"] = r
    # sort the correlation array and the alternate variable list
    idx = numpy.flipud(numpy.argsort(r))
    altvar_list_sorted = [alternate_var_list[j] for j in list(idx)]
    # return the name of the alternate variable that has the highest correlation with the tower data
    return altvar_list_sorted

def gfalternate_getalternatevarlist(ds_alternate,label):
    alternate_var_list = [item for item in ds_alternate.series.keys() if label in item]
    # remove any extraneous Fn labels (alternate has Fn_lw and Fn_sw)
    if label=="Fn":
        alternate_var_list = [item for item in alternate_var_list if "lw" not in item]
        alternate_var_list = [item for item in alternate_var_list if "sw" not in item]
    # check the series in the alternate data
    if len(alternate_var_list)==0:
        print "gfalternate_getalternatevarlist: series "+label+" not in alternate data file"
    return alternate_var_list

def gfalternate_getlag(data_alternate,data_tower,ai):
    if ai["lag"].lower()=="yes":
        lags,corr = qcts.get_laggedcorrelation(data_tower,data_alternate,maxlags=ai["maxlags"],minpoints=ai["min_points"])
        nLags = numpy.argmax(corr)-ai["maxlags"]
    else:
        nLags = 0
    return nLags

def gfalternate_getmrevcorrected(data_plot):
    """
    Fit alternate data to tower data by replacing means and equalising variance.
    """
    results = {}
    # local copies of the data
    data_tower = numpy.ma.copy(data_plot["data_tower"])
    data_alternate = numpy.ma.copy(data_plot["data_alternate"])
    data_twr_hravg = numpy.ma.copy(data_plot["data_tower_hourlyavg"])
    data_alt_hravg = numpy.ma.copy(data_plot["data_alternate_lagcorr_hourlyavg"])
    # calculate the means
    mean_tower = numpy.ma.mean(data_tower)
    mean_alternate = numpy.ma.mean(data_alternate)
    # calculate the variances
    var_twr_hravg = numpy.ma.var(data_twr_hravg)
    var_alt_hravg = numpy.ma.var(data_alt_hravg)
    var_ratio = var_twr_hravg/var_alt_hravg
    # correct the alternate data
    results["fit"] = ((data_alternate - mean_alternate)*var_ratio) + mean_tower
    results["eqnstr"] = "Mean replaced, equal variance"
    results["slope"] = float(0); results["offset"] = float(0)
    return results
    
def gfalternate_getolscorrecteddata(x_in,y_in,alternate_info):
    """
    Calculate the ordinary least squares fit between 2 1D arrays.
    """
    # create results dictionary
    results = {}
    # check that both inputs are either ndarrays or masked arrays
    if numpy.ma.isMA(x_in)!=numpy.ma.isMA(y_in):
        log.error('qcts.getolscorrecteddata: one of x or y is a masked array, the other is not')
        results["eqnstr"] = "X & Y different array types"
        results["slope"] = float(0); results["offset"] = float(0)
        results["fit"] = numpy.ma.copy(x_in)
        return results
    # condition or copy the data
    if numpy.ma.isMA(x_in) and numpy.ma.isMA(y_in):
        mask = numpy.ma.mask_or(x_in.mask,y_in.mask)
        x = numpy.ma.compressed(numpy.ma.array(x_in,mask=mask,copy=True))
        y = numpy.ma.compressed(numpy.ma.array(y_in,mask=mask,copy=True))
    else:
        x = numpy.copy(x_in)
        y = numpy.copy(y_in)
    # get the array lengths
    nx = len(x); ny = len(y)
    # check the input array lengths are equal
    if nx!=ny:
        log.error('qcts.getolscorrecteddata: x and y must be equal length')
        results["eqnstr"] = "X & Y unequal lengths"
        results["slope"] = float(0); results["offset"] = float(0)
        results["fit"] = numpy.ma.copy(x_in)
        return results
    # attempt an OLS fit
    min_points = int(nx*alternate_info["min_percent"]/100)
    if nx>=min_points:
        if alternate_info["thru0"].lower()=="yes":
            resols = sm.OLS(y,x).fit()
            results["fit"] = resols.params[0]*x_in
            results["eqnstr"] = 'y = %.3fx'%(resols.params[0])
            results["slope"] = resols.params[0]; results["offset"] = float(0)
        else:
            resols = sm.OLS(y,sm.add_constant(x,prepend=False)).fit()
            if resols.params.shape[0]==2:
                results["fit"] = resols.params[0]*x_in+resols.params[1]
                results["eqnstr"] = 'y = %.3fx + %.3f'%(resols.params[0],resols.params[1])
                results["slope"] = resols.params[0]; results["offset"] = resols.params[1]
            else:
                results["slope"] = float(0); results["offset"] = float(0)
                results["eqnstr"] = "OLS error, replaced"
                results["fit"] = numpy.ma.copy(x_in)
    else:
        results["slope"] = float(0); results["offset"] = float(0)
        results["eqnstr"] = "Too few points, replaced"
        results["fit"] = numpy.ma.copy(x_in)
    return results

def gfalternate_getreplacedata(x_in):
    results = {}
    results["slope"] = float(1); results["offset"] = float(0)
    results["eqnstr"] = "No OLS, replaced"
    results["fit"] = numpy.ma.copy(x_in)
    return results

def gfalternate_getdataas2d(ldt,odt,data,inds,alternate_info):
    si_wholedays = odt.index(ldt[inds["si"]])
    ei_wholedays = odt.index(ldt[inds["ei"]])
    nperday = alternate_info["nperday"]
    data_wholedays = data[si_wholedays:ei_wholedays+1]
    ndays = len(data_wholedays)/nperday
    return numpy.ma.reshape(data_wholedays,[ndays,nperday])

def gfalternate_updatedict(cf,ds_tower,ds_alt):
    """
    Update the ds.alternate dictionary.  This is done after reading the control file so
    that the user can make changes to the control file while the gap fill GUI is still
    displayed and the re-run the gap filling.  This gives a measure of interactive-like
    behaviour to the gap filling process.
    """
    if "alternate" not in dir(ds_tower): return
    section = "Drivers"
    series_list = cf[section].keys()
    for series in series_list:
        if "GapFillFromAlternate" not in cf[section][series].keys(): continue
        # name of alternate output series in ds
        output_list = cf[section][series]["GapFillFromAlternate"].keys()
        # loop over the outputs listed in the control file
        for output in output_list:
            if output not in ds_tower.alternate.keys(): ds_tower.alternate[output] = {}
            ds_tower.alternate[output]["label_tower"] = series
            # source name
            ds_tower.alternate[output]["source"] = cf[section][series]["GapFillFromAlternate"][output]["source"]
            # site name
            ds_tower.alternate[output]["site_name"] = ds_tower.globalattributes["site_name"]
            # alternate data file name
            ds_tower.alternate[output]["file_name"] = cf[section][series]["GapFillFromAlternate"][output]["file_name"]
            # if the file has not already been read, do it now
            if ds_tower.alternate[output]["file_name"] not in ds_alt:
                ds_alt[ds_tower.alternate[output]["file_name"]] = qcio.nc_read_series(ds_tower.alternate[output]["file_name"])
            # get the type of fit
            ds_tower.alternate[output]["fit_type"] = "OLS"
            if "fit" in cf[section][series]["GapFillFromAlternate"][output]:
                if cf[section][series]["GapFillFromAlternate"][output]["fit"].lower() in ["ols","mrev","replace"]:
                    ds_tower.alternate[output]["fit_type"] = cf[section][series]["GapFillFromAlternate"][output]["fit"]
                else:
                    log.info("gfAlternate: unrecognised fit option for series "+output)
            # force the fit through the origin
            ds_tower.alternate[output]["thru0"] = "no"
            if "thru0" in cf[section][series]["GapFillFromAlternate"][output]:
                if cf[section][series]["GapFillFromAlternate"][output]["thru0"].lower() in ["yes","true"]:
                    ds_tower.alternate[output]["thru0"] = "yes"
                else:
                    log.info("gfAlternate: unrecognised thru0 option for series "+output)
            # correct for lag?
            ds_tower.alternate[output]["lag"] = "yes"
            if "lag" in cf[section][series]["GapFillFromAlternate"][output]:
                if cf[section][series]["GapFillFromAlternate"][output]["lag"].lower() in ["no","false"]:
                    ds_tower.alternate[output]["lag"] = "no"
                else:
                    log.info("gfAlternate: unrecognised lag option for series "+output)
            # alternate data variable name if different from name used in control file
            if "alternate_name" in cf[section][series]["GapFillFromAlternate"][output]:
                ds_tower.alternate[output]["alternate_name"] = cf[section][series]["GapFillFromAlternate"][output]["alternate_name"]
            else:
                ds_tower.alternate[output]["alternate_name"] = series
            # results of best fit for plotting later on
            if "results" not in ds_tower.alternate[output].keys():
                ds_tower.alternate[output]["results"] = {"startdate":[],"enddate":[],"No. points":[],"r":[],
                                                         "Bias":[],"RMSE":[],"Frac Bias":[],"NMSE":[],
                                                         "Avg (tower)":[],"Avg (alternate)":[],
                                                         "Var (tower)":[],"Var (alternate)":[],"Var ratio":[],
                                                         "Lag (uncorrected)":[],"Lag (corrected)":[],
                                                         "Slope":[],"Offset":[]}
            # create an empty series in ds if the alternate output series doesn't exist yet
            if output not in ds_tower.series.keys():
                data,flag,attr = qcutils.MakeEmptySeries(ds_tower,output)
                qcutils.CreateSeries(ds_tower,output,data,Flag=flag,Attr=attr)

def gfalternate_main(ds_tower,ds_alt,alternate_info):
    '''
    This is the main routine for using alternate data to gap fill drivers.
    '''
    startdate = alternate_info["startdate"]
    enddate = alternate_info["enddate"]
    log.info(" Gap filling using alternate data: "+startdate+" to "+enddate)
    # close any open plot windows
    if len(plt.get_fignums())!=0:
        for i in plt.get_fignums(): plt.close(i)
    # read the control file again, this allows the contents of the control file to
    # be changed with the alternate GUI still displayed
    cfname = ds_tower.globalattributes["controlfile_name"]
    cf = qcio.get_controlfilecontents(cfname,mode="quiet")
    # do any QC checks
    qcck.do_qcchecks(cf,ds_tower,mode="quiet")
    # update the ds.alternate dictionary
    gfalternate_updatedict(cf,ds_tower,ds_alt)
    # get local pointer to the datetime series
    ldt_tower = ds_tower.series["DateTime"]["Data"]
    xldt_tower = ds_tower.series["xlDateTime"]["Data"]
    # now loop over the variables to be gap filled using the alternate data
    fig_num = 0
    for output in ds_tower.alternate.keys():
        label_tower = ds_tower.alternate[output]["label_tower"]
        # get a local pointer to the alternate data structure
        ds_alternate = ds_alt[ds_tower.alternate[output]["file_name"]]
        ldt_alternate = ds_alternate.series["DateTime"]["Data"]
        # get the indices of the start and end datetimes
        tower_exact,alternate_exact = gfalternate_getdateindices(ldt_tower,ldt_alternate,alternate_info,"exact","exact")
        alternate_info["tower"]["exact"] = tower_exact
        alternate_info["alternate"]["exact"] = alternate_exact
        si = alternate_info["tower"]["exact"]["si"]
        ei = alternate_info["tower"]["exact"]["ei"]
        # get local pointers for the datetime series for the overlap period
        odt_tower = ldt_tower[tower_exact["si"]:tower_exact["ei"]+1]
        odt_alternate = ldt_alternate[alternate_exact["si"]:alternate_exact["ei"]+1]
        # get the tower data
        data_tower,flag,attr = qcutils.GetSeriesasMA(ds_tower,label_tower,si=tower_exact["si"],ei=tower_exact["ei"])
        alternate_info["min_points"] = int(len(data_tower)*alternate_info["min_percent"]/100)
        # get a copy of the tower data so we can track which gaps have been filled
        data_tower_filled = numpy.ma.copy(data_tower)
        # skip tower variable if there are not enough points
        if numpy.ma.count(data_tower)<alternate_info["min_points"] and ds_tower.alternate[output]["fit_type"]!="replace":
            msg = " Less than "+str(alternate_info["min_percent"])+" % data in tower series "+label_tower+", skipping ..."
            log.info(msg)
            continue
        # save the start and end datetimes for later output
        ds_tower.alternate[output]["results"]["startdate"].append(xldt_tower[tower_exact["si"]])
        ds_tower.alternate[output]["results"]["enddate"].append(xldt_tower[tower_exact["ei"]])
        # put the source for this series in alternate_info for use in plotting
        alternate_info["source"] = ds_tower.alternate[output]["source"]
        alternate_info["fit_type"] = ds_tower.alternate[output]["fit_type"]
        alternate_info["lag"] = ds_tower.alternate[output]["lag"]
        alternate_info["thru0"] = ds_tower.alternate[output]["thru0"]
        alternate_info["label_tower"] = ds_tower.alternate[output]["label_tower"]
        alternate_info["label_alternate"] = output
        # get a list of alternate variables for this tower variable
        alternate_name = ds_tower.alternate[output]["alternate_name"]
        alternate_var_list = gfalternate_getalternatevarlist(ds_alternate,alternate_name)
        # get the alternate series that has the highest correlation with the tower data
        alternate_var_list = gfalternate_getalternatevaratmaxr(alternate_var_list,data_tower,ds_alternate,alternate_info)
        # loop over alternate variables
        if alternate_info["source"].lower()=="access": alternate_var_list = alternate_var_list[0:1]
        for label_alternate in alternate_var_list:
            # get the raw alternate data
            data_alternate,flag,attr = qcutils.GetSeriesasMA(ds_alternate,label_alternate,
                                                             si=alternate_exact["si"],ei=alternate_exact["ei"])
            # skip this alternate variable if there are not enough points
            if numpy.ma.count(data_alternate)<alternate_info["min_points"]: continue
            # check to see if this alternate series will fill any existing gaps
            # indices of gaps in filled tower data
            ind_tower = numpy.ma.where(numpy.ma.getmaskarray(data_tower_filled)==True)[0]
            # indices of good data in alternate series
            ind_alternate = numpy.ma.where(numpy.ma.getmaskarray(data_alternate)==False)[0]
            # boolean array with True where an element of ind_tower appears in ind_alternate
            ind_both = numpy.ma.in1d(ind_tower,ind_alternate)
            if not numpy.ma.any(ind_both): continue
            # clear any existing bom_id in the alternate_info dictionary
            alternate_info.pop("bom_id",None)
            # if this variable has a bom_id attribute, put it in alternate_info
            if "bom_id" in attr.keys(): alternate_info["bom_id"] = attr["bom_id"]
            # load up the data_plot dictionary
            data_plot = {"odt_tower":odt_tower,"data_tower":data_tower,"units_tower":attr["units"],
                         "odt_alternate":odt_alternate,"data_alternate":data_alternate,"units_alternate":attr["units"]}
            # correct for lag in the alternate data if required
            nLags = gfalternate_getlag(data_alternate, data_tower, alternate_info)
            ds_tower.alternate[output]["results"]["Lag (uncorrected)"].append(numpy.float64(nLags*alternate_info["time_step"]))
            si_alternate_lagcorr = alternate_info["alternate"]["exact"]["si"] - nLags
            ei_alternate_lagcorr = alternate_info["alternate"]["exact"]["ei"] - nLags
            data_alternate_lagcorr,f,a = qcutils.GetSeriesasMA(ds_alternate,label_alternate,si=si_alternate_lagcorr,ei=ei_alternate_lagcorr,mode="pad")
            nLags = gfalternate_getlag(data_alternate_lagcorr, data_tower, alternate_info)
            ds_tower.alternate[output]["results"]["Lag (corrected)"].append(numpy.float64(nLags*alternate_info["time_step"]))
            # put the data into the data_plot dictionary
            data_plot["data_alternate_lagcorr"] = data_alternate_lagcorr
            # get daily and diurnal averages of data so far
            tower_wholedays,alternate_wholedays = gfalternate_getdateindices(ldt_tower,ldt_alternate,alternate_info,"startnextday","endpreviousday")
            data_tower_2d = gfalternate_getdataas2d(ldt_tower,odt_tower,data_tower,tower_wholedays,alternate_info)
            data_alternate_2d = gfalternate_getdataas2d(ldt_alternate,odt_alternate,data_alternate,alternate_wholedays,alternate_info)
            data_plot["data_tower_dailyavg"] = numpy.ma.average(data_tower_2d,axis=1)
            data_plot["data_tower_hourlyavg"] = numpy.ma.average(data_tower_2d,axis=0)
            data_plot["data_alternate_dailyavg"] = numpy.ma.average(data_alternate_2d,axis=1)
            data_plot["data_alternate_hourlyavg"] = numpy.ma.average(data_alternate_2d,axis=0)
            data_alternate_lagcorr_2d = gfalternate_getdataas2d(ldt_alternate,odt_alternate,data_alternate_lagcorr,alternate_wholedays,alternate_info)
            data_plot["data_alternate_lagcorr_dailyavg"] = numpy.ma.average(data_alternate_lagcorr_2d,axis=1)
            data_plot["data_alternate_lagcorr_hourlyavg"] = numpy.ma.average(data_alternate_lagcorr_2d,axis=0)
            # best fit to tower data using OLS, MREV or none
            results = gfalternate_getfitcorr(data_plot,alternate_info)
            data_plot["data_alternate_lagfitcorr"] = results["fit"]
            ds_tower.alternate[output]["results"]["Slope"].append(numpy.float64(results["slope"]))
            ds_tower.alternate[output]["results"]["Offset"].append(numpy.float64(results["offset"]))
            # get the daily and diurnal averages of the fitted data
            data_alternate_lagfitcorr_2d = gfalternate_getdataas2d(ldt_alternate,odt_alternate,data_plot["data_alternate_lagfitcorr"],alternate_wholedays,alternate_info)
            data_plot["data_alternate_lagfitcorr_dailyavg"] = numpy.ma.average(data_alternate_lagfitcorr_2d,axis=1)
            data_plot["data_alternate_lagfitcorr_hourlyavg"] = numpy.ma.average(data_alternate_lagfitcorr_2d,axis=0)
            # get the comparison statistics
            gfalternate_getstatistics(ds_tower.alternate[output]["results"],data_tower,data_alternate,alternate_info)
            data_plot["results"] = ds_tower.alternate[output]["results"]
            # plot the data for this period
            pd = gfalternate_initplot()
            gfalternate_plotdetailed(fig_num,label_tower,data_plot,alternate_info,pd)
            fig_num = fig_num + 1
            # put the ordinary least-squares adjusted alternate data into the output series
            if alternate_info["overwrite"]==False:
                ind = numpy.where((abs(ds_tower.series[output]["Data"][si:ei+1]-float(c.missing_value))<c.eps)|
                                  (numpy.ma.getmaskarray(data_plot["data_alternate_lagfitcorr"])==False))[0]
                data_tower_filled[ind] = data_plot["data_alternate_lagfitcorr"][ind]
                ds_tower.series[output]["Data"][si:ei+1][ind] = numpy.ma.filled(data_plot["data_alternate_lagfitcorr"][ind],c.missing_value)
                ds_tower.series[output]["Flag"][si:ei+1][ind] = numpy.int32(20)
            else:
                # only do the overwrite where an alternate exists (don't overwrite with missing data)
                ind = numpy.ma.where(numpy.ma.getmaskarray(data_plot["data_alternate_lagfitcorr"])==False)[0]
                data_tower_filled[ind] = data_plot["data_alternate_lagfitcorr"][ind]
                ds_tower.series[output]["Data"][si:ei+1][ind] = numpy.ma.filled(data_plot["data_alternate_lagfitcorr"][ind],c.missing_value)
                ds_tower.series[output]["Flag"][si:ei+1][ind] = numpy.int32(20)
            # check to see if we have alternate data for this whole period, if so there is no reason to continue
            ind_tower = numpy.where(abs(ds_tower.series[output]["Data"][si:ei+1]-float(c.missing_value))<c.eps)[0]
            if len(ind_tower)==0: break
    # make sure this processing step gets written to the global attribute "Functions"
    if "GapFillFromalternate" not in ds_tower.globalattributes["Functions"]:
        ds_tower.globalattributes["Functions"] = ds_tower.globalattributes["Functions"]+", GapFillFromalternate"

def gfalternate_getfitcorr(data_plot,alternate_info):
    """
    Wrapper for the various methods of fitting the alternate data to the tower data.
    """
    if alternate_info["fit_type"].lower() not in ["ols","mrev","replace"]:
        msg = " Unrecognised fit option for "+alternate_info["label_tower"]+", using OLS ..."
        log.info(msg)
        alternate_info["fit_type"] = "ols"
    if alternate_info["fit_type"].lower()=="ols":
        data_tower = data_plot["data_tower"]
        data_alternate = data_plot["data_alternate_lagcorr"]
        results = gfalternate_getolscorrecteddata(data_alternate,data_tower,alternate_info)
    if alternate_info["fit_type"].lower()=="mrev":
        results = gfalternate_getmrevcorrected(data_plot,alternate_info)
    if alternate_info["fit_type"].lower()=="replace":
        data_tower = data_plot["data_tower"]
        data_alternate = data_plot["data_alternate_lagcorr"]
        results = gfalternate_getreplacedata(data_alternate)
    return results

def gfalternate_initplot(**kwargs):
    pd = {"margin_bottom":0.05,"margin_top":0.05,"margin_left":0.075,"margin_right":0.05,
          "xy_height":0.25,"xy_width":0.20,"xyts_space":0.05,"xyxy_space":0.05,"ts_width":0.9,
          "text_left":0.675,"num_left":0.825,"row_bottom":0.35,"row_space":0.030}
   # calculate bottom of the first time series and the height of the time series plots
    pd["ts_bottom"] = pd["margin_bottom"]+pd["xy_height"]+pd["xyxy_space"]+pd["xy_height"]+pd["xyts_space"]
    pd["ts_height"] = (1.0 - pd["margin_top"] - pd["ts_bottom"])
    for key, value in kwargs.iteritems():
        pd[key] = value
    return pd

def gfalternate_getdatapairasnonmasked(data_tower,data_alternate):
    # get local copies of the data so we do not modify the originals
    tower = numpy.ma.copy(data_tower)
    alternate = numpy.ma.copy(data_alternate)
    # mask both series when either one is missing
    alternate.mask = numpy.ma.mask_or(alternate.mask,tower.mask)
    tower.mask = numpy.ma.mask_or(alternate.mask,tower.mask)
    # get non-masked versions of the data, these are used with the robust statistics module
    alternate_nm = numpy.ma.compressed(alternate)
    tower_nm = numpy.ma.compressed(tower)
    return tower_nm,alternate_nm

def gfalternate_getstatistics(results,tower,alternate,alternate_info):
    npts = numpy.ma.count(tower)
    results["No. points"].append(numpy.float64(npts))
    if npts>=alternate_info["min_points"]:
        results["r"].append(numpy.ma.maximum(alternate_info["r"]))
        diff = tower-alternate
        results["Bias"].append(numpy.ma.average(diff))
        results["Frac Bias"].append(numpy.ma.average((diff)/(0.5*(tower+alternate))))
        rmse = numpy.ma.sqrt(numpy.ma.average(diff*diff))
        results["RMSE"].append(rmse)
        results["NMSE"].append(rmse/(numpy.ma.maximum(tower)-numpy.ma.minimum(tower)))
        var_tow = numpy.ma.var(tower)
        results["Var (tower)"].append(var_tow)
        var_acc = numpy.ma.var(alternate)
        results["Var (alternate)"].append(var_acc)
        if var_acc!=0:
            results["Var ratio"].append(var_tow/var_acc)
        else:
            results["Var ratio"].append(c.missing_value)
        results["Avg (tower)"].append(numpy.ma.average(tower))
        results["Avg (alternate)"].append(numpy.ma.average(alternate))
    else:
        results_list = results.keys()
        for item in ["startdate","enddate","No. points"]:
            if item in results_list: results_list.remove(item)
        for item in results_list:
            results[item].append(float(c.missing_value))

def gfalternate_plotdetailed(nfig,label,data_plot,alternate_info,pd):
    # set up some local pointers
    source = alternate_info["source"]
    if "bom_id" in alternate_info.keys(): source = source+" ("+str(alternate_info["bom_id"])+")"
    # turn on interactive plotting
    if alternate_info["show_plots"]: plt.ion()
    # create the figure canvas
    fig = plt.figure(nfig,figsize=(13,9))
    fig.canvas.set_window_title(label)
    # get the plot title string
    title = alternate_info["site_name"]+" : Comparison of tower and "+source+" data for "+label
    plt.figtext(0.5,0.96,title,ha='center',size=16)
    # top row of XY plots: correlation coefficients
    rect1 = [0.10,pd["margin_bottom"]+pd["margin_bottom"]+pd["xy_height"],pd["xy_width"],pd["xy_height"]]
    ax1 = plt.axes(rect1)
    ind=numpy.arange(len(alternate_info["r"]))
    ax1.bar(ind,alternate_info["r"],0.35)
    ax1.set_ylabel("r")
    ax1.set_xlabel("Grid")
    # top row of XY plots: lagged correlation
    rect2 = [0.40,pd["margin_bottom"]+pd["margin_bottom"]+pd["xy_height"],pd["xy_width"],pd["xy_height"]]
    ax2 = plt.axes(rect2)
    ax2.set_ylabel("r")
    ax2.set_xlabel("Lags")
    # bottom row of XY plots: scatter plot of 30 minute data
    rect3 = [0.10,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax3 = plt.axes(rect3)
    ax3.set_ylabel("Tower ("+data_plot["units_tower"]+")")
    ax3.set_xlabel(source+" ("+data_plot["units_alternate"]+")")
    ax3.text(0.6,0.075,"30 minutes",fontsize=10,horizontalalignment="left",transform=ax3.transAxes)
    # bottom row of XY plots: scatter plot of daily averages
    rect4 = [0.40,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax4 = plt.axes(rect4)
    ax4.set_ylabel("Tower ("+data_plot["units_tower"]+")")
    ax4.set_xlabel(source+" ("+data_plot["units_alternate"]+")")
    ax4.text(0.6,0.075,"Daily average",fontsize=10,horizontalalignment="left",transform=ax4.transAxes)
    # bottom row of XY plots: diurnal average plot
    rect5 = [0.70,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax5 = plt.axes(rect5)
    #ind = numpy.arange(len(alternate_hourly_avg))*float(pd["ts"])/float(60)
    ind = numpy.arange(alternate_info["nperday"])/float(alternate_info["nperhr"])
    ax5.plot(ind,data_plot["data_alternate_hourlyavg"],'b-',label=source)
    ax5.set_ylabel(label+' ('+data_plot["units_tower"]+')')
    ax5.set_xlim(0,24)
    ax5.xaxis.set_ticks([0,6,12,18,24])
    ax5.set_xlabel('Hour')
    # top row: time series
    rect6 = [pd["margin_left"],pd["ts_bottom"],pd["ts_width"],pd["ts_height"]]
    ax6 = plt.axes(rect6)
    ax6.plot(data_plot["odt_alternate"],data_plot["data_alternate"],'b-',label=source)
    ax6.set_ylabel(label+' ('+data_plot["units_tower"]+')')
    # now plot the data if there are more than the minimum number of points
    if numpy.ma.count(data_plot["data_tower"])>=alternate_info["min_points"]:
        # lagged correlations
        tower_nm,alternate_nm = gfalternate_getdatapairasnonmasked(data_plot["data_tower"],data_plot["data_alternate"])
        lag_uncorr = ax2.xcorr(tower_nm,alternate_nm,maxlags=alternate_info["maxlags"],color="b")
        tower_nm,alternate_lagcorr_nm = gfalternate_getdatapairasnonmasked(data_plot["data_tower"],data_plot["data_alternate_lagcorr"])
        lag_corr = ax2.xcorr(tower_nm,alternate_lagcorr_nm,maxlags=alternate_info["maxlags"]/4,color="r")
        # scatter plot of 30 minte data
        ax3.plot(data_plot["data_alternate_lagcorr"],data_plot["data_tower"],'b.')
        results = gfalternate_getolscorrecteddata(data_plot["data_alternate_lagcorr"],data_plot["data_tower"],alternate_info)
        ax3.plot(data_plot["data_alternate_lagcorr"],results["fit"],'g--',linewidth=3)
        ax3.text(0.5,0.9,results["eqnstr"],fontsize=8,horizontalalignment='center',transform=ax3.transAxes,color='green')
        # daily average plot
        ax4.plot(data_plot["data_alternate_lagcorr_dailyavg"],data_plot["data_tower_dailyavg"],'b.')
        results = gfalternate_getolscorrecteddata(data_plot["data_alternate_lagcorr_dailyavg"],data_plot["data_tower_dailyavg"],alternate_info)
        ax4.plot(data_plot["data_alternate_lagcorr_dailyavg"],results["fit"],'g-')
        ax4.text(0.5,0.9,results["eqnstr"],fontsize=8,horizontalalignment='center',transform=ax4.transAxes,color='green')
        # hourly plot
        ax5.plot(ind,data_plot["data_tower_hourlyavg"],'ro',label="Tower")
        ax5.plot(ind,data_plot["data_alternate_lagfitcorr_hourlyavg"],'g-',label=source+" (OLS)")
        # time series plot
        ax6.plot(data_plot["odt_tower"],data_plot["data_tower"],'ro',label="Tower")
        ax6.plot(data_plot["odt_alternate"],data_plot["data_alternate_lagfitcorr"],'g-',label=source+" (OLS)")
    #else:
        #log.error("gfalternate: Less than 100 points available for series "+label+" ...")
    # put up the legends
    ax5.legend(loc='upper right',frameon=False,prop={'size':8})
    ax6.legend(loc='upper right',frameon=False,prop={'size':8})
    # write the comparison statistics
    output_list = ["Lag (corrected)","Lag (uncorrected)","Var (alternate)","Var (tower)","RMSE","Bias","r","No. points"]
    for n,item in enumerate(output_list):
        row_posn = pd["row_bottom"] + n*pd["row_space"]
        plt.figtext(pd["text_left"],row_posn,item)
        plt.figtext(pd["num_left"],row_posn,'%.4g'%(data_plot["results"][item][-1]))
    # save a hard copy of the plot
    sdt = data_plot["odt_tower"][0].strftime("%Y%m%d")
    edt = data_plot["odt_tower"][-1].strftime("%Y%m%d")
    figname = "plots/"+alternate_info["site_name"].replace(" ","")+"_"+source+"_"+label
    figname = figname+"_"+sdt+"_"+edt+'.png'
    fig.savefig(figname,format='png')
    # draw the plot on the screen
    if alternate_info["show_plots"]: plt.draw()
    # turn off interactive plotting
    if alternate_info["show_plots"]: plt.ioff()

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
    # set the default return code
    dsb.returncodes["solo"] = "normal"
    if "solo" not in dir(dsb): return
    # local pointer to the datetime series
    ldt = dsb.series["DateTime"]["Data"]
    startdate = ldt[0]
    enddate = ldt[-1]
    solo_info = {"file_startdate":startdate.strftime("%Y-%m-%d %H:%M"),
                 "file_enddate":enddate.strftime("%Y-%m-%d %H:%M")}
    # set up the GUI
    solo_gui = Tkinter.Toplevel()
    solo_gui.wm_title("SOLO GUI (Fluxes)")
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
    # seventh row
    nrow = nrow + 1
    solo_gui.peropt = Tkinter.IntVar()
    solo_gui.peropt.set(1)
    solo_gui.manualperiod = Tkinter.Radiobutton(solo_gui,text="Manual",variable=solo_gui.peropt,value=1)
    solo_gui.manualperiod.grid(row=nrow,column=0,columnspan=1,sticky="W")
    solo_gui.minptsLabel = Tkinter.Label(solo_gui,text="Min. pts (%)")
    solo_gui.minptsLabel.grid(row=nrow,column=1,columnspan=2,sticky="E")
    solo_gui.minptsEntry = Tkinter.Entry(solo_gui,width=5)
    solo_gui.minptsEntry.grid(row=nrow,column=3,columnspan=2,sticky="W")
    solo_gui.minptsEntry.insert(0,"50")
    # eigth row
    nrow = nrow + 1
    solo_gui.automonthly = Tkinter.Radiobutton(solo_gui,text="Monthly",variable=solo_gui.peropt,value=2)
    solo_gui.automonthly.grid(row=nrow,column=0,columnspan=1,sticky="W")
    solo_gui.daysLabel = Tkinter.Radiobutton(solo_gui,text="No. days",variable=solo_gui.peropt,value=3)
    solo_gui.daysLabel.grid(row=nrow,column=1,columnspan=2,sticky="E")
    solo_gui.daysEntry = Tkinter.Entry(solo_gui,width=5)
    solo_gui.daysEntry.grid(row=nrow,column=3,columnspan=2,sticky="W")
    solo_gui.daysEntry.insert(0,"120")
    # ninth row
    nrow = nrow + 1
    solo_gui.pltopt = Tkinter.IntVar()
    solo_gui.pltopt.set(1)
    solo_gui.showplots = Tkinter.Checkbutton(solo_gui, text="Show plots", variable=solo_gui.pltopt)
    solo_gui.showplots.grid(row=nrow,column=0,columnspan=3,sticky="w")
    solo_gui.owopt = Tkinter.IntVar()
    solo_gui.owopt.set(1)
    solo_gui.overwrite = Tkinter.Checkbutton(solo_gui, text="Overwrite", variable=solo_gui.owopt)
    solo_gui.overwrite.grid(row=nrow,column=3,columnspan=3,sticky="w")
    # tenth row
    nrow = nrow + 1
    solo_gui.doneButton = Tkinter.Button (solo_gui, text="Done",command=lambda:gfSOLO_done(dsb,solo_gui))
    solo_gui.doneButton.grid(row=nrow,column=0,columnspan=2)
    solo_gui.runButton = Tkinter.Button (solo_gui, text="Run",command=lambda:gfSOLO_run(dsa,dsb,solo_gui,solo_info))
    solo_gui.runButton.grid(row=nrow,column=2,columnspan=2)
    solo_gui.quitButton = Tkinter.Button (solo_gui, text="Quit",command=lambda:gfSOLO_quit(dsb,solo_gui))
    solo_gui.quitButton.grid(row=nrow,column=4,columnspan=2)
    # eleventh row
    nrow = nrow + 1
    solo_gui.progress_row = nrow
    solo_gui.progress = Tkinter.Label(solo_gui, text='Waiting for input ...')
    solo_gui.progress.grid(row=nrow,column=0,columnspan=6,sticky="W")

    solo_gui.wait_window(solo_gui)

def gfSOLO_done(ds,solo_gui):
    # plot the summary statistics if required
    if solo_gui.peropt.get()==1: gfSOLO_plotsummary(ds)
    # destroy the SOLO GUI
    solo_gui.destroy()
    # write Excel spreadsheet with fit statistics
    qcio.xl_write_SOLOStats(ds)
    # remove the solo dictionary from the data structure
    ds.returncodes["solo"] = "normal"

def gfSOLO_quit(ds,solo_gui):
    # destroy the GUI
    solo_gui.destroy()
    # put the return code in ds.returncodes
    ds.returncodes["solo"] = "quit"
    
def gfSOLO_getserieslist(cf):
    series_list = []
    if "Drivers" in cf.keys():
        for series in cf["Drivers"].keys():
            if "GapFillUsingSOLO" in cf["Drivers"][series]: series_list.append(series)
    if "Fluxes" in cf.keys():
        for series in cf["Fluxes"].keys():
            if "GapFillUsingSOLO" in cf["Fluxes"][series]: series_list.append(series)
    return series_list

def gfSOLO_main(dsa,dsb,solo_gui,solo_info):
    '''
    This is the main routine for running SOLO, an artifical neural network for gap filling fluxes.
    '''
    startdate = solo_info["startdate"]
    enddate = solo_info["enddate"]
    log.info(" Gap filling using SOLO: "+startdate+" to "+enddate)
    # read the control file again, this allows the contents of the control file to
    # be changed with the SOLO GUI still displayed
    cfname = dsb.globalattributes["controlfile_name"]
    cf = qcio.get_controlfilecontents(cfname,mode="quiet")
    # put the control file object in the solo_info dictionary
    dsb.cf = cf.copy()
    # get the list of series to gap fill using SOLO
    solo_series = gfSOLO_getserieslist(cf)
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
    solo_info["min_points"] = int(nRecs*solo_info["min_percent"]/100)
    # close any open plot windows
    if len(plt.get_fignums())!=0:
        for i in plt.get_fignums(): plt.close(i)
    fig_num = 0
    for series in solo_series:
        # clean up the target series if required
        qcck.do_qcchecks_oneseries(dsb.cf,dsa,series=series)
        dsb.solo[series]["results"]["startdate"].append(xldt[si])
        dsb.solo[series]["results"]["enddate"].append(xldt[ei])
        d,f,a = qcutils.GetSeriesasMA(dsb,series,si=si,ei=ei)
        if numpy.ma.count(d)<solo_info["min_points"]:
            log.error("gfSOLO: Less than "+str(solo_info["min_points"])+" points available for series "+series+" ...")
            dsb.solo[series]["results"]["No. points"].append(float(0))
            results_list = dsb.solo[series]["results"].keys()
            for item in ["startdate","enddate","No. points"]:
                if item in results_list: results_list.remove(item)
            for item in results_list:
                dsb.solo[series]["results"][item].append(float(c.missing_value))
            continue
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

def gfSOLO_progress(solo_gui,text):
    """
        Update progress message in SOLO GUI
        """
    solo_gui.progress.destroy()
    solo_gui.progress = Tkinter.Label(solo_gui, text=text)
    solo_gui.progress.grid(row=solo_gui.progress_row,column=0,columnspan=6,sticky="W")
    solo_gui.update()

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
    f.write(str(c.missing_value)+'\n')
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
    f.write('Line 22: missing data value, default value is c.missing_value.0\n')
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
        index = numpy.where(abs(driver-float(c.missing_value)<c.eps))[0]
        if len(index)!=0:
            log.error(' GapFillUsingSOLO: c.missing_value found in driver '+TheseOnes+' at lines '+str(index))
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
    if platform.system()=="Windows":
        subprocess.call(['./solo/bin/sofm.exe','solo/inf/sofm.inf'],stdout=sofmlogfile)
    else:
        subprocess.call(['./solo/bin/sofm','solo/inf/sofm.inf'],stdout=sofmlogfile)
    sofmlogfile.close()
    # check to see if the sofm output file exists, this is used to indicate that sofm ran correctly
    if os.path.exists('solo/output/sofm_4.out'):
        # write out the details of this run
        #gfSOLO_savethissofmrun(sofmoutname,sofminfname)
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
        index = numpy.where(soloinputdata[:,i]==c.missing_value)[0]
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
    if platform.system()=="Windows":
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

def gfSOLO_run(dsa,dsb,solo_gui,solo_info):
    # populate the solo_info dictionary with things that will be useful
    solo_info["peropt"] = solo_gui.peropt.get()
    solo_info["overwrite"] = True
    if solo_gui.owopt.get()==0: solo_info["overwrite"] = False
    solo_info["show_plots"] = True
    if solo_gui.pltopt.get()==0: solo_info["show_plots"] = False
    solo_info["min_percent"] = int(solo_gui.minptsEntry.get())
    solo_info["site_name"] = dsb.globalattributes["site_name"]
    solo_info["time_step"] = int(dsb.globalattributes["time_step"])
    solo_info["nperhr"] = int(float(60)/solo_info["time_step"]+0.5)
    solo_info["nperday"] = int(float(24)*solo_info["nperhr"]+0.5)
    solo_info["maxlags"] = int(float(12)*solo_info["nperhr"]+0.5)
    solo_info["tower"] = {}
    solo_info["alternate"] = {}
    log.info(" Gap filling "+str(dsb.solo.keys())+" using SOLO")
    if solo_gui.peropt.get()==1:
        gfSOLO_progress(solo_gui,"Starting manual run ...")
        # get the start and end datetimes entered in the SOLO GUI
        solo_info["startdate"] = solo_gui.startEntry.get()
        if len(solo_info["startdate"])==0: solo_info["startdate"] = solo_info["file_startdate"]
        solo_info["enddate"] = solo_gui.endEntry.get()
        if len(solo_info["enddate"])==0: solo_info["enddate"] = solo_info["file_enddate"]
        gfSOLO_main(dsa,dsb,solo_gui,solo_info)
        gfSOLO_progress(solo_gui,"Finished manual run ...")
    elif solo_gui.peropt.get()==2:
        gfSOLO_progress(solo_gui,"Starting auto (monthly) run ...")
        # get the start datetime entered in the SOLO GUI
        solo_info["startdate"] = solo_gui.startEntry.get()
        if len(solo_info["startdate"])==0: solo_info["startdate"] = solo_info["file_startdate"]
        startdate = dateutil.parser.parse(solo_info["startdate"])
        file_startdate = dateutil.parser.parse(solo_info["file_startdate"])
        file_enddate = dateutil.parser.parse(solo_info["file_enddate"])
        enddate = startdate+dateutil.relativedelta.relativedelta(months=1)
        enddate = min([file_enddate,enddate])
        solo_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<file_enddate:
            gfSOLO_main(dsa,dsb,solo_gui,solo_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(months=1)
            solo_info["startdate"] = startdate.strftime("%Y-%m-%d")
            solo_info["enddate"] = enddate.strftime("%Y-%m-%d")
        # plot the summary statistics
        gfSOLO_plotsummary(dsb)
        gfSOLO_progress(solo_gui,"Finished auto (monthly) run ...")
    elif solo_gui.peropt.get()==3:
        gfSOLO_progress(solo_gui,"Starting auto (days) run ...")
        # get the start datetime entered in the SOLO GUI
        solo_info["startdate"] = solo_gui.startEntry.get()
        if len(solo_info["startdate"])==0: solo_info["startdate"] = solo_info["file_startdate"]
        startdate = dateutil.parser.parse(solo_info["startdate"])
        file_startdate = dateutil.parser.parse(solo_info["file_startdate"])
        file_enddate = dateutil.parser.parse(solo_info["file_enddate"])
        nDays = int(solo_gui.daysentry.get())
        enddate = startdate+dateutil.relativedelta.relativedelta(days=nDays)
        enddate = min([file_enddate,enddate])
        solo_info["enddate"] = datetime.datetime.strftime(enddate,"%Y-%m-%d")
        while startdate<file_enddate:
            gfSOLO_main(dsa,dsb,solo_gui,solo_info)
            startdate = enddate
            enddate = startdate+dateutil.relativedelta.relativedelta(days=nDays)
            solo_info["startdate"] = startdate.strftime("%Y-%m-%d")
            solo_info["enddate"] = enddate.strftime("%Y-%m-%d")
        # plot the summary statistics
        gfSOLO_plotsummary(dsb)
        gfSOLO_progress(solo_gui,"Finished auto (days) run ...")
    elif solo_gui.peropt.get()==4:
        pass

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
    # now strip out the bad data
    cind = numpy.zeros(nRecs)
    iind = numpy.arange(nRecs)
    # do only the drivers not the target
    for i in range(ndrivers):
        index = numpy.where(seqsoloinputdata[:,i]==c.missing_value)[0]
        if len(index!=0): cind[index] = 1
    # index of good data
    index = numpy.where(cind==0)[0]
    nRecs_good = len(index)
    gooddata = numpy.zeros((nRecs_good,ndrivers+1))
    for i in range(ndrivers+1):
        gooddata[:,i] = seqsoloinputdata[:,i][index]
    # keep track of the good data indices
    goodindex = iind[index]
    # and then write the seqsolo input file
    seqsolofile = open('solo/input/seqsolo_input.csv','wb')
    wr = csv.writer(seqsolofile,delimiter=',')
    for i in range(gooddata.shape[0]):
        wr.writerow(gooddata[i,0:ndrivers+1])
    seqsolofile.close()
    # if the output file from a previous run exists, delete it
    if os.path.exists('solo/output/seqOut2.out'): os.remove('solo/output/seqOut2.out')
    # now run SEQSOLO
    #log.info(' GapFillUsingSOLO: running SEQSOLO')
    seqsolologfile = open('solo/log/seqsolo.log','wb')
    if platform.system()=="Windows":
        subprocess.call(['./solo/bin/seqsolo.exe','solo/inf/seqsolo.inf'],stdout=seqsolologfile)
    else:
        subprocess.call(['./solo/bin/seqsolo','solo/inf/seqsolo.inf'],stdout=seqsolologfile)
    seqsolologfile.close()
    # check to see if the solo output file exists, this is used to indicate that solo ran correctly
    if os.path.exists('solo/output/seqOut2.out'):
        # now read in the seqsolo results, use the seqOut2 file so that the learning capability of
        # seqsolo can be used via the "learning rate" and "Iterations" GUI options
        seqdata = numpy.genfromtxt('solo/output/seqOut2.out')
        # put the SOLO modelled data back into the data series
        if ei==-1:
            dsb.series[outputlabel]['Data'][si:][goodindex] = seqdata[:,1]
            dsb.series[outputlabel]['Flag'][si:][goodindex] = numpy.int32(30)
        else:
            dsb.series[outputlabel]['Data'][si:ei+1][goodindex] = seqdata[:,1]
            dsb.series[outputlabel]['Flag'][si:ei+1][goodindex] = numpy.int32(30)
        # set the attributes
        for attr in dsa.series[targetlabel]["Attr"].keys():
            dsb.series[outputlabel]["Attr"][attr] = dsa.series[targetlabel]["Attr"][attr]
        dsb.series[outputlabel]["Attr"]["long_name"] = dsb.series[outputlabel]["Attr"]["long_name"]+", modeled by SOLO"
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
    plt.ion()
    fig = plt.figure(pd["fig_num"],figsize=(13,9))
    fig.clf()
    fig.canvas.set_window_title(targetlabel)
    plt.figtext(0.5,0.95,pd["title"],ha='center',size=16)
    # XY plot of the diurnal variation
    rect1 = [0.10,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
    ax1 = plt.axes(rect1)
    # get the diurnal stats of the observations
    mask = numpy.ma.mask_or(obs.mask,mod.mask)
    obs_mor = numpy.ma.array(obs,mask=mask)
    Hr1,Av1,Sd1,Mx1,Mn1 = gf_getdiurnalstats(Hdh,obs_mor,ts)
    ax1.plot(Hr1,Av1,'b-',label="Obs")
    # get the diurnal stats of all SOLO predictions
    mod_mor = numpy.ma.array(mod,mask=mask)
    Hr2,Av2,Sd2,Mx2,Mn2 = gf_getdiurnalstats(Hdh,mod_mor,ts)
    ax1.plot(Hr2,Av2,'r-',label="SOLO(all)")
    if numpy.ma.count_masked(obs)!=0:
        index = numpy.ma.where(numpy.ma.getmaskarray(obs)==False)[0]
        # get the diurnal stats of SOLO predictions when observations are present
        Hr3,Av3,Sd3,Mx3,Mn3=gf_getdiurnalstats(Hdh[index],mod_mor[index],ts)
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
    coefs = numpy.ma.polyfit(numpy.ma.copy(mod),numpy.ma.copy(obs),1)
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
    dsb.solo[targetlabel]["results"]["Bias"].append(bias)
    rmse = numpy.ma.sqrt(numpy.ma.mean((obs-mod)*(obs-mod)))
    plt.figtext(0.65,0.225,'No. points')
    plt.figtext(0.75,0.225,str(numpoints))
    dsb.solo[targetlabel]["results"]["No. points"].append(numpoints)
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
    dsb.solo[targetlabel]["results"]["r"].append(r[0][1])
    plt.figtext(0.815,0.125,'RMSE')
    plt.figtext(0.915,0.125,str(qcutils.round2sig(rmse,sig=4)))
    dsb.solo[targetlabel]["results"]["RMSE"].append(rmse)
    var_obs = numpy.ma.var(obs)
    dsb.solo[targetlabel]["results"]["Var (obs)"].append(var_obs)
    var_mod = numpy.ma.var(mod)
    dsb.solo[targetlabel]["results"]["Var (SOLO)"].append(var_mod)
    dsb.solo[targetlabel]["results"]["Var ratio"].append(var_obs/var_mod)
    dsb.solo[targetlabel]["results"]["Avg (obs)"].append(numpy.ma.average(obs))
    dsb.solo[targetlabel]["results"]["Avg (SOLO)"].append(numpy.ma.average(mod))    
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
    for ThisOne,i in zip(driverlist,range(1,pd["nDrivers"]+1)):
        this_bottom = pd["ts_bottom"] + i*pd["ts_height"]
        rect = [pd["margin_left"],this_bottom,pd["ts_width"],pd["ts_height"]]
        ts_axes.append(plt.axes(rect,sharex=ts_axes[0]))
        data,flag,attr = qcutils.GetSeriesasMA(dsb,ThisOne,si=si,ei=ei)
        ts_axes[i].plot(xdt,data)
        plt.setp(ts_axes[i].get_xticklabels(),visible=False)
        TextStr = ThisOne+'('+dsb.series[ThisOne]['Attr']['units']+')'
        ts_axes[i].text(0.05,0.85,TextStr,color='b',horizontalalignment='left',transform=ts_axes[i].transAxes)
    # save a hard copy of the plot
    sdt = xdt[0].strftime("%Y%m%d")
    edt = xdt[-1].strftime("%Y%m%d")
    figname = "plots/"+pd["site_name"].replace(" ","")+"_SOLO_"+pd["label"]
    figname = figname+"_"+sdt+"_"+edt+'.png'
    fig.savefig(figname,format='png')
    # draw the plot on the screen
    plt.draw()
    # turn off interactive plotting
    plt.ioff()

def gfSOLO_plotsummary(ds):
    """ Plot single pages of summary results for groups of variables. """
    # get a list of variables for which SOLO data was available
    label_list = ds.solo.keys()
    if len(ds.solo[label_list[0]]["results"]["startdate"])==0:
        log.info("gfSOLO: no summary data to plot")
        return
    # get the Excel datemode, needed to convert the Excel datetime to Python datetimes
    datemode = int(ds.globalattributes['xl_datemode'])
    # site name for titles
    site_name = ds.globalattributes["site_name"]
    # datetimes are stored in ds.alternate as Excel datetimes, here we convert to Python datetimes
    # for ease of handling and plotting.
    # start datetimes of the periods compared first
    basedate = datetime.datetime(1899, 12, 30)
    dt_start = []
    for xldt in ds.solo[label_list[0]]["results"]["startdate"]:
        dt_start.append(basedate+datetime.timedelta(days=xldt+1462*datemode))
    startdate = min(dt_start)
    # and then the end datetimes
    dt_end = []
    for xldt in ds.solo[label_list[0]]["results"]["enddate"]:
        dt_end.append(basedate+datetime.timedelta(days=xldt+1462*datemode))
    enddate = max(dt_end)
    # get the major tick locator and label format
    MTLoc = mdt.AutoDateLocator(minticks=3,maxticks=5)
    MTFmt = mdt.DateFormatter('%b')
    # group lists of the resuts to be plotted
    result_list = ["r","Bias","RMSE","Var ratio","m_ols","b_ols"]
    ylabel_list = ["r","Bias","RMSE","Var ratio","Slope","Offset"]
    # turn on interactive plotting
    plt.ion()
    # now loop over the group lists
    for nFig in ds.cf["SOLO_Summary"].keys():
        plot_title = ds.cf["SOLO_Summary"][str(nFig)]["Title"]
        var_list = ast.literal_eval(ds.cf["SOLO_Summary"][str(nFig)]["Variables"])
        # set up the subplots on the page
        fig,axs = plt.subplots(len(result_list),len(var_list),figsize=(13,9))
        fig.canvas.set_window_title("SOLO summary: "+plot_title)
        # make a title string for the plot and render it
        title_str = "SOLO: "+plot_title+"; "+site_name+" "+datetime.datetime.strftime(startdate,"%Y-%m-%d")
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
                result = numpy.ma.masked_equal(ds.solo[label]["results"][rlabel],float(c.missing_value))
                # put the data into the right order to be plotted
                dt,data = gfSOLO_plotsummary_getdata(dt_start,dt_end,result)
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
        figname = "plots/"+site_name.replace(" ","")+"_SOLO_FitStatistics_"+figlab
        figname = figname+"_"+sdt+"_"+edt+".png"
        fig.savefig(figname,format="png")

def gfSOLO_plotsummary_getdata(dt_start,dt_end,result):
    dt = []
    data = []
    for s,e,r in zip(dt_start,dt_end,result):
        dt.append(s)
        data.append(r)
        dt.append(e)
        data.append(r)
    return dt,data

def gf_getdiurnalstats(DecHour,Data,dt):
    nInts = 24*int((60/dt)+0.5)
    Hr = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Av = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Sd = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Mx = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Mn = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    for i in range(nInts):
        Hr[i] = float(i)*dt/60.
        li = numpy.where((abs(DecHour-Hr[i])<c.eps)&(abs(Data-float(c.missing_value))>c.eps))
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

def ImportSeries(cf,ds):
    # check to see if there is an Imports section
    if "Imports" not in cf.keys(): return
    # number of records
    nRecs = int(ds.globalattributes["nc_nrecs"])
    # get the start and end datetime
    ldt = ds.series["DateTime"]["Data"]
    start_date = ldt[0]
    end_date = ldt[-1]
    # loop over the series in the Imports section
    for label in cf["Imports"].keys():
        import_filename = cf["Imports"][label]["file_name"]
        var_name = cf["Imports"][label]["var_name"]
        ds_import = qcio.nc_read_series(import_filename)
        ts_import = ds_import.globalattributes["time_step"]
        ldt_import = ds_import.series["DateTime"]["Data"]
        si = qcutils.GetDateIndex(ldt_import,str(start_date),ts=ts_import,default=0,match="exact")
        ei = qcutils.GetDateIndex(ldt_import,str(end_date),ts=ts_import,default=-1,match="exact")
        data = numpy.ma.ones(nRecs)*float(c.missing_value)
        flag = numpy.ma.ones(nRecs)
        data_import,flag_import,attr_import = qcutils.GetSeriesasMA(ds_import,var_name,si=si,ei=ei)
        ldt_import = ldt_import[si:ei+1]
        index = qcutils.find_indices(ldt,ldt_import)
        data[index] = data_import
        flag[index] = flag_import
        qcutils.CreateSeries(ds,label,data,Flag=flag,Attr=attr_import)