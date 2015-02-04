import ast
import calendar
import constants as c
import datetime
import logging
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy
import time
import qcio
import qcplot
import qcts
import qcutils
import sys
import xlwt

log = logging.getLogger('qc.clim')

def do_2dinterpolation(array_2d):
    """
    Takes a 2d array as input and;
     1) tiles this into a 3 x 3 space (9 repeats of the original 2d array in 3 columns and 3 rows)
     2) removes the missing data (c.missing_value) from the tiled array
     3) does a bi-linear interpolation to replace the the missing data
     4) returns the central tile
     The effect is to replace missing data in the original 2d array with data from a bi-linear
     interpolation, the tiling repeats the original array along its boundaries to avoid problems
     at the array edges.
    """
    WasMA = False
    if numpy.ma.isMA(array_2d):
        WasMA = True
        array_2d = numpy.ma.filled(array_2d,float(c.missing_value))
    # tile the 2d array into a 3 by 3 array
    array_2d_3x3=numpy.tile(array_2d,(3,3))
    # get the dimensions of the tiled array
    nmn=numpy.shape(array_2d_3x3)[1]
    mni=numpy.arange(0,nmn)
    nhr=numpy.shape(array_2d_3x3)[0]
    hri=numpy.arange(0,nhr)
    mn,hr=numpy.meshgrid(mni,hri)
    array_2d_3x3_1d=numpy.reshape(array_2d_3x3,numpy.shape(array_2d_3x3)[0]*numpy.shape(array_2d_3x3)[1])
    mn_1d=numpy.reshape(mn,numpy.shape(mn)[0]*numpy.shape(mn)[1])
    hr_1d=numpy.reshape(hr,numpy.shape(hr)[0]*numpy.shape(hr)[1])
    index=numpy.where(array_2d_3x3_1d!=c.missing_value)
    array_2d_3x3i=griddata(mn_1d[index],hr_1d[index],array_2d_3x3_1d[index],mni,hri)
    array_2di=array_2d_3x3i[nhr/3:2*nhr/3,nmn/3:2*nmn/3]
    #array_2di=numpy.ma.filled(array_2d_3x3i[nhr/3:2*nhr/3,nmn/3:2*nmn/3],0)
    if WasMA:
        array_2di = numpy.ma.masked_where(abs(array_2di-numpy.float64(c.missing_value))<c.eps,array_2di)
        array_2d = numpy.ma.masked_where(abs(array_2d-numpy.float64(c.missing_value))<c.eps,array_2d)
    return array_2di

def write_data_1columnpermonth(xlSheet, data, ts, format_string=''):
    xlCol = 0
    # write the data to the xl file
    nrows = numpy.shape(data)[0]
    ncols = numpy.shape(data)[1]
    xlSheet.write(1,xlCol,'Hour')
    for j in range(nrows+1):
        xlSheet.write(j+2,xlCol,float(j)*ts/60)
    xlCol = xlCol + 1
    if len(format_string)!=0:
        d_xf = xlwt.easyxf(num_format_str=format_string)
    else:
        d_xf = xlwt.easyxf()
    for m in range(1,ncols+1):
        xlSheet.write(0,xlCol,calendar.month_abbr[m])
        xlSheet.write(1,xlCol,'Av')
        for j in range(nrows):
            xlSheet.write(j+2,xlCol,data[j,m-1],d_xf)
        xlCol = xlCol + 1

def write_data_1columnpertimestep(xlSheet, data, ts, year=None, format_string=''):
    tmp = data.copy()
    if numpy.ma.isMA(tmp): tmp = numpy.ma.filled(tmp,float(c.missing_value))
    xlCol = 0
    # write the data to the xl file
    xlSheet.write(1,xlCol,'Day')
    nrows = numpy.shape(tmp)[0]
    ncols = numpy.shape(tmp)[1]
    if year is None:
        for j in range(nrows+1):
            xlSheet.write(j+2,xlCol,j)
    else:
        d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy')
        for j in range(nrows):
            d = datetime.datetime(year, 1, 1) + datetime.timedelta(days=j)
            xlSheet.write(j+2,xlCol,d,d_xf)
    xlCol = xlCol + 1
    if len(format_string)!=0:
        d_xf = xlwt.easyxf(num_format_str=format_string)
    else:
        d_xf = xlwt.easyxf()
    for m in range(1,ncols+1):
        xlSheet.write(1,xlCol,float(m)*ts/60)
        for j in range(nrows):
            xlSheet.write(j+2,xlCol,tmp[j,m-1],d_xf)
        xlCol = xlCol + 1

def do_diurnalstats(Month, Hdh, data, xlSheet, format_string='',ts=30):
    xlCol = 0
    nInts = 24*int((60/ts)+0.5)
    Av_all = numpy.ma.zeros([nInts,12]) + float(c.missing_value)
    if len(format_string)!=0:
        d_xf = xlwt.easyxf(num_format_str=format_string)
    else:
        d_xf = xlwt.easyxf()
    for m in range(1,13):
        mi = numpy.where(Month==m)[0]
        Num,Hr,Av,Sd,Mx,Mn = get_diurnalstats(Hdh[mi],data[mi],ts)
        Av_all[:,m-1] = Av[:]
        Num = numpy.ma.filled(Num,float(c.missing_value))
        Hr = numpy.ma.filled(Hr,float(c.missing_value))
        Av = numpy.ma.filled(Av,float(c.missing_value))
        Sd = numpy.ma.filled(Sd,float(c.missing_value))
        Mx = numpy.ma.filled(Mx,float(c.missing_value))
        Mn = numpy.ma.filled(Mn,float(c.missing_value))
        if m==1:
            xlSheet.write(1,xlCol,'Hour')
            for j in range(len(Hr)):
                xlSheet.write(j+2,xlCol,Hr[j])
            xlCol = xlCol + 1
        xlSheet.write(0,xlCol,calendar.month_abbr[m])
        xlSheet.write(1,xlCol,'Num')
        xlSheet.write(1,xlCol+1,'Av')
        xlSheet.write(1,xlCol+2,'Sd')
        xlSheet.write(1,xlCol+3,'Mx')
        xlSheet.write(1,xlCol+4,'Mn')
        for j in range(len(Hr)):
            xlSheet.write(j+2,xlCol,Num[j])
            xlSheet.write(j+2,xlCol+1,Av[j],d_xf)
            xlSheet.write(j+2,xlCol+2,Sd[j],d_xf)
            xlSheet.write(j+2,xlCol+3,Mx[j],d_xf)
            xlSheet.write(j+2,xlCol+4,Mn[j],d_xf)
        xlCol = xlCol + 5
    return Av_all

def get_diurnalstats(DecHour,Data,ts):
    nInts = 24*int((60/ts)+0.5)
    Num = numpy.ma.zeros(nInts,dtype=int)
    Hr = numpy.ma.zeros(nInts,dtype=float)
    for i in range(nInts):
        Hr[i] = float(i)*ts/60.
    Av = numpy.ma.masked_all(nInts)
    Sd = numpy.ma.masked_all(nInts)
    Mx = numpy.ma.masked_all(nInts)
    Mn = numpy.ma.masked_all(nInts)
    if numpy.size(Data)!=0:
        for i in range(nInts):
            li = numpy.ma.where((abs(DecHour-Hr[i])<c.eps)&(abs(Data-float(c.missing_value))>c.eps))
            Num[i] = numpy.size(li)
            if Num[i]!=0:
                Av[i] = numpy.ma.mean(Data[li])
                Sd[i] = numpy.ma.std(Data[li])
                Mx[i] = numpy.ma.maximum(Data[li])
                Mn[i] = numpy.ma.minimum(Data[li])
    return Num, Hr, Av, Sd, Mx, Mn

def get_rangecheck_limit(cf,label,upr_def=1E10,lwr_def=-1E10):
    upper = float(upr_def)
    lower = float(lwr_def)
    for section in ['Variables']:
        if label in cf[section].keys():
            if 'RangeCheck' in cf[section][label].keys():
                upper = float(cf[section][label]['RangeCheck']['Upper'])
                lower = float(cf[section][label]['RangeCheck']['Lower'])
    return upper,lower

def get_formatstring(cf,label,fmt_def=''):
    fmt_str = fmt_def
    for section in ['Variables']:
        if label in cf[section].keys():
            if 'Format' in cf[section][label].keys():
                fmt_str = str(cf[section][label]['Format'])
    return fmt_str

def climatology(cf):
    nc_filename = qcio.get_infilename_from_cf(cf)
    if not qcutils.file_exists(nc_filename): return
    xl_filename = nc_filename.replace(".nc","_Climatology.xls")
    xlFile = xlwt.Workbook()
    ds = qcio.nc_read_series(nc_filename)
    # calculate Fa if it is not in the data structure
    if "Fa" not in ds.series.keys():
        if "Fn" in ds.series.keys() and "Fg" in ds.series.keys():
            qcts.CalculateAvailableEnergy(ds,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg')
        else:
            log.error("Climatology: Fn or Fg not in data struicture")
    # get the time step
    ts = int(ds.globalattributes['time_step'])
    # get the site name
    SiteName = ds.globalattributes['site_name']
    # get the datetime series
    DateTime = ds.series['DateTime']['Data']
    Hdh = ds.series['Hdh']['Data']
    Month = ds.series['Month']['Data']
    # get the initial start and end dates
    StartDate = str(DateTime[0])
    EndDate = str(DateTime[-1])
    # find the start index of the first whole day (time=00:30)
    si = qcutils.GetDateIndex(DateTime,StartDate,ts=ts,default=0,match='startnextday')
    # find the end index of the last whole day (time=00:00)
    ei = qcutils.GetDateIndex(DateTime,EndDate,ts=ts,default=-1,match='endpreviousday')
    DateTime = DateTime[si:ei+1]
    #print time.strftime('%X')+' Start date; '+str(DateTime[0])+' End date; '+str(DateTime[-1])
    Hdh = Hdh[si:ei+1]
    Month = Month[si:ei+1]
    
    ntsInDay = int(24.0*60.0/float(ts))
    nDays = int(len(DateTime))/ntsInDay
    
    for ThisOne in cf['Variables'].keys():
        if "AltVarName" in cf['Variables'][ThisOne].keys(): ThisOne = cf['Variables'][ThisOne]["AltVarName"]
        if ThisOne in ds.series.keys():
            log.info("Doing climatology for "+ThisOne)
            data,f,a = qcutils.GetSeriesasMA(ds,ThisOne,si=si,ei=ei)
            if numpy.ma.count(data)==0:
                log.error(" No data for "+ThisOne+", skipping ...")
                continue
            fmt_str = get_formatstring(cf,ThisOne,fmt_def='')
            xlSheet = xlFile.add_sheet(ThisOne)
            Av_all = do_diurnalstats(Month,Hdh,data,xlSheet,format_string=fmt_str,ts=ts)
            # now do it for each day
            data_daily = data.reshape(nDays,ntsInDay)
            year = DateTime[0].year
            xlSheet = xlFile.add_sheet(ThisOne+'(day)')
            write_data_1columnpertimestep(xlSheet, data_daily, ts, year=year, format_string=fmt_str)
            data_daily_i = do_2dinterpolation(data_daily)
            xlSheet = xlFile.add_sheet(ThisOne+'i(day)')
            write_data_1columnpertimestep(xlSheet, data_daily_i, ts, year=year, format_string=fmt_str)
        elif ThisOne=="EF":
            log.info("Doing evaporative fraction")
            EF = numpy.ma.zeros([48,12]) + float(c.missing_value)
            Hdh,f,a = qcutils.GetSeriesasMA(ds,'Hdh',si=si,ei=ei)
            Fa,f,a = qcutils.GetSeriesasMA(ds,'Fa',si=si,ei=ei)
            Fe,f,a = qcutils.GetSeriesasMA(ds,'Fe',si=si,ei=ei)
            for m in range(1,13):
                mi = numpy.where(Month==m)[0]
                Fa_Num,Hr,Fa_Av,Sd,Mx,Mn = get_diurnalstats(Hdh[mi],Fa[mi],ts)
                Fe_Num,Hr,Fe_Av,Sd,Mx,Mn = get_diurnalstats(Hdh[mi],Fe[mi],ts)
                index = numpy.ma.where((Fa_Num>4)&(Fe_Num>4))
                EF[:,m-1][index] = Fe_Av[index]/Fa_Av[index]
            # reject EF values greater than upper limit or less than lower limit
            upr, lwr = get_rangecheck_limit(cf,'EF')
            EF = numpy.ma.filled(numpy.ma.masked_where((EF>upr)|(EF<lwr),EF),float(c.missing_value))
            # write the EF to the Excel file
            xlSheet = xlFile.add_sheet('EF')
            write_data_1columnpermonth(xlSheet, EF, ts, format_string='0.00')
            # do the 2D interpolation to fill missing EF values
            EFi = do_2dinterpolation(EF)
            xlSheet = xlFile.add_sheet('EFi')
            write_data_1columnpermonth(xlSheet, EFi, ts, format_string='0.00')
            # now do EF for each day
            Fa,f,a = qcutils.GetSeriesasMA(ds,'Fa',si=si,ei=ei)
            Fe,f,a = qcutils.GetSeriesasMA(ds,'Fe',si=si,ei=ei)
            EF = Fe/Fa
            EF = numpy.ma.filled(numpy.ma.masked_where((EF>upr)|(EF<lwr),EF),float(c.missing_value))
            EF_daily = EF.reshape(nDays,ntsInDay)
            year = DateTime[0].year
            xlSheet = xlFile.add_sheet('EF(day)')
            write_data_1columnpertimestep(xlSheet, EF_daily, ts, year=year, format_string='0.00')
            EFi = do_2dinterpolation(EF_daily)
            xlSheet = xlFile.add_sheet('EFi(day)')
            write_data_1columnpertimestep(xlSheet, EFi, ts, year=year, format_string='0.00')
        elif ThisOne=="BR":
            log.info("Doing Bowen ratio")
            BR = numpy.ma.zeros([48,12]) + float(c.missing_value)
            Fe,f,a = qcutils.GetSeriesasMA(ds,'Fe',si=si,ei=ei)
            Fh,f,a = qcutils.GetSeriesasMA(ds,'Fh',si=si,ei=ei)
            for m in range(1,13):
                mi = numpy.where(Month==m)[0]
                Fh_Num,Hr,Fh_Av,Sd,Mx,Mn = get_diurnalstats(Hdh[mi],Fh[mi],ts)
                Fe_Num,Hr,Fe_Av,Sd,Mx,Mn = get_diurnalstats(Hdh[mi],Fe[mi],ts)
                index = numpy.ma.where((Fh_Num>4)&(Fe_Num>4))
                BR[:,m-1][index] = Fh_Av[index]/Fe_Av[index]
            # reject BR values greater than upper limit or less than lower limit
            upr,lwr = get_rangecheck_limit(cf,'BR')
            BR = numpy.ma.filled(numpy.ma.masked_where((BR>upr)|(BR<lwr),BR),float(c.missing_value))
            # write the BR to the Excel file
            xlSheet = xlFile.add_sheet('BR')
            write_data_1columnpermonth(xlSheet, BR, ts, format_string='0.00')
            # do the 2D interpolation to fill missing EF values
            BRi = do_2dinterpolation(BR)
            xlSheet = xlFile.add_sheet('BRi')
            write_data_1columnpermonth(xlSheet, BRi, ts, format_string='0.00')
            # now do BR for each day ...
            Fe,f,a = qcutils.GetSeriesasMA(ds,'Fe',si=si,ei=ei)
            Fh,f,a = qcutils.GetSeriesasMA(ds,'Fh',si=si,ei=ei)
            BR = Fh/Fe
            BR = numpy.ma.filled(numpy.ma.masked_where((BR>upr)|(BR<lwr),BR),float(c.missing_value))
            BR_daily = BR.reshape(nDays,ntsInDay)
            year = DateTime[0].year
            xlSheet = xlFile.add_sheet('BR(day)')
            write_data_1columnpertimestep(xlSheet, BR_daily, ts, year=year, format_string='0.00')
            BRi = do_2dinterpolation(BR_daily)
            xlSheet = xlFile.add_sheet('BRi(day)')
            write_data_1columnpertimestep(xlSheet, BRi, ts, year=year, format_string='0.00')
        elif ThisOne=="WUE":
            log.info("Doing ecosystem WUE")
            WUE = numpy.ma.zeros([48,12]) + float(c.missing_value)
            Fe,f,a = qcutils.GetSeriesasMA(ds,'Fe',si=si,ei=ei)
            Fc,f,a = qcutils.GetSeriesasMA(ds,'Fc',si=si,ei=ei)
            for m in range(1,13):
                mi = numpy.where(Month==m)[0]
                Fc_Num,Hr,Fc_Av,Sd,Mx,Mn = get_diurnalstats(Hdh[mi],Fc[mi],ts)
                Fe_Num,Hr,Fe_Av,Sd,Mx,Mn = get_diurnalstats(Hdh[mi],Fe[mi],ts)
                index = numpy.ma.where((Fc_Num>4)&(Fe_Num>4))
                WUE[:,m-1][index] = Fc_Av[index]/Fe_Av[index]
            # reject WUE values greater than upper limit or less than lower limit
            upr,lwr = get_rangecheck_limit(cf,'WUE')
            WUE = numpy.ma.filled(numpy.ma.masked_where((WUE>upr)|(WUE<lwr),WUE),float(c.missing_value))
            # write the WUE to the Excel file
            xlSheet = xlFile.add_sheet('WUE')
            write_data_1columnpermonth(xlSheet, WUE, ts, format_string='0.00000')
            # do the 2D interpolation to fill missing EF values
            WUEi = do_2dinterpolation(WUE)
            xlSheet = xlFile.add_sheet('WUEi')
            write_data_1columnpermonth(xlSheet, WUEi, ts, format_string='0.00000')
            # now do WUE for each day ...
            Fe,f,a = qcutils.GetSeriesasMA(ds,'Fe',si=si,ei=ei)
            Fc,f,a = qcutils.GetSeriesasMA(ds,'Fc',si=si,ei=ei)
            WUE = Fc/Fe
            WUE = numpy.ma.filled(numpy.ma.masked_where((WUE>upr)|(WUE<lwr),WUE),float(c.missing_value))
            WUE_daily = WUE.reshape(nDays,ntsInDay)
            year = DateTime[0].year
            xlSheet = xlFile.add_sheet('WUE(day)')
            write_data_1columnpertimestep(xlSheet, WUE_daily, ts, year=year, format_string='0.00000')
            WUEi = do_2dinterpolation(WUE_daily)
            xlSheet = xlFile.add_sheet('WUEi(day)')
            write_data_1columnpertimestep(xlSheet, WUEi, ts, year=year, format_string='0.00000')
        else:
            log.error("qcclim.climatology: requested variable "+ThisOne+" not in data structure")
            continue
    log.info("Saving Excel file "+xl_filename)
    xlFile.save(xl_filename)

def compare_eddypro():
    epname = qcio.get_filename_dialog(title='Choose an EddyPro full output file')
    ofname = qcio.get_filename_dialog(title='Choose an L3 output file')
    
    ds_ep = qcio.read_eddypro_full(epname)
    ds_of = qcio.nc_read_series(ofname)
    
    dt_ep = ds_ep.series['DateTime']['Data']
    dt_of = ds_of.series['DateTime']['Data']
    
    si = dt_of.index(dt_ep[0])
    ei = dt_of.index(dt_ep[-1])
    
    us_of,f,a = qcutils.GetSeriesasMA(ds_of,'ustar',si=si,ei=ei)
    us_ep,f,a = qcutils.GetSeriesasMA(ds_ep,'ustar')
    Fh_of,f,a = qcutils.GetSeriesasMA(ds_of,'Fh',si=si,ei=ei)
    Fh_ep,f,a = qcutils.GetSeriesasMA(ds_ep,'Fh')
    Fe_of,f,a = qcutils.GetSeriesasMA(ds_of,'Fe',si=si,ei=ei)
    Fe_ep,f,a = qcutils.GetSeriesasMA(ds_ep,'Fe')
    Fc_of,f,a = qcutils.GetSeriesasMA(ds_of,'Fc',si=si,ei=ei)
    Fc_ep,f,a = qcutils.GetSeriesasMA(ds_ep,'Fc')
    
    us_of.mask = numpy.ma.mask_or(us_of.mask,us_ep.mask)
    us_ep.mask = numpy.ma.mask_or(us_of.mask,us_ep.mask)
    Fh_of.mask = numpy.ma.mask_or(Fh_of.mask,Fh_ep.mask)
    Fh_ep.mask = numpy.ma.mask_or(Fh_of.mask,Fh_ep.mask)
    Fe_of.mask = numpy.ma.mask_or(Fe_of.mask,Fe_ep.mask)
    Fe_ep.mask = numpy.ma.mask_or(Fe_of.mask,Fe_ep.mask)
    Fc_of.mask = numpy.ma.mask_or(Fc_of.mask,Fc_ep.mask)
    Fc_ep.mask = numpy.ma.mask_or(Fc_of.mask,Fc_ep.mask)
    
    plt.ion()
    fig = plt.figure(1,figsize=(8,8))
    qcplot.xyplot(us_ep,us_of,sub=[2,2,1],regr=1,xlabel='u*_EP (m/s)',ylabel='u*_OF (m/s)')
    qcplot.xyplot(Fh_ep,Fh_of,sub=[2,2,2],regr=1,xlabel='Fh_EP (W/m2)',ylabel='Fh_OF (W/m2)')
    qcplot.xyplot(Fe_ep,Fe_of,sub=[2,2,3],regr=1,xlabel='Fe_EP (W/m2)',ylabel='Fe_OF (W/m2)')
    qcplot.xyplot(Fc_ep,Fc_of,sub=[2,2,4],regr=1,xlabel='Fc_EP (umol/m2/s)',ylabel='Fc_OF (umol/m2/s)')
    plt.tight_layout()
    plt.draw()
    plt.ioff()
