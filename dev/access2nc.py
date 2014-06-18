import os
import sys
# check the scripts directory is present
if not os.path.exists("../scripts/"):
    print "compare_access: the scripts directory is missing"
    sys.exit()
# since the scripts directory is there, try importing the modules
sys.path.append('../scripts')
import constants as c
import datetime
import matplotlib.pyplot as plt
import meteorologicalfunctions as mf
from netCDF4 import MFDataset
import numpy
import pytz
import qcio
import qcutils
from scipy.interpolate import interp1d

def perdelta(start,end,delta):
    curr = start
    while curr <= end:
        yield curr
        curr += delta

# get the control file
cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0: sys.exit()
infilename = qcio.get_infilename_from_cf(cf)
outfilename = qcio.get_outfilename_from_cf(cf)

# create an instance of the data structure
ds_60minutes = qcio.DataStructure()
ds_30minutes = qcio.DataStructure()
# read the netCDF files
f = MFDataset(infilename)
site_tz = pytz.timezone(cf["Global"]["site_timezone"])

valid_date = f.variables["valid_date"]
valid_time = f.variables["valid_time"]
nRecs = len(valid_date)

# set some global attributres
ds_60minutes.globalattributes["nc_nrecs"] = nRecs
ds_60minutes.globalattributes["time_step"] = 60
ds_60minutes.globalattributes["time_zone"] = cf["Global"]["site_timezone"]
# map the ACCESS file global attributes to the OzFluxQC file
for attr in f.ncattrs():
    ds_60minutes.globalattributes[attr] = getattr(f,attr)
    
#dt=[datetime.datetime.strptime(str(valid_date[i]*10000+valid_time[i]),"%Y%m%d%H%M") for i in range(0,nRecs)]
dt_utc_60minutes=[datetime.datetime.strptime(str(valid_date[i]*10000+valid_time[i]),"%Y%m%d%H%M") for i in range(0,len(valid_date))]
# make utc_dt timezone aware
dt_utc_60minutes=[x.replace(tzinfo=pytz.utc) for x in dt_utc_60minutes]
# get local time from UTC
# NOTE: will have to disable daylight saving at some stage, towers stay on Standard Time
dt_loc_60minutes=[x.astimezone(site_tz) for x in dt_utc_60minutes]
# make local time timezone naive to match datetimes in OzFluxQC
dt_loc_60minutes=[x.replace(tzinfo=None) for x in dt_loc_60minutes]
ds_60minutes.series["DateTime"] = {}
ds_60minutes.series["DateTime"]["Data"] = dt_loc_60minutes
# get the year, month etc from the datetime
flag_60minutes = numpy.zeros(nRecs)
xl_date = qcutils.get_xldate_from_datetime(dt_loc_60minutes)
attr = qcutils.MakeAttributeDictionary(long_name="Date/time in Excel format",units="days since 1899-12-31 00:00:00")
qcutils.CreateSeries(ds_60minutes,"xlDateTime",xl_date,Flag=flag_60minutes,Attr=attr)
qcutils.get_ymdhms_from_datetime(ds_60minutes)

# load the data into the data structure
varlist = cf["Variables"].keys()
for label in varlist:
    access_name = cf["Variables"][label]["access_name"]
    attr = {}
    for this_attr in f.variables[access_name].ncattrs():
        attr[this_attr] = getattr(f.variables[access_name],this_attr)
    # loop over all ACCESS grids and give them standard OzFlux names with the grid idices appended
    for i in range(0,3):
        for j in range(0,3):
            if len(f.variables[access_name].shape)==3:
                series = f.variables[access_name][:,i,j]
                label_ij = label+'_'+str(i)+str(j)
                qcutils.CreateSeries(ds_60minutes,label_ij,series,Flag=flag_60minutes,Attr=attr)
            elif len(f.variables[access_name].shape)==4:
                series = f.variables[access_name][:,0,i,j]
                label_ij = label+'_'+str(i)+str(j)
                qcutils.CreateSeries(ds_60minutes,label_ij,series,Flag=flag_60minutes,Attr=attr)
            else:
                print "Unrecognised variable ("+label+") dimension in ACCESS file"
                #sys.exit()

# get derived quantities and adjust units
# air temperature from K to C
attr = qcutils.GetAttributeDictionary(ds_60minutes,"Ta_00")
if attr["units"] == "K":
    for i in range(0,3):
        for j in range(0,3):
            label = "Ta_"+str(i)+str(j)
            Ta,f = qcutils.GetSeriesasMA(ds_60minutes,label)
            Ta = Ta - c.C2K
            attr["units"] = "C"
            qcutils.CreateSeries(ds_60minutes,label,Ta,Flag=flag_60minutes,Attr=attr)
# soil temperature from K to C
attr = qcutils.GetAttributeDictionary(ds_60minutes,"Ts_00")
if attr["units"] == "K":
    for i in range(0,3):
        for j in range(0,3):
            label = "Ts_"+str(i)+str(j)
            Ts,f = qcutils.GetSeriesasMA(ds_60minutes,label)
            Ts = Ts - c.C2K
            attr["units"] = "C"
            qcutils.CreateSeries(ds_60minutes,label,Ts,Flag=flag_60minutes,Attr=attr)
# pressure from Pa to kPa
attr = qcutils.GetAttributeDictionary(ds_60minutes,"ps_00")
if attr["units"] == "Pa":
    for i in range(0,3):
        for j in range(0,3):
            label = "ps_"+str(i)+str(j)
            ps,f = qcutils.GetSeriesasMA(ds_60minutes,label)
            ps = ps/float(1000)
            attr["units"] = "kPa"
            qcutils.CreateSeries(ds_60minutes,label,ps,Flag=flag_60minutes,Attr=attr)
# wind speed from components
for i in range(0,3):
    for j in range(0,3):
        u_label = "u_"+str(i)+str(j)
        v_label = "v_"+str(i)+str(j)
        Ws_label = "Ws_"+str(i)+str(j)
        u,f = qcutils.GetSeriesasMA(ds_60minutes,u_label)
        v,f = qcutils.GetSeriesasMA(ds_60minutes,v_label)
        Ws = numpy.sqrt(u*u+v*v)
        attr = qcutils.MakeAttributeDictionary(long_name="Wind speed",units="m/s",height="10m")
        qcutils.CreateSeries(ds_60minutes,Ws_label,Ws,Flag=f,Attr=attr)
# wind direction from components
for i in range(0,3):
    for j in range(0,3):
        u_label = "u_"+str(i)+str(j)
        v_label = "v_"+str(i)+str(j)
        Wd_label = "Wd_"+str(i)+str(j)
        u,f = qcutils.GetSeriesasMA(ds_60minutes,u_label)
        v,f = qcutils.GetSeriesasMA(ds_60minutes,v_label)
        Wd = float(270) - numpy.ma.arctan2(v,u)*float(180)/numpy.pi
        attr = qcutils.MakeAttributeDictionary(long_name="Wind direction",units="degrees",height="10m")
        qcutils.CreateSeries(ds_60minutes,Wd_label,Wd,Flag=f,Attr=attr)
# relative humidity from temperature, specific humidity and pressure
for i in range(0,3):
    for j in range(0,3):
        q_label = "q_"+str(i)+str(j)
        Ta_label = "Ta_"+str(i)+str(j)
        ps_label = "ps_"+str(i)+str(j)
        RH_label = "RH_"+str(i)+str(j)
        q,f = qcutils.GetSeriesasMA(ds_60minutes,q_label)
        Ta,f = qcutils.GetSeriesasMA(ds_60minutes,Ta_label)
        ps,f = qcutils.GetSeriesasMA(ds_60minutes,ps_label)
        RH = mf.RHfromspecifichumidity(q, Ta, ps)
        attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='not defined')
        qcutils.CreateSeries(ds_60minutes,RH_label,RH,Flag=f,Attr=attr)
# absolute humidity from temperature and relative humidity
for i in range(0,3):
    for j in range(0,3):
        Ta_label = "Ta_"+str(i)+str(j)
        RH_label = "RH_"+str(i)+str(j)
        Ah_label = "Ah_"+str(i)+str(j)
        Ta,f = qcutils.GetSeriesasMA(ds_60minutes,Ta_label)
        RH,f = qcutils.GetSeriesasMA(ds_60minutes,RH_label)
        Ah = mf.absolutehumidityfromRH(Ta, RH)
        attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity',units='g/m3',standard_name='not defined')
        qcutils.CreateSeries(ds_60minutes,Ah_label,Ah,Flag=f,Attr=attr)
# soil moisture from kg/m2 to m3/m3
attr = qcutils.GetAttributeDictionary(ds_60minutes,"Sws_00")
for i in range(0,3):
    for j in range(0,3):
        label = "Sws_"+str(i)+str(j)
        Sws,f = qcutils.GetSeriesasMA(ds_60minutes,label)
        Sws = Sws/float(100)
        attr["units"] = "frac"
        qcutils.CreateSeries(ds_60minutes,label,Sws,Flag=flag_60minutes,Attr=attr)
# net radiation and upwelling short and long wave radiation
for i in range(0,3):
    for j in range(0,3):
        label_Fn = "Fn_"+str(i)+str(j)
        label_Fsd = "Fsd_"+str(i)+str(j)
        label_Fld = "Fld_"+str(i)+str(j)
        label_Fsu = "Fsu_"+str(i)+str(j)
        label_Flu = "Flu_"+str(i)+str(j)
        label_Fn_sw = "Fn_sw_"+str(i)+str(j)
        label_Fn_lw = "Fn_lw_"+str(i)+str(j)
        Fsd,f = qcutils.GetSeriesasMA(ds_60minutes,label_Fsd)
        Fld,f = qcutils.GetSeriesasMA(ds_60minutes,label_Fld)
        Fn_sw,f = qcutils.GetSeriesasMA(ds_60minutes,label_Fn_sw)
        Fn_lw,f = qcutils.GetSeriesasMA(ds_60minutes,label_Fn_lw)
        Fsu = Fsd - Fn_sw
        Flu = Fld - Fn_lw
        Fn = (Fsd-Fsu)+(Fld-Flu)
        attr = qcutils.MakeAttributeDictionary(long_name='Up-welling long wave',
                             standard_name='surface_upwelling_longwave_flux_in_air',units='W/m2')
        qcutils.CreateSeries(ds_60minutes,label_Flu,Flu,Flag=f,Attr=attr)
        attr = qcutils.MakeAttributeDictionary(long_name='Up-welling short wave',
                             standard_name='surface_upwelling_shortwave_flux_in_air',units='W/m2')
        qcutils.CreateSeries(ds_60minutes,label_Fsu,Fsu,Flag=f,Attr=attr)
        attr = qcutils.MakeAttributeDictionary(long_name='Calculated net radiation',
                             standard_name='surface_net_allwave_radiation',units='W/m2')
        qcutils.CreateSeries(ds_60minutes,label_Fn,Fn,Flag=f,Attr=attr)
# ground heat flux as residual
for i in range(0,3):
    for j in range(0,3):
        label_Fg = "Fg_"+str(i)+str(j)
        label_Fn = "Fn_"+str(i)+str(j)
        label_Fh = "Fh_"+str(i)+str(j)
        label_Fe = "Fe_"+str(i)+str(j)
        Fn,f = qcutils.GetSeriesasMA(ds_60minutes,label_Fn)
        Fh,f = qcutils.GetSeriesasMA(ds_60minutes,label_Fh)
        Fe,f = qcutils.GetSeriesasMA(ds_60minutes,label_Fe)
        Fg = Fn - Fh - Fe
        attr = qcutils.MakeAttributeDictionary(long_name='Calculated ground heat flux',
                             standard_name='downward_heat_flux_in_soil',units='W/m2')
        qcutils.CreateSeries(ds_60minutes,label_Fg,Fg,Flag=f,Attr=attr)
# Available energy
for i in range(0,3):
    for j in range(0,3):
        label_Fg = "Fg_"+str(i)+str(j)
        label_Fn = "Fn_"+str(i)+str(j)
        label_Fa = "Fa_"+str(i)+str(j)
        Fn,f = qcutils.GetSeriesasMA(ds_60minutes,label_Fn)
        Fg,f = qcutils.GetSeriesasMA(ds_60minutes,label_Fg)
        Fa = Fn - Fg
        attr = qcutils.MakeAttributeDictionary(long_name='Calculated available energy',
                             standard_name='not defined',units='W/m2')
        qcutils.CreateSeries(ds_60minutes,label_Fa,Fa,Flag=f,Attr=attr)

# interpolate from 60 to 30 minutes if requested
if qcutils.cfoptionskey(cf,"Interpolate"):
    # copy the global attributes
    for this_attr in ds_60minutes.globalattributes.keys():
        ds_30minutes.globalattributes[this_attr] = ds_60minutes.globalattributes[this_attr]
    # update the global attribute "time_step"
    ds_30minutes.globalattributes["time_step"] = 30
    # generate the 30 minute datetime series
    dt_loc_30minutes = [x for x in perdelta(dt_loc_60minutes[0],dt_loc_60minutes[-1],datetime.timedelta(minutes=30))]
    nRecs_30minutes = len(dt_loc_30minutes)
    # update the global attribute "nc_nrecs"
    ds_30minutes.globalattributes['nc_nrecs'] = nRecs_30minutes
    flag_30minutes = numpy.zeros(nRecs_30minutes)
    ds_30minutes.series["DateTime"] = {}
    ds_30minutes.series["DateTime"]["Data"] = dt_loc_30minutes
    # get the year, month etc from the datetime
    xl_date = qcutils.get_xldate_from_datetime(dt_loc_30minutes)
    attr = qcutils.MakeAttributeDictionary(long_name="Date/time in Excel format",units="days since 1899-12-31 00:00:00")
    qcutils.CreateSeries(ds_30minutes,"xlDateTime",xl_date,Flag=flag_30minutes,Attr=attr)
    qcutils.get_ymdhms_from_datetime(ds_30minutes)
    # interpolate to 30 minutes
    x_60minutes = numpy.arange(0,nRecs,1)
    x_30minutes = numpy.arange(0,nRecs-0.5,0.5)
    varlist = ds_60minutes.series.keys()
    # strip out the date and time variables already done
    for this_one in ["DateTime","xlDateTime","Year","Month","Day","Hour","Minute","Second"]:
        if this_one in varlist: varlist.remove(this_one)
    for label in varlist:
        series_60minutes,f = qcutils.GetSeriesasMA(ds_60minutes,label)
        f = interp1d(x_60minutes,series_60minutes)
        series_30minutes = f(x_30minutes)
        attr = qcutils.GetAttributeDictionary(ds_60minutes,label)
        qcutils.CreateSeries(ds_30minutes,label,series_30minutes,Flag=flag_30minutes,Attr=attr)
    # now write out the ACCESS data interpolated to 30 minutes
    ncfile = qcio.nc_open_write(outfilename)
    qcio.nc_write_series(ncfile, ds_30minutes)
else:
    # write out the ACCESS data
    ncfile = qcio.nc_open_write(outfilename)
    qcio.nc_write_series(ncfile, ds_60minutes)
