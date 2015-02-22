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
import glob
import logging
import matplotlib.pyplot as plt
import meteorologicalfunctions as mf
from netCDF4 import MFDataset, Dataset
import numpy
import pytz
import qcio
import qcutils
from scipy.interpolate import interp1d

class ACCESSData(object):
    def __init__(self):
        self.globalattr = {}
        self.variables = {}
        self.varattr = {}

def access_read_mfiles(file_list=[],var_list=[]):
    #ds = qcio.DataStructure()
    f = ACCESSData()
    #print file_list[0]
    ncfile = Dataset(file_list[0])
    # get the global attributes
    for gattr in ncfile.ncattrs():
        f.globalattr[gattr] = getattr(ncfile,gattr)
    dims = ncfile.dimensions
    shape = (len(dims["time"]),len(dims["lat"]),len(dims["lon"]))
    if len(var_list)==0: var_list=ncfile.variables.keys()
    # load the data from the first file into the data structure
    for var in var_list:
        access_name = cf["Variables"][var]["access_name"]
        if access_name in ncfile.variables.keys():
            f.variables[access_name] = ncfile.variables[access_name][:]
            f.varattr[access_name] = {}
            for this_attr in ncfile.variables[access_name].ncattrs():
                f.varattr[access_name][this_attr] = getattr(ncfile.variables[access_name],this_attr)
        else:
            print "ACCESS variable "+access_name+" not found in "+file_list[0]
            f.variables[access_name] = makedummyseries(shape)
    ncfile.close()
    if len(file_list)>1:
        # loop over the remaining files and append the data
        for filename in file_list[1:]:
            #print filename
            ncfile = Dataset(filename)
            for var in var_list:
                access_name = cf["Variables"][var]["access_name"]
                #print var,ncfile.variables[access_name].shape
                if access_name in ncfile.variables.keys():
                    f.variables[access_name] = numpy.concatenate((f.variables[access_name],ncfile.variables[access_name][:]),axis=0)
                    if access_name not in f.varattr.keys():
                        f.varattr[access_name] = {}
                        for this_attr in ncfile.variables[access_name].ncattrs():
                            f.varattr[access_name][this_attr] = getattr(ncfile.variables[access_name],this_attr) 
                else:
                    print "ACCESS variable "+access_name+" not found in "+filename
                    empty_series = makedummyseries(shape)
                    f.variables[access_name] = numpy.concatenate((f.variables[access_name],empty_series),axis=0)
            ncfile.close()
    return f

def makedummyseries(shape):
    return numpy.ma.masked_all(shape)

def perdelta(start,end,delta):
    curr = start
    while curr <= end:
        yield curr
        curr += delta

def find_indices(a,b):
    len_a = len(a)
    len_b = len(b)
    indices = []
    idx = -1
    if len_a>=len_b:
        for item in b:
            idx = a.index(item,idx+1)
            indices.append(idx)
    else:
        for item in a:
            idx = b.index(item,idx+1)
            indices.append(idx)
    return indices

# open the logging file
#logging.basicConfig()
log = qcutils.startlog('access2nc','../logfiles/access2nc.log')

# get the control file
cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0: sys.exit()
var_list = cf["Variables"].keys()
site_list = cf["Sites"].keys()
for site in site_list:
    var_list = cf["Variables"].keys()
    log.info("Starting site: "+site)
    # get the input file mask
    infilename = cf["Sites"][site]["in_filepath"]+cf["Sites"][site]["in_filename"]
    # get a sorted list of files that match the mask in the control file
    file_list = sorted(glob.glob(infilename))
    # number of files to process
    nFiles = len(file_list)
    log.info("Processing "+str(nFiles)+" files")
    # output file name
    outfilename = cf["Sites"][site]["out_filepath"]+cf["Sites"][site]["out_filename"]
    # interpolate to 30 minutes or not
    interpolate = True
    if not cf["Sites"][site].as_bool("interpolate"): interpolate = False
    # get the site timezone
    site_timezone = cf["Sites"][site]["site_timezone"]
    site_tz = pytz.timezone(site_timezone)
    # create an instance of the data structure
    ds_60minutes = qcio.DataStructure()
    ds_30minutes = qcio.DataStructure()
    # read the netCDF files
    #f = MFDataset(infilename)
    f = access_read_mfiles(file_list=file_list,var_list=var_list)
    # get the date and time
    valid_date = f.variables["valid_date"][:]
    valid_time = f.variables["valid_time"][:]
    # get the number of records in the file
    nRecs = len(valid_date)
    # check the files have the expected number of records
    if nRecs!=25*nFiles:
        log.error("Number of records expected was "+str(25*nFiles)+", got "+str(nRecs))
        continue
    # check that every 25th record is at midnight
    index = range(25,nRecs-25,25)
    max_vt=max(valid_time[index])
    min_vt=min(valid_time[index])
    if (max_vt!=0) or (min_vt!=0):
        log.error("Got unexpected time on n*25th record")
        continue
    # seems safe to proceed
    # map the ACCESS file global attributes to the OzFluxQC file
    for attr in f.globalattr:
        ds_60minutes.globalattributes[attr] = f.globalattr[attr]
    # add a value for the Excel datemode, assume Windows date system (1900)
    ds_60minutes.globalattributes['xl_datemode'] = str(0)
    # number of records
    ds_60minutes.globalattributes["nc_nrecs"] = nRecs
    # put the ACCESS data into the 60 minute data structure ds_60minutes
    # make a QC flag with a value of 0
    flag_60minutes = numpy.zeros(nRecs)
    # loop over the variables defined in the control file
    for item in ["valid_date","valid_time"]:
        if item in var_list: var_list.remove(item)
    for label in var_list:
        # get the name of the ACCESS variable
        access_name = cf["Variables"][label]["access_name"]
        if access_name not in f.variables.keys():
            log.error("Requested variable "+access_name+" not found in ACCESS data")
            continue
        attr = {}
        for this_attr in f.varattr[access_name].keys():
            attr[this_attr] = f.varattr[access_name][this_attr]
        attr["missing_value"] = c.missing_value
        # loop over all ACCESS grids and give them standard OzFlux names with the grid idices appended
        for i in range(0,3):
            for j in range(0,3):
                if len(f.variables[access_name].shape)==3:
                    label_ij = label+'_'+str(i)+str(j)
                    series = f.variables[access_name][:,i,j]
                    qcutils.CreateSeries(ds_60minutes,label_ij,series,Flag=flag_60minutes,Attr=attr)
                elif len(f.variables[access_name].shape)==4:
                    label_ij = label+'_'+str(i)+str(j)
                    series = f.variables[access_name][:,0,i,j]
                    qcutils.CreateSeries(ds_60minutes,label_ij,series,Flag=flag_60minutes,Attr=attr)
                else:
                    print "Unrecognised variable ("+label+") dimension in ACCESS file"
                    #sys.exit()
    # trap valid_date==0 occurrences, these happened in some of the files produced
    # in the second batch while the accum_prcp was being sorted out
    index = numpy.where(valid_date==0)[0]
    # if there are some valid_date==0 then replace them with the correct date
    if len(index)!=0:
        for i in index:
            dt = datetime.datetime.strptime(str(int(valid_date[i-1])),"%Y%m%d")
            dt = dt+datetime.timedelta(days=1)
            valid_date[i] = datetime.datetime.strftime(dt,"%Y%m%d")
    # replace missing days with blocks of missing data
    
    # copy the precipitation from the n*24th record to the n*25th record
    # The daily files of ACCESS data provided by Lawrie start at 0000 UTC, have 25 records and
    # end at 0000 UTC on the next day.  The 25th record in the file is included to get the
    # accumulated precipitation ("accum_prcp") between 2300 and 0000.  However, all other data
    # in the 25th record is unreliable and must be discarded.
    # So, after concatenating the daily files (done by MFDataset) we have 2 records for each
    # midnight, the first with the correct precipitation data and the second with everything
    # else.  We handle this by copying the precipitation data from the n*24th record to the
    # n*25th record (thereby putting the correct precipitation data with correct everything
    # else) and then remove the n*24th records.
    ind_mn1=range(24,nRecs,25)
    ind_mn2=range(25,nRecs,25)
    series_list=ds_60minutes.series.keys()
    precip_list = [x for x in series_list if "Precip" in x]
    for precip in precip_list:
        data,flag,attr = qcutils.GetSeriesasMA(ds_60minutes,precip)
        data[ind_mn2]=data[ind_mn1]
        flag[ind_mn2]=flag[ind_mn1]
        qcutils.CreateSeries(ds_60minutes,precip,data,Flag=flag,Attr=attr)
    # get the index to select only the first 24 records per day
    base = range(0,24)
    index = []
    for i in range(0,nFiles):
        index.extend([i*25+j for j in base])
    # loop over the data series in ds_60minutes and extract only the first 24 records for each day
    valid_date = valid_date[index]
    valid_time = valid_time[index]
    for label in ds_60minutes.series.keys():
        data,flag,attr = qcutils.GetSeriesasMA(ds_60minutes,label)
        data = data[index]
        flag = flag[index]
        qcutils.CreateSeries(ds_60minutes,label,data,Flag=flag,Attr=attr)
    nRecs = len(index)
    # set some global attributres
    ds_60minutes.globalattributes["nc_nrecs"] = nRecs
    ds_60minutes.globalattributes["time_step"] = 60
    ds_60minutes.globalattributes["time_zone"] = site_timezone
    ds_60minutes.globalattributes["site_name"] = cf["Sites"][site]["site_name"]
    # now getthe UTC and local datetimes from the ACCESS valid_date and valid_time fields
    #dt=[datetime.datetime.strptime(str(valid_date[i]*10000+valid_time[i]),"%Y%m%d%H%M") for i in range(0,nRecs)]
    dt_utc_60minutes=[datetime.datetime.strptime(str(int(valid_date[i])*10000+int(valid_time[i])),"%Y%m%d%H%M") for i in range(0,len(valid_date))]
    # make utc_dt timezone aware
    dt_utc_60minutes=[x.replace(tzinfo=pytz.utc) for x in dt_utc_60minutes]
    # get local time from UTC
    dt_loc_60minutes=[x.astimezone(site_tz) for x in dt_utc_60minutes]
    # NOTE: will have to disable daylight saving at some stage, towers stay on Standard Time
    # PRI hopes that the following line will do this ...
    dt_loc_60minutes=[x-x.dst() for x in dt_loc_60minutes]
    # make local time timezone naive to match datetimes in OzFluxQC
    dt_loc_60minutes=[x.replace(tzinfo=None) for x in dt_loc_60minutes]
    ds_60minutes.series["DateTime"] = {}
    ds_60minutes.series["DateTime"]["Data"] = dt_loc_60minutes
    ds_60minutes.series["DateTime_UTC"] = {}
    ds_60minutes.series["DateTime_UTC"]["Data"] = dt_utc_60minutes
    # get the year, month etc from the datetime
    flag_60minutes = numpy.zeros(nRecs,dtype=numpy.int32)
    ds_60minutes.series["DateTime"]["Flag"] = flag_60minutes
    qcutils.get_xldatefromdatetime(ds_60minutes)
    qcutils.get_ymdhmsfromdatetime(ds_60minutes)
    # get derived quantities and adjust units
    # air temperature from K to C
    attr = qcutils.GetAttributeDictionary(ds_60minutes,"Ta_00")
    if attr["units"] == "K":
        for i in range(0,3):
            for j in range(0,3):
                label = "Ta_"+str(i)+str(j)
                Ta,f,a = qcutils.GetSeriesasMA(ds_60minutes,label)
                Ta = Ta - c.C2K
                attr["units"] = "C"
                qcutils.CreateSeries(ds_60minutes,label,Ta,Flag=flag_60minutes,Attr=attr)
    # soil temperature from K to C
    attr = qcutils.GetAttributeDictionary(ds_60minutes,"Ts_00")
    if attr["units"] == "K":
        for i in range(0,3):
            for j in range(0,3):
                label = "Ts_"+str(i)+str(j)
                Ts,f,a = qcutils.GetSeriesasMA(ds_60minutes,label)
                Ts = Ts - c.C2K
                attr["units"] = "C"
                qcutils.CreateSeries(ds_60minutes,label,Ts,Flag=flag_60minutes,Attr=attr)
    # pressure from Pa to kPa
    attr = qcutils.GetAttributeDictionary(ds_60minutes,"ps_00")
    if attr["units"] == "Pa":
        for i in range(0,3):
            for j in range(0,3):
                label = "ps_"+str(i)+str(j)
                ps,f,a = qcutils.GetSeriesasMA(ds_60minutes,label)
                ps = ps/float(1000)
                attr["units"] = "kPa"
                qcutils.CreateSeries(ds_60minutes,label,ps,Flag=flag_60minutes,Attr=attr)
    # wind speed from components
    for i in range(0,3):
        for j in range(0,3):
            u_label = "u_"+str(i)+str(j)
            v_label = "v_"+str(i)+str(j)
            Ws_label = "Ws_"+str(i)+str(j)
            u,f,a = qcutils.GetSeriesasMA(ds_60minutes,u_label)
            v,f,a = qcutils.GetSeriesasMA(ds_60minutes,v_label)
            Ws = numpy.sqrt(u*u+v*v)
            attr = qcutils.MakeAttributeDictionary(long_name="Wind speed",units="m/s",height="10m")
            qcutils.CreateSeries(ds_60minutes,Ws_label,Ws,Flag=f,Attr=attr)
    # wind direction from components
    for i in range(0,3):
        for j in range(0,3):
            u_label = "u_"+str(i)+str(j)
            v_label = "v_"+str(i)+str(j)
            Wd_label = "Wd_"+str(i)+str(j)
            u,f,a = qcutils.GetSeriesasMA(ds_60minutes,u_label)
            v,f,a = qcutils.GetSeriesasMA(ds_60minutes,v_label)
            Wd = float(270) - numpy.ma.arctan2(v,u)*float(180)/numpy.pi
            index = numpy.ma.where(Wd>360)[0]
            if len(index)>0: Wd[index] = Wd[index] - float(360)
            attr = qcutils.MakeAttributeDictionary(long_name="Wind direction",units="degrees",height="10m")
            qcutils.CreateSeries(ds_60minutes,Wd_label,Wd,Flag=f,Attr=attr)
    # relative humidity from temperature, specific humidity and pressure
    for i in range(0,3):
        for j in range(0,3):
            q_label = "q_"+str(i)+str(j)
            Ta_label = "Ta_"+str(i)+str(j)
            ps_label = "ps_"+str(i)+str(j)
            RH_label = "RH_"+str(i)+str(j)
            q,f,a = qcutils.GetSeriesasMA(ds_60minutes,q_label)
            Ta,f,a = qcutils.GetSeriesasMA(ds_60minutes,Ta_label)
            ps,f,a = qcutils.GetSeriesasMA(ds_60minutes,ps_label)
            RH = mf.RHfromspecifichumidity(q, Ta, ps)
            attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='not defined')
            qcutils.CreateSeries(ds_60minutes,RH_label,RH,Flag=f,Attr=attr)
    # absolute humidity from temperature and relative humidity
    for i in range(0,3):
        for j in range(0,3):
            Ta_label = "Ta_"+str(i)+str(j)
            RH_label = "RH_"+str(i)+str(j)
            Ah_label = "Ah_"+str(i)+str(j)
            Ta,f,a = qcutils.GetSeriesasMA(ds_60minutes,Ta_label)
            RH,f,a = qcutils.GetSeriesasMA(ds_60minutes,RH_label)
            Ah = mf.absolutehumidityfromRH(Ta, RH)
            attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity',units='g/m3',standard_name='not defined')
            qcutils.CreateSeries(ds_60minutes,Ah_label,Ah,Flag=f,Attr=attr)
    # soil moisture from kg/m2 to m3/m3
    attr = qcutils.GetAttributeDictionary(ds_60minutes,"Sws_00")
    for i in range(0,3):
        for j in range(0,3):
            label = "Sws_"+str(i)+str(j)
            Sws,f,a = qcutils.GetSeriesasMA(ds_60minutes,label)
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
            Fsd,f,a = qcutils.GetSeriesasMA(ds_60minutes,label_Fsd)
            Fld,f,a = qcutils.GetSeriesasMA(ds_60minutes,label_Fld)
            Fn_sw,f,a = qcutils.GetSeriesasMA(ds_60minutes,label_Fn_sw)
            Fn_lw,f,a = qcutils.GetSeriesasMA(ds_60minutes,label_Fn_lw)
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
            Fn,f,a = qcutils.GetSeriesasMA(ds_60minutes,label_Fn)
            Fh,f,a = qcutils.GetSeriesasMA(ds_60minutes,label_Fh)
            Fe,f,a = qcutils.GetSeriesasMA(ds_60minutes,label_Fe)
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
            Fn,f,a = qcutils.GetSeriesasMA(ds_60minutes,label_Fn)
            Fg,f,a = qcutils.GetSeriesasMA(ds_60minutes,label_Fg)
            Fa = Fn - Fg
            attr = qcutils.MakeAttributeDictionary(long_name='Calculated available energy',
                                 standard_name='not defined',units='W/m2')
            qcutils.CreateSeries(ds_60minutes,label_Fa,Fa,Flag=f,Attr=attr)

    # dump to an Excel file so we can see what is going on
    xlfullname= outfilename.replace('.nc','.xls')
    log.info("Writing to file: "+xlfullname)
    qcio.xl_write_series(ds_60minutes, xlfullname, outputlist=None)

    # check for time gaps in the file
    log.info("Checking for time gaps")
    has_gaps = qcutils.CheckTimeStep(ds_60minutes,mode="fix")
    
    # dump corrected data to an Excel file
    #xlfullname= xlfullname.replace('.xls','_nogaps.xls')
    #log.info("Writing to file: "+xlfullname)
    #qcio.xl_write_series(ds_60minutes, xlfullname, outputlist=None)

    # interpolate from 60 to 30 minutes if requested
    if interpolate:
        log.info("Interpolating site "+site+" to 30 minute time step")
        # copy the global attributes
        for this_attr in ds_60minutes.globalattributes.keys():
            ds_30minutes.globalattributes[this_attr] = ds_60minutes.globalattributes[this_attr]
        # update the global attribute "time_step"
        ds_30minutes.globalattributes["time_step"] = 30
        # generate the 30 minute datetime series
        dt_loc_30minutes = [x for x in perdelta(dt_loc_60minutes[0],dt_loc_60minutes[-1],datetime.timedelta(minutes=30))]
        dt_utc_30minutes = [x for x in perdelta(dt_utc_60minutes[0],dt_utc_60minutes[-1],datetime.timedelta(minutes=30))]
        nRecs_30minutes = len(dt_loc_30minutes)
        # update the global attribute "nc_nrecs"
        ds_30minutes.globalattributes['nc_nrecs'] = nRecs_30minutes
        flag_30minutes = numpy.zeros(nRecs_30minutes)
        ds_30minutes.series["DateTime"] = {}
        ds_30minutes.series["DateTime"]["Data"] = dt_loc_30minutes
        ds_30minutes.series["DateTime_UTC"] = {}
        ds_30minutes.series["DateTime_UTC"]["Data"] = dt_utc_30minutes
        # get the year, month etc from the datetime
        qcutils.get_xldatefromdatetime(ds_30minutes)
        qcutils.get_ymdhmsfromdatetime(ds_30minutes)
        # interpolate to 30 minutes
        nRecs_60 = len(ds_60minutes.series["DateTime"]["Data"])
        nRecs_30 = len(ds_30minutes.series["DateTime"]["Data"])
        x_60minutes = numpy.arange(0,nRecs_60,1)
        x_30minutes = numpy.arange(0,nRecs_60-0.5,0.5)
        varlist_60 = ds_60minutes.series.keys()
        # strip out the date and time variables already done
        for item in ["DateTime","DateTime_UTC","xlDateTime","Year","Month","Day","Hour","Minute","Second","Hdh"]:
            if item in varlist_60: varlist_60.remove(item)
        # now do the interpolation (its OK to interpolate accumulated precipitation)
        for label in varlist_60:
            series_60minutes,flag,attr = qcutils.GetSeries(ds_60minutes,label)
            ci_60minutes = numpy.zeros(len(series_60minutes))
            index = numpy.where(abs(series_60minutes-float(c.missing_value))<c.eps)[0]
            ci_60minutes[index] = float(1)
            int_fn = interp1d(x_60minutes,series_60minutes)
            series_30minutes = int_fn(x_30minutes)
            int_fn = interp1d(x_60minutes,ci_60minutes)
            ci_30minutes = int_fn(x_30minutes)
            index = numpy.where(abs(ci_30minutes-float(0))>c.eps)[0]
            series_30minutes[index] = numpy.float64(c.missing_value)
            qcutils.CreateSeries(ds_30minutes,label,series_30minutes,Flag=flag_30minutes,Attr=attr)
        # now get precipitation per time step from the interpolated precipitation accumulated over the day
        for i in range(0,3):
            for j in range(0,3):
                label = "Precip_"+str(i)+str(j)
                accum_24hr,flag,attr = qcutils.GetSeriesasMA(ds_30minutes,label)
                index = numpy.ma.where(accum_24hr<0.001)[0]
                accum_24hr[index] = float(0)
                precip = numpy.ma.ediff1d(accum_24hr,to_begin=0)
                index = [x for x in range(len(dt_utc_30minutes)) if (dt_utc_30minutes[x].hour==0) and (dt_utc_30minutes[x].minute==30)]
                precip[index] = float(0)
                index = [x for x in range(len(dt_utc_30minutes)) if (dt_utc_30minutes[x].hour==1) and (dt_utc_30minutes[x].minute==0)]
                precip[index] = accum_24hr[index]
                attr["long_name"] = "Precipitation total over time step"
                attr["units"] = "mm/30 minutes"
                qcutils.CreateSeries(ds_30minutes,label,precip,Flag=flag_30minutes,Attr=attr)
        # now write out the ACCESS data interpolated to 30 minutes
        ncfile = qcio.nc_open_write(outfilename)
        qcio.nc_write_series(ncfile, ds_30minutes,ndims=1)
        log.info("Finished site : "+site)
    else:
        # now get precipitation per time step from the precipitation accumulated over the day
        dt_utc_60minutes=ds_60minutes.series["DateTime_UTC"]["Data"]
        idx_0100 = [x for x in range(len(dt_utc_60minutes)) if (dt_utc_60minutes[x].hour==1) and (dt_utc_60minutes[x].minute==0)]
        for i in range(0,3):
            for j in range(0,3):
                label = "Precip_"+str(i)+str(j)
                accum_24hr,flag,attr = qcutils.GetSeriesasMA(ds_60minutes,label)
                index = numpy.ma.where(accum_24hr<0.001)[0]
                accum_24hr[index] = float(0)
                precip = numpy.ma.ediff1d(accum_24hr,to_begin=0)
                precip[idx_0100] = accum_24hr[idx_0100]
                attr["long_name"] = "Precipitation total over time step"
                attr["units"] = "mm/hr"
                qcutils.CreateSeries(ds_60minutes,label,precip,Flag=flag,Attr=attr)
        # write out the ACCESS data
        ncfile = qcio.nc_open_write(outfilename)
        qcio.nc_write_series(ncfile, ds_60minutes,ndims=1)
        log.info("Finished site : "+site)
