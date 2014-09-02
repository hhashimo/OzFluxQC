import sys
sys.path.append('../scripts')
import constants as c
import glob
import meteorologicalfunctions as mf
import netCDF4
import numpy
import qcio
import qcutils

# open the logging file
log = qcutils.startlog('bios2nc','../logfiles/bios2nc.log')

# get the control file
cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0: sys.exit()
var_list = cf["Variables"].keys()
site_list = cf["Sites"].keys()
for site in site_list:
    log.info("Starting site: "+site)
    # get a data structure
    ds = qcio.DataStructure()
    # get the input file mask
    infilename = cf["Sites"][site]["in_filepath"]+cf["Sites"][site]["in_filename"]
    outfilename = cf["Sites"][site]["out_filepath"]+cf["Sites"][site]["out_filename"]
    # interpolate to 30 minutes or not
    interpolate = True
    if not cf["Sites"][site].as_bool("interpolate"): interpolate = False
    # get the site time zone
    site_timezone = cf["Sites"][site]["site_timezone"]
    # read the BIOS file
    bios_ncfile = netCDF4.Dataset(infilename)
    time = bios_ncfile.variables["time"][:]
    nRecs = len(time)
    flag = numpy.zeros(nRecs)
    # set some global attributes
    ds.globalattributes["time_step"] = 30
    ds.globalattributes["time_zone"] = site_timezone
    ds.globalattributes["nc_nrecs"] = nRecs
    time_units = getattr(bios_ncfile.variables["time"],"units")
    qcutils.get_datetimefromnctime(ds,time,time_units)
    qcutils.get_ymdhms_from_datetime(ds)
    xl_date_loc = qcutils.get_xldate_from_datetime(ds.series["DateTime"]["Data"])
    attr = qcutils.MakeAttributeDictionary(long_name="Date/time (local) in Excel format",units="days since 1899-12-31 00:00:00")
    qcutils.CreateSeries(ds,"xlDateTime",xl_date_loc,Flag=flag,Attr=attr)
    # get the data
    for label in var_list:
        bios_name = cf["Variables"][label]["bios_name"]
        if len(bios_ncfile.variables[bios_name].shape)==1:
            #print label+" has 1 dimension"
            data = bios_ncfile.variables[bios_name][:]
        elif len(bios_ncfile.variables[bios_name].shape)==2:
            #print label+" has 2 dimensions"
            data = bios_ncfile.variables[bios_name][:,0]
        elif len(bios_ncfile.variables[bios_name].shape)==3:
            #print label+" has 3 dimensions"
            data = bios_ncfile.variables[bios_name][:,0,0]
        attr = {}
        for this_attr in bios_ncfile.variables[bios_name].ncattrs():
            attr[this_attr] = getattr(bios_ncfile.variables[bios_name],this_attr)
        attr["missing_value"] = c.missing_value
        qcutils.CreateSeries(ds,label,data,Flag=flag,Attr=attr)
    # close the netCDF file
    bios_ncfile.close()
    # convert precipitation from kg/m2/s to mm/30 minutes
    precip,flag,attr = qcutils.GetSeriesasMA(ds,"Precip")
    precip = float(1800)*precip
    attr["units"] = "mm"
    qcutils.CreateSeries(ds,"Precip",precip,Flag=flag,Attr=attr)
    # convert Ta from K to C
    Ta,flag,attr = qcutils.GetSeriesasMA(ds,"Ta")
    Ta = Ta - c.C2K
    attr["units"] = "C"
    qcutils.CreateSeries(ds,"Ta",Ta,Flag=flag,Attr=attr)
    # convert Ts from K to C
    Ts,flag,attr = qcutils.GetSeriesasMA(ds,"Ts")
    Ts = Ts - c.C2K
    attr["units"] = "C"
    qcutils.CreateSeries(ds,"Ts",Ts,Flag=flag,Attr=attr)
    # convert ps from hPa to kPa
    ps,flag,attr = qcutils.GetSeriesasMA(ds,"ps")
    ps = ps/float(10)
    attr["units"] = "kPa"
    qcutils.CreateSeries(ds,"ps",ps,Flag=flag,Attr=attr)
    # calculate relative humidity
    q,f,a = qcutils.GetSeriesasMA(ds,"q")
    Ta,f,a = qcutils.GetSeriesasMA(ds,"Ta")
    ps,f,a = qcutils.GetSeriesasMA(ds,"ps")
    RH = mf.RHfromspecifichumidity(q, Ta, ps)
    attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='not defined')
    qcutils.CreateSeries(ds,"RH",RH,Flag=flag,Attr=attr)
    # calculate absolute humidity
    Ta,f,a = qcutils.GetSeriesasMA(ds,"Ta")
    RH,f,a = qcutils.GetSeriesasMA(ds,"RH")
    Ah = mf.absolutehumidityfromRH(Ta, RH)
    attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity',units='g/m3',standard_name='not defined')
    qcutils.CreateSeries(ds,"Ah",Ah,Flag=flag,Attr=attr)
    # calculate net radiation
    Fsd,f,a = qcutils.GetSeriesasMA(ds,"Fsd")
    Fld,f,a = qcutils.GetSeriesasMA(ds,"Fld")
    Fn_sw,f,a = qcutils.GetSeriesasMA(ds,"Fn_sw")
    Fn_lw,f,a = qcutils.GetSeriesasMA(ds,"Fn_lw")
    Fsu = Fsd - Fn_sw
    Flu = Fld - Fn_lw
    Fn = (Fsd-Fsu)+(Fld-Flu)
    attr = qcutils.MakeAttributeDictionary(long_name='Up-welling long wave',
                         standard_name='surface_upwelling_longwave_flux_in_air',units='W/m2')
    qcutils.CreateSeries(ds,"Flu",Flu,Flag=flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Up-welling short wave',
                         standard_name='surface_upwelling_shortwave_flux_in_air',units='W/m2')
    qcutils.CreateSeries(ds,"Fsu",Fsu,Flag=flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Calculated net radiation',
                         standard_name='surface_net_allwave_radiation',units='W/m2')
    qcutils.CreateSeries(ds,"Fn",Fn,Flag=flag,Attr=attr)
    # calculate available energy
    Fn,f,a = qcutils.GetSeriesasMA(ds,"Fn")
    Fg,f,a = qcutils.GetSeriesasMA(ds,"Fg")
    Fa = Fn - Fg
    attr = qcutils.MakeAttributeDictionary(long_name='Calculated available energy',
                         standard_name='not defined',units='W/m2')
    qcutils.CreateSeries(ds,"Fa",Fa,Flag=flag,Attr=attr)
    # write the output file
    ncfile = qcio.nc_open_write(outfilename)
    qcio.nc_write_series(ncfile,ds,ndims=1)
    log.info("Finished site: "+site)

print "All done"