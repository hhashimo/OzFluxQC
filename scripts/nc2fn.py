import copy
import csv
import datetime
import dateutil
import meteorologicalfunctions as mf
import numpy
import qcio
import qcutils
import sys

def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta

# open the logging file
log = qcutils.startlog('ncconcat','../logfiles/ncconcat.log')

# get the control file contents
# was there an argument on the command line?
if len(sys.argv)>1:
    try:
        # was it a control file name?
        cf = qcio.get_controlfilecontents(sys.argv[1])
    except:
        # oh, well, let's do it the hard way
        cf = qcio.load_controlfile(path='../controlfiles')
else:
    cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0: sys.exit()
# get the file names
ncFileName = qcio.get_infilename_from_cf(cf)
csvFileName = qcio.get_outfilename_from_cf(cf)
# open the csv file
csvfile = open(csvFileName,'wb')
writer = csv.writer(csvfile)
# read the netCDF file
ds = qcio.nc_read_series(ncFileName)
# Tumbarumba doesn't have RH in the netCDF files
if "RH" not in ds.series.keys():
    Ah,f,a = qcutils.GetSeriesasMA(ds,'Ah')
    Ta,f,a = qcutils.GetSeriesasMA(ds,'Ta')
    RH = mf.RHfromabsolutehumidity(Ah, Ta)
    attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='not defined')
    qcutils.CreateSeries(ds,"RH",RH,FList=['Ta','Ah'],Attr=attr)
ts = int(ds.globalattributes["time_step"])
ts_delta = datetime.timedelta(minutes=ts)
# get the datetime series
dt = ds.series["DateTime"]["Data"]
# check the start datetime of the series and adjust if necessary
start_datetime = dateutil.parser.parse(str(cf["General"]["start_datetime"]))
if dt[0]<start_datetime:
    # requested start_datetime is after the start of the file
    print "nc2fn: truncating start of file"
    si = qcutils.GetDateIndex(dt,str(start_datetime),ts=ts,match="exact")
    for thisone in ds.series.keys():
        ds.series[thisone]["Data"] = ds.series[thisone]["Data"][si:]
        ds.series[thisone]["Flag"] = ds.series[thisone]["Flag"][si:]
    ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
elif dt[0]>start_datetime:
    # requested start_datetime is before the start of the file
    print "nc2fn: padding start of file"
    dt_patched = [ldt for ldt in perdelta(start_datetime, dt[0], ts_delta)]
    data_patched = numpy.ones(len(dt_patched))*float(-9999)
    flag_patched = numpy.ones(len(dt_patched))
    # list of series in the data structure
    series_list = ds.series.keys()
    # ds.series["DateTime"]["Data"] is a list not a numpy array so we must treat it differently
    ds.series["DateTime"]["Data"] = dt_patched+ds.series["DateTime"]["Data"]
    ds.series["DateTime"]["Flag"] = numpy.concatenate((flag_patched,ds.series["DateTime"]["Flag"]))
    series_list.remove("DateTime")
    for thisone in series_list:
        ds.series[thisone]["Data"] = numpy.concatenate((data_patched,ds.series[thisone]["Data"]))
        ds.series[thisone]["Flag"] = numpy.concatenate((flag_patched,ds.series[thisone]["Flag"]))
    ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
    # refresh the year, month, day etc arrays now that we have padded the datetime series
    qcutils.get_ymdhms_from_datetime(ds)
# now check the end datetime of the file
end_datetime = dateutil.parser.parse(str(cf["General"]["end_datetime"]))
if dt[-1]>end_datetime:
    # requested end_datetime is before the end of the file
    print "nc2fn: truncating end of file",dt[-1],end_datetime
    ei = qcutils.GetDateIndex(dt,str(end_datetime),ts=ts,match="exact")
    for thisone in ds.series.keys():
        ds.series[thisone]["Data"] = ds.series[thisone]["Data"][:ei+1]
        ds.series[thisone]["Flag"] = ds.series[thisone]["Flag"][:ei+1]
    ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
elif dt[-1]<end_datetime:
    # requested start_datetime is before the start of the file
    print "nc2fn: padding end of file",dt[-1],end_datetime
    dt_patched = [ldt for ldt in perdelta(dt[-1]+ts_delta, end_datetime+ts_delta, ts_delta)]
    data_patched = numpy.ones(len(dt_patched))*float(-9999)
    flag_patched = numpy.ones(len(dt_patched))
    # list of series in the data structure
    series_list = ds.series.keys()
    # ds.series["DateTime"]["Data"] is a list not a numpy array so we must treat it differently
    ds.series["DateTime"]["Data"] = ds.series["DateTime"]["Data"]+dt_patched
    ds.series["DateTime"]["Flag"] = numpy.concatenate((ds.series["DateTime"]["Flag"],flag_patched))
    series_list.remove("DateTime")
    for thisone in series_list:
        ds.series[thisone]["Data"] = numpy.concatenate((ds.series[thisone]["Data"],data_patched))
        ds.series[thisone]["Flag"] = numpy.concatenate((ds.series[thisone]["Flag"],flag_patched))
    ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
    # refresh the year, month, day etc arrays now that we have padded the datetime series
    qcutils.get_ymdhms_from_datetime(ds)
if ts==30:
    nRecs_year = 17520
    nRecs_leapyear = 17568
elif ts==60:
    nRecs_year = 8760
    nRecs_leapyear = 8784
else:
    print "nc2fn: unrecognised time step ("+str(ts)+")"
    sys.exit()
if (int(ds.globalattributes["nc_nrecs"])!=nRecs_year) & (int(ds.globalattributes["nc_nrecs"])!=nRecs_leapyear):
    print "nc2fn: number of records in file does not equal "+str(nRecs_year)+" or "+str(nRecs_leapyear)
    print len(ds.series["DateTime"]["Data"]),ds.series["DateTime"]["Data"][0],ds.series["DateTime"]["Data"][-1]
    sys.exit()
# get the date and time data
Day,flag,attr = qcutils.GetSeries(ds,'Day')
Month,flag,attr = qcutils.GetSeries(ds,'Month')
Year,flag,attr = qcutils.GetSeries(ds,'Year')
Hour,flag,attr = qcutils.GetSeries(ds,'Hour')
Minute,flag,attr = qcutils.GetSeries(ds,'Minute')
# get the data
data = {}
series_list = cf["Variables"].keys()
for series in series_list:
    ncname = cf["Variables"][series]["ncname"]
    if ncname not in ds.series.keys():
        log.error("Series "+ncname+" not in netCDF file, skipping ...")
        series_list.remove(series)
        continue
    data[series] = ds.series[ncname]
    fmt = cf["Variables"][series]["format"]
    if "." in fmt:
        numdec = len(fmt) - (fmt.index(".") + 1)
        strfmt = "{0:."+str(numdec)+"f}"
    else:
        strfmt = "{0:d}"
    data[series]["fmt"] = strfmt
#adjust units if required
for series in series_list:
    if series=="FC" and data[series]["Attr"]["units"]=='mg/m2/s':
        data[series]["Data"] = mf.Fc_umolpm2psfrommgpm2ps(data[series]["Data"])
        data[series]["Attr"]["units"] = "umol/m2/s"
    if series=="CO2" and data[series]["Attr"]["units"]=='mg/m3':
        CO2 = data["CO2"]["Data"]
        TA = data["TA"]["Data"]
        PA = data["PA"]["Data"]
        data[series]["Data"] = mf.co2_ppmfrommgpm3(CO2,TA,PA)
        data[series]["Attr"]["units"] = "umol/mol"
    if series=="H2O" and data[series]["Attr"]["units"]=='g/m3':
        H2O = data["H2O"]["Data"]
        TA = data["TA"]["Data"]
        PA = data["PA"]["Data"]
        data[series]["Data"] = mf.h2o_mmolpmolfromgpm3(H2O,TA,PA)
        data[series]["Attr"]["units"] = "mmol/mol"
    if series=="RH" and data[series]["Attr"]["units"] in ["fraction","frac"]:
        data[series]["Data"] = float(100)*data[series]["Data"]
        data[series]["Attr"]["units"] = "%"
# write the general information to csv file
for item in cf["General"]:
    writer.writerow([item,str(cf['General'][item])])
# write the variable names to the csv file
row_list = ['DateTime','Year','Month','Day','HHMM']
for item in series_list:
    row_list.append(item)
writer.writerow(row_list)
# write the units line to the csv file
units_list = ["-","-","-","-","-"]
for item in series_list:
    units_list.append(data[item]["Attr"]["units"])
writer.writerow(units_list)
# now write the data
for i in range(len(Year)):
    # get the datetime string
    dtstr = '%02d/%02d/%d %02d:%02d'%(Day[i],Month[i],Year[i],Hour[i],Minute[i])
    hrmn = '%02d%02d'%(Hour[i],Minute[i])
    dttup = datetime.datetime(Year[i],Month[i],Day[i],Hour[i],Minute[i]).timetuple()
    doy = float(dttup.tm_yday) + float(dttup.tm_hour)/24 + float(dttup.tm_min)/1440
    data_list = [dtstr,'%d'%(Year[i]),'%02d'%(Month[i]),'%02d'%(Day[i]),hrmn]
    for series in series_list:
        strfmt = data[series]["fmt"]
        if "d" in strfmt:
            data_list.append(strfmt.format(int(round(data[series]["Data"][i]))))
        else:
            data_list.append(strfmt.format(data[series]["Data"][i]))
    writer.writerow(data_list)
# close the csv file
csvfile.close()

print 'nc2fn: All Done'
