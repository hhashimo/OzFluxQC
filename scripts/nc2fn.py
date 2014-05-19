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
Day,flag = qcutils.GetSeries(ds,'Day')
Month,flag = qcutils.GetSeries(ds,'Month')
Year,flag = qcutils.GetSeries(ds,'Year')
Hour,flag = qcutils.GetSeries(ds,'Hour')
Minute,flag = qcutils.GetSeries(ds,'Minute')
# get the data
ust = ds.series[cf['Variables']['ust']['ncname']]
FC = ds.series[cf['Variables']['FC']['ncname']]
CO2 = ds.series[cf['Variables']['CO2_top']['ncname']]
TA = ds.series[cf['Variables']['TA']['ncname']]
RGin = ds.series[cf['Variables']['RG_in']['ncname']]
H2O = ds.series[cf['Variables']['H2O']['ncname']]
LE = ds.series[cf['Variables']['LE']['ncname']]
H = ds.series[cf['Variables']['H']['ncname']]
G1 = ds.series[cf['Variables']['G1']['ncname']]
PRECIP = ds.series[cf['Variables']['PRECIP']['ncname']]
SWC1 = ds.series[cf['Variables']['SWC1']['ncname']]
TS1 = ds.series[cf['Variables']['TS1']['ncname']]
RNET = ds.series[cf['Variables']['RNET']['ncname']]
SWin = ds.series[cf['Variables']['SWin']['ncname']]
SWout = ds.series[cf['Variables']['SWout']['ncname']]
LWin = ds.series[cf['Variables']['LWin']['ncname']]
LWout = ds.series[cf['Variables']['LWout']['ncname']]
PA = ds.series[cf['Variables']['PA']['ncname']]
WD = ds.series[cf['Variables']['WD']['ncname']]
WS = ds.series[cf['Variables']['WS']['ncname']]
if cf['Variables']['RH']['ncname'] not in ds.series.keys():
    Ah,f = qcutils.GetSeriesasMA(ds,'Ah')
    Ta,f = qcutils.GetSeriesasMA(ds,'Ta')
    RH = mf.RHfromabsolutehumidity(Ah, Ta)
    attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='not defined')
    qcutils.CreateSeries(ds,cf['Variables']['RH']['ncname'],RH,FList=['Ta','Ah'],Attr=attr)
RH = ds.series[cf['Variables']['RH']['ncname']]
# adjust units if required
if FC['Attr']['units']=='mg/m2/s':
    FC['Data'] = mf.Fc_umolpm2psfrommgpm2ps(FC['Data'])
    FC['Attr']['units'] = 'umol/m2/s'
if CO2['Attr']['units']=='mg/m3':
    CO2['Data'] = mf.co2_ppmfrommgpm3(CO2['Data'],TA['Data'],PA['Data'])
    CO2['Attr']['units'] = 'umol/mol'
if H2O['Attr']['units']=='g/m3':
    H2O['Data'] = mf.h2o_mmolpmolfromgpm3(H2O['Data'],TA['Data'],PA['Data'])
    H2O['Attr']['units'] = 'mmol/mol'
if RH['Attr']['units'] in ["fraction","frac"]:
    RH['Data'] = float(100)*RH['Data']
    RH['Attr']['units'] = '%'
# write the general information to csv file
for item in cf["General"]:
    writer.writerow([item,str(cf['General'][item])])
# write the variable names to the csv file
writer.writerow(['DateTime','Year','Month','Day','HHMM',
                 'FC','CO2','UST','RGin','TA','H2O',
                 'LE','H','G1','PRECIP','SWC1','Ts1',
                 'RNET','SWin','SWout','LWin','LWout',
                 'PA','WD','WS','RH'])
# write the units line to the csv file
writer.writerow(['-','-','-','-','-',
                 FC['Attr']['units'],
                 CO2['Attr']['units'],
                 ust['Attr']['units'],
                 RGin['Attr']['units'],
                 TA['Attr']['units'],
                 H2O['Attr']['units'],
                 LE['Attr']['units'],
                 H['Attr']['units'],
                 G1['Attr']['units'],
                 PRECIP['Attr']['units'],
                 SWC1['Attr']['units'],
                 TS1['Attr']['units'],
                 RNET['Attr']['units'],
                 SWin['Attr']['units'],
                 SWout['Attr']['units'],
                 LWin['Attr']['units'],
                 LWout['Attr']['units'],
                 PA['Attr']['units'],
                 WD['Attr']['units'],
                 WS['Attr']['units'],
                 RH['Attr']['units']])
for i in range(len(Year)):
    # get the datetime string
    dtstr = '%02d/%02d/%d %02d:%02d'%(Day[i],Month[i],Year[i],Hour[i],Minute[i])
    hrmn = '%02d%02d'%(Hour[i],Minute[i])
    dttup = datetime.datetime(Year[i],Month[i],Day[i],Hour[i],Minute[i]).timetuple()
    doy = float(dttup.tm_yday) + float(dttup.tm_hour)/24 + float(dttup.tm_min)/1440
    writer.writerow([dtstr,'%d'%(Year[i]),'%02d'%(Month[i]),'%02d'%(Day[i]),hrmn,
                     '%.2f'%(FC['Data'][i]),
                     '%.1f'%(CO2['Data'][i]),
                     '%.2f'%(ust['Data'][i]),
                     '%d'%(RGin['Data'][i]),
                     '%.2f'%(TA['Data'][i]),
                     '%.2f'%(H2O['Data'][i]),
                     '%d'%(LE['Data'][i]),
                     '%d'%(H['Data'][i]),
                     '%d'%(G1['Data'][i]),
                     '%.2f'%(PRECIP['Data'][i]),
                     '%.2f'%(SWC1['Data'][i]),
                     '%.2f'%(TS1['Data'][i]),
                     '%d'%(RNET['Data'][i]),
                     '%d'%(SWin['Data'][i]),
                     '%d'%(SWout['Data'][i]),
                     '%d'%(LWin['Data'][i]),
                     '%d'%(LWout['Data'][i]),
                     '%.2f'%(PA['Data'][i]),
                     '%d'%(WD['Data'][i]),
                     '%.2f'%(WS['Data'][i]),
                     '%d'%(RH['Data'][i])])
# close the csv file
csvfile.close()

print 'nc2fn: All Done'
