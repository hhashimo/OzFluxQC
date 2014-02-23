import constants as c
import csv
import datetime
import logging
import meteorologicalfunctions as mf
import numpy
import os
import qcio
import qcts
import qcutils
import sys
import time
import xlrd

# start the log file
log = qcutils.startlog('fn2nc','../logfiles/fn2nc.log')
# get the control dile contents
cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0:
    print 'fn2nc: no control file specified, exiting ...'
    sys.exit()
# set platform specific features
if 'General' in cf:
    if 'Platform' in cf['General'].keys():
        if cf['General']['Platform']=='PC':
            datemode = 0
        else:
            datemode = 1
    else:
        print 'fn2nc: Platform not given in General section of control file'
        sys.exit()
else:
    print 'fn2nc: No General section in control file'
    sys.exit()
# get the processing level
level = cf['Files'].keys()[0]
# get the file names
csvFileName = cf['Files'][level]['csvFilePath'] + cf['Files'][level]['csvFileName']
# get a list of the series to write to the netCDF file
outSeries = cf['Variables'].keys()
# get a list of the series names to read from the csv file
inSeries = []
for ThisOne in outSeries: inSeries.append(cf['Variables'][ThisOne]['csv']['Name'])
# open the csv file
csvfile = open(csvFileName,'rb')
# read lines to the header line
for i in range(int(cf['Files'][level]['HeaderRow'])-1): csvfile.readline()
# read the header line
headerLine = csvfile.readline().strip().split(',')
headerList = [l.strip() for l in headerLine]
#print headerList
# go back to the start of the csv file
csvfile.seek(0)
# read to the first data row
for i in range(int(cf['Files'][level]['FirstDataRow'])-1): csvfile.readline()
# get the csv reader object, read file as dictionary
csvreader = csv.DictReader(csvfile, headerList, delimiter=',')
# transfer the contents of the csv file into the dictionary
dd = {}
for item in csvreader:
    for name in headerList:
        dd.setdefault(name,[]).append(item[name])
# close the csv file
csvfile.close()
# get the dictionary contents as numpy arrays
# instance the data structure
ds = qcio.DataStructure()
# set the attributes
ds.globalattributes['csvFullName'] = csvFileName
s = os.stat(csvFileName)
t = time.localtime(s.st_mtime)
ds.globalattributes['csvModDateTime'] = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
# global attribute Functions needed by CalculateMeteorologicalVariables
ds.globalattributes['Functions'] = ''
if 'Global' in cf.keys():
    for gattr in cf['Global'].keys():
        ds.globalattributes[gattr] = cf['Global'][gattr]
if 'nc_level' not in ds.globalattributes:
    ds.globalattributes['nc_level'] = str(level)
# set up the data series in ds
# NOTE: QC flags are done separately after all data series have been defined
# do xlDateTime first, exit routine if not found in control file
if 'xlDateTime' not in outSeries:
    log.error('fn2nc: xlDateTime not found in control file, exiting ...')
    sys.exit()
else:
    ds.series['xlDateTime'] = {}
    s = dd[cf['Variables']['xlDateTime']['csv']['Name']]
    xldtlist = []
    for item in s:
        dt_obj = datetime.datetime.strptime(item,'%d/%m/%Y %H:%M')
        xldt = xlrd.xldate.xldate_from_datetime_tuple((dt_obj.year,dt_obj.month,dt_obj.day,
                                                      dt_obj.hour,dt_obj.minute,dt_obj.second),datemode)
        xldtlist.append(xldt)
    ds.series['xlDateTime']['Data'] = numpy.array(xldtlist,dtype=numpy.float64)
    qcutils.get_ymdhmsfromxldate(ds)
    qcutils.get_datetimefromymdhms(ds)
    outSeries.remove('xlDateTime')
nRecs = len(ds.series['xlDateTime']['Data'])
ds.globalattributes['nc_nrecs'] = str(nRecs)
# now do the rest of the data series
for ThisOne in outSeries:
    print ThisOne
    ds.series[ThisOne] = {}
    ds.series[ThisOne]['Data'] = numpy.array(dd[cf['Variables'][ThisOne]['csv']['Name']],dtype=numpy.float64)
for ThisOne in outSeries:
    if 'Attr' in cf['Variables'][ThisOne].keys():
        ds.series[ThisOne]['Attr'] = {}
        for attr in cf['Variables'][ThisOne]['Attr'].keys():
            ds.series[ThisOne]['Attr'][attr] = cf['Variables'][ThisOne]['Attr'][attr]
# get the QC flags
qcts.get_qcflag(ds)
# get meteorological variables
if 'Ah' not in ds.series.keys():
    if 'RH' in ds.series.keys() and 'Ta' in ds.series.keys():
        Ta,Ta_flag = qcutils.GetSeriesasMA(ds,'Ta')
        RH,RH_flag = qcutils.GetSeriesasMA(ds,'RH')
        Ah = mf.absolutehumidityfromRH(Ta,RH)
        attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity',units='g/m3',standard_name='not defined')
        qcutils.CreateSeries(ds,'Ah',Ah,Flag=RH_flag,Attr=attr)
    else:
        log.error('fn2nc: Neither Ah nor RH found, exiting ...')
        sys.exit()
qcts.CalculateMeteorologicalVariables(ds)
# get available energy if it is not in the file
if 'Fa' not in ds.series.keys():
    if 'Fn' in ds.series.keys() and 'Fg' in ds.series.keys():
        qcts.CalculateAvailableEnergy(ds,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg')
# write the netCDF file
outfilename = qcio.get_outfilename_from_cf(cf)
ncFile = qcio.nc_open_write(outfilename)
qcio.nc_write_series(ncFile,ds)
