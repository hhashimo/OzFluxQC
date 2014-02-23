import datetime
import logging
import numpy
import os
import qcio
import qcutils
import sys

# open the logging file
log = qcutils.startlog('ncconcat','../logfiles/ncconcat.log')

# initialise logicals
TimeGap = False

ds = qcio.DataStructure()

# get the control file
cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0: sys.exit()

InFile_list = cf['Files']['In'].keys()
# read in the first file
ncFileName = cf['Files']['In'][InFile_list[0]]
log.info(' Reading data from '+ncFileName)
ds_n = qcio.nc_read_series(ncFileName)
if len(ds_n.series.keys())==0:
    log.error(' ncconcat: an error occurred reading the netCDF file')
    sys.exit()
# fill the global attributes
for ThisOne in ds_n.globalattributes.keys():
    ds.globalattributes[ThisOne] = ds_n.globalattributes[ThisOne]
# fill the variables
for ThisOne in ds_n.series.keys():
    ds.series[ThisOne] = {}
    ds.series[ThisOne]['Data'] = ds_n.series[ThisOne]['Data']
    ds.series[ThisOne]['Flag'] = ds_n.series[ThisOne]['Flag']
    ds.series[ThisOne]['Attr'] = {}
    for attr in ds_n.series[ThisOne]['Attr'].keys():
        ds.series[ThisOne]['Attr'][attr] = ds_n.series[ThisOne]['Attr'][attr]
ts = int(ds.globalattributes['time_step'])
# loop over the remaining files given in the control file
for n in InFile_list[1:]:
    ncFileName = cf['Files']['In'][InFile_list[int(n)]]
    log.info(' Reading data from '+ncFileName)
    #print 'ncconcat: reading data from '+ncFileName
    ds_n = qcio.nc_read_series(ncFileName)
    if len(ds.series.keys())==0:
        log.error(' ncconcat: an error occurred reading the netCDF file')
        sys.exit()
    dt_n = ds_n.series['DateTime']['Data']
    dt = ds.series['DateTime']['Data']
    nRecs_n = len(ds_n.series['xlDateTime']['Data'])
    nRecs = len(ds.series['xlDateTime']['Data'])
    #print ds.series['DateTime']['Data'][-1],ds_n.series['DateTime']['Data'][-1]
    #print dt[-1],dt[-1]+datetime.timedelta(minutes=ts),dt_n[0]
    if dt_n[0]<dt[-1]+datetime.timedelta(minutes=ts):
        log.info(' Overlapping times detected in consecutive files')
        si = qcutils.GetDateIndex(dt_n,str(dt[-1]))+1
        ei = -1
    if dt_n[0]==dt[-1]+datetime.timedelta(minutes=ts):
        log.info(' Start and end times OK in consecutive files')
        si = 0; ei = -1
    if dt_n[0]>dt[-1]+datetime.timedelta(minutes=ts):
        log.info(' Gap between start and end times in consecutive files')
        si = 0; ei = -1
        TimeGap = True
    # loop over the data series in the concatenated file
    for ThisOne in ds.series.keys():
        # does this series exist in the file being added to the concatenated file
        if ThisOne in ds_n.series.keys():
            # if so, then append this series to the concatenated series
            ds.series[ThisOne]['Data'] = numpy.append(ds.series[ThisOne]['Data'],ds_n.series[ThisOne]['Data'][si:ei])
            ds.series[ThisOne]['Flag'] = numpy.append(ds.series[ThisOne]['Flag'],ds_n.series[ThisOne]['Flag'][si:ei])
        else:
            # if not, then create a dummy series and concatenate that
            ds_n.series[ThisOne] = {}
            ds_n.series[ThisOne]['Data'] = numpy.array([-9999]*nRecs_n,dtype=numpy.float64)
            ds_n.series[ThisOne]['Flag'] = numpy.array([1]*nRecs_n,dtype=numpy.int32)
            ds.series[ThisOne]['Data'] = numpy.append(ds.series[ThisOne]['Data'],ds_n.series[ThisOne]['Data'][si:ei])
            ds.series[ThisOne]['Flag'] = numpy.append(ds.series[ThisOne]['Flag'],ds_n.series[ThisOne]['Flag'][si:ei])
    # and now loop over the series in the file being concatenated
    for ThisOne in ds_n.series.keys():
        # does this series exist in the concatenated data
        if ThisOne not in ds.series.keys():
            # if not then add it
            ds.series[ThisOne] = {}
            ds.series[ThisOne]['Data'] = numpy.array([-9999]*nRecs,dtype=numpy.float64)
            ds.series[ThisOne]['Flag'] = numpy.array([1]*nRecs,dtype=numpy.int32)
            ds.series[ThisOne]['Data'] = numpy.append(ds.series[ThisOne]['Data'],ds_n.series[ThisOne]['Data'][si:ei])
            ds.series[ThisOne]['Flag'] = numpy.append(ds.series[ThisOne]['Flag'],ds_n.series[ThisOne]['Flag'][si:ei])
            ds.series[ThisOne]['Attr'] = {}
            for attr in ds_n.series[ThisOne]['Attr'].keys():
                ds.series[ThisOne]['Attr'][attr] = ds_n.series[ThisOne]['Attr'][attr]
    # now sort out any time gaps
    if TimeGap:
        qcutils.FixTimeGaps(ds)
        TimeGap = False

ds.globalattributes['nc_nrecs'] = str(len(ds.series['xlDateTime']['Data']))

# write the netCDF file
ncFileName = qcio.get_keyvalue_from_cf(cf['Files']['Out'],'ncFileName')
log.info(' Writing data to '+ncFileName)
ncFile = qcio.nc_open_write(ncFileName)
qcio.nc_write_series(ncFile,ds)
