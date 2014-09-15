import sys
sys.path.append('../scripts')
import constants as c
import csv
import matplotlib.pyplot as plt
import meteorologicalfunctions as mf
import netCDF4
import numpy
import os
import qcio
import qcts
import qcutils

# open the logging file
log = qcutils.startlog('bios2nc','../logfiles/bios2nc.log')

# load a control file
cf = qcio.load_controlfile(path='../controlfiles')
# get the variable list
var_list = cf["Variables"].keys()
# get the site list
site_list = cf["Sites"].keys()
# get the variable names
xldt_label = cf["Variables"]["xlDateTime"]["dingo_name"]
evii_label = cf["Variables"]["EVI_MODIS_interp"]["dingo_name"]
evis_label = cf["Variables"]["EVI_MODIS_smooth"]["dingo_name"]
for site in site_list:
    # get the input file mask
    infilename = cf["Sites"][site]["in_filepath"]+cf["Sites"][site]["in_filename"]
    if not os.path.isfile(infilename):
        log.error("CSV file "+infilename+" not found, skipping ...")
        continue
    log.info("Starting site: "+site)
    csvfile = open(infilename,'rb')
    csvreader = csv.reader(csvfile)
    header = csvreader.next()
    dt_col=header.index(xldt_label)
    eviint_col = header.index(evii_label)
    evismo_col = header.index(evis_label)
    csvfile.close()

    evi=numpy.genfromtxt(infilename,skip_header=1,delimiter=",",usecols=(dt_col,eviint_col,evismo_col))
    
    ds=qcio.DataStructure()
    ds.globalattributes["xl_datemode"] = 0
    ds.globalattributes["time_step"] = cf["Sites"][site]["time_step"]
    ds.globalattributes["time_zone"] = cf["Sites"][site]["time_zone"]
    ds.globalattributes["latitude"] = cf["Sites"][site]["latitude"]
    ds.globalattributes["longitude"] = cf["Sites"][site]["longitude"]
    
    ds.series["xlDateTime"] = {}
    ds.series["xlDateTime"]["Data"] = numpy.array(evi[:,0]+1,dtype=numpy.float64)
    nRecs = len(ds.series["xlDateTime"]["Data"])
    ds.globalattributes["nc_nrecs"] = nRecs
    flag = numpy.zeros(nRecs,dtype=numpy.int32)
    ds.series["xlDateTime"]["Flag"] = flag
    ds.series["xlDateTime"]["Attr"] = {}
    ds.series["xlDateTime"]["Attr"]["long_name"] = "Date/time in Excel format"
    ds.series["xlDateTime"]["Attr"]["units"] = "days since 1899-12-31 00:00:00"
    ds.series["xlDateTime"]["Attr"]["standard_name"] = "not defined"
    ds.series["xlDateTime"]["Attr"]["cf_role"] = "timeseries_id"
    
    qcutils.get_datetimefromxldate(ds)
    qcutils.get_ymdhms_from_datetime(ds)
    
    attr=qcutils.MakeAttributeDictionary(long_name="MODIS EVI, 250m, interpolated",units="none")
    qcutils.CreateSeries(ds,"EVI_MODIS_interp",evi[:,1],Flag=flag,Attr=attr)
    attr=qcutils.MakeAttributeDictionary(long_name="MODIS EVI, 250m, interpolated and smooothed",units="none")
    qcutils.CreateSeries(ds,"EVI_MODIS_smooth",evi[:,2],Flag=flag,Attr=attr)
    
    outfilename = cf["Sites"][site]["out_filepath"]+cf["Sites"][site]["out_filename"]
    ncfile=qcio.nc_open_write(outfilename)
    qcio.nc_write_series(ncfile,ds)