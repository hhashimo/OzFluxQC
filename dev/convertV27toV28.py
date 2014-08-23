# the following line needed for unicode character in convert_anglestring
# -*- coding: latin-1 -*-
import datetime
import netCDF4
import numpy
import os
import sys
import time
import Tkinter,tkFileDialog,tkSimpleDialog

class DataStructure(object):
    def __init__(self):
        self.series = {}
        self.globalattributes = {}
        self.dimensions = {}
        self.mergeserieslist = []
        self.averageserieslist = []
        self.soloserieslist = []
        self.climatologyserieslist = []

def convert_anglestring(anglestring):
    """
    Purpose:
     Attempt to convert an angle string to a float.
    Usage:
     a = qcutils.convert_anglestring(astr)
     Acceptable input formats:
      astr = '''34 12' 24" S'''
      astr = '''34 12 24S'''
      astr = '''34 12'24.123"S''
      astr = '''34.123 S'''
      astr = '''-34.123'''
    """
    quadlist=["N","E","S","W"]
    direction = {'N':1, 'S':-1, 'E': 1, 'W':-1}
    try:
        # simple casting may work, who knows?
        return float(anglestring)
    except ValueError:
        # replace the degrees, minutes and seconds symbols with spaces
        new = anglestring.replace(u'°',' ').replace('\'',' ').replace('"',' ')
        # check there is a space between the quadrant letter (assumed to be one of N, E, W or S)
        # and the next character to the left
        # find out which of N, E, S, or W is in the string
        for item in quadlist:
            if item in new: quadletter=item
        # now get the index of this character in the string
        i=new.index(quadletter)
        # check that the next character to the left is a space character
        if new[i-1] != " ": new = new[0:i]+" "+new[i:]
        # now split the string on space characters
        new = new.split()
        # get the quadrant letter
        new_dir = new.pop()
        # make sure we have 3 parts
        new.extend([0,0,0])
        # return with the string converted to a float
        return (float(new[0])+float(new[1])/60.0+float(new[2])/3600.0) * direction[new_dir]    

def get_datetimefromymdhms(ds):
    ''' Creates a series of Python datetime objects from the year, month,
    day, hour, minute and second series stored in the netCDF file.'''
    SeriesList = ds.series.keys()
    if 'Year' not in SeriesList or 'Month' not in SeriesList or 'Day' not in SeriesList or 'Hour' not in SeriesList or 'Minute' not in SeriesList or 'Second' not in SeriesList:
        log.info(' get_datetimefromymdhms: unable to find all datetime fields required')
        return
    print ' Getting the date and time series'
    nRecs = int(ds.globalattributes["nc_nrecs"])
    ds.series[unicode('DateTime')] = {}
    ds.series['DateTime']['Data'] = [None]*nRecs
    for i in range(nRecs):
        ds.series['DateTime']['Data'][i] = datetime.datetime(int(ds.series['Year']['Data'][i]),
                                                       int(ds.series['Month']['Data'][i]),
                                                       int(ds.series['Day']['Data'][i]),
                                                       int(ds.series['Hour']['Data'][i]),
                                                       int(ds.series['Minute']['Data'][i]),
                                                       int(ds.series['Second']['Data'][i]))
    ds.series['DateTime']['Flag'] = numpy.zeros(nRecs)
    ds.series['DateTime']['Attr'] = {}
    ds.series['DateTime']['Attr']['long_name'] = 'Date-time object'
    ds.series['DateTime']['Attr']['units'] = 'None'

def get_filename_dialog(path='.',title='Choose a file'):
    '''
    Put up a file open dialog.
    USEAGE:
     fname = qcio.get_filename_dialog(path=<path_to_file>,title=<tile>)
    INPUT:
     path  - the path to the file location, optional
     title - the title for the file open dialog, optional
    RETURNS:
     fname - the full file name including path, string
    '''
    root = Tkinter.Tk(); root.withdraw()
    FileName = tkFileDialog.askopenfilename(parent=root,initialdir=path,title=title)
    root.destroy()
    return str(FileName)

def get_ncdtype(Series):
    sd = Series.dtype.name
    dt = 'f'
    if sd=='float64': dt = 'd'
    if sd=='int32': dt = 'i'
    if sd=='int64': dt = 'l'
    return dt

def nc_read_series(ncFullName):
    ''' Read a netCDF file and put the data and meta-data into a DataStructure'''
    print ' Reading netCDF file '+ncFullName
    netCDF4.default_encoding = 'latin-1'
    ds = DataStructure()
    # check to see if the requested file exists, return empty ds if it doesn't
    if not os.path.exists(ncFullName):
        raise Exception("nc_read_series: file "+ncFullName+" not found")
    # file probably exists, so let's read it
    ncFile = netCDF4.Dataset(ncFullName,'r')
    # now deal with the global attributes
    gattrlist = ncFile.ncattrs()
    if len(gattrlist)!=0:
        for gattr in gattrlist:
            ds.globalattributes[gattr] = getattr(ncFile,gattr)
    # get a list of the variables in the netCDF file (not their QC flags)
    varlist = [x for x in ncFile.variables.keys() if "_QCFlag" not in x]
    for ThisOne in varlist:
        # skip variables that do not have time as a dimension
        dimlist = [x.lower() for x in ncFile.variables[ThisOne].dimensions]
        if "time" not in dimlist: continue
        # create the series in the data structure
        ds.series[unicode(ThisOne)] = {}
        # get the data and the QC flag
        data,flag,attr = nc_read_var(ncFile,ThisOne)
        ds.series[ThisOne]["Data"] = data
        ds.series[ThisOne]["Flag"] = flag
        ds.series[ThisOne]["Attr"] = attr
    ncFile.close()
    # get a series of Python datetime objects
    get_datetimefromymdhms(ds)
    return ds

def nc_read_var(ncFile,ThisOne):
    """ Reads a variable from a netCDF file and returns the data, the QC flag and the variable
        attribute dictionary.
    """
    # check the number of dimensions
    nDims = len(ncFile.variables[ThisOne].shape)
    if nDims not in [1,3]:
        msg = "nc_read_var: unrecognised number of dimensions ("+str(nDims)
        msg = msg+") for netCDF variable "+ ThisOne
        raise Exception(msg)
    if nDims==1:
        # single dimension
        data = ncFile.variables[ThisOne][:]
        # netCDF4 returns a masked array if the "missing_variable" attribute has been set
        # for the variable, here we trap this and force the array in ds.series to be ndarray
        if numpy.ma.isMA(data): data = numpy.ma.filled(data,float(-9999))
        # check for a QC flag
        if ThisOne+'_QCFlag' in ncFile.variables.keys():
            # load it from the netCDF file
            flag = ncFile.variables[ThisOne+'_QCFlag'][:]
        else:
            # create an empty flag series if it does not exist
            nRecs = numpy.size(data)
            flag = numpy.zeros(nRecs,dtype=numpy.int32)
    elif nDims==3:
        # 3 dimensions
        data = ncFile.variables[ThisOne][:,0,0]
        # netCDF4 returns a masked array if the "missing_variable" attribute has been set
        # for the variable, here we trap this and force the array in ds.series to be ndarray
        if numpy.ma.isMA(data): data = numpy.ma.filled(data,float(-9999))
        # check for a QC flag
        if ThisOne+'_QCFlag' in ncFile.variables.keys():
            # load it from the netCDF file
            flag = ncFile.variables[ThisOne+'_QCFlag'][:,0,0]
        else:
            # create an empty flag series if it does not exist
            nRecs = numpy.size(data)
            flag = numpy.zeros(nRecs,dtype=numpy.int32)
    # get the variable attributes
    vattrlist = ncFile.variables[ThisOne].ncattrs()
    attr = {}
    if len(vattrlist)!=0:
        for vattr in vattrlist:
            attr[vattr] = getattr(ncFile.variables[ThisOne],vattr)
    return data,flag,attr

def nc_open_write(ncFullName,nctype='NETCDF4'):
    print ' Opening netCDF file '+ncFullName+' for writing'
    try:
        ncFile = netCDF4.Dataset(ncFullName,'w',format=nctype)
    except:
        print ' Unable to open netCDF file '+ncFullName+' for writing'
        ncFile = ''
    return ncFile

def nc_write_series(ncFile,ds,outputlist=None,ndims=3):
    ldt = ds.series["DateTime"]["Data"]
    ds.globalattributes['QC_version'] = "OzFluxQC V2.8.1"
    for ThisOne in ds.globalattributes.keys():
        setattr(ncFile,ThisOne,ds.globalattributes[ThisOne])
    t = time.localtime()
    rundatetime = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    setattr(ncFile,'nc_rundatetime',rundatetime)
    # we specify the size of the Time dimension because netCDF4 is slow to write files
    # when the Time dimension is unlimited
    nRecs = int(ds.globalattributes['nc_nrecs'])
    ncFile.createDimension("time",nRecs)
    if ndims==3:
        ncFile.createDimension("latitude",1)
        ncFile.createDimension("longitude",1)
        dims = ("time","latitude","longitude")
    else:
        dims = ("time",)
    if outputlist==None:
        outputlist = ds.series.keys()
    else:
        for ThisOne in outputlist:
            if ThisOne not in ds.series.keys():
                log.info(' Requested series '+ThisOne+' not found in data structure')
                outputlist.remove(ThisOne)
        if len(outputlist)==0: outputlist = ds.series.keys()
    # can't write an array of Python datetime objects to a netCDF file
    # actually, this could be written as characters
    for ThisOne in ["DateTime","DateTime_UTC"]:
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    # write the time variable
    if "time" not in outputlist:
        nc_time = netCDF4.date2num(ldt,"days since 1800-01-01 00:00:00.0",calendar="gregorian")
        ncVar = ncFile.createVariable("time","d",("time",))
        ncVar[:] = nc_time
        setattr(ncVar,"long_name","time")
        setattr(ncVar,"standard_name","time")
        setattr(ncVar,"units","days since 1800-01-01 00:00:00.0")
        setattr(ncVar,"calendar","gregorian")
    # now write the latitude and longitude variables
    if ndims==3:
        if "latitude" not in outputlist:
            ncVar = ncFile.createVariable("latitude","d",("latitude",))
            ncVar[:] = convert_anglestring(str(ds.globalattributes["latitude"]))
            setattr(ncVar,'long_name','latitude')
            setattr(ncVar,'standard_name','latitude')
            setattr(ncVar,'units','degrees north')
        if "longitude" not in outputlist:
            ncVar = ncFile.createVariable("longitude","d",("longitude",))
            ncVar[:] = convert_anglestring(str(ds.globalattributes["longitude"]))
            setattr(ncVar,'long_name','longitude')
            setattr(ncVar,'standard_name','longitude')
            setattr(ncVar,'units','degrees east')
    # now make sure the date and time series are in outputlist
    datetimelist = ['xlDateTime','xlDateTime_UTC','Year','Month','Day','Hour','Minute','Second','Hdh']
    # and write them to the netCDF file
    for ThisOne in sorted(datetimelist):
        if ThisOne in ds.series.keys(): nc_write_var(ncFile,ds,ThisOne,dims)
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    # write everything else to the netCDF file
    for ThisOne in sorted(outputlist):
        nc_write_var(ncFile,ds,ThisOne,dims)
    # write the coordinate reference system (crs) variable
    if "crs" not in outputlist:
        ncVar = ncFile.createVariable("crs","i",())
        setattr(ncVar,"grid_mapping_name","latitude_longitude")
        setattr(ncVar,"long_name","WGS 1984 datum")
        setattr(ncVar,"longitude_of_prime_meridian","0.0")
        setattr(ncVar,"semi_major_axis","6378137.0")
        setattr(ncVar,"inverse_flattening","298.257223563")
    ncFile.close()

def nc_write_var(ncFile,ds,ThisOne,dim):
    """
    PURPOSE:
     Function to write data from a series in the data structure to a netCDF variable.
    USAGE:
     nc_write_var(ncFile,ds,ThisOne,("time","latitude","longitude"))
      where ncFile is a netCDF file object
            ds is the data structure
            ThisOne is the label of a series in ds
            ("time","latitude","longitude") is the dimension tuple
    AUTHOR: PRI
    DATE: August 2014
    """
    # get the data type of the series in ds
    dt = get_ncdtype(ds.series[ThisOne]['Data'])
    # create the netCDF variable
    ncVar = ncFile.createVariable(ThisOne,dt,dim)
    # different writes to the variable depending on whether it is 1D or 3D
    if len(dim)==1: ncVar[:] = ds.series[ThisOne]['Data'].tolist()
    if len(dim)==3: ncVar[:,0,0] = ds.series[ThisOne]['Data'].tolist()
    # write the attributes
    for attr in ds.series[ThisOne]['Attr']:
        setattr(ncVar,attr,ds.series[ThisOne]['Attr'][attr])
    # get the data type of the QC flag
    dt = get_ncdtype(ds.series[ThisOne]['Flag'])
    # create the variable
    ncVar = ncFile.createVariable(ThisOne+'_QCFlag',dt,dim)
    # write 1D or 3D
    if len(dim)==1: ncVar[:] = ds.series[ThisOne]['Flag'].tolist()
    if len(dim)==3: ncVar[:,0,0] = ds.series[ThisOne]['Flag'].tolist()
    # set the attributes
    setattr(ncVar,'long_name',ThisOne+'QC flag')
    setattr(ncVar,'units','none')

def get_timezone(site_name):
    """ Return the time zone based on the site name."""
    tz_dict = {"adelaideriver":"Australia/Darwin",
               "alicespringsmulga":"Australia/Darwin",
               "arcturus":"Australia/Brisbane",
               "calperum":"Australia/Adelaide",
               "capetribulation":"Australia/Brisbane",
               "cumberlandplains":"Australia/Sydney",
               "cup_ec":"Australia/Sydney",
               "daintree":"Australia/Brisbane",
               "dalypasture":"Australia/Darwin",
               "dalyregrowth":"Australia/Darwin",
               "dalyuncleared":"Australia/Darwin",
               "dargo":"Australia/Melbourne",
               "dryriver":"Australia/Darwin",
               "foggdam":"Australia/Darwin",
               "gingin":"Australia/Perth",
               "greatwestern":"Australia/Perth",
               "howardsprings":"Australia/Darwin",
               "litchfield":"Australia/Darwin",
               "nimmo":"Australia/Sydney",
               "reddirt":"Australia/Darwin",
               "riggs":"Australia/Melbourne",
               "robson":"Australia/Brisbane",
               "samford":"Australia/Brisbane",
               "sturtplains":"Australia/Darwin",
               "titreeeast":"Australia/Darwin",
               "tumbarumba":"Australia/Canberra",
               "wallaby":"Australia/Melbourne",
               "warra":"Australia/Hobart",
               "whroo":"Australia/Melbourne",
               "wombat":"Australia/Melbourne",
               "yanco_jaxa":"Australia/Sydney"}
    # strip out spaces and commas from the site name
    site_name = site_name.replace(" ","").replace(",","")
    found_tz = False
    for item in tz_dict.keys():
        if item in site_name.lower():
            time_zone = tz_dict[item]
            found_tz = True
            break
    if not found_tz:
        # cant find the site in the dictionary so ask the user
        print item
        print site_name.lower()
        root = Tkinter.Tk(); root.withdraw()
        time_zone = tkSimpleDialog.askstring("Time zone","Enter time zone eg Australia/Melbourne")
        root.destroy()
    return time_zone

""" Convert V2.7 (1D) netCDF files to V2.8 (3D). """
ncV27name_list = []
# use the following lines of code to have the script prompt the user for a file name
# put up an open file dialog to get the file name
# !!! start open file dialog block !!!
#ncV27name = get_filename_dialog(path="../../Sites",title="Choose a V27 netCDF file ...")
#if len(ncV27name)==0: sys.exit()
#ncV27name_list.append(ncV27name)
# !!! end open file dialog block !!!
# use the following line of code to read a file name from the command line
# get the file name from the command line
# !!! start command line entry block !!!
#ncV27name_list = []
#ncV27name = sys.argv[1]
#ncV27name_list.append(ncV27name)
# !!! end command line entry block !!!
# use the following lines of code to read a control file containing file names
# !!! start control file block !!!
filename_listfile = get_filename_dialog(path="../controlfiles",title="Choose a control file ...")
if len(filename_listfile)==0: sys.exit()
f = open(filename_listfile)
while 1:
    line = f.readline()
    if not line: break
    ncV27name_list.append(line.strip("\n"))
# !!! end control file block !!!

for ncV27name in ncV27name_list:
    if not os.path.exists(ncV27name):
        print "File "+ncV27name+" not found ..."
        continue
    print "Processing file: "+ncV27name
    # read the V2.7 file
    ds = nc_read_series(ncV27name)
    # add the "time_zone" global attribute if it is not present
    if "time_zone" not in ds.globalattributes.keys():
        for gattr in ["site_name","SiteName"]:
            if gattr in ds.globalattributes.keys():
                time_zone = get_timezone(ds.globalattributes[gattr])
        ds.globalattributes["time_zone"] = time_zone
    # add the "missing_value" variable attribute if it is not present
    for ThisOne in ds.series.keys():
        if "missing_value" not in ds.series[ThisOne]["Attr"].keys():
            ds.series[ThisOne]["Attr"]["missing_value"] = str(-9999)
    # write the V2.8 file
    ncV28name = ncV27name.replace(".nc","_V28.nc")
    ncFile = nc_open_write(ncV28name,nctype='NETCDF4')
    nc_write_series(ncFile, ds)
