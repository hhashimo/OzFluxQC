import sys
sys.path.append('../scripts')
import csv
import datetime
import constants as c
import glob
import logging
import meteorologicalfunctions as mf
import netCDF4
import numpy
import qcio
import qcts
import qcutils
import xlrd

logging.basicConfig()
# open the logging file
log = qcutils.startlog('aws2nc','../logfiles/aws2nc.log')

# get the site information and the AWS stations to use
xlname = "../../BoM/Locations/AWS_Locations.xls"
wb = xlrd.open_workbook(xlname)
sheet = wb.sheet_by_name("OzFlux")
xl_row = 10
xl_col = 0
bom_sites_info = {}
for n in range(xl_row,sheet.nrows):
    xlrow = sheet.row_values(n)
    bom_sites_info[str(xlrow[0])] = {}
    bom_sites_info[xlrow[0]]["latitude"] = xlrow[1]
    bom_sites_info[xlrow[0]]["longitude"] = xlrow[2]
    bom_sites_info[xlrow[0]]["elevation"] = xlrow[3]
    for i in [4,10,16,22]:
        if xlrow[i]!="":
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))] = {}
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))]["site_name"] = xlrow[i]
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))]["latitude"] = xlrow[i+2]
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))]["longitude"] = xlrow[i+3]
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))]["elevation"] = xlrow[i+4]
            bom_sites_info[str(xlrow[0])][str(int(xlrow[i+1]))]["distance"] = xlrow[i+5]

in_path = "../../BoM/AWS/Current/"
in_filename = in_path+"HM01X_Data*.csv"
file_list = sorted(glob.glob(in_filename))

site_list = bom_sites_info.keys()
for site_name in sorted(site_list):
    log.info("Starting site: "+site_name)
    sname = site_name.replace(" ","")
    ncname = in_path+sname+"_AWS.nc"
    site_latitude = bom_sites_info[site_name]["latitude"]
    site_longitude = bom_sites_info[site_name]["longitude"]
    site_elevation = bom_sites_info[site_name]["elevation"]
    site_number_list = bom_sites_info[site_name].keys()
    for item in ["latitude","longitude","elevation"]:
        if item in site_number_list: site_number_list.remove(item)
    # read the CSV files and put the contents into data_dict
    data_dict = {}
    for idx,sn in enumerate(site_number_list):
        # get a list of file names that contain the relevent station numbers
        csvname = [fn for fn in file_list if str(sn) in fn]
        # continue to next station if this station not in file_list
        if len(csvname)==0: continue
        log.info("Reading CSV file: "+str(csvname[0]))
        # columns are:
        # file data content
        #  1    0    station number
        #  7    1    year, local standard time
        #  8    2    month, local standard time
        #  9    3    day, local standard time
        #  10   4    hour, local standard time
        #  11   5    minute, local standard time
        #  12   6    precip since 0900
        #  14   7    air temperature, C
        #  16   8    dew point temperature, C
        #  18   9    relative humidity, %
        #  20   10   wind speed, m/s
        #  22   11   wind direction, degT
        #  24   12   gust in last 10 minutes, m/s
        #  26   13   station pressure, hPa
        data=numpy.genfromtxt(csvname[0],skip_header=1,delimiter=",",usecols=(1,7,8,9,10,11,12,14,16,18,20,22,24,26),
                              missing_values=-9999,filling_values=-9999)
        data = numpy.ma.masked_equal(data,float(-9999),copy=True)
        data_dict[sn] = data
    # now pull the data out and put it in separate data structures, one per station, all
    # of which are held in a data structure dictionary
    ds_dict = {}
    for bom_id in data_dict.keys():
        log.info("Processing BoM station: "+str(bom_id))
        ds=qcio.DataStructure()
        # put the year, month, day, hour and minute into the data structure
        nRecs = data_dict[bom_id].shape[0]
        ds.globalattributes["nc_nrecs"] = nRecs
        ds.globalattributes["time_step"] = 30
        ds.globalattributes["latitude"] = bom_sites_info[site_name][str(bom_id)]["latitude"]
        ds.globalattributes["longitude"] = bom_sites_info[site_name][str(bom_id)]["longitude"]
        flag = numpy.zeros(nRecs,dtype=numpy.int32)
        Seconds = numpy.zeros(nRecs,dtype=numpy.float64)
        qcutils.CreateSeries(ds,'Year',data_dict[bom_id][:,1],Flag=flag,Attr=qcutils.MakeAttributeDictionary(long_name='Year',units='none'))
        qcutils.CreateSeries(ds,'Month',data_dict[bom_id][:,2],Flag=flag,Attr=qcutils.MakeAttributeDictionary(long_name='Month',units='none'))
        qcutils.CreateSeries(ds,'Day',data_dict[bom_id][:,3],Flag=flag,Attr=qcutils.MakeAttributeDictionary(long_name='Day',units='none'))
        qcutils.CreateSeries(ds,'Hour',data_dict[bom_id][:,4],Flag=flag,Attr=qcutils.MakeAttributeDictionary(long_name='Hour',units='none'))
        qcutils.CreateSeries(ds,'Minute',data_dict[bom_id][:,5],Flag=flag,Attr=qcutils.MakeAttributeDictionary(long_name='Minute',units='none'))
        qcutils.CreateSeries(ds,'Second',Seconds,Flag=flag,Attr=qcutils.MakeAttributeDictionary(long_name='Second',units='none'))
        # now get the Python datetime
        qcutils.get_datetimefromymdhms(ds)
        # fix any time stamp issues
        if qcutils.CheckTimeStep(ds):
            qcutils.FixTimeStep(ds)
            # update the Year, Month, Day etc from the Python datetime
            qcutils.get_ymdhmsfromdatetime(ds)
        ldt = ds.series["DateTime"]["Data"]
        year = ds.series["Year"]["Data"]
        month = ds.series["Month"]["Data"]
        day = ds.series["Day"]["Data"]
        hour = ds.series["Hour"]["Data"]
        minute = ds.series["Minute"]["Data"]
        #print bom_id,ldt[-1],year[-1],month[-1],day[-1],hour[-1],minute[-1]
        # now put the data into the data structure
        attr=qcutils.MakeAttributeDictionary(long_name='Precipitation since 0900',units='mm',
                                             bom_id=str(bom_id),bom_name=bom_sites_info[site_name][str(bom_id)]["site_name"],
                                             bom_dist=bom_sites_info[site_name][str(bom_id)]["distance"])
        qcutils.CreateSeries(ds,'Precip',data_dict[bom_id][:,6],Flag=flag,Attr=attr)
        attr=qcutils.MakeAttributeDictionary(long_name='Air temperature',units='C',
                                             bom_id=str(bom_id),bom_name=bom_sites_info[site_name][str(bom_id)]["site_name"],
                                             bom_dist=bom_sites_info[site_name][str(bom_id)]["distance"])
        qcutils.CreateSeries(ds,'Ta',data_dict[bom_id][:,7],Flag=flag,Attr=attr)
        attr=qcutils.MakeAttributeDictionary(long_name='Dew point temperature',units='C',
                                             bom_id=str(bom_id),bom_name=bom_sites_info[site_name][str(bom_id)]["site_name"],
                                             bom_dist=bom_sites_info[site_name][str(bom_id)]["distance"])
        qcutils.CreateSeries(ds,'Td',data_dict[bom_id][:,8],Flag=flag,Attr=attr)
        attr=qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',
                                             bom_id=str(bom_id),bom_name=bom_sites_info[site_name][str(bom_id)]["site_name"],
                                             bom_dist=bom_sites_info[site_name][str(bom_id)]["distance"])
        qcutils.CreateSeries(ds,'RH',data_dict[bom_id][:,9],Flag=flag,Attr=attr)
        attr=qcutils.MakeAttributeDictionary(long_name='Wind speed',units='m/s',
                                             bom_id=str(bom_id),bom_name=bom_sites_info[site_name][str(bom_id)]["site_name"],
                                             bom_dist=bom_sites_info[site_name][str(bom_id)]["distance"])
        qcutils.CreateSeries(ds,'Ws',data_dict[bom_id][:,10],Flag=flag,Attr=attr)
        attr=qcutils.MakeAttributeDictionary(long_name='Wind direction',units='degT',
                                             bom_id=str(bom_id),bom_name=bom_sites_info[site_name][str(bom_id)]["site_name"],
                                             bom_dist=bom_sites_info[site_name][str(bom_id)]["distance"])
        qcutils.CreateSeries(ds,'Wd',data_dict[bom_id][:,11],Flag=flag,Attr=attr)
        attr=qcutils.MakeAttributeDictionary(long_name='Wind gust',units='m/s',
                                             bom_id=str(bom_id),bom_name=bom_sites_info[site_name][str(bom_id)]["site_name"],
                                             bom_dist=bom_sites_info[site_name][str(bom_id)]["distance"])
        qcutils.CreateSeries(ds,'Wd',data_dict[bom_id][:,12],Flag=flag,Attr=attr)
        data_dict[bom_id][:,13] = data_dict[bom_id][:,13]/float(10)
        attr=qcutils.MakeAttributeDictionary(long_name='Air Pressure',units='kPa',
                                             bom_id=str(bom_id),bom_name=bom_sites_info[site_name][str(bom_id)]["site_name"],
                                             bom_dist=bom_sites_info[site_name][str(bom_id)]["distance"])
        qcutils.CreateSeries(ds,'ps',data_dict[bom_id][:,13],Flag=flag,Attr=attr)
        # now interpolate
        for label in ["Precip","Ta","Td","RH","Ws","Wd","ps"]:
            qcts.InterpolateOverMissing(ds,series=label,maxlen=2)
        # put this stations data into the data structure dictionary
        ds_dict[bom_id] = ds

    # get the earliest start datetime and the latest end datetime
    log.info("Finding the start and end dates")
    bom_id_list = ds_dict.keys()
    ds0 = ds_dict[bom_id_list[0]]
    ldt = ds0.series["DateTime"]["Data"]
    #print bom_id_list[0],":",ldt[0],ldt[-1]
    start_date = ldt[0]
    end_date = ldt[-1]
    bom_id_list.remove(bom_id_list[0])
    for bom_id in bom_id_list:
        dsn = ds_dict[bom_id]
        ldtn = dsn.series["DateTime"]["Data"]
        #print bom_id,":",ldtn[0],ldtn[-1]
        start_date = min([start_date,ldtn[0]])
        end_date = max([end_date,ldtn[-1]])
    #print start_date,end_date

    # merge the individual data structures into a single one
    log.info("Merging file contents")
    ds_all = qcio.DataStructure()
    ds_all.globalattributes["time_step"] = 30
    ds_all.globalattributes["xl_datemode"] = 0
    ds_all.globalattributes["site_name"] = site_name
    ds_all.globalattributes["latitude"] = site_latitude
    ds_all.globalattributes["longitude"] = site_longitude
    ds_all.globalattributes["elevation"] = site_elevation
    ts = int(ds_all.globalattributes["time_step"])
    ldt_all = [result for result in qcutils.perdelta(start_date,end_date,datetime.timedelta(minutes=ts))]
    nRecs = len(ldt_all)
    ds_all.globalattributes["nc_nrecs"] = nRecs
    ds_all.series["DateTime"] = {}
    ds_all.series["DateTime"]["Data"] = ldt_all
    flag = numpy.zeros(nRecs,dtype=numpy.int32)
    ds_all.series["DateTime"]["Flag"] = flag
    ds_all.series["DateTime"]["Attr"] = {}
    ds_all.series['DateTime']["Attr"]["long_name"] = "Date-time object"
    ds_all.series['DateTime']["Attr"]["units"] = "None"
    # get the year, month, day, hour, minute and seconds from the Python datetime
    qcutils.get_ymdhmsfromdatetime(ds_all)
    # get the xlDateTime from the 
    xlDateTime = qcutils.get_xldatefromdatetime(ds_all)
    attr = qcutils.MakeAttributeDictionary(long_name="Date/time in Excel format",units="days since 1899-12-31 00:00:00")
    qcutils.CreateSeries(ds_all,"xlDateTime",xlDateTime,Flag=flag,Attr=attr)
    # loop over the stations
    for idx,bom_id in enumerate(ds_dict.keys()):
        log.info("Merging BoM site: "+str(bom_id))
        ds = ds_dict[bom_id]
        ldt = ds.series["DateTime"]["Data"]
        index = qcutils.find_indices(ldt_all,ldt)
        # loop over the variables
        for label in ["Precip","Ta","Td","RH","Ws","Wd","ps"]:
            data_all = numpy.ma.ones(nRecs,dtype=numpy.float64)*float(c.missing_value)
            flag_all = numpy.zeros(nRecs,dtype=numpy.int32)
            data,flag,attr = qcutils.GetSeriesasMA(ds,label)
            data_all[index] = data
            flag_all[index] = flag
            output_label = label+"_"+str(idx)
            attr["bom_id"] = str(bom_id)
            qcutils.CreateSeries(ds_all,output_label,data_all,Flag=flag_all,Attr=attr)
    # get precipitation per time step
    # now get precipitation per time step from the interpolated precipitation accumulated over the day
    precip_list = [x for x in ds_all.series.keys() if ("Precip" in x) and ("_QCFlag" not in x)]
    #print precip_list
    log.info("Converting 24 hour accumulated precipitation")
    for output_label in precip_list:
        accum_24hr,flag,attr = qcutils.GetSeriesasMA(ds_all,output_label)
        index = numpy.ma.where(accum_24hr<0.001)[0]
        accum_24hr[index] = float(0)
        precip = numpy.ma.ediff1d(accum_24hr,to_begin=0)
        index = [x for x in range(len(ldt_all)) if (ldt_all[x].hour==8) and (ldt_all[x].minute==30)]
        precip[index] = float(0)
        index = [x for x in range(len(ldt_all)) if (ldt_all[x].hour==9) and (ldt_all[x].minute==0)]
        precip[index] = accum_24hr[index]
        attr["long_name"] = "Precipitation total over time step"
        attr["units"] = "mm/30 minutes"
        qcutils.CreateSeries(ds_all,output_label,precip,Flag=flag,Attr=attr)

    # now write the data structure to file
    ncfile = qcio.nc_open_write(ncname)
    qcio.nc_write_series(ncfile,ds_all,ndims=1)
    log.info("Finished site: "+site_name)

print "aws2nc: All done"