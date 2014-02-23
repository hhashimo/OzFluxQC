import csv
import datetime
import meteorologicalfunctions as mf
import qcio
import qcutils
import sys

def is_yearstart(ds):
    dt = ds.series["DateTime"]["Data"]
    ts = int(ds.globalattributes["time_step"])
    if ts==60:
        start_hour = 1
        start_minute = 0
    elif ts==30:
        start_hour = 0
        start_minute = 30
    else:
        print "is_yearstart: unrecognised time step"
    result = (dt[0].month==1)&(dt[0].day==1)&\
             (dt[0].hour==start_hour)&(dt[0].minute==start_minute)
    return result

def is_yearend(ds):
    dt = ds.series["DateTime"]["Data"]
    result = (dt[-1].month==1)&(dt[-1].day==1)&\
             (dt[-1].hour==0)&(dt[-1].minute==0)
    return result

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
# get the datetime series
dt = ds.series["DateTime"]["Data"]
start_dt - dt[0]
end_dt = dt[-1]
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
