import constants as c
import datetime
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import meteorologicalfunctions as mf
import numpy
import qcck
import qcio
import qcts
import qcutils
import qcplot
import sys

def plot_fingerprint(data,xlabel=None,sub=[1,1,1],extent=None,ticks=None):
    loc,fmt = qcplot.get_ticks(datetime.datetime.fromordinal(sd),datetime.datetime.fromordinal(ed))
    ax = plt.subplot(sub[0],sub[1],sub[2])
    plt.imshow(data,extent=extent,aspect='auto',origin='lower')
    ax.yaxis.set_major_locator(loc)
    ax.yaxis.set_major_formatter(fmt)
    plt.colorbar(orientation='horizontal',fraction=0.02,pad=0.075,ticks=ticks)
    plt.xticks([0,6,12,18,24])
    if xlabel != None: plt.xlabel(xlabel)
    if sub[2] != 1: plt.setp(ax.get_yticklabels(), visible=False)

# open the logging file
log = qcutils.startlog('qc','../logfiles/fingerprint.log')
# get the control file
cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0: sys.exit()
infilename = qcio.get_infilename_from_cf(cf)
# get the default plot width and height
PlotWidth = float(cf['General']['PlotWidth'])
PlotHeight = float(cf['General']['PlotHeight'])
# read the netCDF file and return the data structure "ds"
ds = qcio.nc_read_series(infilename)
if len(ds.series.keys())==0: print time.strftime('%X')+' netCDF file '+infilename+' not found'; sys.exit()
StartDate = str(ds.series['DateTime']['Data'][0])
EndDate = str(ds.series['DateTime']['Data'][-1])
# check to see if the required fields are in ds.globalattributes, if they aren't, get them from the control file
# get the time step
if 'time_step' in ds.globalattributes.keys():
    ts = int(qcutils.GetGlobalAttributeValue(cf,ds,'time_step'))
elif 'time_step' in cf['General'].keys():
    ts = int(cf['General']['time_step'])
    ds.globalattributes['time_step'] = ts
else:
    print 'fingerprint: cant find a value for the time step'
    sys.exit()
# get the site name
if 'site_name' in ds.globalattributes.keys():
    site_name = str(qcutils.GetGlobalAttributeValue(cf,ds,'site_name'))
elif 'site_name' in cf['General'].keys():
    site_name = str(cf['General']['site_name'])
    ds.globalattributes['site_name'] = site_name
else:
    print 'fingerprint: cant find a value for the site name'
    sys.exit()
# get the data level
if 'nc_level' in ds.globalattributes.keys():
    level = str(qcutils.GetGlobalAttributeValue(cf,ds,'nc_level'))
elif 'nc_level' in cf['General'].keys():
    level = str(cf['General']['nc_level'])
    ds.globalattributes['nc_level'] = level
else:
    print 'fingerprint: cant find a value for the nc_level'
    sys.exit()
if 'nc_nrecs' in ds.globalattributes.keys():
    nRecsInFile = int(qcutils.GetGlobalAttributeValue(cf,ds,'nc_nrecs'))
else:
    nRecsInFile = len(ds.series['DateTime']['Data'])
    ds.globalattributes['nc_nrecs'] = nRecsInFile
# check for time gaps in the file
has_gaps = qcutils.CheckTimeStep(ds)
if has_gaps:
    qcutils.FixTimeGaps(ds)
TitleStr = site_name+' '+level
nPerHr = int(float(60)/ts+0.5)
nPerDay = int(float(24)*nPerHr+0.5)

# get the datetime series
DateTime = ds.series['DateTime']['Data']

# find the start index of the first whole day (time=00:30)
si = qcutils.GetDateIndex(DateTime,StartDate,ts=ts,default=0,match='startnextday')
# find the end index of the last whole day (time=00:00)
ei = qcutils.GetDateIndex(DateTime,EndDate,ts=ts,default=-1,match='endpreviousday')

DateTime = DateTime[si:ei+1]
sd = datetime.datetime.toordinal(DateTime[0])
ed = datetime.datetime.toordinal(DateTime[-1])
TitleStr = TitleStr+' from '+str(DateTime[0])+' to '+str(DateTime[-1])
nDays = len(DateTime)/nPerDay

for nFig in cf['Plots'].keys():
    n = 0
    fig = plt.figure(nFig,figsize=[15,10])
    plt.figtext(0.5,0.95,TitleStr,horizontalalignment='center')
    SeriesList = qcutils.GetPlotVariableNamesFromCF(cf,nFig)
    nPlots = len(SeriesList)
    for ThisOne in SeriesList:
        n += 1
        VarName = qcutils.GetAltNameFromCF(cf,ThisOne)
        ticks = qcutils.GetcbTicksFromCF(cf,ThisOne)
        lower, upper = qcutils.GetRangesFromCF(cf,ThisOne)
        data_30min,flag = qcutils.GetSeriesasMA(ds,VarName,si=si,ei=ei)
        data_30min = qcck.cliptorange(data_30min, lower, upper)
        data_daily = data_30min.reshape(nDays,nPerDay)
        if 'Units' in ds.series[VarName]['Attr']: units = str(ds.series[VarName]['Attr']['Units'])
        if 'units' in ds.series[VarName]['Attr']: units = str(ds.series[VarName]['Attr']['units'])
        label = VarName + ' (' + units + ')'
        plot_fingerprint(data_daily,xlabel=label,sub=[1,nPlots,n],extent=[0,24,sd,ed],ticks=ticks)
    pngname = '../plots/'+site_name.replace(' ','')+'_'+level+'_'
    pngname = pngname+qcutils.GetPlotTitleFromCF(cf,nFig).replace(' ','_')+'.png'
    fig.savefig(pngname,format='png')

plt.draw()
plt.show()
