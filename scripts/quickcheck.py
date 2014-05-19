import constants as c
import sys
import logging
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdt
import meteorologicalfunctions as mf
import numpy
import scipy.ndimage as ndimage
import qcio
from qcutils import GetAltName, GetDateIndex, GetSeriesasMA, GetSeries, startlog, GetUnitsFromds, SetUnitsInds, startlog
from qcplot import tsplot

def xyplot(x,y,sub=[1,1,1],regr=0,thru0=0,title=None,xlabel=None,ylabel=None,fname=None):
    '''Generic XY scatter plot routine'''
    wspace = 0.0
    hspace = 0.0
    plt.subplot(sub[0],sub[1],sub[2])
    plt.plot(x,y,'b.')
    ax = plt.gca()
    if xlabel!=None:
        plt.xlabel(xlabel)
    if ylabel!=None:
        plt.ylabel(ylabel)
        wspace = 0.3
    if title!=None:
        plt.title(title)
        hspace = 0.3
    if regr!=0:
        coefs = numpy.ma.polyfit(x,y,1)
        xfit = numpy.ma.array([numpy.ma.minimum(x),numpy.ma.maximum(x)])
        yfit = numpy.polyval(coefs,xfit)
        r = numpy.ma.corrcoef(x,y)
        eqnstr = 'y = %.3fx + %.3f, r = %.3f'%(coefs[0],coefs[1],r[0][1])
        plt.plot(xfit,yfit,'r--',linewidth=3)
        plt.text(0.5,0.925,eqnstr,fontsize=8,horizontalalignment='center',transform=ax.transAxes)
    if thru0!=0:
        x = x[:,numpy.newaxis]
        a, _, _, _ = numpy.linalg.lstsq(x, y)
        eqnstr = 'y = %.3fx'%(a)
        plt.text(0.5,0.875,eqnstr,fontsize=8,horizontalalignment='center',transform=ax.transAxes)
    plt.subplots_adjust(wspace=wspace,hspace=hspace)

def hrplot(x,y,sub=[1,1,1],title=None,xlabel=None,ylabel=None,colours=None):
    plt.subplot(sub[0],sub[1],sub[2])
    if (y.all() is numpy.ma.masked):
        y = numpy.ma.zeros(len(y))
    if colours!=None:
        plt.scatter(x,y,c=colours)
    else:
        plt.scatter(x,y)
    plt.xlim(0,24)
    plt.xticks([0,6,12,18,24])
    if title!=None:
        plt.title(title)
    if ylabel!=None:
        plt.ylabel(ylabel)
    if xlabel!=None:
        plt.xlabel(xlabel)

nFig = 0
# start the log file
log = startlog('quickcheck','../logfiles/quickcheck.log')
# get the control file
cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0: sys.exit()
# get the netCDF filename
ncfilename = qcio.get_infilename_from_cf(cf)
# get the plot width and height
PlotWidth_landscape = float(cf['General']['PlotWidth_landscape'])
PlotHeight_landscape = float(cf['General']['PlotHeight_landscape'])
PlotWidth_portrait = float(cf['General']['PlotWidth_portrait'])
PlotHeight_portrait = float(cf['General']['PlotHeight_portrait'])
# read the netCDF file and return the data structure "ds"
log.info(' Opening and reading netCDF file '+ncfilename)
ds = qcio.nc_read_series(ncfilename)
if len(ds.series.keys())==0: log.error(' netCDF file '+ncfilename+' not found'); sys.exit()
# get the time step
ts = int(ds.globalattributes['time_step'])
# get the site name
SiteName = ds.globalattributes['site_name']
# get the datetime series
DateTime = ds.series['DateTime']['Data']
# get the initial start and end dates
StartDate = str(DateTime[0])
EndDate = str(DateTime[-1])
# find the start index of the first whole day (time=00:30)
si = GetDateIndex(DateTime,StartDate,ts=ts,default=0,match='startnextday')
# find the end index of the last whole day (time=00:00)
ei = GetDateIndex(DateTime,EndDate,ts=ts,default=-1,match='endpreviousday')
DateTime = DateTime[si:ei+1]
PlotTitle = SiteName + ': ' + str(DateTime[0]) + ' to ' + str(DateTime[-1])
# get the final start and end dates
StartDate = str(DateTime[0])
EndDate = str(DateTime[-1])
# get the 30 minute data from the data structure
log.info(' Getting data from data structure ')
#  radiation first ...
Mnth_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Month'),si=si,ei=ei)
Hour_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Hour'),si=si,ei=ei)
Mnit_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Minute'),si=si,ei=ei)
Fsd,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Fsd'),si=si,ei=ei)
if 'Fsd_syn' in ds.series.keys():
    Fsd_syn,flag = GetSeriesasMA(ds,'Fsd_syn',si=si,ei=ei)
    index = numpy.ma.where(Fsd.mask==True)[0]
    Fsd[index] = Fsd_syn[index]
night_mask = (Fsd<10)
day_mask = (Fsd>=10)
Fsd_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Fsd'),si=si,ei=ei)
Fsu_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Fsu'),si=si,ei=ei)
Fld_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Fld'),si=si,ei=ei)
Flu_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Flu'),si=si,ei=ei)
Fn_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Fn'),si=si,ei=ei)
#  then fluxes ...
Fg_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Fg'),si=si,ei=ei)
Fa2_30min = Fn_30min - Fg_30min
Fa_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Fa'),si=si,ei=ei)
index = numpy.where((Fa_30min.mask==True)&(Fa2_30min.mask==False))[0]
Fa_30min[index] = Fa2_30min[index]
Fe_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Fe'),si=si,ei=ei)
Fh_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Fh'),si=si,ei=ei)
Fc_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Fc'),si=si,ei=ei)
Fc_units = ds.series[GetAltName(cf,ds,'Fc')]['Attr']['units']
us_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'ustar'),si=si,ei=ei)
#  then meteorology ...
Ta_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Ta'),si=si,ei=ei)
H2O_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'H2O'),si=si,ei=ei)
H2O_units = ds.series[GetAltName(cf,ds,'H2O')]['Attr']['units']
CO2_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'CO2'),si=si,ei=ei)
CO2_units = ds.series[GetAltName(cf,ds,'CO2')]['Attr']['units']
Rain_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Precip'),si=si,ei=ei)
Ws_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Ws'),si=si,ei=ei)
#  then soil ...
Sws_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Sws'),si=si,ei=ei)
Ts_30min,flag = GetSeriesasMA(ds,GetAltName(cf,ds,'Ts'),si=si,ei=ei)

# get the number of days in the data set
ntsInDay = float(24.0*60.0/float(ts))
if math.modf(ntsInDay)[0]!=0:
    print 'quickcheck: Time step is not a sub-multiple of 60 minutes ', ts
    sys.exit
ntsInDay = int(ntsInDay)
nDays = float(len(DateTime))/ntsInDay
if math.modf(nDays)[0]!=0:
    print 'quickcheck: Not a whole number of days ', nDays
    sys.exit
nDays = int(nDays)

# *** start of section based on 30 minute data ***
# scatter plot of (Fh+Fe) versys Fa, all data
log.info(' Doing surface energy balance plots ')
mask = numpy.ma.mask_or(Fa_30min.mask,Fe_30min.mask)
mask = numpy.ma.mask_or(mask,Fh_30min.mask)
Fa_SEB = numpy.ma.array(Fa_30min,mask=mask)     # apply the mask
FhpFe_SEB = numpy.ma.array(Fh_30min,mask=mask) + numpy.ma.array(Fe_30min,mask=mask)
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(8,8))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
xyplot(Fa_SEB,FhpFe_SEB,sub=[2,2,1],regr=1,title="All hours",xlabel='Fa (W/m2)',ylabel='Fh+Fe (W/m2)')
# scatter plot of (Fh+Fe) versus Fa, 24 hour averages
Fa_daily = Fa_30min.reshape(nDays,ntsInDay)
Fe_daily = Fe_30min.reshape(nDays,ntsInDay)
Fh_daily = Fh_30min.reshape(nDays,ntsInDay)
mask = numpy.ma.mask_or(Fa_daily.mask,Fe_daily.mask)
mask = numpy.ma.mask_or(mask,Fh_daily.mask)
Fa_daily = numpy.ma.array(Fa_daily,mask=mask)         # apply the mask
Fe_daily = numpy.ma.array(Fe_daily,mask=mask)
Fh_daily = numpy.ma.array(Fh_daily,mask=mask)
Fa_daily_avg = numpy.ma.average(Fa_daily,axis=1)      # get the daily average
Fe_daily_avg = numpy.ma.average(Fe_daily,axis=1)
Fh_daily_avg = numpy.ma.average(Fh_daily,axis=1)
FhpFe_daily_avg = Fh_daily_avg + Fe_daily_avg
xyplot(Fa_daily_avg,FhpFe_daily_avg,sub=[2,2,2],regr=1,thru0=1,title="Daily Average",xlabel='Fa (W/m2)',ylabel='Fh+Fe (W/m2)')
# scatter plot of (Fh+Fe) versus Fa, day time
Fa_day = numpy.ma.masked_where(day_mask==False,Fa_30min)
Fe_day = numpy.ma.masked_where(day_mask==False,Fe_30min)
Fh_day = numpy.ma.masked_where(day_mask==False,Fh_30min)
mask = numpy.ma.mask_or(Fa_day.mask,Fe_day.mask)
mask = numpy.ma.mask_or(mask,Fh_day.mask)
Fa_day = numpy.ma.array(Fa_day,mask=mask)         # apply the mask
Fe_day = numpy.ma.array(Fe_day,mask=mask)
Fh_day = numpy.ma.array(Fh_day,mask=mask)
FhpFe_day = Fh_day + Fe_day
xyplot(Fa_day,FhpFe_day,sub=[2,2,3],regr=1,title="Day",xlabel='Fa (W/m2)',ylabel='Fh+Fe (W/m2)')
# scatter plot of (Fh+Fe) versus Fa, night time
Fa_night = numpy.ma.masked_where(night_mask==False,Fa_30min)
Fe_night = numpy.ma.masked_where(night_mask==False,Fe_30min)
Fh_night = numpy.ma.masked_where(night_mask==False,Fh_30min)
mask = numpy.ma.mask_or(Fa_night.mask,Fe_night.mask)
mask = numpy.ma.mask_or(mask,Fh_night.mask)
Fa_night = numpy.ma.array(Fa_night,mask=mask)         # apply the mask
Fe_night = numpy.ma.array(Fe_night,mask=mask)
Fh_night = numpy.ma.array(Fh_night,mask=mask)
FhpFe_night = Fh_night + Fe_night
xyplot(Fa_night,FhpFe_night,sub=[2,2,4],regr=1,title="Night",xlabel='Fa (W/m2)',ylabel='Fh+Fe (W/m2)')
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'SEB_30minutes.png'
fig.savefig(figname,format='png')

# *** start of section based on daily averages ***
log.info(' Getting daily averages from 30 minute data ')
MTFmt = mdt.DateFormatter('%m/%Y')
# reshape the 1D array of 30 minute data into a 2D array of (nDays,ntsInDay)
DT_daily = DateTime[0::ntsInDay]
Mnth_daily = Mnth_30min.reshape(nDays,ntsInDay)
Hour_daily = Hour_30min.reshape(nDays,ntsInDay)
Mnit_daily = Mnit_30min.reshape(nDays,ntsInDay)
dm_daily = day_mask.reshape(nDays,ntsInDay)
nm_daily = night_mask.reshape(nDays,ntsInDay)
Fn_daily = Fn_30min.reshape(nDays,ntsInDay)
Fa_daily = Fa_30min.reshape(nDays,ntsInDay)
Fe_daily = Fe_30min.reshape(nDays,ntsInDay)
Fh_daily = Fh_30min.reshape(nDays,ntsInDay)
Fc_daily = Fc_30min.reshape(nDays,ntsInDay)
Rain_daily = Rain_30min.reshape(nDays,ntsInDay)
Sws_daily = Sws_30min.reshape(nDays,ntsInDay)
Ts_daily = Ts_30min.reshape(nDays,ntsInDay)
us_daily = us_30min.reshape(nDays,ntsInDay)

# get the SEB ratio
# get the daytime data, defined by Fsd>10 W/m2
Fa_day = numpy.ma.masked_where(nm_daily==True,Fa_daily)
Fe_day = numpy.ma.masked_where(nm_daily==True,Fe_daily)
Fh_day = numpy.ma.masked_where(nm_daily==True,Fh_daily)
mask = numpy.ma.mask_or(Fa_day.mask,Fe_day.mask)  # mask based on dependencies, set all to missing if any missing
mask = numpy.ma.mask_or(mask,Fh_day.mask)
Fa_day = numpy.ma.array(Fa_day,mask=mask)         # apply the mask
Fe_day = numpy.ma.array(Fe_day,mask=mask)
Fh_day = numpy.ma.array(Fh_day,mask=mask)
Fa_day_avg = numpy.ma.average(Fa_day,axis=1)      # get the daily average
Fe_day_avg = numpy.ma.average(Fe_day,axis=1)
Fh_day_avg = numpy.ma.average(Fh_day,axis=1)      # get the number of values in the daily average
SEB_day_num = numpy.ma.count(Fh_day,axis=1)       # get the SEB ratio
SEB_day_avg = (Fe_day_avg+Fh_day_avg)/Fa_day_avg
SEB_day_avg = numpy.ma.masked_where(SEB_day_num<=5,SEB_day_avg)
index = numpy.ma.where(SEB_day_avg.mask==True)
SEB_day_num[index] = 0

# get the EF
# get the daytime data, defined by Fsd>10 W/m2
Fa_day = numpy.ma.masked_where(nm_daily==True,Fa_daily)
Fe_day = numpy.ma.masked_where(nm_daily==True,Fe_daily)
mask = numpy.ma.mask_or(Fa_day.mask,Fe_day.mask)  # mask based on dependencies, set all to missing if any missing
Fa_day = numpy.ma.array(Fa_day,mask=mask)         # apply the mask
Fe_day = numpy.ma.array(Fe_day,mask=mask)
Fa_day_avg = numpy.ma.average(Fa_day,axis=1)      # get the daily average
Fe_day_avg = numpy.ma.average(Fe_day,axis=1)
EF_day_num = numpy.ma.count(Fe_day,axis=1)        # get the number of values in the daily average
EF_day_avg = Fe_day_avg/Fa_day_avg                # get the EF ratio
EF_day_avg = numpy.ma.masked_where(EF_day_num<=5,EF_day_avg)
index = numpy.ma.where(EF_day_avg.mask==True)
EF_day_num[index] = 0

# get the BR
# get the daytime data, defined by Fsd>10 W/m2
Fe_day = numpy.ma.masked_where(nm_daily==True,Fe_daily)
Fh_day = numpy.ma.masked_where(nm_daily==True,Fh_daily)
mask = numpy.ma.mask_or(Fe_day.mask,Fh_day.mask)  # mask based on dependencies, set all to missing if any missing
Fe_day = numpy.ma.array(Fe_day,mask=mask)         # apply the mask
Fh_day = numpy.ma.array(Fh_day,mask=mask)
Fe_day_avg = numpy.ma.average(Fe_day,axis=1)      # get the daily average
Fh_day_avg = numpy.ma.average(Fh_day,axis=1)
BR_day_num = numpy.ma.count(Fh_day,axis=1)        # get the number of values in the daily average
BR_day_avg = Fh_day_avg/Fe_day_avg                # get the BR ratio
BR_day_avg = numpy.ma.masked_where(BR_day_num<=5,BR_day_avg)
index = numpy.ma.where(BR_day_avg.mask==True)
BR_day_num[index] = 0

# get the Wue
# get the daytime data, defined by Fsd>10 W/m2
Fe_day = numpy.ma.masked_where(nm_daily==True,Fe_daily)
Fc_day = numpy.ma.masked_where(nm_daily==True,Fc_daily)
mask = numpy.ma.mask_or(Fe_day.mask,Fc_day.mask)  # mask based on dependencies, set all to missing if any missing
Fe_day = numpy.ma.array(Fe_day,mask=mask)         # apply the mask
Fc_day = numpy.ma.array(Fc_day,mask=mask)
Fe_day_avg = numpy.ma.average(Fe_day,axis=1)      # get the daily average
Fc_day_avg = numpy.ma.average(Fc_day,axis=1)
WUE_day_num = numpy.ma.count(Fc_day,axis=1)       # get the number of values in the daily average
WUE_day_avg = Fc_day_avg/Fe_day_avg
WUE_day_avg = numpy.ma.masked_where(WUE_day_num<=5,WUE_day_avg)
index = numpy.ma.where(WUE_day_avg.mask==True)
WUE_day_num[index] = 0
# get the soil moisture
Sws_daily_avg = numpy.ma.average(Sws_daily,axis=1)
Sws_daily_num = numpy.ma.count(Sws_daily,axis=1)
# get the rainfall
Rain_daily_sum = numpy.ma.sum(Rain_daily,axis=1)
Rain_daily_num = numpy.ma.count(Rain_daily,axis=1)
# plot the SEB, EF and Wue
log.info(' Doing the daily ratios plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_landscape,PlotHeight_landscape))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
tsplot(DT_daily,SEB_day_avg,sub=[6,1,1],colours=SEB_day_num,ylabel='(Fh+Fe)/Fa',lineat=1)
tsplot(DT_daily,EF_day_avg,sub=[6,1,2],colours=EF_day_num,ylabel='EF=Fe/Fa')
tsplot(DT_daily,BR_day_avg,sub=[6,1,3],colours=BR_day_num,ylabel='BR=Fh/Fe')
tsplot(DT_daily,WUE_day_avg,sub=[6,1,4],colours=WUE_day_num,ylabel='WUE=Fc/Fe',lineat=0)
tsplot(DT_daily,Sws_daily_avg,sub=[6,1,5],colours=Sws_daily_num,ylabel='Sws')
tsplot(DT_daily,Rain_daily_sum,sub=[6,1,6],colours=Rain_daily_num,ylabel='Rain')
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DailyRatios.png'
fig.savefig(figname,format='png')

# now we do the daily averages of the fluxes and the meteorology
# get the 1D array of 30 minute data into a 2D array with a dimension for
#  the day number and a dimension for the time of day
Fsd_daily = Fsd_30min.reshape(nDays,ntsInDay)
Fsu_daily = Fsu_30min.reshape(nDays,ntsInDay)
Fld_daily = Fld_30min.reshape(nDays,ntsInDay)
Flu_daily = Flu_30min.reshape(nDays,ntsInDay)
Fn_daily = Fn_30min.reshape(nDays,ntsInDay)
Fg_daily = Fg_30min.reshape(nDays,ntsInDay)
Fsd_day = numpy.ma.masked_where(nm_daily==True,Fsd_daily)
Fsu_day = numpy.ma.masked_where(nm_daily==True,Fsu_daily)
Fld_day = numpy.ma.masked_where(nm_daily==True,Fld_daily)
Flu_day = numpy.ma.masked_where(nm_daily==True,Flu_daily)
Fn_day = numpy.ma.masked_where(nm_daily==True,Fn_daily)
Fg_day = numpy.ma.masked_where(nm_daily==True,Fg_daily)
Fsd_day_avg = numpy.ma.average(Fsd_day,axis=1)
Fsu_day_avg = numpy.ma.average(Fsu_day,axis=1)
Fld_day_avg = numpy.ma.average(Fld_day,axis=1)
Flu_day_avg = numpy.ma.average(Flu_day,axis=1)
Fn_day_avg = numpy.ma.average(Fn_day,axis=1)
Fg_day_avg = numpy.ma.average(Fg_day,axis=1)
Fsd_day_num = numpy.ma.count(Fsd_day,axis=1)
Fsu_day_num = numpy.ma.count(Fsu_day,axis=1)
Fld_day_num = numpy.ma.count(Fld_day,axis=1)
Flu_day_num = numpy.ma.count(Flu_day,axis=1)
Fn_day_num = numpy.ma.count(Fn_day,axis=1)
Fg_day_num = numpy.ma.count(Fg_day,axis=1)
log.info(' Doing the daily radiation plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_landscape,PlotHeight_landscape))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
tsplot(DT_daily,Fsd_day_avg,sub=[6,1,1],colours=Fsd_day_num,ylabel='Fsd (W/m2)')
tsplot(DT_daily,Fsu_day_avg,sub=[6,1,2],colours=Fsu_day_num,ylabel='Fsu (W/m2)')
tsplot(DT_daily,Fld_day_avg,sub=[6,1,3],colours=Fld_day_num,ylabel='Fld (W/m2)')
tsplot(DT_daily,Flu_day_avg,sub=[6,1,4],colours=Flu_day_num,ylabel='Flu (W/m2)')
tsplot(DT_daily,Fn_day_avg,sub=[6,1,5],colours=Fn_day_num,ylabel='Fn (W/m2)')
tsplot(DT_daily,Fg_day_avg,sub=[6,1,6],colours=Fg_day_num,ylabel='Fg (W/m2)')
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DailyRadn.png'
fig.savefig(figname,format='png')

Fsd_daily = Fsd_30min.reshape(nDays,ntsInDay)
Fa_daily = Fa_30min.reshape(nDays,ntsInDay)
Fe_daily = Fe_30min.reshape(nDays,ntsInDay)
Fh_daily = Fh_30min.reshape(nDays,ntsInDay)
Fc_daily = Fc_30min.reshape(nDays,ntsInDay)
# ... then get the day time values only (defined by Fsd>10 W/m2)
Fsd_day = numpy.ma.masked_where(nm_daily==True,Fsd_daily)
Fa_day = numpy.ma.masked_where(nm_daily==True,Fa_daily)
Fe_day = numpy.ma.masked_where(nm_daily==True,Fe_daily)
Fh_day = numpy.ma.masked_where(nm_daily==True,Fh_daily)
Fc_day = numpy.ma.masked_where(nm_daily==True,Fc_daily)
Fc_night = numpy.ma.masked_where(nm_daily==True,Fc_daily)
# ... then get the daily averages
Fsd_day_avg = numpy.ma.average(Fsd_day,axis=1)      # get the daily average
Fa_day_avg = numpy.ma.average(Fa_day,axis=1)      # get the daily average
Fe_day_avg = numpy.ma.average(Fe_day,axis=1)      # get the daily average
Fh_day_avg = numpy.ma.average(Fh_day,axis=1)      # get the daily average
Fc_day_avg = numpy.ma.average(Fc_day,axis=1)      # get the daily average
Fc_night_avg = numpy.ma.average(Fc_night,axis=1)      # get the daily average
# ... then the number of values in each day time block
Fsd_day_num = numpy.ma.count(Fsd_day,axis=1)
Fa_day_num = numpy.ma.count(Fa_day,axis=1)
Fe_day_num = numpy.ma.count(Fe_day,axis=1)
Fh_day_num = numpy.ma.count(Fh_day,axis=1)
Fc_day_num = numpy.ma.count(Fc_day,axis=1)
Fc_night_num = numpy.ma.count(Fc_night,axis=1)
# ... now plot the day time averages with the colour of the points controlled
#     by the number of values used to get the average
log.info(' Doing the daily fluxes plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_landscape,PlotHeight_landscape))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
tsplot(DT_daily,Fsd_day_avg,sub=[5,1,1],colours=Fsd_day_num,ylabel='Fsd (W/m2)')
tsplot(DT_daily,Fa_day_avg,sub=[5,1,2],colours=Fa_day_num,ylabel='Fa (W/m2)')
tsplot(DT_daily,Fe_day_avg,sub=[5,1,3],colours=Fe_day_num,ylabel='Fe (W/m2)')
tsplot(DT_daily,Fh_day_avg,sub=[5,1,4],colours=Fh_day_num,ylabel='Fh (W/m2)')
tsplot(DT_daily,Fc_day_avg,sub=[5,1,5],colours=Fc_day_num,ylabel='Fc ('+Fc_units+')',lineat=0)
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DailyFluxes.png'
fig.savefig(figname,format='png')

Ta_daily = Ta_30min.reshape(nDays,ntsInDay)
H2O_daily = H2O_30min.reshape(nDays,ntsInDay)
CO2_daily = CO2_30min.reshape(nDays,ntsInDay)
Ws_daily = Ws_30min.reshape(nDays,ntsInDay)
CO2_day = numpy.ma.masked_where(nm_daily==True,CO2_daily)
Ta_daily_avg = numpy.ma.average(Ta_daily,axis=1)      # get the daily average
Ta_daily_num = numpy.ma.count(Ta_daily,axis=1)
H2O_daily_avg = numpy.ma.average(H2O_daily,axis=1)      # get the daily average
H2O_daily_num = numpy.ma.count(H2O_daily,axis=1)
CO2_day_avg = numpy.ma.average(CO2_day,axis=1)          # get the daily average
CO2_day_num = numpy.ma.count(CO2_day,axis=1)
Ws_daily_avg = numpy.ma.average(Ws_daily,axis=1)      # get the daily average
Ws_daily_num = numpy.ma.count(Ws_daily,axis=1)
log.info(' Doing the daily meteorology plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_landscape,PlotHeight_landscape))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
tsplot(DT_daily,Ta_daily_avg,sub=[5,1,1],colours=Ta_daily_num,ylabel='Ta (C)')
tsplot(DT_daily,H2O_daily_avg,sub=[5,1,2],colours=H2O_daily_num,ylabel='H2O ('+H2O_units+')')
tsplot(DT_daily,CO2_day_avg,sub=[5,1,3],colours=CO2_day_num,ylabel='CO2 ('+CO2_units+')')
tsplot(DT_daily,Ws_daily_avg,sub=[5,1,4],colours=Ws_daily_num,ylabel='WS (m/s)')
tsplot(DT_daily,Rain_daily_sum,sub=[5,1,5],colours=Rain_daily_num,ylabel='Rain (mm)')
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DailyMet.png'
fig.savefig(figname,format='png')

Ta_daily = Ta_30min.reshape(nDays,ntsInDay)
Ts_daily = Ts_30min.reshape(nDays,ntsInDay)
Sws_daily = Sws_30min.reshape(nDays,ntsInDay)
Fg_daily = Fg_30min.reshape(nDays,ntsInDay)
Rain_daily = Rain_30min.reshape(nDays,ntsInDay)
Ta_daily_avg = numpy.ma.average(Ta_daily,axis=1)      # get the daily average
Ta_daily_num = numpy.ma.count(Ta_daily,axis=1)
Ts_daily_avg = numpy.ma.average(Ts_daily,axis=1)      # get the daily average
Ts_daily_num = numpy.ma.count(Ts_daily,axis=1)
Sws_day_avg = numpy.ma.average(Sws_daily,axis=1)          # get the daily average
Sws_day_num = numpy.ma.count(Sws_daily,axis=1)
Fg_daily_avg = numpy.ma.average(Fg_daily,axis=1)      # get the daily average
Fg_daily_num = numpy.ma.count(Fg_daily,axis=1)
Fg_daily_avg = numpy.ma.average(Fg_daily,axis=1)      # get the daily average
Rain_daily_sum = numpy.ma.sum(Rain_daily,axis=1)
Rain_daily_num = numpy.ma.count(Rain_daily,axis=1)
log.info(' Doing the daily soil data plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_landscape,PlotHeight_landscape))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
tsplot(DT_daily,Ta_daily_avg,sub=[5,1,1],colours=Ta_daily_num,ylabel='Ta (C)')
tsplot(DT_daily,Ts_daily_avg,sub=[5,1,2],colours=Ts_daily_num,ylabel='Ts (C)')
tsplot(DT_daily,Sws_day_avg,sub=[5,1,3],colours=Sws_day_num,ylabel='Sws (frac)')
tsplot(DT_daily,Fg_daily_avg,sub=[5,1,4],colours=Fg_daily_num,ylabel='Fg (W/m2)')
tsplot(DT_daily,Rain_daily_sum,sub=[5,1,5],colours=Rain_daily_num,ylabel='Rain (mm)')
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DailySoil.png'
fig.savefig(figname,format='png')

MnthList = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
# plot Fsd
log.info(' Doing the diurnal Fsd by month plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_portrait,PlotHeight_portrait))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
j = 0
for i in [12,1,2,3,4,5,6,7,8,9,10,11]:
    j = j + 1
    index = numpy.where(Mnth_daily==i)[0]
    if len(index)!=0:
        hr = Hour_daily[index]+Mnit_daily[index]/float(60)
        Fsd_hr_avg = numpy.ma.average(Fsd_daily[index],axis=0)
        Fsd_hr_num = numpy.ma.count(Fsd_daily[index],axis=0)
        if j in [1,2,3,4,5,6,7,8,9]:
            xlabel = None
        else:
            xlabel = 'Hour'
        if j in [2,3,5,6,8,9,11,12]:
            ylabel = None
        else:
            ylabel = 'Fsd (W/m2)'
        hrplot(hr[0],Fsd_hr_avg,sub=[4,3,j],
               title=MnthList[i-1],xlabel=xlabel,ylabel=ylabel,
               colours=Fsd_hr_num)
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DiurnalFsdByMonth.png'
fig.savefig(figname,format='png')

# plot Fa
log.info(' Doing the diurnal Fa by month plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_portrait,PlotHeight_portrait))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
j = 0
for i in [12,1,2,3,4,5,6,7,8,9,10,11]:
    j = j + 1
    index = numpy.where(Mnth_daily==i)[0]
    if len(index)!=0:
        hr = Hour_daily[index]+Mnit_daily[index]/float(60)
        Fa_hr_avg = numpy.ma.average(Fa_daily[index],axis=0)
        Fa_hr_num = numpy.ma.count(Fa_daily[index],axis=0)
        if j in [1,2,3,4,5,6,7,8,9]:
            xlabel = None
        else:
            xlabel = 'Hour'
        if j in [2,3,5,6,8,9,11,12]:
            ylabel = None
        else:
            ylabel = 'Fa (W/m2)'
        hrplot(hr[0],Fa_hr_avg,sub=[4,3,j],
               title=MnthList[i-1],xlabel=xlabel,ylabel=ylabel,
               colours=Fa_hr_num)
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DiurnalFaByMonth.png'
fig.savefig(figname,format='png')

# plot Fn
log.info(' Doing the diurnal Fn by month plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_portrait,PlotHeight_portrait))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
j = 0
for i in [12,1,2,3,4,5,6,7,8,9,10,11]:
    j = j + 1
    index = numpy.where(Mnth_daily==i)[0]
    if len(index)!=0:
        hr = Hour_daily[index]+Mnit_daily[index]/float(60)
        Fn_hr_avg = numpy.ma.average(Fn_daily[index],axis=0)
        Fn_hr_num = numpy.ma.count(Fn_daily[index],axis=0)
        if j in [1,2,3,4,5,6,7,8,9]:
            xlabel = None
        else:
            xlabel = 'Hour'
        if j in [2,3,5,6,8,9,11,12]:
            ylabel = None
        else:
            ylabel = 'Fn (W/m2)'
        hrplot(hr[0],Fn_hr_avg,sub=[4,3,j],
               title=MnthList[i-1],xlabel=xlabel,ylabel=ylabel,
               colours=Fn_hr_num)
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DiurnalFnByMonth.png'
fig.savefig(figname,format='png')

# plot Fg
log.info(' Doing the diurnal Fg by month plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_portrait,PlotHeight_portrait))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
j = 0
for i in [12,1,2,3,4,5,6,7,8,9,10,11]:
    j = j + 1
    index = numpy.where(Mnth_daily==i)[0]
    if len(index)!=0:
        hr = Hour_daily[index]+Mnit_daily[index]/float(60)
        Fg_hr_avg = numpy.ma.average(Fg_daily[index],axis=0)
        Fg_hr_num = numpy.ma.count(Fg_daily[index],axis=0)
        if j in [1,2,3,4,5,6,7,8,9]:
            xlabel = None
        else:
            xlabel = 'Hour'
        if j in [2,3,5,6,8,9,11,12]:
            ylabel = None
        else:
            ylabel = 'Fg (W/m2)'
        hrplot(hr[0],Fg_hr_avg,sub=[4,3,j],
               title=MnthList[i-1],xlabel=xlabel,ylabel=ylabel,
               colours=Fg_hr_num)
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DiurnalFgByMonth.png'
fig.savefig(figname,format='png')

# plot Ts
log.info(' Doing the diurnal Ts by month plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_portrait,PlotHeight_portrait))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
j = 0
for i in [12,1,2,3,4,5,6,7,8,9,10,11]:
    j = j + 1
    index = numpy.where(Mnth_daily==i)[0]
    if len(index)!=0:
        hr = Hour_daily[index]+Mnit_daily[index]/float(60)
        Ts_hr_avg = numpy.ma.average(Ts_daily[index],axis=0)
        Ts_hr_num = numpy.ma.count(Ts_daily[index],axis=0)
        if j in [1,2,3,4,5,6,7,8,9]:
            xlabel = None
        else:
            xlabel = 'Hour'
        if j in [2,3,5,6,8,9,11,12]:
            ylabel = None
        else:
            ylabel = 'Ts (C)'
        hrplot(hr[0],Ts_hr_avg,sub=[4,3,j],
               title=MnthList[i-1],xlabel=xlabel,ylabel=ylabel,
               colours=Fg_hr_num)
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DiurnalTsByMonth.png'
fig.savefig(figname,format='png')

# plot Fh
log.info(' Doing the diurnal Fh by month plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_portrait,PlotHeight_portrait))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
j = 0
for i in [12,1,2,3,4,5,6,7,8,9,10,11]:
    j = j + 1
    index = numpy.where(Mnth_daily==i)[0]
    if len(index)!=0:
        hr = Hour_daily[index]+Mnit_daily[index]/float(60)
        Fh_hr_avg = numpy.ma.average(Fh_daily[index],axis=0)
        Fh_hr_num = numpy.ma.count(Fh_daily[index],axis=0)
        if j in [1,2,3,4,5,6,7,8,9]:
            xlabel = None
        else:
            xlabel = 'Hour'
        if j in [2,3,5,6,8,9,11,12]:
            ylabel = None
        else:
            ylabel = 'Fh (W/m2)'
        hrplot(hr[0],Fh_hr_avg,sub=[4,3,j],
               title=MnthList[i-1],xlabel=xlabel,ylabel=ylabel,
               colours=Fh_hr_num)
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DiurnalFhByMonth.png'
fig.savefig(figname,format='png')

# plot Fe
log.info(' Doing the diurnal Fe by month plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_portrait,PlotHeight_portrait))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
j = 0
for i in [12,1,2,3,4,5,6,7,8,9,10,11]:
    j = j + 1
    index = numpy.where(Mnth_daily==i)[0]
    if len(index)!=0:
        hr = Hour_daily[index]+Mnit_daily[index]/float(60)
        Fe_hr_avg = numpy.ma.average(Fe_daily[index],axis=0)
        Fe_hr_num = numpy.ma.count(Fe_daily[index],axis=0)
        if j in [1,2,3,4,5,6,7,8,9]:
            xlabel = None
        else:
            xlabel = 'Hour'
        if j in [2,3,5,6,8,9,11,12]:
            ylabel = None
        else:
            ylabel = 'Fe (W/m2)'
        hrplot(hr[0],Fe_hr_avg,sub=[4,3,j],
               title=MnthList[i-1],xlabel=xlabel,ylabel=ylabel,
               colours=Fe_hr_num)
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DiurnalFeByMonth.png'
fig.savefig(figname,format='png')

# plot Fc
log.info(' Doing the diurnal Fc by month plot ')
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(PlotWidth_portrait,PlotHeight_portrait))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
j = 0
for i in [12,1,2,3,4,5,6,7,8,9,10,11]:
    j = j + 1
    index = numpy.where(Mnth_daily==i)[0]
    if len(index)!=0:
        hr = Hour_daily[index]+Mnit_daily[index]/float(60)
        Fc_hr_avg = numpy.ma.average(Fc_daily[index],axis=0)
        Fc_hr_num = numpy.ma.count(Fc_daily[index],axis=0)
        if j in [1,2,3,4,5,6,7,8,9]:
            xlabel = None
        else:
            xlabel = 'Hour'
        if j in [2,3,5,6,8,9,11,12]:
            ylabel = None
        else:
            ylabel = 'Fc ('+Fc_units+')'
        hrplot(hr[0],Fc_hr_avg,sub=[4,3,j],
               title=MnthList[i-1],xlabel=xlabel,ylabel=ylabel,
               colours=Fc_hr_num)
#fig.show()
figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_QC_'+'DiurnalFcByMonth.png'
fig.savefig(figname,format='png')

plt.show()