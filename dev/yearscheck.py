import matplotlib.pyplot as plt
import numpy
import sys
sys.path.append("../scripts/")
import qcio
import qcutils

def get_daynightmask(ds,Fsd_threshold=10,si=0,ei=-1):
    Fsd,f = qcutils.GetSeriesasMA(ds,'Fsd',si=si,ei=ei)
    if "Fsd_syn" in ds.series.keys():
        Fsd_syn,f = qcutils.GetSeriesasMA(ds,'Fsd_syn',si=si,ei=ei)
        index = numpy.ma.where(Fsd.mask==True)[0]
        Fsd[index] = Fsd_syn[index]
    night_mask = (Fsd<Fsd_threshold)
    day_mask = (Fsd>=Fsd_threshold)
    return day_mask,night_mask

def get_dailyaverage(series,ts,mask=None):
    ntsInDay = int(float(24.0*60.0/float(ts)))
    nDays = float(len(series))/ntsInDay
    series_daily = series.reshape(nDays,ntsInDay)
    if mask!=None:
        mask_daily = mask.reshape(nDays,ntsInDay)
        series_daily = numpy.ma.masked_where(mask_daily==False,series_daily)
    series_daily_avg = numpy.ma.average(series_daily,axis=1)
    series_daily_num = numpy.ma.count(series_daily,axis=1)
    return series_daily_avg,series_daily_num
    
def get_sieiforwholeyear(dt,yr,ts):
    sd = str(yr)+"-01-01 00:30"
    ed = str(yr+1)+"-01-01 00:00"
    si = qcutils.GetDateIndex(dt,sd,ts=ts,default=0,match="startnextday")
    ei = qcutils.GetDateIndex(dt,ed,ts=ts,default=-1,match="endpreviousday")
    return si,ei

def match_masks(s1,s2):
    mask = numpy.ma.mask_or(s1.mask,s2.mask)
    s1m = numpy.ma.array(s1,mask=mask)
    s2m = numpy.ma.array(s2,mask=mask)
    return s1m,s2m

log = qcutils.startlog('yearscheck','../logfiles/yearscheck.log')

nFig = 0

#duname=qcio.get_filename_dialog()
duname = "../../Sites/DalyUncleared/Data/Processed/all/DalyUncleared_2008_to_2013_L3.nc"
ds = qcio.nc_read_series(duname)
ldt_all = ds.series["DateTime"]["Data"]
# get the site name, time step etc
site_name = ds.globalattributes["site_name"]
ts = ds.globalattributes["time_step"]
# get a list of years in the file
start_year = ldt_all[0].year
end_year = ldt_all[-1].year
years = range(start_year,end_year+1)
# do the daily ratios plot
log.info(' Doing the daily ratios plot ')
plot_title = site_name + ': ' + str(ldt_all[0]) + ' to ' + str(ldt_all[-1])
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(10.9,7.5))
plt.figtext(0.5,0.95,plot_title,horizontalalignment='center',size=16)
# SEB ratio
SEB_ax = plt.subplot(4,1,1)
EF_ax = plt.subplot(4,1,2)
BR_ax = plt.subplot(4,1,3)
WUE_ax = plt.subplot(4,1,4)
for yr in years:
    si,ei = get_sieiforwholeyear(ldt_all,yr,ts)
    # get the day/night mask
    day_mask,night_mask = get_daynightmask(ds,si=si,ei=ei)
    # get the data
    Fa,f = qcutils.GetSeriesasMA(ds,'Fa',si=si,ei=ei)
    Fe,f = qcutils.GetSeriesasMA(ds,'Fe',si=si,ei=ei)
    Fh,f = qcutils.GetSeriesasMA(ds,'Fh',si=si,ei=ei)
    Fc,f = qcutils.GetSeriesasMA(ds,'Fc',si=si,ei=ei)
    # match the missing data masks
    Fa_SEB,Fe_SEB = match_masks(Fa,Fe)
    Fa_SEB,Fh_SEB = match_masks(Fa_SEB,Fh)
    # get the daily average
    Fa_daily_avg,Fa_daily_num = get_dailyaverage(Fa_SEB,ts,mask=day_mask)
    Fe_daily_avg,Fe_daily_num = get_dailyaverage(Fe_SEB,ts,mask=day_mask)
    Fh_daily_avg,Fe_daily_num = get_dailyaverage(Fh_SEB,ts,mask=day_mask)
    # get the SEB ratio
    SEB_daily_avg = (Fe_daily_avg+Fh_daily_avg)/Fa_daily_avg
    SEB_daily_avg = numpy.ma.masked_where(Fa_daily_num<=5,SEB_daily_avg)
    SEB_daily_avg = numpy.ma.masked_where(SEB_daily_avg>2,SEB_daily_avg)
    SEB_doy = numpy.float64(range(len(SEB_daily_avg)))
    SEB_ax.plot(SEB_doy,SEB_daily_avg,"o",label=str(yr))
    SEB_ax.legend(loc="upper right",frameon=False,prop={'size':8})
    SEB_ax.set_ylabel('(Fh+Fe)/Fa')
    # evaporative fraction
    # match the missing data masks
    Fa_EF,Fe_EF = match_masks(Fa,Fe)
    # get the daily average
    Fa_daily_avg,Fa_daily_num = get_dailyaverage(Fa_EF,ts,mask=day_mask)
    Fe_daily_avg,Fe_daily_num = get_dailyaverage(Fe_EF,ts,mask=day_mask)
    # get the evaporative fraction
    EF_daily_avg = Fe_daily_avg/Fa_daily_avg
    EF_daily_avg = numpy.ma.masked_where(Fa_daily_num<=5,EF_daily_avg)
    EF_daily_avg = numpy.ma.masked_where(EF_daily_avg>1.5,EF_daily_avg)
    EF_doy = numpy.float64(range(len(EF_daily_avg)))
    EF_ax.plot(EF_doy,EF_daily_avg,"o",label=str(yr))
    EF_ax.legend(loc="upper right",frameon=False,prop={'size':8})
    EF_ax.set_ylabel('EF=Fe/Fa')
    # Bowen ratio
    # match the missing data masks
    Fh_BR,Fe_BR = match_masks(Fh,Fe)
    # get the daily average
    Fh_daily_avg,Fh_daily_num = get_dailyaverage(Fh_BR,ts,mask=day_mask)
    Fe_daily_avg,Fe_daily_num = get_dailyaverage(Fe_BR,ts,mask=day_mask)
    # get the Bowen ratio
    BR_daily_avg = Fh_daily_avg/Fe_daily_avg
    BR_daily_avg = numpy.ma.masked_where(Fh_daily_num<=5,BR_daily_avg)
    BR_daily_avg = numpy.ma.masked_where(BR_daily_avg>5,BR_daily_avg)
    BR_doy = numpy.float64(range(len(BR_daily_avg)))
    BR_ax.plot(BR_doy,BR_daily_avg,"o",label=str(yr))
    BR_ax.legend(loc="upper right",frameon=False,prop={'size':8})
    BR_ax.set_ylabel('BR=Fh/Fe')
    # ecosystem water use efficiency
    # match the missing data masks
    Fc_WUE,Fe_WUE = match_masks(Fc,Fe)
    # get the daily average
    Fc_daily_avg,Fc_daily_num = get_dailyaverage(Fc_WUE,ts,mask=day_mask)
    Fe_daily_avg,Fe_daily_num = get_dailyaverage(Fe_WUE,ts,mask=day_mask)
    # get the ecosystem water use efficiency
    WUE_daily_avg = Fc_daily_avg/Fe_daily_avg
    WUE_daily_avg = numpy.ma.masked_where(Fc_daily_num<=5,WUE_daily_avg)
    WUE_daily_avg = numpy.ma.masked_where((WUE_daily_avg>0.02)|(WUE_daily_avg<-0.2),WUE_daily_avg)
    WUE_doy = numpy.float64(range(len(WUE_daily_avg)))
    WUE_ax.plot(WUE_doy,WUE_daily_avg,"o",label=str(yr))
    WUE_ax.legend(loc="upper right",frameon=False,prop={'size':8})
    WUE_ax.set_ylabel('WUE=Fc/Fe')

figname='../plots/'+ds.globalattributes['site_name'].replace(' ','')+'_'+ds.globalattributes['nc_level']+'_YC_'+'DailyRatios.png'
fig.savefig(figname,format='png')
plt.show()
