# coding: utf-8
duname=qcio.get_filename_dialog()
ds=qcio.nc_read_series(duname)
ldt=ds.series["DateTime"]["Data"]
ts=ds.globalattributes["time_step"]
si_2009=qcutils.GetDateIndex(ldt,"2009-01-01 00:30",ts=ts,match="startnextday")
ei_2009=qcutils.GetDateIndex(ldt,"2010-01-01 00:00",ts=ts,match="endpreviousday")
Fa_2009,f=qcutils.GetSeriesasMA(ds,'Fa',si=si_2009,ei=ei_2009)
Fe_2009,f=qcutils.GetSeriesasMA(ds,'Fe',si=si_2009,ei=ei_2009)
Fsd_2009,f=qcutils.GetSeriesasMA(ds,'Fsd',si=si_2009,ei=ei_2009)
Fsd_syn,f=qcutils.GetSeriesasMA(ds,'Fsd_syn',si=si_2009,ei=ei_2009)
index=numpy.ma.where(Fsd_2009.mask==True)[0]
Fsd_2009[index]=Fsd_syn[index]
night_mask=(Fsd_2009<10)
day_mask=(Fsd_2009>=10)
Fsd_2009,f=qcutils.GetSeriesasMA(ds,'Fsd',si=si_2009,ei=ei_2009)
ldt_2009=ldt[si_2009:ei_2009+1]
nDays = float(len(ldt_2009))/ntsInDay
Fa_2009_daily=Fa_2009.reshape(nDays,ntsInDay)
Fe_2009_daily=Fe_2009.reshape(nDays,ntsInDay)
dm_daily = day_mask.reshape(nDays,ntsInDay)
Fa_2009_day=numpy.ma.masked_where(dm_daily==False,Fa_2009_daily)
Fe_2009_day=numpy.ma.masked_where(dm_daily==False,Fe_2009_daily)
mask=numpy.ma.mask_or(Fa_2009_day.mask,Fe_2009_day.mask)
Fa_2009_day=numpy.ma.array(Fa_2009_day,mask=mask)
Fe_2009_day=numpy.ma.array(Fe_2009_day,mask=mask)
Fa_2009_day_avg=numpy.ma.average(Fa_2009_day,axis=1)
Fe_2009_day_avg=numpy.ma.average(Fe_2009_day,axis=1)
EF_2009_day_avg=Fe_2009_day_avg/Fa_2009_day_avg
EF_2009_day_num = numpy.ma.count(Fe_2009_day,axis=1)
EF_2008_day_avg = numpy.ma.masked_where(EF_2008_day_num<=5,EF_2008_day_avg)
index = numpy.ma.where(EF_2008_day_avg.mask==True)
index = numpy.ma.where(EF_2009_day_avg.mask==True)
EF_2009_day_avg = numpy.ma.masked_where(EF_2009_day_num<=5,EF_2009_day_avg)
index = numpy.ma.where(EF_2009_day_avg.mask==True)
EF_2009_day_num[index] = 0
plot(EF_2009_day_avg,'r+')