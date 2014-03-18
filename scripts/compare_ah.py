import constants as c
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy
import os
import qcio
import qcplot
import qcutils
import statsmodels.api as sm
import sys
import time

nfig = 0
plotwidth = 10.9
plotheight = 7.5
# load the control file
cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0: sys.exit()
min_n = int(cf["General"]["minimum_number"])
min_r = float(cf["General"]["minimum_correlation"])
# get the input file name
fname = qcio.get_infilename_from_cf(cf)
if not os.path.exists(fname):
    print " compare_ah: Input netCDF file "+fname+" doesn't exist"
    sys.exit()
# read the input file and return the data structure
ds = qcio.nc_read_series(fname)
if len(ds.series.keys())==0: print time.strftime('%X')+' netCDF file '+fname+' not found'; sys.exit()
# get the site name
SiteName = ds.globalattributes['site_name']
# get the time step
ts = int(ds.globalattributes['time_step'])
# get the datetime series
DateTime = ds.series['DateTime']['Data']
# get the initial start and end dates
StartDate = str(DateTime[0])
EndDate = str(DateTime[-1])
# find the start index of the first whole day (time=00:30)
si = qcutils.GetDateIndex(DateTime,StartDate,ts=ts,default=0,match='startnextday')
# find the end index of the last whole day (time=00:00)
ei = qcutils.GetDateIndex(DateTime,EndDate,ts=ts,default=-1,match='endpreviousday')
DateTime = DateTime[si:ei+1]
print time.strftime('%X')+' Start date; '+str(DateTime[0])+' End date; '+str(DateTime[-1])
Hdh = ds.series['Hdh']['Data'][si:ei+1]
Month = ds.series['Month']['Data'][si:ei+1]

nrecs = len(DateTime)
nperhr = int(float(60)/ts+0.5)
nperday = int(float(24)*nperhr+0.5)
ndays = nrecs/nperday
nrecs=ndays*nperday

Ah_7500_name = str(cf['Variables']['Ah_7500'])
Ah_HMP_name = str(cf['Variables']['Ah_HMP'])
ah_7500_30min_1d,flag = qcutils.GetSeriesasMA(ds,Ah_7500_name,si=si,ei=ei)
ah_HMP1_30min_1d,flag = qcutils.GetSeriesasMA(ds,Ah_HMP_name,si=si,ei=ei)
month_30min_1d,flag = qcutils.GetSeriesasMA(ds,'Month',si=si,ei=ei)
ah_7500_30min_2d = numpy.ma.reshape(ah_7500_30min_1d,[ndays,nperday])
ah_HMP1_30min_2d = numpy.ma.reshape(ah_HMP1_30min_1d,[ndays,nperday])
month_30min_2d = numpy.ma.reshape(month_30min_1d,[ndays,nperday])

mask = numpy.ma.mask_or(ah_7500_30min_2d.mask,ah_HMP1_30min_2d.mask)  # mask based on dependencies, set all to missing if any missing
ah_7500_30min_2d = numpy.ma.array(ah_7500_30min_2d,mask=mask)         # apply the mask
ah_HMP1_30min_2d = numpy.ma.array(ah_HMP1_30min_2d,mask=mask)
month_30min_2d = numpy.ma.array(month_30min_2d,mask=mask)

month_daily_avg = numpy.ma.average(month_30min_2d,axis=1)
ah_7500_daily_avg = numpy.ma.average(ah_7500_30min_2d,axis=1)
ah_HMP1_daily_avg = numpy.ma.average(ah_HMP1_30min_2d,axis=1)
ah_7500_daily_std = numpy.ma.std(ah_7500_30min_2d,axis=1)
ah_HMP1_daily_std = numpy.ma.std(ah_HMP1_30min_2d,axis=1)
ah_7500_daily_max = numpy.ma.max(ah_7500_30min_2d,axis=1)
ah_HMP1_daily_max = numpy.ma.max(ah_HMP1_30min_2d,axis=1)
ah_7500_daily_min = numpy.ma.min(ah_7500_30min_2d,axis=1)
ah_HMP1_daily_min = numpy.ma.min(ah_HMP1_30min_2d,axis=1)

ah_avgdiff_daily = ah_7500_daily_avg - ah_HMP1_daily_avg
ah_stdratio_daily = ah_HMP1_daily_std/ah_7500_daily_std
ah_7500range_daily = ah_7500_daily_max - ah_7500_daily_min
ah_HMP1range_daily = ah_HMP1_daily_max - ah_HMP1_daily_min
ah_rangeratio_daily = (ah_HMP1_daily_max - ah_HMP1_daily_min)/(ah_7500_daily_max - ah_7500_daily_min)
DT_daily = DateTime[0:nrecs:nperday]
# time series plot of daily averaged absolute humidities and differencves
nfig = nfig + 1
fig = plt.figure(nfig,figsize=(plotwidth,plotheight))
plt.figtext(0.5,0.95,SiteName,horizontalalignment='center',size=16)
qcplot.tsplot(DT_daily,ah_7500_daily_avg,sub=[3,1,1],ylabel='Ah_7500')
qcplot.tsplot(DT_daily,ah_HMP1_daily_avg,sub=[3,1,2],ylabel='Ah_HMP_01')
qcplot.tsplot(DT_daily,ah_avgdiff_daily,sub=[3,1,3],ylabel='7500-HMP')
# scatter plots of absolute humidities by month
MnthList = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
nfig = nfig + 1
fig = plt.figure(nfig,figsize=(plotwidth,plotheight))
plt.figtext(0.5,0.95,SiteName,horizontalalignment='center',size=16)
j = 0
for i in [1,2,3,4,5,6,7,8,9,10,11,12]:
    j = j + 1
    index = numpy.where(month_30min_1d==i)[0]
    if len(index)!=0:
        y = ah_HMP1_30min_1d[index]
        x = ah_7500_30min_1d[index]
        if j in [1,2,3,4,5,6,7,8,9]:
            xlabel = None
        else:
            xlabel = '7500 (g/m3)'
        if j in [2,3,5,6,8,9,11,12]:
            ylabel = None
        else:
            ylabel = 'HMP (g/m3)'
        qcplot.xyplot(x,y,sub=[4,3,j],regr=2,title=MnthList[i-1],xlabel=xlabel,ylabel=ylabel)
plt.tight_layout()
# daily regressions
slope = numpy.ones(ndays)
offset = numpy.zeros(ndays)
correl = numpy.ones(ndays)
number = numpy.zeros(ndays)
for i in range(0,ndays-1):
    x = ah_7500_30min_2d[i,:]
    y = ah_HMP1_30min_2d[i,:]
    x_nm = numpy.ma.compressed(x)
    x_nm = sm.add_constant(x_nm)
    y_nm = numpy.ma.compressed(y)
    if len(y_nm)>1:
        resrlm = sm.RLM(y_nm,x_nm,M=sm.robust.norms.TukeyBiweight()).fit()
        coefs = resrlm.params
        r = numpy.ma.corrcoef(x,y)
        number[i] = numpy.ma.count(x)
        slope[i] = coefs[0]
        offset[i] = coefs[1]
        correl[i] = r[0][1]
correl2 = numpy.ma.masked_where((correl<min_r)|(number<min_n),correl)
number2 = numpy.ma.masked_where((correl<min_r)|(number<min_n),number)
slope2 = numpy.ma.masked_where((correl<min_r)|(number<min_n),slope)
offset2 = numpy.ma.masked_where((correl<min_r)|(number<min_n),offset)
sdratio2 = numpy.ma.masked_where((correl<min_r)|(number<min_n),ah_stdratio_daily)
nfig = nfig + 1
figts = plt.figure(nfig,figsize=(plotwidth,plotheight))
plt.figtext(0.5,0.95,SiteName,horizontalalignment='center',size=16)
qcplot.tsplot(DT_daily,correl2,sub=[5,1,1],ylabel='Correl',colours=number)
qcplot.tsplot(DT_daily,number2,sub=[5,1,2],ylabel='Number',colours=correl)
qcplot.tsplot(DT_daily,slope2,sub=[5,1,3],ylabel='Slope',colours=correl)
qcplot.tsplot(DT_daily,offset2,sub=[5,1,4],ylabel='Offset',colours=correl)
qcplot.tsplot(DT_daily,sdratio2,sub=[5,1,5],ylabel='Sd(HMP)/Sd(7500)',colours=correl)

class PointBrowser:
    def __init__(self):
        self.ind_day = 0
        self.start_ind_day = 0
        self.ind_30min = 0
        self.start_ind_30min = 0
        self.end_ind = 0
        self.nfig = nfig
        self.slope = []
        self.offset = []
        self.correl = []
        self.start_date = []
        self.end_date = []
        self.stdratio = []
        self.rangeratio = []
        self.last_index = []

    def onpress(self, event):
        if self.ind_day is None: return
        if event.key=='n': self.new()
        if event.key=='f': self.forward()
        if event.key=='b': self.backward()
        if event.key=='q': self.quitprog()
        if event.key not in ('n', 'f', 'b', 'q'): return

    def new(self):
        print 'Creating new XY plot ...'
        # save the summary results from the last period
        if self.ind_day!=0:
            self.start_date.append(DT_daily[self.start_ind_day])
            self.end_date.append(DT_daily[self.ind_day])
            self.slope.append(self.coefs[0])
            self.offset.append(self.coefs[1])
            self.correl.append(self.r[0][1])
            self.stdratio.append(numpy.ma.std(self.y_obs)/numpy.ma.std(self.x_obs))
            self.rangeratio.append((numpy.ma.max(self.y_obs)-numpy.ma.min(self.y_obs))/(numpy.ma.max(self.x_obs)-numpy.ma.min(self.x_obs)))
            self.start_ind_day = self.ind_day
            self.ind_30min = DateTime.index(DT_daily[self.ind_day])
            self.start_ind_30min = self.ind_30min
        # initialise some arrays
        self.x_obs = numpy.ma.array([])
        self.y_obs = numpy.ma.array([])
        self.last_index = []
        # put up the new XY plot
        self.nfig += 1
        self.figxy = plt.figure(self.nfig,figsize=(5,4))
        self.figxy.subplots_adjust(bottom=0.15,left=0.15)
        self.axxy = self.figxy.add_subplot(111)
        self.axxy.set_xlabel('Ah_7500 (g/m3)')
        self.axxy.set_ylabel('Ah_HMP (g/m3)')
        plt.show()

    def forward(self):
        self.ind_day += 1
        self.ind_day = numpy.clip(self.ind_day,0,len(DT_daily)-1)
        self.ind_30min = DateTime.index(DT_daily[self.ind_day])
        self.ind_start = DateTime.index(DT_daily[self.ind_day-1])
        x = ah_7500_30min_1d[self.ind_start:self.ind_30min]
        y = ah_HMP1_30min_1d[self.ind_start:self.ind_30min]
        mask = (x.mask)|(y.mask)
        x.mask = mask
        y.mask = mask
        print DateTime[self.ind_start],DateTime[self.ind_30min],numpy.ma.count(x),numpy.ma.corrcoef(x,y)[0][1]
        if numpy.ma.count(x)>min_n:
            if (numpy.ma.corrcoef(x,y)[0][1]>min_r):
                x_nm = numpy.ma.compressed(x)
                y_nm = numpy.ma.compressed(y)
                self.x_obs = numpy.ma.concatenate([self.x_obs,x_nm])
                self.y_obs = numpy.ma.concatenate([self.y_obs,y_nm])
                self.r = numpy.ma.corrcoef(self.x_obs,self.y_obs)
                self.last_index.append(len(self.x_obs))
                resrlm = sm.RLM(self.y_obs,sm.add_constant(self.x_obs),M=sm.robust.norms.TukeyBiweight()).fit()
                self.coefs = resrlm.params
        self.update()

    def backward(self):
        self.ind_day += -1
        self.ind_day = numpy.clip(self.ind_day,0,len(DT_daily)-1)
        self.ind_30min = DateTime.index(DT_daily[self.ind_day])
        self.ind_start = DateTime.index(DT_daily[self.ind_day-1])
        self.y_obs = self.y_obs[:self.last_index[-1]]
        self.x_obs = self.x_obs[:self.last_index[-1]]
        self.last_index = self.last_index[:-1]
        resrlm = sm.RLM(self.y_obs,sm.add_constant(self.x_obs),M=sm.robust.norms.TukeyBiweight()).fit()
        self.coefs = resrlm.params
        self.update()

    def update(self):
        self.axxy.cla()
        self.axxy.plot(self.x_obs,self.y_obs,'b.')
        self.axxy.set_xlabel('Ah_7500 (g/m3)')
        self.axxy.set_ylabel('Ah_HMP (g/m3)')
        if len(self.x_obs)!=0:
            self.axxy.plot(self.x_obs,self.coefs[0]*self.x_obs+self.coefs[1],'r--',linewidth=3)
            eqnstr = 'y = %.3fx + %.3f, r = %.3f'%(self.coefs[0],self.coefs[1],self.r[0][1])
            self.axxy.text(0.5,0.875,eqnstr,fontsize=8,horizontalalignment='center',transform=self.axxy.transAxes)
        dtstr = str(DT_daily[self.start_ind_day]) + ' to ' + str(DT_daily[self.ind_day])
        self.axxy.text(0.5,0.925,dtstr,fontsize=8,horizontalalignment='center',transform=self.axxy.transAxes)
        self.figxy.canvas.draw()

    def quitprog(self):
        self.start_date.append(DT_daily[self.start_ind_day])
        self.end_date.append(DT_daily[self.ind_day])
        self.slope.append(self.coefs[0])
        self.offset.append(self.coefs[1])
        self.correl.append(self.r[0][1])
        self.stdratio.append(numpy.ma.average(ah_stdratio_daily[self.start_ind_day:self.ind_day]))
        self.rangeratio.append(numpy.ma.average(ah_rangeratio_daily[self.start_ind_day:self.ind_day]))
        for i in range(len(self.slope)):
            eqnstr = '%.3f, %.3f, %.3f, %.3f, %.3f'%(self.slope[i],self.offset[i],self.correl[i],self.stdratio[i],self.rangeratio[i])
            print self.start_date[i], self.end_date[i], eqnstr
        plt.close('all')

browser = PointBrowser()

figts.canvas.mpl_connect('key_press_event', browser.onpress)

plt.show()
