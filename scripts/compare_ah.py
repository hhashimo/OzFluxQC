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
fname = qcio.get_infilenamefromcf(cf)
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
# find the start index of the first whole day (time=00:30)
si = qcutils.GetDateIndex(DateTime,str(DateTime[0]),ts=ts,default=0,match='startnextday')
# find the end index of the last whole day (time=00:00)
ei = qcutils.GetDateIndex(DateTime,str(DateTime[-1]),ts=ts,default=-1,match='endpreviousday')
# clip the datetime series to a whole number of days
DateTime = DateTime[si:ei+1]
StartDate = DateTime[0]
EndDate = DateTime[-1]
print time.strftime('%X')+' Start date; '+str(StartDate)+' End date; '+str(EndDate)
Hdh = ds.series['Hdh']['Data'][si:ei+1]
Month = ds.series['Month']['Data'][si:ei+1]

nrecs = len(DateTime)
nperhr = int(float(60)/ts+0.5)
nperday = int(float(24)*nperhr+0.5)
ndays = nrecs/nperday
nrecs=ndays*nperday

Ah_7500_name = str(cf['Variables']['Ah_7500'])
Ah_HMP_name = str(cf['Variables']['Ah_HMP'])
# get local data series from the data structure
ah_7500_30min_1d,flag,attr = qcutils.GetSeriesasMA(ds,Ah_7500_name,si=si,ei=ei)
ah_HMP1_30min_1d,flag,attr = qcutils.GetSeriesasMA(ds,Ah_HMP_name,si=si,ei=ei)
month_30min_1d,flag,attr = qcutils.GetSeriesasMA(ds,'Month',si=si,ei=ei)
# mask data points unless both 7500 and HMP present
mask = numpy.ma.mask_or(ah_7500_30min_1d.mask,ah_HMP1_30min_1d.mask)
ah_7500_30min_1d = numpy.ma.array(ah_7500_30min_1d,mask=mask)
ah_HMP1_30min_1d = numpy.ma.array(ah_HMP1_30min_1d,mask=mask)
month_30min_1d = numpy.ma.array(month_30min_1d,mask=mask)
# reshape the 1D time series into 2D arrays
ah_7500_30min_2d = numpy.ma.reshape(ah_7500_30min_1d,[ndays,nperday])
ah_HMP1_30min_2d = numpy.ma.reshape(ah_HMP1_30min_1d,[ndays,nperday])
month_30min_2d = numpy.ma.reshape(month_30min_1d,[ndays,nperday])
# get the daily statistics
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
    x_nm = sm.add_constant(x_nm,prepend=False)
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

for i in range(0,ndays-1):
    x = ah_7500_30min_2d[i,:]
    y = ah_HMP1_30min_2d[i,:]
    x_nm = numpy.ma.compressed(x)
    y_nm = numpy.ma.compressed(y)
    nx = numpy.ma.count(x_nm)
    ny = numpy.ma.count(y_nm)
    r = numpy.ma.corrcoef(x_nm,y_nm)
    if (nx<min_n) or (r[0][1]<min_r):
        ah_7500_30min_2d[i,:].mask = True
        ah_HMP1_30min_2d[i,:].mask = True

class PointBrowser:
    def __init__(self):
        self.si = 0
        self.ei = 0
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
        #if self.ind_day is None: return
        if event.key=='n': self.new()
        if event.key=='f': self.forward()
        if event.key=='b': self.backward()
        if event.key=='q': self.quitprog()
        if event.key not in ('n', 'f', 'b', 'q'): return

    def new(self):
        print 'Creating new XY plot ...'
        # save the summary results from the last period
        if self.ei!=0:
            self.start_date.append(DT_daily[self.si])
            self.end_date.append(DT_daily[self.ei])
            self.slope.append(self.coefs[0])
            self.offset.append(self.coefs[1])
            self.correl.append(self.r[0][1])
            self.stdratio.append(self.sd)
            self.rangeratio.append(self.rr)
            self.si = self.ei
        # put up the new XY plot
        self.nfig += 1
        self.figxy = plt.figure(self.nfig,figsize=(5,4))
        self.figxy.subplots_adjust(bottom=0.15,left=0.15)
        self.axxy = self.figxy.add_subplot(111)
        self.axxy.set_xlabel('Ah_7500 (g/m3)')
        self.axxy.set_ylabel('Ah_HMP (g/m3)')
        plt.show()

    def forward(self):
        self.ei += 1
        self.update()

    def backward(self):
        self.ei += -1
        self.update()

    def update(self):
        self.ei = numpy.clip(self.ei,self.si,len(DT_daily)-1)
        x = ah_7500_30min_2d[self.ei,:]
        y = ah_HMP1_30min_2d[self.ei,:]
        if min([numpy.ma.count(x),numpy.ma.count(y)])<=0:
            print DT_daily[self.ei],'%g'%(numpy.ma.count(x))
        else:
            print DT_daily[self.ei],'%g %.3f %.3f %.3f'%(numpy.ma.count(x),numpy.ma.corrcoef(x,y)[0][1],
                                                         numpy.ma.polyfit(numpy.ma.copy(x),numpy.ma.copy(y),1)[0],
                                                         numpy.ma.polyfit(numpy.ma.copy(x),numpy.ma.copy(y),1)[1])
        x = ah_7500_30min_2d[self.si:self.ei+1,:]
        y = ah_HMP1_30min_2d[self.si:self.ei+1,:]
        x_nm = numpy.ma.compressed(x)
        y_nm = numpy.ma.compressed(y)
        if len(x_nm)!=0:
            self.r = numpy.corrcoef(x_nm,y_nm)
            self.sd = numpy.std(y_nm)/numpy.std(x_nm)
            self.rr = (numpy.max(y_nm)-numpy.min(y_nm))/(numpy.max(x_nm)-numpy.min(x_nm))
            resrlm = sm.RLM(y_nm,sm.add_constant(x_nm,prepend=False),M=sm.robust.norms.TukeyBiweight()).fit()
            self.coefs = resrlm.params
            m = self.coefs[0]; b = self.coefs[1]
            self.axxy.cla()
            self.axxy.plot(x_nm,y_nm,'b.')
            self.axxy.set_xlabel('Ah_7500 (g/m3)')
            self.axxy.set_ylabel('Ah_HMP (g/m3)')
            self.axxy.plot(x_nm,self.coefs[0]*x_nm+self.coefs[1],'r--',linewidth=3)
            eqnstr = 'y = %.3fx + %.3f, r = %.3f'%(self.coefs[0],self.coefs[1],self.r[0][1])
            self.axxy.text(0.5,0.875,eqnstr,fontsize=8,horizontalalignment='center',transform=self.axxy.transAxes)
            #print DT_daily[self.ei],'%g %.3f %.3f %.3f'%(numpy.ma.count(x),numpy.ma.corrcoef(x,y)[0][1],m,b)
        #else:
            #print DT_daily[self.ei],numpy.ma.count(x),numpy.ma.corrcoef(x,y)[0][1]
            #print str(DT_daily[self.ei])+'%g %.3f'%(numpy.ma.count(x),numpy.ma.corrcoef(x,y)[0][1])
            #m = numpy.ma.zeros(1); m.mask = True; m=m[0]
            #b = numpy.ma.zeros(1); b.mask = True; b=b[0]
        dtstr = str(DT_daily[self.si]) + ' to ' + str(DT_daily[self.ei])
        self.axxy.text(0.5,0.925,dtstr,fontsize=8,horizontalalignment='center',transform=self.axxy.transAxes)
        self.figxy.canvas.draw()

    def quitprog(self):
        self.start_date.append(DT_daily[self.si])
        self.end_date.append(DT_daily[self.ei])
        self.slope.append(self.coefs[0])
        self.offset.append(self.coefs[1])
        self.correl.append(self.r[0][1])
        self.stdratio.append(self.sd)
        self.rangeratio.append(self.rr)
        # print everything
        print '*** all results ***'
        for i in range(len(self.slope)):
            eqnstr = '%.3f, %.3f, %.3f, %.3f, %.3f'%(self.slope[i],self.offset[i],self.correl[i],self.stdratio[i],self.rangeratio[i])
            print self.start_date[i], self.end_date[i], eqnstr
        # print the linear fit for correcting Ah_7500_Av
        print '*** corrections for Ah_7500_Av ***'
        for i in range(len(self.slope)):
            eqnstr = '%.3f,%.3f'%(self.slope[i],self.offset[i])
            print str(i)+'='+'"['+"'"+self.start_date[i].strftime('%Y-%m-%d %H:%M')+"'"+','\
                      +"'"+self.end_date[i].strftime('%Y-%m-%d %H:%M')+"'"+","+eqnstr+']"'
        # print the ratio of the standard deviations for correcting the covariances
        print '*** corrections for covariances UxA, UyA and UzA'
        for i in range(len(self.slope)):
            eqnstr = '%.3f'%(self.stdratio[i])
            print str(i)+'='+'"['+"'"+self.start_date[i].strftime('%Y-%m-%d %H:%M')+"'"+','\
                      +"'"+self.end_date[i].strftime('%Y-%m-%d %H:%M')+"'"+","+eqnstr\
                      +',0.0]"'
        plt.close('all')

browser = PointBrowser()

figts.canvas.mpl_connect('key_press_event', browser.onpress)

plt.show()
