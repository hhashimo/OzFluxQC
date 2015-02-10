import ast
import sys
sys.path.append('../scripts')
import constants as c
import datetime
import matplotlib.pyplot as plt
import numpy
import os
import qcgf
import qcio
import qcutils
from scipy.optimize import curve_fit
import scipy.ndimage as ndimage
import sys

def Reco_LloydTaylor_2param(T,rb,E0):
    t1 = 1/(c.Tref-c.T0)
    t2 = 1/(T-c.T0)
    return rb*numpy.exp(E0*(t1-t2))

def Reco_LloydTaylor_1param(T,rb):
    t1 = 1/(c.Tref-c.T0)
    t2 = 1/(T-c.T0)
    return rb*numpy.exp(c.E0_long*(t1-t2))

def NEE_Lasslop_pri(data,alpha,beta0,gamma):
    Fsd = numpy.array(data[0])
    D = numpy.array(data[1])
    beta = beta0/(1 + (D-c.D0)/c.D0)
    index = numpy.where(D<=c.D0)
    beta[index] = beta0
    NEE = (alpha*beta*Fsd/(alpha*Fsd + beta)) + gamma
    return NEE

def NEE_Lasslop_org(data,alpha,beta0,k,gamma):
    Fsd = numpy.array(data[0])
    D = numpy.array(data[1])
    beta = beta0*exp(-k*(D-c.D0))
    index = numpy.where(D<=c.D0)
    beta[index] = beta0
    NEE = (alpha*beta*Fsd/(alpha*Fsd + beta)) + gamma
    return NEE

def xyplot(x,y,sub=[1,1,1],regr=0,thru0=0,title=None,xlabel=None,ylabel=None,colours=None,fname=None):
    '''Generic XY scatter plot routine'''
    wspace = 0.0
    hspace = 0.0
    plt.subplot(sub[0],sub[1],sub[2])
    if colours is None:
        plt.plot(x,y,'b.')
    else:
        plt.scatter(x,y,c=colours)
    ax = plt.gca()
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
        wspace = 0.3
    if title is not None:
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

# get the control file
cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0: sys.exit()
# get the netCDF filename
ncfilename = qcio.get_infilenamefromcf(cf)
# get the Fsdand ustar thresholds
Fsd_lower = float(cf['Params']['Fsd_lower'])
Fsd_upper = float(cf['Params']['Fsd_upper'])
ustar_threshold = float(cf['Params']['ustar_threshold'])
# read the netCDF file
ds3 = qcio.nc_read_series(ncfilename)
if len(ds3.series.keys())==0: print time.strftime('%X')+' netCDF file '+ncfilename+' not found'; sys.exit()
SiteName = ds3.globalattributes['site_name']
DateTime = ds3.series['DateTime']['Data']
PlotTitle = SiteName + ': ' + str(DateTime[0]) + ' to ' + str(DateTime[-1])

# first figure is general plots of Fc as a function of ustar, Ts and Sws
# get the data as masked arrays
Fc,f,a=qcutils.GetSeriesasMA(ds3,'Fc')
Fc_units = ds3.series['Fc']['Attr']['units']
us,f,a=qcutils.GetSeriesasMA(ds3,'ustar')
Fsd,f,a=qcutils.GetSeriesasMA(ds3,'Fsd')
nFig = 1
fig = plt.figure(nFig,figsize=(8,8))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
# scatter plot of Fc versus ustar, night time
Fc_night = numpy.ma.masked_where(Fsd>Fsd_lower,Fc)
us_night = numpy.ma.masked_where(Fsd>Fsd_lower,us)
mask = numpy.ma.mask_or(Fc_night.mask,us_night.mask)
Fc_night = numpy.ma.array(Fc_night,mask=mask)         # apply the mask
us_night = numpy.ma.array(us_night,mask=mask)
xyplot(us_night,Fc_night,sub=[2,2,1],title="Night",xlabel='u* (m/s)',ylabel='Fc ('+Fc_units+')')
# scatter plot of Fc binned on ustar
bin_size = 0.05
# get the data as numpy array (not masked)
Fc,f,a = qcutils.GetSeries(ds3,'Fc')
us,f,a = qcutils.GetSeries(ds3,'ustar')
Fsd,f,a = qcutils.GetSeries(ds3,'Fsd')
Ts,f,a = qcutils.GetSeries(ds3,'Ts')
Sws,f,a = qcutils.GetSeries(ds3,'Sws')
month,f,a = qcutils.GetSeries(ds3,'Month')
# make sure all data is missing if one series is missing
#  - there are more more efficient ways to do this
for s1 in [Fc,us,Fsd,Ts,Sws]:
    i = numpy.where(abs(s1-float(-9999))<c.eps)[0]
    Fc[i] = float(-9999)
    us[i] = float(-9999)
    Fsd[i] = float(-9999)
    Ts[i] = float(-9999)
    Sws[i] = float(-9999)
# get night time data
i=numpy.where(Fsd>Fsd_lower)[0]
for s in [Fc,us,Fsd,Ts,Sws]:
    s[i] = -9999
labels = (us/bin_size) + 1
labels = labels.astype(int)
index = numpy.arange(1,numpy.max(labels)+1)
bin_upper = index*bin_size
Fc_bin_avg = ndimage.mean(Fc, labels=labels, index=index)
xyplot(bin_upper,Fc_bin_avg,sub=[2,2,2],title="Fc binned on u*",xlabel='u* (m/s)',ylabel='Fc ('+Fc_units+')')
# apply the ustar filter
i=numpy.where(us<ustar_threshold)[0]
for s in [Fc,us,Fsd,Ts,Sws]:
    s[i] = -9999
# scatter plot of Fc versus soil temperature coloured by soil moisture
idx=numpy.where((Fc!=float(-9999))&(Ts!=float(-9999))&(Sws!=float(-9999)))[0]
Ts_night = Ts[idx]
Fc_night = Fc[idx]
Sws_night = Sws[idx]
xyplot(Ts_night,Fc_night,sub=[2,2,3],title="Fc vs Ts (colour is Sws)",xlabel='Ts (C)',ylabel='Fc ('+Fc_units+')',colours=Sws_night)
# scatter plot of Fc versus soil moisture coloured by soil temperature
xyplot(Sws_night,Fc_night,sub=[2,2,4],title="Fc vs Sws (colour is Ts)",xlabel='Sws (frac)',ylabel='Fc ('+Fc_units+')',colours=Ts_night)
figname='../plots/'+ds3.globalattributes['site_name'].replace(' ','')+'_'+ds3.globalattributes['nc_level']+'_'+'FcvsustarTsSws.png'
fig.savefig(figname,format='png')

# second figure plots Fc as a fucntion of Ts for 4 ranges of soil moisture
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(8,8))
plt.figtext(0.5,0.95,PlotTitle,horizontalalignment='center',size=16)
# scatter plot of Fc binned on temperature for 4 ranges of soil moisture
Sws_min = numpy.min(Sws_night)
Sws_max = numpy.max(Sws_night)
nSwsranges = 4
Swsbin_size = (Sws_max - Sws_min)/float(nSwsranges)
# bin Fc on Ts for the 4 soil moisture bins
bin_means = []
bin_stds = []
mid_pts = []
leg_str = []
for j in range(nSwsranges):
    lower = Sws_min + j * Swsbin_size
    upper = Sws_min + (j+1) * Swsbin_size
    idx = numpy.where((Sws_night>lower)&(Sws_night<upper))[0]
    Ts_filt = Ts_night[idx]
    Fc_filt = Fc_night[idx]
    Ts_bins = numpy.linspace(Ts_filt.min(),Ts_filt.max(),7)
    digitized = numpy.digitize(Ts_filt,Ts_bins)
    bin_means.append(numpy.array([Fc_filt[digitized==i].mean() for i in range(1,len(Ts_bins))]))
    bin_stds.append(numpy.array([Fc_filt[digitized==i].std() for i in range(1,len(Ts_bins))]))
    firstmid = (Ts_bins[0]+Ts_bins[1])/2
    lastmid = (Ts_bins[-2]+Ts_bins[-1])/2
    mid_pts.append(numpy.linspace(firstmid,lastmid,len(Ts_bins)-1))
    ax = plt.subplot(3,2,j+1)
    plt.errorbar(mid_pts[j],bin_means[j],yerr=bin_stds[j],fmt='o')
    plt.xlabel('Ts (C)')
    plt.xlim(20,40)
    plt.ylim(0,10)
    plt.ylabel('Fc ('+Fc_units+')')
    plt.title('Fc vs Ts')
    leg_str.append('%.2f < Sws < %.2f'%(lower,upper))
    ax.text(0.5,0.9,leg_str[j],horizontalalignment='center',transform=ax.transAxes)

ax = plt.subplot(3,2,5)
for x,y,e in zip(mid_pts,bin_means,bin_stds):
    plt.errorbar(x,y,yerr=e,fmt='o')
plt.xlabel('Ts (C)')
plt.xlim(20,40)
plt.ylim(0,10)
plt.ylabel('Fc ('+Fc_units+')')
plt.title('Fc vs Ts')
plt.legend([s for s in leg_str], loc='center left', bbox_to_anchor=(1, 0.5))
plt.subplots_adjust(wspace=0.3,hspace=0.3)
figname='../plots/'+ds3.globalattributes['site_name'].replace(' ','')+'_'+ds3.globalattributes['nc_level']+'_'+'FcvsTs.png'
fig.savefig(figname,format='png')

# third figure is monthly plots of Fc vs Ts with Lloyd-Taylor fit
MnthList = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
# get the activation energy over the whole data set (Reichstein et al 2005)
# NB: this may not be a sensible thing to do for Australian ecosystems, allowing
#     E0 to vary as a function of Sws might make more sense
index = numpy.where((Fsd!=float(-9999))&(us!=float(-9999))&(Fc!=float(-9999))&(Ts!=float(-9999)))
Xfit_long = numpy.linspace(numpy.min(Ts[index]),numpy.max(Ts[index]), 100)
popt_long, pcov_long = curve_fit(Reco_LloydTaylor_2param, Ts[index], Fc[index],p0=(10,100))
Yfit_long = Reco_LloydTaylor_2param(Xfit_long, *popt_long)
c.E0_long = popt_long[1]
print popt_long[0],popt_long[1]
nFig = nFig + 1
fig = plt.figure(nFig,figsize=(7.5,10.9))
plt.figtext(0.5,0.95,'Respiration: Lloyd-Taylor',horizontalalignment='center',size=16)
j = 0
for i in [12,1,2,3,4,5,6,7,8,9,10,11]:
    j = j + 1
    index = numpy.where((Fsd<Fsd_lower)&(us>ustar_threshold)&(month==i)&(Fc!=float(-9999))&(Ts!=float(-9999)))[0]
#    index = numpy.where((Fsd<Fsd_lower)&(us>ustar_threshold)&(month==i)&(Fc!=float(-9999))&(Ts!=float(-9999))&(Fc>float(0)))[0]
    if len(index)>10:
        if j in [1,2,3,4,5,6,7,8,9]:
            xlabel = None
        else:
            xlabel = 'Ts (C)'
        if j in [2,3,5,6,8,9,11,12]:
            ylabel = None
        else:
            ylabel = 'Fc (umol/m2/s)'
        ax = plt.subplot(4,3,j)
        Xfit = numpy.linspace(numpy.min(Ts[index]),numpy.max(Ts[index]), 100)
        popt, pcov = curve_fit(Reco_LloydTaylor_2param, Ts[index], Fc[index],p0=(10,100))
        Yfit = Reco_LloydTaylor_2param(Xfit, *popt)
        textstr = '(%.1f,%.1f)'%(popt[0],popt[1])
        ax.text(0.1,0.9,textstr,horizontalalignment='left',transform=ax.transAxes,fontsize=10)
        popt, pcov = curve_fit(Reco_LloydTaylor_1param, Ts[index], Fc[index],p0=(10))
        Yfit_long = Reco_LloydTaylor_1param(Xfit, *popt)
        textstr = '(%.1f,%.1f)'%(popt[0],c.E0_long)
        ax.text(0.1,0.8,textstr,horizontalalignment='left',transform=ax.transAxes,fontsize=10)
        plt.scatter(Ts[index],Fc[index],c=Sws[index])
        plt.plot(Xfit,Yfit, 'r-')
        plt.plot(Xfit,Yfit_long, 'r--',linewidth=3.0)
        plt.title(MnthList[i-1])
        if xlabel is not None: plt.xlabel(xlabel)
        if ylabel is not None: plt.ylabel(ylabel)
        plt.subplots_adjust(wspace=0.275,hspace=0.275)
figname='../plots/'+ds3.globalattributes['site_name'].replace(' ','')+'_'+ds3.globalattributes['nc_level']+'_'+'FcLloydTaylor_monthly.png'
fig.savefig(figname,format='png')

## fourth figure is thumbnails of Fc vs Fsd for day time conditions with
## the Lasslop et al 2010 LUE fitted
#window_size = 10    # number of days in the data window
#nWindows = 36
#nrows = 6
#ncols = 6
#ldt = ds3.series['DateTime']['Data']
#StartDate = ldt[0]
#nFig = nFig + 1
##fig = plt.figure(nFig,figsize=(7.5,10.9))
#fig,ax = plt.subplots(6,6,sharex=True,sharey=True)
#fig.subplots_adjust(hspace=0)
#fig.subplots_adjust(wspace=0)
#plt.figtext(0.5,0.95,'Respiration: LUE intercept',horizontalalignment='center',size=16)
#for n in range(nWindows):
    #EndDate = StartDate + datetime.timedelta(days=window_size)
    #si = qcutils.GetDateIndex(ldt,str(StartDate),ts=30,default=0,match='exact')
    #ei = qcutils.GetDateIndex(ldt,str(EndDate),ts=30,default=-1,match='exact')
    #Fsd,f,a = qcutils.GetSeries(ds3,'Fsd',si=si,ei=ei)
    #Fc,f,a = qcutils.GetSeries(ds3,'Fc',si=si,ei=ei)
    #VPD,f,a = qcutils.GetSeries(ds3,'VPD',si=si,ei=ei)
    #us,f,a = qcutils.GetSeries(ds3,'ustar',si=si,ei=ei)
    ## get the day time data
    #index = numpy.where((Fsd>Fsd_lower)&(Fsd<Fsd_upper)&(Fsd!=float(-9999))&(Fc!=float(-9999))&(VPD!=float(-9999))&(us>0.10))[0]
    #Fsd_day = Fsd[index]
    #Fc_day = Fc[index]
    #VPD_day = VPD[index]
    #VPD_day = float(10) * VPD_day
    #ax = plt.subplot(6,6,n)
    #plt.plot(Fsd_day,Fc_day,'b.')
    #plt.xlim(0,1200)
    #plt.ylim(-50,20)
    #plt.gca().invert_yaxis()
    #ax.text(0.075,0.85,StartDate.strftime('%Y-%m-%d'),horizontalalignment='left',transform=ax.transAxes,fontsize=10)
    #if n not in [1,7,13,19,25,31]: plt.setp(ax.get_yticklabels(),visible=False)
    #if n not in [31,32,33,34,35,0]: plt.setp(ax.get_xticklabels(),visible=False)
    ## now do the curve fit
    #try:
        #data = []
        #data.append(Fsd_day)
        #data.append(VPD_day)
        #popt, pcov = curve_fit(NEE_Lasslop_org, data, Fc_day,p0=(-0.1,-100,0,10))
        #Yfit = NEE_Lasslop_org(data, *popt)
        #plt.plot(Fsd_day,Yfit,'r.')
        #textstr = '(%.3f,%.1f,%.1f)'%(popt[0],popt[1],popt[2])
        #ax.text(0.9,0.075,textstr,horizontalalignment='right',transform=ax.transAxes,fontsize=10)
    #except:
        #pass
    #StartDate = EndDate
#figname='../plots/'+ds3.globalattributes['site_name'].replace(' ','')+'_'+ds3.globalattributes['nc_level']+'_'+'LUE_10day.png'
#fig.savefig(figname,format='png')

## let's see how well SOFM/SOLO does
#driverlist = ast.literal_eval(cf['SOLO']['drivers'])
## get the output variable name
#outlabel = str(cf['SOLO']['output'])
## get the ustar threshold from the control file
#Fsd_lower = cf['Params'].as_float('Fsd_lower')
#ustar_threshold = cf['Params'].as_float('ustar_threshold')
## filter out the day time data
#index = numpy.where(ds3.series['Fsd']['Data']>Fsd_lower)[0]
#ds3.series['Fc']['Data'][index]=numpy.float64(-9999)
#ds3.series['Fc']['Flag'][index]=numpy.int32(18)
## apply the ustar filter
#index = numpy.where(ds3.series['ustar']['Data']<ustar_threshold)[0]
#ds3.series['Fc']['Data'][index]=numpy.float64(-9999)
#ds3.series['Fc']['Flag'][index]=numpy.int32(18)
## run SOFM
#nRecs = int(ds3.globalattributes['nc_nrecs'])
#cwd = os.getcwd()
#os.chdir('../')
#qcgf.gfSOLO_runsofm(cf,ds3,driverlist,'Fc',nRecs)
## run SOLO
## note that we pass in <series> but runsolo will use <series>_L3
#qcgf.gfSOLO_runsolo(cf,ds3,driverlist,'Fc',nRecs)
## run seqsolo and put the solo_modelled data into the ds series
## note that we pass in <series> but runseqsolo will use <series>_L3
#qcgf.gfSOLO_runseqsolo(cf,ds3,driverlist,'Fc',outlabel,nRecs)
## range check the seqsolo results?
## plot the results
##qcgf.gfSOLO_plotresults(cf,ds3,driverlist,'Fc',outlabel)

#os.chdir(cwd)

## fifth figure is monthly plots of Fc vs Ts for data and SOLO
#Fsd,f = qcutils.GetSeriesasMA(ds3,'Fsd')
#us,f = qcutils.GetSeriesasMA(ds3,'ustar')
#month,f = qcutils.GetSeriesasMA(ds3,'Month')
#Fc,f = qcutils.GetSeriesasMA(ds3,'Fc')
#Ts,f = qcutils.GetSeriesasMA(ds3,'Ts')
#Sws,f = qcutils.GetSeriesasMA(ds3,'Sws')
#Reco_SOLO,f = qcutils.GetSeriesasMA(ds3,'Reco_SOLO')
#MnthList = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#nFig = nFig + 1
#fig = plt.figure(nFig,figsize=(7.5,10.9))
#plt.figtext(0.5,0.95,'Respiration: SOLO',horizontalalignment='center',size=16)
#j = 0
#for i in [12,1,2,3,4,5,6,7,8,9,10,11]:
    #j = j + 1
    #index = numpy.where((Fsd<Fsd_lower)&(us>ustar_threshold)&(month==i)&(Fc!=float(-9999))&(Ts!=float(-9999))&(Fc>float(0)))[0]
    #if len(index)>10:
        #if j in [1,2,3,4,5,6,7,8,9]:
            #xlabel = None
        #else:
            #xlabel = 'Ts (C)'
        #if j in [2,3,5,6,8,9,11,12]:
            #ylabel = None
        #else:
            #ylabel = 'Fc (umol/m2/s)'
        #ax = plt.subplot(4,3,j)
        #plt.scatter(Ts[index],Fc[index],c=Sws[index])
        #plt.scatter(Ts[index],Reco_SOLO[index],color='r')
        #plt.title(MnthList[i-1])
        #if xlabel is not None: plt.xlabel(xlabel)
        #if ylabel is not None: plt.ylabel(ylabel)
        #plt.subplots_adjust(wspace=0.275,hspace=0.275)
#figname='../plots/'+ds3.globalattributes['site_name'].replace(' ','')+'_'+ds3.globalattributes['nc_level']+'_'+'Reco_SOLO_monthly.png'
#fig.savefig(figname,format='png')

plt.show()
