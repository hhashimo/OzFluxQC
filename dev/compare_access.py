import os
import sys
# check the scripts directory is present
if not os.path.exists("../scripts/"):
    print "compare_access: the scripts directory is missing"
    sys.exit()
# since the scripts directory is there, try importing the modules
sys.path.append('../scripts')
import numpy
import matplotlib.pyplot as plt
import meteorologicalfunctions as mf
import statsmodels.api as sm
import qcio
import qcutils

# load the control file
cf = qcio.load_controlfile(path='../controlfiles')
if len(cf)==0: sys.exit()
tow_name = cf["Files"]["tower_filename"]
acc_name = cf["Files"]["access_filename"]
# read the data series
ds_tow = qcio.nc_read_series(tow_name)
ds_acc = qcio.nc_read_series(acc_name)
# get the time step and the site name
ts = int(ds_tow.globalattributes["time_step"])
site_name = str(ds_tow.globalattributes["site_name"])
# get the start and end indices for the first and last whole days
dt_acc = ds_acc.series["DateTime"]["Data"]
si_acc = qcutils.GetDateIndex(dt_acc,str(dt_acc[0]),ts=ts,match="startnextday")
ei_acc = qcutils.GetDateIndex(dt_acc,str(dt_acc[0]),ts=ts,match="endpreviousday")
dt_acc = dt_acc[si_acc:ei_acc+1]

nrecs = len(dt_acc)
nperhr = int(float(60)/ts+0.5)
nperday = int(float(24)*nperhr+0.5)
ndays = nrecs/nperday
nrecs=ndays*nperday

dt_tow = ds_tow.series["DateTime"]["Data"]
if "Ws" not in ds_tow.series.keys() and "Ws_CSAT" in ds_tow.series.keys():
    Ws,f = qcutils.GetSeriesasMA(ds_tow,"Ws_CSAT")
    attr = qcutils.GetAttributeDictionary(ds_tow,"Ws_CSAT")
    qcutils.CreateSeries(ds_tow,"Ws",Ws,Flag=f,Attr=attr)
if "Wd" not in ds_tow.series.keys() and "Wd_CSAT" in ds_tow.series.keys():
    Wd,f = qcutils.GetSeriesasMA(ds_tow,"Wd_CSAT")
    attr = qcutils.GetAttributeDictionary(ds_tow,"Wd_CSAT")
    qcutils.CreateSeries(ds_tow,"Wd",Wd,Flag=f,Attr=attr)
if "q" not in ds_tow.series.keys():
    if "RH" not in ds_tow.series.keys():
        Ah,f = qcutils.GetSeriesasMA(ds_tow,"Ah")
        Ta,f = qcutils.GetSeriesasMA(ds_tow,"Ta")
        RH = mf.RHfromabsolutehumidity(Ah,Ta)
        attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='not defined')
        qcutils.CreateSeries(ds_tow,"RH",RH,FList=["Ta","Ah"],Attr=attr)
    RH,f = qcutils.GetSeriesasMA(ds_tow,"RH")
    Ta,f = qcutils.GetSeriesasMA(ds_tow,"Ta")
    ps,f = qcutils.GetSeriesasMA(ds_tow,"ps")
    q = mf.qfromrh(RH, Ta, ps)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific humidity',units='kg/kg',standard_name='specific_humidity')
    qcutils.CreateSeries(ds_tow,'q',q,FList=["Ta","ps","RH"],Attr=attr)

start_dt = str(dt_acc[0])
end_dt = str(dt_acc[-1])

si_tow = qcutils.GetDateIndex(dt_tow,start_dt,ts=ts)
ei_tow = qcutils.GetDateIndex(dt_tow,end_dt,ts=ts)
dt_tow = dt_tow[si_tow:ei_tow+1]

# margins
margin_bottom = 0.05
margin_top = 0.05
margin_left = 0.075
margin_right = 0.05
# plot heights, plot widths and spaces between plots
xy_height = 0.25
xy_width = 0.20
xyts_space = 0.05
xyxy_space = 0.05
ts_width = 0.9
# calculate bottom of the first time series and the height of the time series plots
ts_bottom = margin_bottom + xy_height + xyxy_space + xy_height + xyts_space
#ts_height = (1.0 - margin_top - ts_bottom)/float(nDrivers+1)
ts_height = (1.0 - margin_top - ts_bottom)
fig_num = 0
for label in ["Fsd","Fld","Fn","Ta","q","Ws","Ts","Sws","ps","ustar","Fh","Fe"]:
    print "compare_access: doing "+label
    fig_num = fig_num + 1
    # get the tower data
    data_tow,f = qcutils.GetSeriesasMA(ds_tow,label,si=si_tow,ei=ei_tow)
    # get the units
    attr = qcutils.GetAttributeDictionary(ds_tow,label)
    units = attr["units"]
    # get the correlation between the tower data and the 3x3 grid cells
    r_array=numpy.zeros((3,3))
    for i in range(0,3):
        for j in range(0,3):
            label_acc=label+"_"+str(i)+str(j)
            data_acc,f=qcutils.GetSeriesasMA(ds_acc,label_acc,si=si_acc,ei=ei_acc)
            r=numpy.ma.corrcoef(data_tow,data_acc)
            r_array[i,j]=r[0][1]
    # find the grid indices of the maximum correlation
    max_ij = numpy.unravel_index(r_array.argmax(), r_array.shape)
    # get data from the grid cell with highest correlation to tower data
    label_acc=label+"_"+str(max_ij[0])+str(max_ij[1])
    data_acc_all,f=qcutils.GetSeriesasMA(ds_acc,label_acc,si=si_acc,ei=ei_acc)
    # flatten the correlation array for plotting
    r_flat=r_array.flatten()
    # mask both series when either one is missing
    data_acc.mask = (data_acc_all.mask)|(data_tow.mask)
    data_tow.mask = (data_acc_all.mask)|(data_tow.mask)
    # get non-masked versions of the data
    data_acc_nm = numpy.ma.compressed(data_acc)
    data_tow_nm = numpy.ma.compressed(data_tow)
    # create the figure canvas
    fig=plt.figure(fig_num,figsize=(13,9))
    plt.figtext(0.5,0.96,site_name+' : Comparison of tower and ACCESS data for '+label,ha='center',size=16)
    # top row of XY plots
    # correlation coefficients
    rect1 = [0.10,margin_bottom+margin_bottom+xy_height,xy_width,xy_height]
    ax1 = plt.axes(rect1)
    ind=numpy.arange(len(r_flat))
    ax1.bar(ind,r_flat,0.35)
    ax1.set_ylabel('r')
    ax1.set_xlabel('Grid')
    # lagged correlation
    rect2 = [0.40,margin_bottom+margin_bottom+xy_height,xy_width,xy_height]
    ax2 = plt.axes(rect2)
    l=ax2.xcorr(data_tow_nm,data_acc_nm,maxlags=nperday)
    ax2.set_ylabel('r')
    ax2.set_xlabel('Lags')
    max_lag = ts*(numpy.argmax(l[1])-nperday)
    # bottom row of XY plots
    # scatter plot of 30 minute data
    rect3 = [0.10,margin_bottom,xy_width,xy_height]
    ax3 = plt.axes(rect3)
    ax3.plot(data_acc,data_tow,'b.')
    ax3.set_ylabel('Tower ('+units+')')
    ax3.set_xlabel('ACCESS ('+units+')')
    resrlm = sm.RLM(data_tow_nm,sm.add_constant(data_acc_nm),M=sm.robust.norms.TukeyBiweight()).fit()
    m_rlm = resrlm.params[0]; b_rlm = resrlm.params[1]
    eqnstr = 'y = %.3fx + %.3f'%(m_rlm,b_rlm)
    ax3.plot(data_acc_nm,resrlm.fittedvalues,'r--',linewidth=3)
    ax3.text(0.5,0.915,eqnstr,fontsize=8,horizontalalignment='center',transform=ax3.transAxes,color='red')
    resols = sm.OLS(data_tow_nm,sm.add_constant(data_acc_nm)).fit()
    m_ols = resols.params[0]; b_ols = resols.params[1]
    eqnstr = 'y = %.3fx + %.3f'%(m_ols,b_ols)
    ax3.plot(data_acc_nm,resrlm.fittedvalues,'g--',linewidth=3)
    ax3.text(0.5,0.85,eqnstr,fontsize=8,horizontalalignment='center',transform=ax3.transAxes,color='green')
    ax3.text(0.6,0.075,'30 minutes',fontsize=10,horizontalalignment='left',transform=ax3.transAxes)
    # scatter plot of daily averages
    rect4 = [0.40,margin_bottom,xy_width,xy_height]
    ax4 = plt.axes(rect4)
    data_tow_2d = numpy.ma.reshape(data_tow,[ndays,nperday])
    data_tow_daily_avg = numpy.ma.average(data_tow_2d,axis=1)
    data_acc_2d = numpy.ma.reshape(data_acc,[ndays,nperday])
    data_acc_daily_avg = numpy.ma.average(data_acc_2d,axis=1)
    ax4.plot(data_acc_daily_avg,data_tow_daily_avg,'b.')
    ax4.set_ylabel('Tower ('+units+')')
    ax4.set_xlabel('ACCESS ('+units+')')
    data_acc_daily_avg.mask = (data_acc_daily_avg.mask)|(data_tow_daily_avg.mask)
    data_tow_daily_avg.mask = (data_acc_daily_avg.mask)|(data_tow_daily_avg.mask)
    data_acc_daily_avg_nm = numpy.ma.compressed(data_acc_daily_avg)
    data_tow_daily_avg_nm = numpy.ma.compressed(data_tow_daily_avg)
    resrlm = sm.RLM(data_tow_daily_avg_nm,sm.add_constant(data_acc_daily_avg_nm),M=sm.robust.norms.TukeyBiweight()).fit()
    eqnstr = 'y = %.3fx + %.3f'%(resrlm.params[0],resrlm.params[1])
    ax4.plot(data_acc_daily_avg_nm,resrlm.fittedvalues,'r--',linewidth=3)
    ax4.text(0.5,0.915,eqnstr,fontsize=8,horizontalalignment='center',transform=ax4.transAxes,color='red')
    resrlm = sm.OLS(data_tow_daily_avg_nm,sm.add_constant(data_acc_daily_avg_nm)).fit()
    eqnstr = 'y = %.3fx + %.3f'%(resrlm.params[0],resrlm.params[1])
    ax4.plot(data_acc_daily_avg_nm,resrlm.fittedvalues,'g--',linewidth=3)
    ax4.text(0.5,0.85,eqnstr,fontsize=8,horizontalalignment='center',transform=ax4.transAxes,color='green')
    ax4.text(0.6,0.075,'Daily average',fontsize=10,horizontalalignment='left',transform=ax4.transAxes)
    # diurnal average plot
    rect5 = [0.70,margin_bottom,xy_width,xy_height]
    ax5 = plt.axes(rect5)
    data_tow_hourly_avg = numpy.ma.average(data_tow_2d,axis=0)
    data_acc_hourly_avg = numpy.ma.average(data_acc_2d,axis=0)
    data_acc_hourly_rlm = numpy.ma.average(m_rlm*data_acc_2d+b_rlm,axis=0)
    data_acc_hourly_ols = numpy.ma.average(m_ols*data_acc_2d+b_ols,axis=0)
    ind = numpy.arange(len(data_tow_hourly_avg))*float(ts)/float(60)
    ax5.plot(ind,data_tow_hourly_avg,'ro',label='Tower')
    ax5.plot(ind,data_acc_hourly_avg,'b-',label='ACCESS-A')
    ax5.plot(ind,data_acc_hourly_rlm,'r-',label='ACCESS-A (RLM)')
    ax5.plot(ind,data_acc_hourly_ols,'g-',label='ACCESS-A (OLS)')
    ax5.set_ylabel(label+' ('+units+')')
    ax5.set_xlim(0,24)
    ax5.xaxis.set_ticks([0,6,12,18,24])
    ax5.set_xlabel('Hour')
    ax5.legend(loc='upper right',frameon=False,prop={'size':8})
    # time series
    rect_ts = [margin_left,ts_bottom,ts_width,ts_height]
    axes_ts = plt.axes(rect_ts)
    axes_ts.plot(dt_tow,data_tow,'ro',label="Tower")
    axes_ts.plot(dt_acc,data_acc_all,'b-',label="ACCESS-A")
    axes_ts.plot(dt_acc,m_rlm*data_acc_all+b_rlm,'r-',label="ACCESS-A (RLM)")
    axes_ts.plot(dt_acc,m_ols*data_acc_all+b_ols,'g-',label="ACCESS-A (OLS)")
    axes_ts.set_ylabel(label+' ('+units+')')
    axes_ts.legend(loc='upper right',frameon=False,prop={'size':8})
    # now get some statistics
    numpoints = numpy.ma.count(data_tow)
    diff = data_tow - data_acc
    bias = numpy.ma.average(diff)
    rmse = numpy.ma.sqrt(numpy.ma.average(diff*diff))
    var_tow = numpy.ma.var(data_tow)
    var_acc = numpy.ma.var(data_acc)
    text_left = 0.70
    num_left = 0.80
    row_bottom = 0.375
    row_space = 0.030
    i = 0
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'Lag (minutes)')
    plt.figtext(num_left,row_posn,'%.4g'%(max_lag))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'Var (ACCESS)')
    plt.figtext(num_left,row_posn,'%.4g'%(var_acc))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'Var (tower)')
    plt.figtext(num_left,row_posn,'%.4g'%(var_tow))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'RMSE')
    plt.figtext(num_left,row_posn,'%.4g'%(rmse))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'Bias')
    plt.figtext(num_left,row_posn,'%.4g'%(bias))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'r')
    plt.figtext(num_left,row_posn,'%.4g at (%.2g,%.2g)'%(numpy.ma.maximum(r_array),max_ij[0],max_ij[1]))
    i = i + 1
    row_posn = row_bottom + i*row_space
    plt.figtext(text_left,row_posn,'No. points')
    plt.figtext(num_left,row_posn,str(numpoints))
    # save a hard copy of the plot
    figname='../plots/'+site_name.replace(' ','')+'_'+'ACCESS'+'_'+label+'.png'
    fig.savefig(figname,format='png')
    # show the plot

plt.show()
