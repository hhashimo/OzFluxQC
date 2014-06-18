import matplotlib.pyplot as plt
import numpy
import qcio
import qcutils

#Things to do:
#1) include storage term in Fc
#2) implement Papale et al (2006) spike removal for Fc
#3) test with synthetic data
#   - done, correctly picks u*_th
#4) implement Eqn 1a (diagnostic version)
#5) implement QA as per Barr et al (2013)
#   - check F values for significance
#   - divide into "deficit" and "excess" classes, discard the less populous
#   - reject site-year if less than 4000 (?) or 20% change points accepted
#   - use values of a1/a2 as diagnostic (Fig 3)
#6) implement annual average of accepted u*_th
#7) implement bootstrap process
#   - Papale et al (2006) use 100 bootstraps, Barr et al (2013)
#     seem to use 1000.
#8) implement robust least-squares from scikits
#9) draw Fmax lines on Fc vs ustar plots
#   - done
#10) implement "seasonal" output
#    - start_date,end_date,u*@Fmax_1,u*max_1,u*@Fmax_2,u*max_2,u*@Fmax_3,u*max_3,u*@Fmax_4,u*max_4,Fmax_1,Fmax_2,Fmax_3,Fmax_4,Ta_min,Ta_12,Ta_23,Ta_34,Ta_max
#11) implement bootstrap output
#    - bootstrap_number,u*@Fmax_1,u*max_1,u*@Fmax_2,u*max_2,u*@Fmax_3,u*max_3,u*@Fmax_4,u*max_4,Fmax_1,Fmax_2,Fmax_3,Fmax_4,Ta_min,Ta_12,Ta_23,Ta_34,Ta_max
#
# Questions to answer:
#1) what is beta_0?
#2) how do you trap cases where Fc continues to increase with us (F above Table 2 values)?

def SSE_2b(fc):
    SSE_2b = {}
    SSE_2b['beta_0'] = numpy.ma.average(fc)
    resid = fc - SSE_2b['beta_0']
    SSE_2b['values'] = numpy.ma.sum(resid*resid)
    #print 'SSE_2b= ',SSE_2b
    return SSE_2b

def SSE_1b(fc,us,c):
    SSE_1b = {}
    # local pointers to data below and above the current guess at the change point
    x_below = us[0:c+1]
    x_above = us[c+1:]
    y_below = fc[0:c+1]
    y_above = fc[c+1:]
    # least-squares best fit to below and above points
    #print 'c is ',c
    #print 'x_below=',x_below
    #print 'y_below=',y_below
    #print 'x_above=',x_above
    #print 'y_above=',y_above
    b = numpy.ma.polyfit(x_below,y_below,1)
    SSE_1b['b'] = b
    #print 'b_below=',b_below
    #b_above = numpy.ma.polyfit(x_above,y_above,1)
    #print 'b_above=',b_above
    # residuals for below and above groups of points
    r_below = y_below - (b[1]+b[0]*x_below)
    r_above = y_above - (b[1]+b[0]*x_below[-1])
    # calculate the SSE 1b term
    SSE_1b_below = numpy.ma.sum(r_below*r_below)
    SSE_1b_above = numpy.ma.sum(r_above*r_above)
    SSE_1b['values'] = SSE_1b_below + SSE_1b_above
    #print 'SSE_1b_below=',SSE_1b_below,' SSE_1b_above=',SSE_1b_above,' SSE_1b=',SSE_1b
    return SSE_1b

def Fc_1b(fc,us,c):
    Fc_1b = {}
    n = numpy.ma.count(us)
    sse_2b = SSE_2b(fc)
    Fc_1b['SSE_2b'] = sse_2b
    sse_1b = SSE_1b(fc,us,c)
    Fc_1b['SSE_1b'] = sse_1b
    Fc_1b['values'] = (sse_2b['values']-sse_1b['values'])/(sse_2b['values']/(n-2))
    #print 'Fc_1b=',Fc_1b
    return Fc_1b

def get_Fscore_1b(fc,us):
    Fscore = {}
    Fscore['values'] = numpy.ma.zeros(len(us))
    Fscore['b'] = {}
    for i in range(1,len(us)):
        fc_1b = Fc_1b(fc,us,i)
        Fscore['values'][i] = fc_1b['values']
        Fscore['b'][str(i)] = fc_1b['SSE_1b']['b']
    Fscore['Fmax'] = float(numpy.ma.maximum(Fscore['values']))
    Fscore['iatFmax'] = int(numpy.ma.where(Fscore['values']==Fscore['Fmax'])[0])
    Fscore['usatFmax'] = float(us[Fscore['iatFmax']])
    return Fscore

# initialise some constants
nFig = 0
nTemp = 4
nustarbins = 50
nustarpointsperbin = 5
npointsperseason = nustarpointsperbin*nustarbins*nTemp
npointsperjump = npointsperseason/2

ncname = '../../Sites/HowardSprings/Data/Processed/2011/HowardSprings_2011_L3.nc'
# read the netCDF file
ds = qcio.nc_read_series(ncname)
nRecs = int(ds.globalattributes['nc_nrecs'])
# get the data from the data structure
Fsd,f = qcutils.GetSeriesasMA(ds,'Fsd')
Ta,f = qcutils.GetSeriesasMA(ds,'Ta')
ustar,f = qcutils.GetSeriesasMA(ds,'ustar')
# uncomment following line to use real data
Fc,f = qcutils.GetSeriesasMA(ds,'Fc')
# uncomment following 3 lines to use synthetic data
#Fc = numpy.ma.ones(nRecs)*float(5)
#index = numpy.ma.where(ustar<0.25)[0]
#Fc[index] = float(20)*ustar[index]
dt = ds.series['DateTime']['Data']
# get the night time values where Fc is not masked
index = numpy.ma.where((Fsd<10)&(Fc.mask==False)&(ustar.mask==False)&(Ta.mask==False))[0]
Fc_n = Fc[index]
us_n = ustar[index]
Ta_n = Ta[index]
dt_n = [dt[i] for i in index]
# loop over the "seasons"
si = 0; ei = 0; n = 0
nPoints = len(index)
while ei<=nPoints-npointsperjump:
#while si==0:
    ei = npointsperseason + n*npointsperjump - 1
    #print n, si, ei
    Fc_ns = Fc_n[si:ei]
    us_ns = us_n[si:ei]
    Ta_ns = Ta_n[si:ei]
    dt_ns = dt_n[si:ei]
    print 'Processing ',str(dt_ns[0]),' to ',str(dt_ns[-1])
    # get the temperature quantiles
    Ta_bins = [min(Ta_ns)]
    for i in range(0,nTemp-1):
        percent = (i+1)*100/nTemp
        Ta_bins.append(numpy.percentile(Ta_ns,percent))
    Ta_bins.append(max(Ta_ns))
    # separate the data into the 4 temperature categories
    Fc_Tclass = {}
    us_Tclass = {}
    for i in range(0,nTemp):
        index = numpy.ma.where((Ta_ns>=Ta_bins[i])&(Ta_ns<Ta_bins[i+1]))
        Fc_Tclass[str(i)] = Fc_ns[index]
        us_Tclass[str(i)] = us_ns[index]
    # now get the average Fc and ustar across "nustarbins" equally sized bins
    nFig = nFig + 1
    fig = plt.figure()
    plt.figtext(0.5,0.94,str(dt_ns[0])+' to '+str(dt_ns[-1]),ha='center',size=16)
    # loop over the temperature classes
    for i in range(0,nTemp):
        # local pointers to the appropriate dictionary entries
        fc = Fc_Tclass[str(i)]
        us = us_Tclass[str(i)]
        # first element of the ustar bin edges is the minimum value
        us_bins = [min(us)]
        # now get the rest of the bin edges, here set so that each bin contains 5 points
        for j in range(0,nustarbins-1):
            percent = (j+1)*100/nustarbins
            us_bins.append(numpy.percentile(us,percent))
        # last element is the maximum value
        us_bins.append(max(us))
        # use numpy.histogram to get the number of points in each ustar bin
        num,_ = numpy.histogram(us,us_bins)
        # use numpy.histogram to get the sum of data in each bin
        sumFc_bins,_ = numpy.histogram(us,us_bins,weights=fc)
        sumus_bins,_ = numpy.histogram(us,us_bins,weights=us)
        # and then get the average in each bin
        avgFc_bins = sumFc_bins/num
        avgus_bins = sumus_bins/num
        # now plot the bin averages
        axl = fig.add_subplot(4,2,2*i+1)
        axl.plot(avgus_bins,avgFc_bins,'b.')
        axl.set_xlabel('u* (m/s)')
        axl.set_ylabel('Fc (umol/m2/s)')
        textstr = '{0:.1f} to {1:.1f}'.format(Ta_bins[i],Ta_bins[i+1])
        axl.text(0.075,0.875,textstr,transform=axl.transAxes,ha='left',size=12)
        # now get the Fc score and plot this as a function of bin average ustar
        Fscore = get_Fscore_1b(avgFc_bins,avgus_bins)
        #print Fscore
        # now plot the F scores
        axr = fig.add_subplot(4,2,2*(i+1))
        axr.plot(avgus_bins[1:-2],Fscore['values'][1:-2],'b.')
        axr.set_xlabel('u* (m/s)')
        axr.set_ylabel('F score')
        textstr = 'Fmax {0:.1f} at u* = {1:.2f} (m/s)'.format(Fscore['Fmax'],Fscore['usatFmax'])
        axr.text(0.075,0.875,textstr,transform=axr.transAxes,ha='left',size=12)
        # now plot the best fit lines on the binned us and Fc data
        iatFmax = Fscore['iatFmax']
        b0 = Fscore['b'][str(iatFmax)][0]
        b1 = Fscore['b'][str(iatFmax)][1]
        axl.plot(avgus_bins[0:iatFmax+1],b1+b0*avgus_bins[0:iatFmax+1],'r-')
        axl.plot(avgus_bins[iatFmax:],[b1+b0*avgus_bins[iatFmax]]*len(avgus_bins[iatFmax:]),'r-')
    si = ei - npointsperjump + 1
    n = n + 1
plt.show()
# write the "seasonal" data to 