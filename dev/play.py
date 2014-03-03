import csv
import matplotlib.pyplot as plt
import numpy
import os
import qcio
import qcutils
import qcts
from scipy.signal import butter,filtfilt
import subprocess
import sys
sys.path.append("/usr/local/lib/python2.7/dist-packages")
from ffnet import ffnet, mlgraph

fname = qcio.get_filename_dialog()
ds = qcio.nc_read_series(fname)

nRecs =int(ds.globalattributes['nc_nrecs'])

Fc,f=qcutils.GetSeriesasMA(ds,'Fc')
Month,f=qcutils.GetSeriesasMA(ds,'Month')
us,f=qcutils.GetSeriesasMA(ds,'ustar')
Fsd,f=qcutils.GetSeriesasMA(ds,'Fsd')
dt = ds.series['DateTime']['Data']

qcts.InterpolateOverMissing(ds,'Sws')
Sws,f=qcutils.GetSeriesasMA(ds,'Sws')
qcts.InterpolateOverMissing(ds,'Ts')
Ts,f=qcutils.GetSeriesasMA(ds,'Ts')

b,a = butter(3,0.1)
Sws_notMA,WasMA = qcutils.MAtoSeries(Sws)
Sws_notMA_filt = filtfilt(b,a,Sws_notMA)
Sws_filt = numpy.ma.masked_where(Sws_notMA_filt<0.0,Sws_notMA_filt)
Sws_filt = numpy.ma.masked_where(Sws_notMA_filt>0.5,Sws_filt)
Sws_filt = Sws

Ts_notMA,WasMA = qcutils.MAtoSeries(Ts)
Ts_notMA_filt = filtfilt(b,a,Ts_notMA)
Ts_filt,WasND = qcutils.SeriestoMA(Ts_notMA_filt)
Ts_filt = Ts

Fc_night=numpy.ma.masked_where(Fsd>10,Fc)
Fc_filt=numpy.ma.masked_where(Fc<0.05,Fc_night)
Fc_turb=numpy.ma.masked_where(us<0.3,Fc_filt)

# plot the drivers and the target
#fig1,(ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,sharex=True)
#ax1.plot(dt,Fc,'b.')
#ax2.plot(dt,Fc_night,'b.')
#ax3.plot(dt,Fc_turb,'b.')
#ax4.plot(dt,Ts,'b.',dt,Ts_filt,'r-')
#ax5.plot(dt,Sws,'b.',dt,Sws_filt,'r-')
#plt.draw()

# change the working directory
cwd = os.getcwd()
os.chdir('../')
# specify the drivers
driverlist = [Ts_filt,Sws_filt]
# synchronise any gaps due to missing data
cind = numpy.zeros(nRecs)
for ThisOne in driverlist:
    index = numpy.ma.where(ThisOne.mask==True)[0]
    cind[index] = 1
index = numpy.where(cind==0)[0]
nRecs_filt = len(index)
nDrivers = len(driverlist)
sofminputdata = numpy.zeros((nRecs_filt,nDrivers))
for i,ThisOne in zip(range(nDrivers),driverlist):
    sofminputdata[:,i] = ThisOne[index]
# write sofm input file
sofmfile = open('solo/input/sofm_input.csv','wb')
wr = csv.writer(sofmfile,delimiter=',')
for i in range(sofminputdata.shape[0]):
    wr.writerow(sofminputdata[i,0:nDrivers])
sofmfile.close()
# run sofm
#subprocess.call(['./solo/bin/sofm','solo/inf/sofm.inf'],stdout=sofmlogfile)
subprocess.call(['./solo/bin/sofm','solo/inf/sofm.inf'])

# write solo input file
# specify the drivers
targetlist = [Fc_turb]
sololist = driverlist+targetlist
# synchronise any gaps due to missing data
cind = numpy.zeros(nRecs)
for ThisOne in sololist:
    index = numpy.ma.where(ThisOne.mask==True)[0]
    cind[index] = 1
index = numpy.where(cind==0)[0]
nRecs_filt = len(index)
nDrivers = len(sololist)
soloinputdata = numpy.zeros((nRecs_filt,nDrivers))
for i,ThisOne in zip(range(nDrivers),sololist):
    soloinputdata[:,i] = ThisOne[index]
# write solo input file
solofile = open('solo/input/solo_input.csv','wb')
wr = csv.writer(solofile,delimiter=',')
for i in range(soloinputdata.shape[0]):
    wr.writerow(soloinputdata[i,0:nDrivers])
solofile.close()
# run solo
subprocess.call(['./solo/bin/solo','solo/inf/solo.inf'])

# write seqsolo input file
# synchronise any gaps due to missing data
cind = numpy.zeros(nRecs)
for ThisOne in driverlist:
    index = numpy.ma.where(ThisOne.mask==True)[0]
    cind[index] = 1
index = numpy.where(cind==0)[0]
nRecs_filt = len(index)
nDrivers = len(driverlist)
seqsoloinputdata = numpy.zeros((nRecs_filt,nDrivers+1))
for i,ThisOne in zip(range(nDrivers),driverlist):
    seqsoloinputdata[:,i] = ThisOne[index]
target,WasMA = qcutils.MAtoSeries(Fc_turb)
seqsoloinputdata[:,nDrivers] = target[index]
# write seqsolo input file
seqsolofile = open('solo/input/seqsolo_input.csv','wb')
wr = csv.writer(seqsolofile,delimiter=',')
for i in range(seqsoloinputdata.shape[0]):
    wr.writerow(seqsoloinputdata[i,0:nDrivers+1])
seqsolofile.close()
# run seqsolo
subprocess.call(['./solo/bin/seqsolo','solo/inf/seqsolo.inf'])
# get the seqsolo results
seqdata = numpy.genfromtxt('solo/output/seqOut0.out')
Reco_SOLO = seqdata[:,1]
Reco_SOLO_filt = filtfilt(b,a,Reco_SOLO)

# plot the drivers, target and the SOLO predictions
fig2,(ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,sharex=True)
ax1.plot(dt,Fc,'b.')
ax2.plot(dt,Fc_night,'b.')
ax3.plot(dt,Fc_turb,'b.',dt,Reco_SOLO,'r-',dt,Reco_SOLO_filt,'g-')
#ax3.plot(dt,Fc_turb,'b.',dt,Reco_SOLO_filt,'r-')
ax4.plot(dt,Ts,'b.',dt,Ts_filt,'r-')
ax5.plot(dt,Sws,'b.',dt,Sws_filt,'r-')
#plt.show()

# let's try ffnet
# get the training and target data
cind = numpy.zeros(nRecs)
for ThisOne in [Fc_turb,Ts_filt,Sws_filt]:
    index = numpy.ma.where(ThisOne.mask==True)[0]
    cind[index] = 1
index = numpy.where(cind==0)[0]
trainingdata = numpy.zeros((len(index),2))
trainingdata[:,0] = Ts_filt[index]
trainingdata[:,1] = Sws_filt[index]
targetdata = numpy.zeros(len(index))
targetdata[:] = Fc_turb[index]
# build the network
conec = mlgraph((2,5,1))
net = ffnet(conec)
# train the network
net.train_tnc(trainingdata,targetdata,maxfun=100)
# get the training statistics
net.test(trainingdata,targetdata)
# get the prediction data
cind = numpy.zeros(nRecs)
for ThisOne in [Ts_filt,Sws_filt]:
    index = numpy.ma.where(ThisOne.mask==True)[0]
    cind[index] = 1
index = numpy.where(cind==0)[0]
predictiondata = numpy.zeros((len(index),2))
predictiondata[:,0] = Ts_filt[index]
predictiondata[:,1] = Sws_filt[index]
# get the predicted Fc
Reco_ffnet = net(predictiondata)
# plot the drivers, target and the SOLO predictions
fig3,(ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,sharex=True)
ax1.plot(dt,Fc,'b.')
ax2.plot(dt,Fc_night,'b.')
ax3.plot(dt,Fc_turb,'b.',dt,Reco_ffnet,'r-')
ax4.plot(dt,Ts,'b.',dt,Ts_filt,'r-')
ax5.plot(dt,Sws,'b.',dt,Sws_filt,'r-')
plt.show()

# clean up
os.chdir(cwd)
