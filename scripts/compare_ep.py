import numpy
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import qcio
import qcplot
import qcutils

epname = qcio.get_filename_dialog(title='Choose an EddyPro full output file')
ofname = qcio.get_filename_dialog(title='Choose an L3 output file')

ds_ep = qcio.read_eddypro_full(epname)
ds_of = qcio.nc_read_series(ofname)

dt_ep = ds_ep.series['DateTime']['Data']
dt_of = ds_of.series['DateTime']['Data']

si = dt_of.index(dt_ep[0])
ei = dt_of.index(dt_ep[-1])

us_of,f,a = qcutils.GetSeriesasMA(ds_of,'ustar',si=si,ei=ei)
us_ep,f,a = qcutils.GetSeriesasMA(ds_ep,'ustar')
Fh_of,f,a = qcutils.GetSeriesasMA(ds_of,'Fh',si=si,ei=ei)
Fh_ep,f,a = qcutils.GetSeriesasMA(ds_ep,'Fh')
Fe_of,f,a = qcutils.GetSeriesasMA(ds_of,'Fe',si=si,ei=ei)
Fe_ep,f,a = qcutils.GetSeriesasMA(ds_ep,'Fe')
Fc_of,f,a = qcutils.GetSeriesasMA(ds_of,'Fc',si=si,ei=ei)
Fc_ep,f,a = qcutils.GetSeriesasMA(ds_ep,'Fc')

us_of.mask = numpy.ma.mask_or(us_of.mask,us_ep.mask)
us_ep.mask = numpy.ma.mask_or(us_of.mask,us_ep.mask)
print "us_of: ",numpy.ma.maximum(us_of),numpy.ma.minimum(us_of)
print "us_ep: ",numpy.ma.maximum(us_ep),numpy.ma.minimum(us_ep)
Fh_of.mask = numpy.ma.mask_or(Fh_of.mask,Fh_ep.mask)
Fh_ep.mask = numpy.ma.mask_or(Fh_of.mask,Fh_ep.mask)
print "Fh_of: ",numpy.ma.maximum(Fh_of),numpy.ma.minimum(Fh_of)
print "Fh_ep: ",numpy.ma.maximum(Fh_ep),numpy.ma.minimum(Fh_ep)
Fe_of.mask = numpy.ma.mask_or(Fe_of.mask,Fe_ep.mask)
Fe_ep.mask = numpy.ma.mask_or(Fe_of.mask,Fe_ep.mask)
print "Fe_of: ",numpy.ma.maximum(Fe_of),numpy.ma.minimum(Fe_of)
print "Fe_ep: ",numpy.ma.maximum(Fe_ep),numpy.ma.minimum(Fe_ep)
Fc_of.mask = numpy.ma.mask_or(Fc_of.mask,Fc_ep.mask)
Fc_ep.mask = numpy.ma.mask_or(Fc_of.mask,Fc_ep.mask)
print "Fc_of: ",numpy.ma.maximum(Fc_of),numpy.ma.minimum(Fc_of)
print "Fc_ep: ",numpy.ma.maximum(Fc_ep),numpy.ma.minimum(Fc_ep)

fig = plt.figure(1,figsize=(8,8))
qcplot.xyplot(us_ep,us_of,sub=[2,2,1],regr=1,xlabel='u*_EP (m/s)',ylabel='u*_OF (m/s)')
qcplot.xyplot(Fh_ep,Fh_of,sub=[2,2,2],regr=1,xlabel='Fh_EP (W/m2)',ylabel='Fh_OF (W/m2)')
qcplot.xyplot(Fe_ep,Fe_of,sub=[2,2,3],regr=1,xlabel='Fe_EP (W/m2)',ylabel='Fe_OF (W/m2)')
qcplot.xyplot(Fc_ep,Fc_of,sub=[2,2,4],regr=1,xlabel='Fc_EP (umol/m2/s)',ylabel='Fc_OF (umol/m2/s)')
plt.tight_layout()
plt.show()