import matplotlib.pyplot as plt
import numpy

x = numpy.arange(0,360*10,1)
y = numpy.sin(x*numpy.pi/float(180))

nDrivers = 3
n_ts = nDrivers + 1
# margins
margin_bottom = 0.05
margin_top = 0.10
margin_left = 0.05
Margin_right = 0.05
# plot heights, plot widths and spaces between plots
xy_height = 0.20
xy_width = 0.20
xyts_space = 0.05
ts_width = 0.9
# calculate bottom of the first time series and the height of the time series plots
ts_bottom = margin_bottom + xy_height + xyts_space
ts_height = (1.0 - margin_top - ts_bottom)/float(nDrivers+1)
# make the figure
fig=plt.figure(1,figsize=(13,9))
rect1 = [0.10,margin_bottom,xy_width,xy_height]
ax1 = plt.axes(rect1)
rect2 = [0.40,margin_bottom,xy_width,xy_height]
ax2 = plt.axes(rect2)
rect3 = [0.70,margin_bottom,xy_width,xy_height]
ax3 = plt.axes(rect3)

ts_axes = []
for i in range(n_ts):
    this_bottom = ts_bottom + i*ts_height
    rect = [margin_left,this_bottom,ts_width,ts_height]
    if i==0:
        ts_axes.append(plt.axes(rect))
    else:
        ts_axes.append(plt.axes(rect,sharex=ts_axes[0]))

plt.show()