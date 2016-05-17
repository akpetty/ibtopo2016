############################################################## 
# Date: 20/01/16
# Name: plot_coast_distance.py
# Author: Alek Petty
# Description: Script to plot coastal proximity
# Input requirements: coastal proximity data
# Output: coastlien proximity map

import matplotlib
matplotlib.use("AGG")
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import mpl_toolkits.basemap.pyproj as pyproj
from pylab import *
import IB_functions as ro
import numpy.ma as ma
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import mark_inset
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
from netCDF4 import Dataset
from matplotlib import rc
from glob import glob
from scipy.interpolate import griddata
import os

#-------------- GET IB Projection ------------------
mplot=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)


rawdatapath='../../../DATA/'
figpath = './Figures/'

my_cmap=ro.perceptual_colormap("Linear_L", rawdatapath+'OTHER/CMAPS/', reverse=1)


coast_file = Dataset(rawdatapath+'/OTHER/distance-to-coast_2m.nc','r')
lat_c = coast_file.variables['lat'][:]
lon_c = coast_file.variables['lon'][:]
z_c = coast_file.variables['z'][:]

#REMOVE SOME OF THE COAST DISTANCE DATA THAT ISN'T NEEDED
lat_c=lat_c[4250:-1]
lon_c=lon_c[::20]
z_c=z_c[4250:-1, ::20]

xpts_c, ypts_c = mplot(*np.meshgrid(lon_c, lat_c))

lonlatBC = [-170., -120., 69., 79.]
lonlatCA = [-150., 10., 81., 90.]
xptsBC, yptsBC = ro.get_box_xy(mplot, lonlatBC)
xptsCA, yptsCA = ro.get_box_xy(mplot, lonlatCA)

aspect = mplot.ymax/mplot.xmax


textwidth=3.
fig = figure(figsize=(textwidth+0.8,(textwidth*aspect)+0.4))

ax1=subplot(1, 1, 1)

vals=[0, 100, 200, 300, 400, 500, 600, 700,800, 1000]
im1 = mplot.contourf(xpts_c , ypts_c, z_c, zorder=1, levels=vals, extend='both', cmap=my_cmap, rasterized=True)
im2 = contour(xpts_c , ypts_c, z_c,0, colors='k', linewidth = 0.25)

#mplot.fillcontinents(color='w',lake_color='grey', zorder=3)
#mplot.drawcoastlines(linewidth=0.25, zorder=5)
mplot.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mplot.drawparallels(np.arange(90,-90,-10), labels=[False,False,True,False], fontsize=8,linewidth = 0.25, zorder=5)
mplot.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mplot.drawmeridians(np.arange(-180.,180.,60.),labels=[False,False,False,True], fontsize=8, linewidth = 0.25, zorder=5)

mplot.plot(xptsCA, yptsCA, '--', linewidth = 1, color='r', zorder=12)
mplot.plot(xptsBC, yptsBC, '--', linewidth = 1, color='b', zorder=12)


cax = fig.add_axes([0.83, 0.05, 0.03, 0.9])
cbar = colorbar(im1,cax=cax, orientation='vertical', extend='both',use_gridspec=True)
#cbar.set_alpha(1)
#cbar.draw_all()
cbar.set_label('Distance (km)', rotation=270, labelpad=9)
#xticks = np.linspace(minval, maxval, 4)
#cbar.set_ticks(xticks)
cbar.solids.set_rasterized(True)
subplots_adjust(bottom=0.05, left = 0.01, right = 0.82, top=0.95)

#savefig(out_path+'/Ridges_6years'+ftype+out_var+'.pdf', dpi=200)
savefig(figpath+'/figureS7.png', dpi=200)









