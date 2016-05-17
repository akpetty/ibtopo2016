############################################################## 
# Date: 20/01/16
# Name: plot_ridges6_levpercent.py
# Author: Alek Petty
# Description: Script to plot indy topography statistics (level percent)
# Input requirements: Topography stats across all years 2009-2014
# Output: plot of level percentile on ice type background for each year

import matplotlib
matplotlib.use("AGG")
# basemap import
from mpl_toolkits.basemap import Basemap, shiftgrid
# Numpy import
#from osgeo import osr, gdal
import numpy as np
import mpl_toolkits.basemap.pyproj as pyproj
from pylab import *
import IB_functions as ro
import numpy.ma as ma
from matplotlib.patches import Polygon
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import mark_inset
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
#import cubehelix
from netCDF4 import Dataset
from matplotlib import rc
import itertools
from glob import glob
from operator import truediv, mul



rcParams['axes.labelsize'] =9
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9
rcParams['legend.fontsize']=9
rcParams['font.size']=9
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})


mplot=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)
#-------------- IB Projection------------------
mib=pyproj.Proj("+init=EPSG:3413")

thresh=20
plot_max=1
fadd = ''
ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd
print 'file path', ftype
print 'plot type:', plot_max
figpath = './Figures/'
datapath='./Data_output/'+ftype+'/'
rawdatapath='../../../DATA/'
my_cmap=ro.perceptual_colormap("Linear_L", rawdatapath+'/OTHER/CMAPS/', reverse=1)

xpts=[]
ypts=[]
levpercent=[]
for year in xrange(2009, 2015):
	#BIG IS IN TERMS OF AREA - WHERE WE REMOVE RIDGES LESS THAN 100m2
	xptsT, yptsT, levpercentT = ro.get_bulk_ridge_stats(mib, mplot, year, datapath, getlevpercent=1)
	xpts.append(xptsT)
	ypts.append(yptsT)
	levpercent.append(levpercentT)

label_str='Level ice elevation percenitle (%)'
out_str='lev_percent'
minval=20
maxval=50


ice_type=[]
for year in xrange(2009, 2015):
	ice_typeT, xpts_type, ypts_type = ro.get_mean_ice_type(mplot, rawdatapath, year, res=1)
	ice_typeT2=np.zeros((shape(ice_typeT)))
	ice_typeT2[where(ice_typeT>0.9)]=0.6
	ice_typeT2[where((ice_typeT<0.9) & (ice_typeT>0.6))]=0.4
	ice_typeT2[where((ice_typeT<0.6) & (ice_typeT>0.4))]=0.2
	ice_type.append(ice_typeT2)

lonlatBC = [-170., -120., 69., 79.]
lonlatCA = [-150., 10., 81., 90.]
xptsBC, yptsBC = ro.get_box_xy(mplot, lonlatBC)
xptsCA, yptsCA = ro.get_box_xy(mplot, lonlatCA)

axesname = ['ax1', 'ax2', 'ax3', 'ax4', 'ax5', 'ax6']

aspect = mplot.ymax/mplot.xmax

textwidth=5.
fig = figure(figsize=(textwidth,(textwidth*(3./2.)*aspect)+0.9))
subplots_adjust(left = 0.01, right = 0.99, bottom=0.15, top = 0.95, wspace = 0.02, hspace=0.01)

for i in xrange(6):
	vars()['ax'+str(i)] = subplot(3, 2, i+1)
	axT=gca()
	im0 = mplot.pcolormesh(xpts_type , ypts_type, ice_type[i], edgecolors='white', vmin=0, vmax=1, cmap=cm.Greys,shading='gouraud', zorder=1, rasterized=True)
	im1 = vars()['ax'+str(i)].hexbin(xpts[i], ypts[i], C=levpercent[i], vmin=minval, vmax=maxval, cmap=my_cmap, zorder=2, rasterized=True)
	mplot.fillcontinents(color='w',lake_color='grey', zorder=3)
	mplot.drawcoastlines(linewidth=0.25, zorder=5)
	mplot.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	mplot.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	if (i<2):
		mplot.drawparallels(np.arange(90,-90,-10),labels=[False,False,True,False], fontsize=8, linewidth = 0.25, zorder=10)
	if (i>3):
		mplot.drawmeridians(np.arange(-180.,180.,60.),labels=[False,False,False,True], fontsize=8,latmax=85, linewidth = 0.25, zorder=10)

	mplot.plot(xptsCA, yptsCA, '--', linewidth = 1, color='r', zorder=12)
	mplot.plot(xptsBC, yptsBC, '--', linewidth = 1, color='b', zorder=12)
	xS, yS = mplot(177, 64.8)
	vars()['ax'+str(i)].text(xS, yS, str(i+2009), zorder = 11)

cax = fig.add_axes([0.25, 0.09, 0.5, 0.03])
cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both',use_gridspec=True)

cbar.set_label(label_str, labelpad=2)
xticks = np.linspace(minval, maxval, 6)
cbar.set_ticks(xticks)
cbar.solids.set_rasterized(True)
#savefig(out_path+'Ridges_6years_bulkvar'+ftype+'.png', dpi=300)
savefig(figpath+'/figureS4.png', dpi=200)




