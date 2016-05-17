############################################################## 
# Date: 20/01/16
# Name: plot_IB_guess.py
# Author: Alek Petty
# Description: Script to plot estimate of thickness using sail heights and specified regression parameter
# Input requirements: indy results and IB 10km averages, produced using calc_sail_thickness_ave.py script


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
import os
import itertools
from scipy.spatial import cKDTree as KDTree
from scipy import stats
import scipy.optimize as optimization

rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9
rcParams['legend.fontsize']=9
rcParams['font.size']=9
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

mplot = Basemap(projection='npstere',boundinglat=66,lon_0=0, resolution='l'  )
mplot2=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)

thresh=20
fadd=''
ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd

datapath = './Data_output/'+ftype+'/'
avefilepath=datapath+'/CORRELATIONS/INDY/'
figpath = './Figures/'
rawdatapath='../../../DATA/'

my_cmap=ro.perceptual_colormap("Linear_L", rawdatapath+'/OTHER/CMAPS/', reverse=1)


year=2015
mean_height_years = load(avefilepath+'/Sailheightave_10km_'+str(year)+'.txt')
mean_thicknessIBR_years = load(avefilepath+'/IBthickness_10km_'+str(year)+'.txt')
#mean_snowIBR_years.append(load(savepath+'/IBsnow_10km_'+str(year)+'.txt'))
xpts_matchT = load(avefilepath+'/xpts_10km_'+str(year)+'.txt')
ypts_matchT = load(avefilepath+'/ypts_10km_'+str(year)+'.txt')

lonT, latT=mplot(xpts_matchT, ypts_matchT,  inverse=True)
xpts_match, ypts_match = mplot2(lonT, latT)



ice_typeT, xpts_type, ypts_type = ro.get_mean_ice_type(mplot2, rawdatapath, year, res=1)
#lon_type, lat_type=m(xpts_type, ypts_type,  inverse=True)
#xpts_type, ypts_type = mplot2(lon_type, lat_type)
ice_type=np.zeros((shape(ice_typeT)))
ice_type[where(ice_typeT>0.9)]=0.6
ice_type[where((ice_typeT<0.9) & (ice_typeT>0.6))]=0.4
ice_type[where((ice_typeT<0.6) & (ice_typeT>0.4))]=0.2




axesname = ['ax1', 'ax2', 'ax3', 'ax4', 'ax5', 'ax6']
textwidth=5.
fig = figure(figsize=(textwidth,textwidth*0.65))
subplots_adjust(left = 0.01, right = 0.99, bottom=0.14, top = 0.99, wspace=0.03, hspace=0.2)


ax1 = subplot(2, 2, 1)
minval=0
maxval=2
im1 = mplot2.pcolormesh(xpts_type , ypts_type, ice_type, edgecolors='white', vmin=0, vmax=1, cmap=cm.Greys,shading='gouraud', zorder=1, rasterized=True)
im11 = ax1.hexbin(xpts_match, ypts_match, C=mean_height_years, gridsize=100, vmin=minval, vmax=maxval, cmap=my_cmap, zorder=2, rasterized=True)
mplot2.fillcontinents(color='w',lake_color='grey', zorder=3)
mplot2.drawcoastlines(linewidth=0.25, zorder=5)
mplot2.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mplot2.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mplot2.drawparallels(np.arange(90,-90,-10),labels=[False,False,True,False], fontsize=7, linewidth = 0.25, zorder=10)

		

xS, yS = mplot2(177, 64.2)


ax2 = subplot(2, 2, 2)
minval=0
maxval=5
im2 = mplot2.pcolormesh(xpts_type , ypts_type, ice_type, edgecolors='white', vmin=0, vmax=1, cmap=cm.Greys,shading='gouraud', zorder=1, rasterized=True)
im21 = ax2.hexbin(xpts_match, ypts_match, C=((mean_height_years)/0.72)**2, gridsize=100, vmin=minval, vmax=maxval, cmap=my_cmap, zorder=2, rasterized=True)
mplot2.fillcontinents(color='w',lake_color='grey', zorder=3)
mplot2.drawcoastlines(linewidth=0.25, zorder=5)
mplot2.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mplot2.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mplot2.drawparallels(np.arange(90,-90,-10),labels=[False,False,True,False], fontsize=7, linewidth = 0.25, zorder=10)

xS, yS = mplot2(177, 64.2)
#ax2.text(xS, yS, str(year), zorder = 11)
xS, yS = mplot2(-150, 60.2)
bbox_props = dict(boxstyle="square,pad=0.3", fc="white")
ax2.annotate('Estimate', xy=(0.04, 0.15), bbox=bbox_props, xycoords='axes fraction', horizontalalignment='left', verticalalignment='top', zorder=10)


ax3 = subplot(2, 2, 3)
minval=0
maxval=5
im3 = mplot2.pcolormesh(xpts_type , ypts_type, ice_type, edgecolors='white', vmin=0, vmax=1, cmap=cm.Greys,shading='gouraud', zorder=1, rasterized=True)
im31 = ax3.hexbin(xpts_match, ypts_match, C=mean_thicknessIBR_years, gridsize=100, vmin=minval, vmax=maxval, cmap=my_cmap, zorder=2, rasterized=True)
mplot2.fillcontinents(color='w',lake_color='grey', zorder=3)
mplot2.drawcoastlines(linewidth=0.25, zorder=5)
mplot2.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mplot2.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mplot2.drawmeridians(np.arange(-180.,180.,60.),labels=[False,False,False,True], fontsize=7,latmax=85, linewidth = 0.25, zorder=10)
#xS, yS = mplot2(177, 64.2)
#ax3.text(xS, yS, str(year), zorder = 11)
ax3.annotate('IB', xy=(0.04, 0.15), bbox=bbox_props, xycoords='axes fraction', horizontalalignment='left', verticalalignment='top', zorder=10)


ax4 = subplot(2, 2, 4)
minval=-2
maxval=2
im4 = mplot2.pcolormesh(xpts_type , ypts_type, ice_type, edgecolors='white', vmin=0, vmax=1, cmap=cm.Greys,shading='gouraud', zorder=1, rasterized=True)
im41 = ax4.hexbin(xpts_match, ypts_match, C=(((mean_height_years)/0.72)**2) - mean_thicknessIBR_years, gridsize=100, vmin=minval, vmax=maxval, cmap=cm.RdBu, zorder=2, rasterized=True)
mplot2.fillcontinents(color='w',lake_color='grey', zorder=3)
mplot2.drawcoastlines(linewidth=0.25, zorder=5)
mplot2.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mplot2.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mplot2.drawmeridians(np.arange(-180.,180.,60.),labels=[False,False,False,True], fontsize=7,latmax=85, linewidth = 0.25, zorder=10)
#xS, yS = mplot2(177, 64.2)
#ax4.text(xS, yS, str(i+start_year), zorder = 11)
ax4.annotate('Est-IB', xy=(0.04, 0.15), bbox=bbox_props, xycoords='axes fraction', horizontalalignment='left', verticalalignment='top', zorder=10)


cax = fig.add_axes([0.1, 0.595, 0.3, 0.03])
cbar = colorbar(im11,cax=cax, orientation='horizontal', extend='both',use_gridspec=True)
xticks = np.linspace(0, 2, 6)
cbar.set_ticks(xticks)
cbar.solids.set_rasterized(True)
#cbar.ax.xaxis.set_ticks_position('top')
#cbar.ax.xaxis.set_label_position('top')
cbar.set_label('Surface feature height (m)', labelpad=2)

cax2 = fig.add_axes([0.6, 0.595, 0.3, 0.03])
cbar2 = colorbar(im21,cax=cax2, orientation='horizontal', extend='both',use_gridspec=True)
xticks2 = np.linspace(0, 5, 6)
cbar2.set_ticks(xticks2)
cbar2.solids.set_rasterized(True)
#cbar2.ax.xaxis.set_ticks_position('top')
#cbar2.ax.xaxis.set_label_position('top')
cbar2.set_label('Ice thickness (m)', labelpad=2)

cax3 = fig.add_axes([0.1, 0.1, 0.3, 0.03])
cbar3 = colorbar(im31,cax=cax3, orientation='horizontal', extend='both',use_gridspec=True)

xticks3 = np.linspace(0, 5, 6)
cbar3.set_ticks(xticks3)
cbar3.solids.set_rasterized(True)
#cbar3.ax.yaxis.set_ticks_position('left')
#cbar3.ax.yaxis.set_label_position('left')
cbar3.set_label('Ice thickness (m)', labelpad=2)

cax4 = fig.add_axes([0.6, 0.1, 0.3, 0.03])
cbar4 = colorbar(im41,cax=cax4, orientation='horizontal', extend='both',use_gridspec=True)
cbar4.set_label('Difference (m)', labelpad=2)
xticks4 = np.linspace(-2, 2, 5)
cbar4.set_ticks(xticks4)
cbar4.solids.set_rasterized(True)

#savefig(figpath+'/IB_thickness_guess'+ftype+'nosnow.png', dpi=300)
savefig(figpath+'/figure13.pdf', dpi=300)


