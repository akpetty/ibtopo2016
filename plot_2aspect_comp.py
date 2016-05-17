############################################################## 
# Date: 20/01/16
# Name: plot_2aspect_comp.py
# Author: Alek Petty
# Description: Script to plot aspect ratio across different years
# Input requirements: bulk distributions (All/MYI/FYI in the CA/BC)
# Output: 6 panels highlighting the distributions

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
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9
rcParams['legend.fontsize']=9
rcParams['font.size']=9
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def return_aspect(year, ftype):
	datapathT=datapath+ftype+'/'

	xptsT, yptsT, heightT, nearmax_heightT, max_heightT, sizeRT = ro.get_indy_mean_max_height(mib, mplot, year, datapathT)
	c00T, c01T, c10T, c11T, sect_numRT= ro.get_indy_covar(year, datapathT)
	region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsT, yptsT), method='nearest')
	mask = where((region_maskR<11)&(max_heightT<10))
	xptsT=xptsT[mask]
	yptsT=yptsT[mask]
	heightT=heightT[mask]
	max_heightT=max_heightT[mask]
	sect_numRT=sect_numRT[mask]
	c00T=c00T[mask]
	c01T=c01T[mask]
	c10T=c10T[mask]
	c11T=c11T[mask]

	#CHECK BODMAS
	l1T=(c00T+c11T-np.sqrt((c00T+c11T)**2.0-4*(c00T*c11T-c01T**2.0)))/2.0
	l2T=(c00T+c11T+np.sqrt((c00T+c11T)**2.0-4*(c00T*c11T-c01T**2.0)))/2.0
	#theta=0.5*180.0/np.pi*np.arctan2(2*c01T,c00T-c11T)
	#b=np.sqrt((dx**2)*sizeRT*np.sqrt(l2T/l1T))

	aspectT = sqrt(l2T/l1T)
	return xptsT, yptsT, aspectT

mplot=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)
#-------------- GET DMS Projection ------------------
mib=pyproj.Proj("+init=EPSG:3413")

thresh=20
plot_max=1
fadd = ''

ftype1='1km_xyres2m_'+str(20)+'cm'+fadd
ftype2='1km_xyres2m_'+str(80)+'cm'+fadd

datapath = './Data_output/'
figpath = './Figures/'
rawdatapath='../../../DATA/'
xpts=[]
ypts=[]
height=[]
aspect=[]
region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)

year=2012

xptsT1, yptsT1, aspectT1 = return_aspect(year, ftype1)
xptsT2, yptsT2, aspectT2 = return_aspect(year, ftype2)
		
ridge_var=aspect
minval = 1
maxval = 4
label_str = 'Aspect ratio'
out_var = 'aspect_ratio'


ice_typeT, xpts_type, ypts_type = ro.get_mean_ice_type(mplot, rawdatapath, year, res=1)
ice_type=np.zeros((shape(ice_typeT)))
ice_type[where(ice_typeT>0.9)]=0.6
ice_type[where((ice_typeT<0.9) & (ice_typeT>0.6))]=0.4
ice_type[where((ice_typeT<0.6) & (ice_typeT>0.4))]=0.2

my_cmap=ro.perceptual_colormap("Linear_L", rawdatapath+'OTHER/CMAPS/', reverse=1)
axesname = ['ax1', 'ax2', 'ax3', 'ax4', 'ax5', 'ax6']

lonlatBC = [-170., -120., 69., 79.]
lonlatCA = [-150., 10., 81., 90.]
xptsBC, yptsBC = ro.get_box_xy(mplot, lonlatBC)
xptsCA, yptsCA = ro.get_box_xy(mplot, lonlatCA)

aspect = mplot.ymax/mplot.xmax

textwidth=5.
fig = figure(figsize=(textwidth,(textwidth*(1./2.)*aspect)+0.8))
subplots_adjust(left = 0.01, right = 0.99, bottom=0.27, top = 0.95, wspace = 0.02, hspace=0.01)

ax1 = subplot(1, 2, 1)
im0 = mplot.pcolormesh(xpts_type , ypts_type, ice_type, edgecolors='white', vmin=0, vmax=1, cmap=cm.Greys,shading='gouraud', zorder=1, rasterized=True)
im1 = ax1.hexbin(xptsT1, yptsT1, C=aspectT1, gridsize=100, vmin=minval, vmax=maxval, cmap=my_cmap, zorder=2, rasterized=True)

mplot.fillcontinents(color='w',lake_color='grey', zorder=3)
mplot.drawcoastlines(linewidth=0.25, zorder=5)
mplot.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mplot.drawparallels(np.arange(90,-90,-10), labels=[False,False,True,False], fontsize=8,linewidth = 0.25, zorder=5)
mplot.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mplot.drawmeridians(np.arange(-180.,180.,60.),labels=[False,False,False,True], fontsize=8, linewidth = 0.25, zorder=5)

mplot.plot(xptsCA, yptsCA, '--', linewidth = 1, color='r', zorder=12)
mplot.plot(xptsBC, yptsBC, '--', linewidth = 1, color='b', zorder=12)

xS, yS = mplot(177, 64.2)
ax1.text(xS, yS, str(20)+' cm', zorder = 11)

ax2 = subplot(1, 2, 2)
im2 = mplot.pcolormesh(xpts_type , ypts_type, ice_type, edgecolors='white', vmin=0, vmax=1, cmap=cm.Greys,shading='gouraud', zorder=1, rasterized=True)
im3 = ax2.hexbin(xptsT2, yptsT2, C=aspectT2, gridsize=100, vmin=minval, vmax=maxval, cmap=my_cmap, zorder=2, rasterized=True)

mplot.fillcontinents(color='w',lake_color='grey', zorder=3)
mplot.drawcoastlines(linewidth=0.25, zorder=5)
mplot.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
mplot.drawparallels(np.arange(90,-90,-10), labels=[False,False,True,False], fontsize=8,linewidth = 0.25, zorder=5)
mplot.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
mplot.drawmeridians(np.arange(-180.,180.,60.),labels=[False,False,False,True], fontsize=8, linewidth = 0.25, zorder=5)

mplot.plot(xptsCA, yptsCA, '--', linewidth = 1, color='r', zorder=12)
mplot.plot(xptsBC, yptsBC, '--', linewidth = 1, color='b', zorder=12)

xS, yS = mplot(177, 64.2)
ax2.text(xS, yS, str(80)+' cm', zorder = 11)

cax = fig.add_axes([0.3, 0.2, 0.4, 0.05])
cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both',use_gridspec=True)
cbar.set_label(label_str, labelpad=2)
xticks = np.linspace(minval, maxval, 4)
cbar.set_ticks(xticks)
cbar.solids.set_rasterized(True)


savefig(figpath+'/figure14.png', dpi=200)





