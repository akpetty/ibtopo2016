############################################################## 
# Date: 20/01/16
# Name: plot_ridges_bulk.py
# Author: Alek Petty
# Description: Script to plot wind and ice type and IB flight lines
# Input requirements: ERA-I wind data, ice type, and IB flines
# Extra info: check the wind/ice_type/IB flightline functions for more info on where to put the data.

import matplotlib
matplotlib.use("AGG")

# basemap import
from mpl_toolkits.basemap import Basemap, shiftgrid
# Numpy import
import numpy as np
from pylab import *
from scipy.io import netcdf
import numpy.ma as ma
import string
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import mark_inset
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
from scipy import stats
import IB_functions as ro
from matplotlib import rc
from netCDF4 import Dataset
from glob import glob


rcParams['axes.labelsize'] =10
rcParams['xtick.labelsize']=10
rcParams['ytick.labelsize']=10
rcParams['legend.fontsize']=10
rcParams['font.size']=10
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

#mpl.rc('text', usetex=True)
m=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)

rawdatapath='../../../DATA/'
figpath='./Figures/'
xptsW, yptsW, xvel, yvel, wind_speed = ro.get_era_winds(m, rawdatapath, 2009, 2014, 0, 3, 75)


xpts_all=[]
ypts_all=[]
for year in xrange(2009, 2014+1, 1):
	print year
	lonsT, latsT = ro.calc_icebridge_flights(year,rawdatapath, 'GR')
	xptsT, yptsT=m(lonsT, latsT)
	xpts_all.append(xptsT)
	ypts_all.append(yptsT)


#xpts_all, ypts_all = ro.calc_icebridge_flights_years(m,2009,2014, 'GR')

#xpts_type, ypts_type, ice_type_mean, ice_type, latsT = get_mean_icetype(m, 1)

ice_type=[]

for year in xrange(2009, 2015):
	ice_typeT, xpts_type, ypts_type = ro.get_mean_ice_type(m, rawdatapath, year, res=1)
	ice_type.append(ice_typeT)


ice_type_meanT=np.mean(ice_type, axis=0)
ice_type_mean=np.copy(ice_type_meanT)
ice_type_mean[where(ice_type_meanT>0.9)]=0.75
ice_type_mean[where((ice_type_meanT<0.9) & (ice_type_meanT>0.6))]=0.5
ice_type_mean[where((ice_type_meanT<0.6) & (ice_type_meanT>0.4))]=0.25

#2011
dms2011a_x, dms2011a_y = m(-146.78,72.96)

dms2010a_x, dms2010a_y = m(-101.85,81.18)

dms2011b_x, dms2011b_y = m(-146.78, 73.02)
dms2013_x, dms2013_y = m(-113.95,85.7)

dms2012_x, dms2012_y = m(-132.68,75.64)

dms2014_x, dms2014_y = m(-38.78,86.01)


minval=0
maxval=5
scale_vec=10
vector_val=2
res=5

aspect = m.ymax/m.xmax
textwidth=5.
fig = figure(figsize=(textwidth,(textwidth*(1.8)*aspect)))
ax1 = subplot(2, 1, 1)

im0 = m.pcolormesh(xpts_type , ypts_type, ice_type_mean, edgecolors='white', vmin=0.25, vmax=0.75, cmap=cm.Greys,shading='gouraud', zorder=1)

im1 = m.plot(xpts_all[0] , ypts_all[0], color = 'b', zorder=4)
im2 = m.plot(xpts_all[1] , ypts_all[1], color = 'g', zorder=4)
im3 = m.plot(xpts_all[2] , ypts_all[2], color = 'y', zorder=4)
im4 = m.plot(xpts_all[3] , ypts_all[3], color = 'm', zorder=4)
im5 = m.plot(xpts_all[4] , ypts_all[4], color = 'c', zorder=4)
im6 = m.plot(xpts_all[5] , ypts_all[5], color = 'r', zorder=4)


im7=m.plot(dms2011a_x, dms2011a_y, 'y', marker='*', markersize=10, zorder=6)
im7=m.plot(dms2010a_x, dms2010a_y, 'g', marker='*', markersize=10, zorder=5)
#im8=m.plot(dms2011b_x, dms2011b_y, 'y', marker='v', markersize=8, zorder=5)
im10=m.plot(dms2012_x, dms2012_y, 'b', marker='*', markersize=10, zorder=5)

#im9=m.plot(dms2012_x, dms2012_y, 'm', marker='*', markersize=10, zorder=5)

#im11=m.plot(dms2014_x, dms2014_y, 'r', marker='*', markersize=10, zorder=5)


plts = im1+im2+im3+im4+im5+im6
varnames=['2009', '2010', '2011', '2012', '2013', '2014']

leg = ax1.legend(plts, varnames, loc=1, ncol=1,columnspacing=0.8, labelspacing=0.5,handletextpad=0.1, borderaxespad=0.05,bbox_to_anchor=(1.185, 0.75), frameon=False)
llines = leg.get_lines()
setp(llines, linewidth=2.0)
leg.set_zorder(20)

m.drawparallels(np.arange(90,-90,-10), labels=[False,False,True,False], fontsize=8,linewidth = 0.25, zorder=5)
m.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)
m.fillcontinents(color='white',lake_color='white', zorder=5)
m.drawcoastlines(linewidth=.5, zorder=10)


cax = fig.add_axes([0.88, 0.93, 0.09, 0.025])
cbar = colorbar(im0,cax=cax, orientation='horizontal', use_gridspec=True)

cbar.set_ticks([0.25, 0.75])
cbar.set_ticklabels(['FY', 'MY'])
cbar.set_clim(0., 1.)

ax2 = subplot(2, 1, 2)

im10 = m.pcolormesh(xptsW , yptsW, wind_speed , cmap=cm.YlOrRd,vmin=minval, vmax=maxval,shading='gouraud', zorder=1)
# LOWER THE SCALE THE LARGER THE ARROW
Q = m.quiver(xptsW[::res, ::res], yptsW[::res, ::res], xvel[::res, ::res], yvel[::res, ::res], linewidths=(0.5,),color='k', edgecolors=('k'), pivot='mid',units='inches',scale=scale_vec, width = 0.01, zorder=7)
xS, yS = m(0, 68.5)
qk = quiverkey(Q, xS, yS, vector_val, str(vector_val)+r' m s$^{-1}$', fontproperties={'size': 'medium'}, coordinates='data', zorder = 11)   
 
m.drawparallels(np.arange(90,-90,-10),linewidth = 0.25, zorder=5)
m.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)
m.drawmeridians(np.arange(-180.,180.,60.),labels=[False,False,False,True], fontsize=8, linewidth = 0.25, zorder=5)


m.fillcontinents(color='white',lake_color='white', zorder=5)
m.drawcoastlines(linewidth=.5, zorder=5)
#ax1.set_ylim(m.ymax*0.1,m.ymax*0.55) 

cax1 = fig.add_axes([0.88, 0.05, 0.035, 0.3])
cbar1 = colorbar(im10,cax=cax1, orientation='vertical', extend='both', use_gridspec=True)
cbar1.set_ticks(np.linspace(minval, maxval, 6))
cbar1.set_label('Wind speed (m/s)')

subplots_adjust( bottom=0.04, top=0.95, left = 0.01, hspace=0.05)
savefig(figpath+'figure2.png', dpi=300)
close(fig)





