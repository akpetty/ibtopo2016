############################################################## 
# Date: 20/01/16
# Name: plot_6distributions.py
# Author: Alek Petty
# Description: Script to plot indy stats distributions
# Input requirements: indy distributions (All/MYI/FYI in the CA/BC)
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
#import cubehelix
from netCDF4 import Dataset
from matplotlib import rc
from glob import glob
import os

rcParams['axes.labelsize'] =9
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9
rcParams['legend.fontsize']=9
rcParams['font.size']=9
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

# Plotting projection
mplot = Basemap(projection='npstere',boundinglat=68,lon_0=0, resolution='l'  )
# IB projection
mib=pyproj.Proj("+init=EPSG:3413")

thresh=20
fadd=''
ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd

figpath = './Figures/'
datapath = './Data_output/'+ftype+'/DISTS'

start_year=2009
end_year=2014
num_years = end_year - start_year + 1
years = np.arange(start_year, end_year+1)
xmax=3
xmin=0
bin_width=0.1

bins = load(datapath+'/binsT_r'+str(0)+'_t'+str(0)+'_'+ftype+str(2009)+str(bin_width)+'.txt') 
bins = bins[:-1]
histALL=np.zeros((2, 3, num_years, size(bins)))
statsALL=np.zeros((2, 3, num_years+1, 4))

h_type = 0

for r in xrange(2):
	for t in xrange(3):
		print r
		statsALL[r, t] = loadtxt(datapath+'/statsALL_r'+str(r)+'_t'+str(t)+'_'+ftype+str(bin_width)+'.txt', skiprows=1, delimiter='&') 
		for year in xrange(start_year, end_year+1):
			
			histT = load(datapath+'/histT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(year)+str(bin_width)+'.txt') 
			#convert to probability density
			histT = 100*histT/(float(sum(histT)))
			histALL[r, t, year-start_year] = histT

colors=['k', 'r', 'g', 'b', 'c', 'y']

#region_strs=[, ]
textwidth=5.

fig = figure(figsize=(textwidth,textwidth*0.9))
plotnum=0

for t in xrange(3):
	for r in xrange(2):
		plotnum+=1
		vars()['ax'+str(plotnum)] = subplot(3, 2, plotnum)
		axT=gca()
		for x in xrange(num_years):
			vars()['p'+str(x+1)] = plot(bins+(bin_width/2.),histALL[r, t, x],color=colors[x],ls='-',lw=1.)
	
			axvline(statsALL[r, t, x, 0],ymin=0, ymax=1,color=colors[x],ls='-', lw=0.5)
			vars()['p1'+str(x+1)] = axvline(statsALL[r,t, x, 2],ymin=0, ymax=1,color=colors[x], ls='--', lw=0.5)
		xlim(xmin,xmax)
		ylim(0, int(ceil(np.amax(histALL))))
		axT.yaxis.set_major_locator(MaxNLocator(7))
		axT.xaxis.set_major_locator(MaxNLocator(6))
		if (plotnum<5):
			axT.set_xticklabels([])
#ylabel('Probablity (%)')
#
ax2.set_yticklabels([])
ax4.set_yticklabels([])
ax6.set_yticklabels([])
#ax.set_xticklabels([])
ax5.set_xlabel('Feature height (m)', labelpad=3)
ax6.set_xlabel('Feature height (m)', labelpad=3)

ax1.set_ylabel('Probablity (%)')
ax3.set_ylabel('Probablity (%)')
ax5.set_ylabel('Probablity (%)')
	

ax1.annotate('Central Arctic (CA)', xy=(0.25, 1.05), xycoords='axes fraction', horizontalalignment='middle')
ax2.annotate('Beaufort/Chukchi (BC)', xy=(0.25, 1.05), xycoords='axes fraction', horizontalalignment='middle')

ax1.annotate('(a) All\nCA', xy=(0.95, 0.75), xycoords='axes fraction', horizontalalignment='right')
ax2.annotate('(b) All\nBC', xy=(0.95, 0.75), xycoords='axes fraction', horizontalalignment='right')
ax3.annotate('(c) FYI\nCA', xy=(0.95, 0.75), xycoords='axes fraction', horizontalalignment='right')
ax4.annotate('(d) FYI\nBC', xy=(0.95, 0.75), xycoords='axes fraction', horizontalalignment='right')
ax5.annotate('(e) MYI\nCA', xy=(0.95, 0.75), xycoords='axes fraction', horizontalalignment='right')
ax6.annotate('(f) MYI\nBC', xy=(0.95, 0.75), xycoords='axes fraction', horizontalalignment='right')


years.astype('str')
plts_net = p1+p2+p3+p4+p5+p6
leg = ax2.legend(plts_net, years, loc=1, ncol=2,columnspacing=0.3, handletextpad=0.0001,bbox_to_anchor=(1.02, 0.8), frameon=False)
leg.get_frame().set_alpha(0.5)
leg.get_frame().set_linewidth(2)

subplots_adjust(left = 0.09, right = 0.97, bottom = 0.08, top = 0.95, hspace=0.09, wspace = 0.1)

savefig(figpath+'/figure6.pdf', dpi=300)









