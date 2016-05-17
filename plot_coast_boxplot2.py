############################################################## 
# Date: 20/01/16
# Name: plot_coast_boxplot2.py
# Author: Alek Petty
# Description: Script to plot indy coastline boxplots
# Input requirements: indy distributions (All/MYI/FYI in the CA/BC)
# Output: 2 panels highlighting the distributions in coastline bins

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
from scipy.interpolate import griddata

rcParams['font.size']=9
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9
rcParams['legend.fontsize']=9

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def get_ridge_dists_precentsINDY(region_str):
	z_c_int_mayears=[]
	max_height_mayears=[]
	z_c_int_maALL=[]
	max_height_maALL=[]

	
	for x in xrange(size(years)):
		z_c_int_maT=load(datapath+'z_c_int_ma'+region_str+str(years[x]+start_year)+'.txt')
		z_c_int_mayears.append(z_c_int_maT)
		z_c_int_maALL.extend(z_c_int_maT)
		max_height_maT=load(datapath+'max_height_ma'+region_str+str(years[x]+start_year)+'.txt') 
		max_height_mayears.append(max_height_maT)
		max_height_maALL.extend(max_height_maT)

	z_c_int_maALL=array(z_c_int_maALL)
	max_height_maALL=array(max_height_maALL)
	sort_indices = np.argsort(z_c_int_maALL, axis=0)
	z_c_int_maALL=z_c_int_maALL[sort_indices]
	max_height_maALL=max_height_maALL[sort_indices]

	max_height_maALL_sections=[]

	num_bins = int(ceil(np.amax(z_c_int_maALL)/bin_width))
	for bin_num in xrange(num_bins):
		dist_idxs=where((z_c_int_maALL>(bin_num*bin_width))&(z_c_int_maALL<((bin_num+1)*bin_width)))
		max_height_maALL_sections.append(max_height_maALL[dist_idxs])


	return max_height_maALL_sections, num_bins

def get_ridge_dists_precentsINDYYEARS(region_str,x):
	z_c_int_maT=load(datapath+'z_c_int_ma'+region_str+str(years[x]+start_year)+'.txt')
	max_height_maT=load(datapath+'max_height_ma'+region_str+str(years[x]+start_year)+'.txt') 
	
	sort_indices = np.argsort(z_c_int_maT, axis=0)
	z_c_int_maT=z_c_int_maT[sort_indices]
	max_height_maT=max_height_maT[sort_indices]

	max_height_maT_sections=[]

	for bin_num in xrange(num_bins):
		dist_idxs=where((z_c_int_maT>(bin_num*bin_width))&(z_c_int_maT<((bin_num+1)*bin_width)))
		max_height_maT_sections.append(max_height_maT[dist_idxs])


	return max_height_maT_sections, num_bins


start_year=2009
end_year=2014
num_years = end_year - start_year + 1
num_bins = 9

thresh = 20
add=''
ftype='1km_xyres2m_'+str(thresh)+'cm'+add


years=[0, 1, 2, 3, 4, 5]

datapath = './Data_output/'+ftype+'/COAST/INDY/'
figpath = './Figures/'
bin_width = 100

height_maALLCA_sections, num_binsCA = get_ridge_dists_precentsINDY('CA')
height_maALLBC_sections, num_binsBC = get_ridge_dists_precentsINDY('BC')

height_maALLCA_sections_years=[]
height_maALLBC_sections_years=[]
for x in xrange(num_years):
	height_maALLCA_sections_year, num_binsCA = get_ridge_dists_precentsINDYYEARS('CA', x)
	height_maALLBC_sections_year, num_binsBC = get_ridge_dists_precentsINDYYEARS('BC', x)

	height_maALLCA_sections_years.append(height_maALLCA_sections_year)
	height_maALLBC_sections_years.append(height_maALLBC_sections_year)

textwidth=3.
colors=['m', 'r', 'g', 'b', 'c', 'y', 'k']

fig = figure(figsize=(textwidth,(textwidth)))

ax1=subplot(2, 1, 1)
bp1=boxplot(height_maALLCA_sections, positions=np.arange(0.5, num_binsCA), widths=0.9, sym='', whis=[5, 95])
for x in xrange(num_years):
	for y in xrange(num_bins):
		if (size(height_maALLCA_sections_years[x][y])>0):
			vars()['p'+str(x+1)]=ax1.axhline(y=median(height_maALLCA_sections_years[x][y]), xmin=(float(y)/num_bins)+((1./9.)*0.1), xmax=((float(y)+1)/num_bins)-((1./9.)*0.1), color=colors[x])
			vars()['p'+str(x+1)]=ax1.axhline(y=np.percentile(height_maALLCA_sections_years[x][y], 95), xmin=(float(y)/num_bins)+0.01, xmax=((float(y)+1)/num_bins), color=colors[x], ls='dashed')
		#vars()['p'+str(x+1)] =ax1.plot((0.1*float(x))+float(y)+0.2, median(height_maALLCA_sections_years[x][y]), marker='x', color=colors[x])
setp(bp1['boxes'], color='black')
setp(bp1['whiskers'], color='black', ls='solid')
setp(bp1['fliers'], color='black')
setp(bp1['medians'], color='black')
ax1.set_xlim(0, 9)
ax1.set_ylim(0, 4)
ax1.set_ylabel('Feature height (m)')
ax1.set_xticks(np.arange(0, 1*(num_binsCA)+1, 1))
ax1.set_xticklabels([])
ax1.yaxis.set_major_locator(MaxNLocator(6))

ax2=subplot(2, 1, 2)
bp2=boxplot(height_maALLBC_sections, positions=np.arange(0.5, num_binsBC), widths=0.9, sym='', whis=[5, 95])
for x in xrange(num_years):
	for y in xrange(num_bins):
		if (size(height_maALLBC_sections_years[x][y])>0):
			vars()['p'+str(x+1)]=ax2.axhline(y=median(height_maALLBC_sections_years[x][y]), xmin=(float(y)/num_bins)+((1./9.)*0.1), xmax=((float(y)+1)/num_bins)-((1./9.)*0.1), color=colors[x])
			vars()['p'+str(x+1)]=ax2.axhline(y=np.percentile(height_maALLBC_sections_years[x][y], 95), xmin=(float(y)/num_bins)+0.01, xmax=((float(y)+1)/num_bins), color=colors[x], ls='dashed')
setp(bp2['boxes'], color='black')
setp(bp2['whiskers'], color='black')
setp(bp2['fliers'], color='black')
setp(bp2['medians'], color='black')
ax2.set_xlim(0, 9)
ax2.set_ylim(0, 4)
ax2.set_ylabel('Feature height (m)')
ax2.set_xticks(np.arange(0, 1*(num_binsCA)+1, 1))
ax2.set_xticklabels(np.arange(0, 1*(num_binsCA)+1, 1))
ax2.yaxis.set_major_locator(MaxNLocator(6))

ax2.yaxis.set_major_locator(MaxNLocator(6))
ax2.set_xlabel('Distance to coast (100 km)', labelpad=2)

ax1.xaxis.grid(True)
ax2.xaxis.grid(True)

ax1.annotate('(a) Central Arctic', xy=(0.95, 0.85), xycoords='axes fraction', horizontalalignment='right')
ax2.annotate('(b) Beaufort/Chukchi', xy=(0.95, 0.85), xycoords='axes fraction', horizontalalignment='right')

years = np.arange(start_year, end_year+1)
year_strs =  list(years.astype('str'))
year_strs.append('ALL')
#plts_net = p1+p2+p3+p4+p5+p6
for x in xrange(num_years+1):
	ax1.text(0.05+(x*0.13), 1.01, year_strs[x],horizontalalignment='center',verticalalignment='bottom',color=colors[x],transform=ax1.transAxes)


print 'saving...'
subplots_adjust(left = 0.15, right = 0.98, bottom = 0.12, top = 0.93, hspace=0.14, wspace=0.05)
savefig(figpath+'figure10.pdf', dpi=300)








