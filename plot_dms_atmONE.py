############################################################## 
# Date: 20/01/16
# Name: plot_dms_atmONE.py
# Author: Alek Petty
# Description: Script to plot a detection case study
# Input requirements: DMS image, posAV and ATM data for specific case studies
# Output: 6 panels highlighting the detection technique

import matplotlib
matplotlib.use("AGG")
import IB_functions as ro
import matplotlib.colors as colors
# basemap import
from mpl_toolkits.basemap import Basemap, shiftgrid
from skimage import measure 
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
from scipy.interpolate import griddata
from matplotlib import rc
from netCDF4 import Dataset
from glob import glob
import os
from osgeo import osr, gdal
import mpl_toolkits.basemap.pyproj as pyproj 
from scipy import ndimage
import matplotlib.tri as tri
import scipy.interpolate
import time
import h5py
from scipy.spatial import cKDTree as KDTree

rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9
rcParams['legend.fontsize']=9
rcParams['font.size']=9
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})


def get_atm_from_dms(res):
	#FIND ATM FILE THAT STARTED BEFORE THE DMS GPS TIME
	#- SHOULD MAINLY WORK BUT COULD BE PROBLEMS DUE TO SLIGHT OFFSET
	atm_files = glob(ATM_path+str(year)+'/'+date+'/*')
	#print size(atm_files)
	for j in xrange(size(atm_files)):
		atm_time = float(atm_files[j][71:77])
		atm_time = float(atm_files[j][71:77])
		#print atm_time
		#print float(dms_time[0:-2])
		if (atm_time>float(dms_time[0:-2])):
			break

	return atm_files[j-1]

def get_dms(image_path):
	geo = gdal.Open(image_path) 
	band1 = geo.GetRasterBand(1)
	band2 = geo.GetRasterBand(2)
	band3 = geo.GetRasterBand(3)
	red = band1.ReadAsArray()
	green = band2.ReadAsArray()
	blue = band3.ReadAsArray()

	dms = (0.299*red + 0.587*green + 0.114*blue)
	dms = ma.masked_where(dms<1, dms)

	trans = geo.GetGeoTransform()
	width = geo.RasterXSize
	height = geo.RasterYSize

	x1 = np.linspace(trans[0], trans[0] + width*trans[1] + height*trans[2], width)
	y1 = np.linspace(trans[3], trans[3] + width*trans[4] + height*trans[5], height)
	xT, yT = meshgrid(x1, y1)
	return xT, yT, dms, geo

def plot_6plot_thresh(dms_plot=0):
	
	minvalD = 0
	maxvalD = 255
	sizex = np.amax(xT) - np.amin(xT)
	sizey = np.amax(yT) - np.amin(yT)
	ratio = sizey/sizex
	minvalL = 0
	if (found_ridges==1):
		elevation2d_ridge_comp = ma.compressed(elevation2d_ridge_maL)
		maxvalL = np.round(np.percentile(elevation2d_ridge_comp, 99), decimals=1)
	else:
		maxvalL=min_ridge_height
	lowerp = 5
	upperp = 99.5
	minval = np.round(np.percentile(elevation_ma, lowerp), decimals=1)
	maxval = np.round(np.percentile(elevation_ma, upperp), decimals=1)
	res = 2
	#dms_plot=0
	textwidth=5.5
	fig = figure(figsize=(textwidth,textwidth*1.25*ratio))

	ax1 = subplot(321)
	ax1.annotate('(a) Raw DMS' , xy=(0.03, 0.9), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')

	if (dms_plot==1):
		im1 = pcolormesh(xT[::res, ::res], yT[::res, ::res], dms[::res, ::res], vmin = minvalD, vmax = maxvalD, cmap = cm.gist_gray, rasterized=True)

	ax2 = subplot(322)
	ax2.annotate('(b) Raw DMS + ATM' , xy=(0.03, 0.9), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')

	if (dms_plot==1):
		im2 = pcolormesh(xT[::res, ::res], yT[::res, ::res], dms[::res, ::res], vmin = minvalD, vmax = maxvalD, cmap = cm.gist_gray, rasterized=True)
	im21 = scatter(xatm, yatm, c = elevation_ma, vmin = minval, vmax = maxval, s=1, lw = 0, cmap = cm.RdYlBu_r, rasterized=True)
	
	ax3 = subplot(323)

	xe = [i for i in xrange(100)]
	ye = [np.percentile(elevation_ma, i) for i in xrange(100)]
	im3 = plot(xe, ye, 'k')
	mid_percent = (pint*min_index)+(pwidth/2)
	low_percent = (pint*min_index)
	up_percent = (pint*min_index)+pwidth
	#axhline(y=level_elev, xmin=0, xmax=mid_percent, color='b')
	axhline(y=level_elev, linestyle='--', color='b')
	axhline(thresh, color='r')
	ax3.annotate('Elevation threshold' , xy=(10, thresh), xycoords='data', color='r', horizontalalignment='left', verticalalignment='bottom')
	ax3.annotate('Level ice ('+str(low_percent)+'-'+str(up_percent)+'%)' , xy=(10, level_elev), xycoords='data', color='b', horizontalalignment='left', verticalalignment='top')
	xlim([0, 100])

	ax3.annotate('(c) Elevation distribution' , xy=(0.03, 0.9), xycoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')

	ax3.annotate(str(min_ridge_height)+' m' , xy=(mid_percent-15, level_elevl + (thresh-level_elevu)*0.5), color='k', horizontalalignment='left', verticalalignment='bottom')
	#ax4.annotate('(d)' , xy=(0.03, 0.93), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
	ax3.set_ylim(np.amin(ye), np.amax(ye))
	norm = ro.MidpointNormalize(midpoint=0)
	ax4 = subplot(324)
	ax4.annotate('(d) Gridded ('+str(xy_res)+' m) ATM' , xy=(0.03, 0.9), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
	
	im4 = pcolormesh(xx2d, yy2d, elevation2d-level_elev, norm=norm,vmin=minvalEL, vmax=maxvalEL, cmap = cm.RdBu_r)
	#
	im4.set_clim(minvalEL, maxvalEL)
	im41 = contour(xx2d, yy2d, level_ice, levels=[np.amin(level_ice), np.amax(level_ice)], colors='k', linewidths=0.5)

	#cs.set_clim(50, 210)
	ax5 = subplot(325)
	ax5.annotate('(e) High topography (>0.2 m)' , xy=(0.03, 0.9), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')

	#im3 = pcolormesh(xx2d, yy2d, level_ice, vmin = minvalL, vmax = maxvalL, cmap = cm.RdYlBu_r)
	if (found_ridges==1):
		im31 = pcolormesh(xx2d, yy2d, elevation2d_ridge_maL, vmin = minvalL, vmax = maxvalL, cmap = cm.RdYlBu_r)
		im51 = contour(xx2d, yy2d, label_im, np.arange(0, num_ridges+1), colors='k', linewidths=1)

	ax6 = subplot(326)
	ax6.annotate(r'(f) Unique features' , xy=(0.03, 0.9), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
	if (found_big_ridge==1):
		minvalLAB=0
		maxvalLAB=np.amax(label_im)
		label_im_ma=ma.masked_where(label_im<0.5, label_im)
		im6 = pcolormesh(xx2d, yy2d, label_im_ma, vmin = minvalLAB, vmax = maxvalLAB, cmap = my_cmap)
		im61 = plot(ridge_stats[:, 0], ridge_stats[:, 1], marker='o', markersize=2, linestyle = 'None', color='k')

		#for i in xrange(ridge_stats.shape[0]):
		#	im62 = plot([ridge_stats[i, 0]-2*ridge_stats[i, 2], ridge_stats[i, 0]+2*ridge_stats[i, 2]], [ridge_stats[i, 1]-2*ridge_stats[i, 3], ridge_stats[i, 1]+2*ridge_stats[i, 3]], marker='None', linestyle = '-', color='k')
	
	ax3.set_xlabel('Percentile (%)', labelpad=1)
	ax3.set_ylabel('Relative elevation (m)', labelpad=1)
	axesname = ['ax1', 'ax2', 'ax3', 'ax4', 'ax5', 'ax6']

	for plotnum in [0, 1, 3, 4, 5]:
		vars()[axesname[plotnum]].set_xlim(np.amin(xT),np.amax(xT))
		vars()[axesname[plotnum]].set_ylim(np.amin(yT),np.amax(yT))
	for plotnum in [0, 4, 5]:
		vars()[axesname[plotnum]].set_xlabel('y (m)', labelpad=1)
		
	for plotnum in [0, 3, 4]:
		vars()[axesname[plotnum]].set_ylabel('x (m)', labelpad=1)
	

	for plotnum in xrange(6):
		vars()[axesname[plotnum]].yaxis.grid(True)
		vars()[axesname[plotnum]].xaxis.grid(True)
	letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
	#for plotnum in xrange(6):
	
	fig.text(0.15, 0.98, 'DMS Date: '+date+'  DMS Time: '+dms_time+'     '+latDMS_str+'N, '+lonDMS_str+'E')

	cax = fig.add_axes([0.58, 0.71, 0.02, 0.08])
	cbar = colorbar(im21,cax=cax, orientation='vertical', extend='both', use_gridspec=True)
	cbar.set_label('Elevation to \n WGS84 (m)', labelpad=28, rotation=0)
	xticks1 = np.linspace(minval, maxval, 3)
	cbar.set_ticks(xticks1)

	cax1 = fig.add_axes([0.58, 0.4, 0.02, 0.08])
	cbar1 = colorbar(im4,cax=cax1, orientation='vertical', extend='max', use_gridspec=True)
	cbar1.set_label('Elevation to \n level ice (m)', labelpad=28, rotation=0)
	xticksL = np.linspace(minvalEL, maxvalEL, 4)
	cbar1.set_ticks(xticksL)

	cax2 = fig.add_axes([0.1, 0.08, 0.02, 0.08])
	cbar2 = colorbar(im31,cax=cax2, orientation='vertical', extend='max', use_gridspec=True)
	cbar2.set_label('Elevation to \n level ice (m)', labelpad=30, rotation=0)
	xticks2 = np.linspace(minvalL, maxvalL, 3)
	cbar2.set_ticks(xticks2)

	if (found_big_ridge==1):
		cax3 = fig.add_axes([0.58, 0.08, 0.02, 0.08])
		cbar3 = colorbar(im6,cax=cax3, orientation='vertical', extend='neither', use_gridspec=True)
		cbar3.set_label('Feature\nidentifier', labelpad=20, rotation=0)
		xticks2 = np.linspace(minvalLAB, maxvalLAB, 2)
		cbar3.set_ticks(xticks2)
		cbar3.solids.set_rasterized(True)

	cbar.solids.set_rasterized(True)
	cbar1.solids.set_rasterized(True)
	cbar2.solids.set_rasterized(True)
	
	im4.cmap.set_under('w')
	#plt.tight_layout()
	print 'Saving figure...'
	subplots_adjust(bottom=0.07, left=0.09, top = 0.96, right=0.98, hspace=0.22)

	savefig(figpath+'atm_dms6plot_'+str(xy_res)+'mxy_'+date+'_'+dms_time+int_method+str(int(min_ridge_height*100))+'_cm2.png', dpi=300)

def plot_probdist():

	fig = figure(figsize=(3,3))
	ax=gca()
	im3 = plot(bins[0:-1], hist, 'k')

	axvline(x=level_elev, color='b')
	axvline(x=level_elevl, linestyle='--', color='r')
	axvline(x=level_elevu, linestyle='--', color='r')
	ax.set_ylabel('Prob Density', labelpad=1)
	ax.set_xlabel('Relative elevation (m)', labelpad=1)
	tight_layout()
	savefig(figpath+'atm_dmsdist_'+str(xy_res)+'mxy_'+date+'_'+dms_time+int_method+str(int(min_ridge_height*100))+'_cm2.png', dpi=300)


datapath='./Data_output/'
rawdatapath = '../../../DATA/ICEBRIDGE/'
ATM_path = rawdatapath+'/ATM/ARCTIC/'
dms_path = rawdatapath+'/DMS/'
posAV_path =rawdatapath+'/POSAV/'
figpath='./Figures/'

my_cmap=ro.perceptual_colormap("cube1", '../../../DATA/OTHER/CMAPS/', reverse=1)


norm = ro.MidpointNormalize(midpoint=0)

xy_res=2
pwidth=20
pint=5
num_points_req=100/(xy_res**2)

dms_image = 7
min_ridge_height = .2

minvalEL=-1
maxvalEL=1
sh=0
if (sh==1):
	print sys.argv[1]
	dms_image = int(sys.argv[1])
	min_ridge_height = float(sys.argv[2])/100.
print 'DMS IMAGE:', dms_image
print 'THRESHOLD:', min_ridge_height

if (dms_image==1):
	#GOOD RIDGE
	date = '20110323'
	dms_time =  '17440152'
if (dms_image==2):
	date = '20100405'
	dms_time = '12253130'
if (dms_image==3):
	#BIG RIDGE
	date = '20110323'
	dms_time = '17445601'
if (dms_image==4):
	date = '20130424'
	dms_time = '13284560'
if (dms_image==5):
	#BIG RIDGE
	date = '20110323'
	dms_time = '19024646'
if (dms_image==6):
	#BIG RIDGE
	date = '20110323'
	dms_time = '18045984'
if (dms_image==7):
	#BIG RIDGE
	date = '20120322'
	dms_time = '14280792'
if (dms_image==8):
	#BIG RIDGE
	date = '20110323'
	dms_time = '17450800'
if (dms_image==9):
	#BIG RIDGE
	date = '20110323'
	dms_time = '17450551'
if (dms_image==10):
	#GOOD RIDGE
	date = '20110323'
	dms_time =  '17440950'
if (dms_image==11):
	date = '20130424'
	dms_time = '13183977'
if (dms_image==12):
	#BIG RIDGE
	date = '20110323'
	dms_time = '18235793'
if (dms_image==13):
	date = '20130424'
	dms_time = '13293109'
if (dms_image==14):
	date = '20140331'
	dms_time = '15374024'
if (dms_image==15):
	date = '20140331'
	dms_time = '15155031'
if (dms_image==16):
	#BIG RIDGE
	date = '20110415'
	dms_time = '14163295'

print date
print dms_time
int_method='linear'

year = int(date[0:4])
image_path = glob(dms_path+str(year)+'/'+date+'/*'+date+'_'+dms_time+'.tif')
xT, yT, dms, geo = get_dms(image_path[0])

#-------------- GET ATM ------------------
print 'Finding atm...'
atm_file = get_atm_from_dms(res=1)

lon, lat, elevation = ro.get_atmqih5(atm_file, year, 1, utc_time=0)
#--------------------------------------------------

#-------------- PROJECT ATM ------------------------
ProjDMS=pyproj.Proj("+init=EPSG:3413")
x, y = ProjDMS(lon, lat)
#--------------------------------------------------

lonDMS, latDMS = ProjDMS(xT[0, 0], yT[0, 0], inverse=True)
lonDMS_str = '%.2f' %lonDMS
latDMS_str = '%.2f' %latDMS


minX = xT[-1, 0]
minY = yT[-1, 0]
xT = xT - minX
yT = yT - minY
x = x - minX
y = y - minY

print lonDMS, latDMS, date[0:4]
#-------------- CHECK BOUNDS OF ATM/DMS ----------------
print 'tiff lims x:', np.amin(xT),np.amax(xT)
print 'tiff lims y:', np.amin(yT),np.amax(yT)
print 'ATM lims x:', np.amin(x),np.amax(x)
print 'ATM lims y:', np.amin(y),np.amax(y)
#--------------------------------------------------
#-------------- MASK ELEVATION AND X/Y  ------------------
elevation_maL = ma.masked_where((x<np.amin(xT)) | (x>np.amax(xT)) | (y<np.amin(yT)) | (y>np.amax(yT)), elevation)
elevation_ma = ma.compressed(elevation_maL)
mask = ma.getmask(elevation_maL)
indexT = np.where(mask==False)
xatm = x[indexT]
yatm = y[indexT]
#--------------------------------------------------

level_elev, level_elevl, level_elevu, min_index, thresh = ro.calc_level_ice(elevation_ma, pint, pwidth, min_ridge_height, lev_bounds=1)
#thresh = level_elev+min_ridge_height

hist, bins = histogram(elevation_ma, bins=100, density=True)
plot_probdist()

#-------------- CREATE 2D ATM GIRD ------------------
print 'GRID ATM...'

xx = np.arange(np.amin(xT),np.amax(xT), xy_res)
yy = np.arange(np.amin(yT),np.amax(yT), xy_res)
xx2d, yy2d = meshgrid(xx, yy)

elevation2d, elevation2d_ridge_ma, ridge_area = ro.grid_elevation(xatm, yatm,elevation_ma, xx2d, yy2d, thresh,kdtree=1)
found_ridges=0
if (ridge_area>0):
	found_ridges=1
	
level_ice = ma.masked_where((elevation2d<level_elevl) | (elevation2d>level_elevu), elevation2d)

elevation2d_ridge_maL =elevation2d_ridge_ma-level_elev
level_ice=level_ice-level_elev
#-------------- LABEL ARRAY ------------------
label_im  = ro.get_labels(elevation2d_ridge_maL, xy_res, min_ridge_height)

found_big_ridge=0
if (np.amax(label_im)>0):
	found_big_ridge=1
	num_ridges=np.amax(label_im)
	#-------------- LABEL STATS ------------------
	ridge_stats, ridge_height_mesh, cov_matrix, index = ro.calc_ridge_stats(elevation2d_ridge_ma, num_ridges, label_im, xx2d, yy2d, level_elev,0, calc_orientation=1)

#-------------- GENERAL PLOTTING VARIABLES  ------------------
print 'Plot ATM/DMS...'
plot_6plot_thresh(dms_plot=1)






