############################################################## 
# Date: 20/01/16
# Name: plot_spacing.py
# Author: Alek Petty
# Description: Script to plot spacing between ATM points
# Input requirements: ATM and DMS for a given section
# Output: 1 km ATM spacing

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
#rc('font', family='sans-serif') 
#rc('font', serif='Helvetica Neue') 
#rc('text', usetex='false')

def plot_section_3():

	sizex = np.amax(xatm_km) - np.amin(xatm_km)
	sizey = np.amax(yatm_km) - np.amin(yatm_km)
	ratio = sizey/sizex
	#minvalL = 0
	elevation2d_ridge_comp = ma.compressed(elevation2d_ridge_ma)
	#maxvalL = np.round(np.percentile(elevation2d_ridge_comp, 99), decimals=1)
	lowerp = 5
	upperp = 99
	minval = np.round(np.percentile(elevation_km, lowerp), decimals=1)
	maxval = np.round(np.percentile(elevation_km, upperp), decimals=1)
	minvalR = 1.0 #round(np.amin(ridge_height_mesh), 1)
	maxvalR = 1.3 #round(np.amax(ridge_height_mesh), 1)

	#dms_plot=0
	textwidth=5.
	fig = figure(figsize=(textwidth,textwidth*0.45*ratio))


	ax1 = subplot(131)

	ax1.annotate('(a) Raw ATM' , xy=(0.03, 1.03), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')

	#dx_in_points = ax1.transData.transform(1, 1)
	#print dx_in_points
	#im1 = plt.Circle(xatm_km, yatm_km, 1, c=elevation_km, edgecolor='none')
	im1 = scatter(xatm_km, yatm_km, c = elevation_km, vmin = minval, vmax = maxval, s=0.3, lw = 0, cmap = cm.RdYlBu_r, rasterized=True)
	#patch = patches.PathPatch(poly_path, facecolor='none', lw=2)
	#ax1.add_patch(patch)
	ax1.set_aspect('equal')
	ax2 = subplot(132)
	ax2.annotate('(b) Gridded ('+str(xy_res)+' m) ATM' , xy=(0.03, 1.03), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')

	elevation_rel = ma.masked_where((elevation2d-level_elev)<-0.1, (elevation2d-level_elev))

	im2 = pcolormesh(xx2d, yy2d, elevation2d-level_elev, vmin=-(maxval-level_elev), vmax = maxval-level_elev, cmap = cm.RdYlBu_r)
	#ax2.add_patch(patch)
	ax2.set_aspect('equal')
	ax3 = subplot(133)
	ax3.annotate('(c) Unique ridges' , xy=(0.03, 1.03), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
	#ax3.add_patch(patch)

	ax3.set_aspect('equal')


	if (np.amax(label_im)>=1):
		minvalLAB=0
		maxvalLAB=np.amax(label_im)
		label_im_ma=ma.masked_where(label_im<0.5, label_im)
		im3 = pcolormesh(xx2d, yy2d, label_im_ma, vmin = minvalLAB, vmax = maxvalLAB, cmap = my_cmap)
		
		#im3 = pcolormesh(xx2d, yy2d, ridge_height_mesh, vmin = minvalR, vmax = maxvalR, cmap = cm.YlOrRd)
		im31 = plot(ridge_stats[:, 0], ridge_stats[:, 1], marker='o', markersize=1, linestyle = 'None', color='k')

		#for i in xrange(ridge_stats.shape[0]):
		#	im32 = plot([ridge_stats[i, 0]-2*ridge_stats[i, 2], ridge_stats[i, 0]+2*ridge_stats[i, 2]], [ridge_stats[i, 1]-2*ridge_stats[i, 3], ridge_stats[i, 1]+2*ridge_stats[i, 3]], marker='None', linewidth=0.5, linestyle = '-', color='k')

	axesname = ['ax1', 'ax2', 'ax3']
	for plotnum in xrange(3):
		vars()[axesname[plotnum]].set_xlim(np.amin(xatm_km),np.amax(xatm_km))
		vars()[axesname[plotnum]].set_ylim(np.amin(yatm_km),np.amax(yatm_km))
	
		vars()[axesname[plotnum]].xaxis.set_major_locator(MaxNLocator(5))
	
		vars()[axesname[plotnum]].yaxis.grid(True)
		vars()[axesname[plotnum]].xaxis.grid(True)

	ax1.set_ylabel('False northing (m)', labelpad=1)
	ax2.set_xlabel('False easting (m)', labelpad=1)
	ax2.set_yticklabels([])
	ax3.set_yticklabels([])
	#fig.text(0.3, 0.98, 'DMS Date: '+date+'  DMS Time: '+dms_time)

	cax = fig.add_axes([0.1, 0.13, 0.2, 0.03])
	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both', use_gridspec=True)
	cbar.set_label('Elevation to WGS84 (m)', labelpad=2)
	xticks1 = np.linspace(minval, maxval, 3)
	cbar.set_ticks(xticks1)

	cax2 = fig.add_axes([0.4, 0.13, 0.2, 0.03])
	cbar2 = colorbar(im2,cax=cax2, orientation='horizontal', extend='both', use_gridspec=True)
	cbar2.set_label('Elevation to level ice (m)', labelpad=2)
	xticks2 = np.linspace(-(maxval-level_elev), maxval-level_elev, 3)
	cbar2.set_ticks(xticks2)


	if (np.amax(label_im)>=1):
		cax3 = fig.add_axes([0.73, 0.13, 0.2, 0.03])
		cbar3 = colorbar(im3,cax=cax3, orientation='horizontal', extend='both', use_gridspec=True)
		cbar3.set_label('Sail number', labelpad=0)
		xticks2 = np.linspace(minvalLAB, maxvalLAB, 2)
		cbar3.set_ticks(xticks2)
		cbar3.solids.set_rasterized(True)

	cbar.solids.set_rasterized(True)
	cbar3.solids.set_rasterized(True)

	ax3.annotate('# sails:'+str(num_ridges)+'\nSail area:'+ridge_areaL+r' m$^{2}$', xy=(0.03, 0.02), textcoords='axes fraction', color='k', horizontalalignment='left', verticalalignment='bottom')
	

	
	
	#plt.tight_layout()
	print 'Saving figure...'
	subplots_adjust(left=0.09, top = 0.98, right=0.98, bottom=0.23, wspace=0.05)

	savefig(figpath+'/3PLOTS/3plot_'+str(year)+str(days)+str(atm_file)+'_'+str(section)+'_'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm'+nstr+'.png', dpi=300)

def plot_spacing():

	sizex = np.amax(xatm_km) - np.amin(xatm_km)
	sizey = np.amax(yatm_km) - np.amin(yatm_km)
	ratio = sizey/sizex

	minval = 0.
	maxval = 4.
	#dms_plot=0
	textwidth = 3.
	fig = figure(figsize=(textwidth,textwidth*ratio))
	subplots_adjust(left=0.14, bottom=0.1, top=0.98)

	ax1 = subplot(111)

	im1 = scatter(xatm_km, yatm_km, c = spacing, vmin = minval, vmax = maxval, s=1.0, lw = 0, cmap = cm.YlOrRd, rasterized=True)

	ax1.set_aspect('equal')

	axesname = ['ax1']
	ax1.set_xlim(np.amin(xatm_km),np.amax(xatm_km))
	ax1.set_ylim(np.amin(yatm_km),np.amax(yatm_km))
	ax1.xaxis.set_major_locator(MaxNLocator(5))

	ax1.set_xlabel('False easting (m)', labelpad=1)

	ax1.set_ylabel('False northing (m)', labelpad=1)

	ax1.annotate('Altitude:'+mean_altSTR+'m\nVelocity:'+mean_velSTR+r' m s$^{-1}$' , xy=(0.5, 0.95), textcoords='axes fraction', color='k', horizontalalignment='left', verticalalignment='top')
	ax1.annotate('Mean:'+meanSSTR+'m\np(99%):'+nearmaxSSTR+'m\nMax:'+maxSSTR+'m' , xy=(0.03, 0.3), textcoords='axes fraction', color='k', horizontalalignment='left', verticalalignment='bottom')

	#fig.text(0.3, 0.98, 'DMS Date: '+date+'  DMS Time: '+dms_time)

	cax = fig.add_axes([0.23, 0.27, 0.3, 0.03])
	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='max', use_gridspec=True)
	cbar.set_label('Shot spacing (m)', labelpad=1)
	xticks1 = np.linspace(minval, maxval, 3)
	cbar.set_ticks(xticks1)
	cbar.solids.set_rasterized(True)


	savefig(figpath+'figureS5.png', dpi=1000)

def calc_bulk_stats(stats_found, num_pts_section):
	
	if (stats_found==1):	
		ice_area = ma.count(elevation2d)*(xy_res**2)
		ridge_area_all = ma.count(elevation2d_ridge_ma)*(xy_res**2)
		mean_ridge_height_all = np.mean(elevation2d_ridge_ma) - level_elev
		mean_ridge_heightL = np.mean(ridge_height_mesh)
		ridge_areaL = ma.count(ridge_height_mesh)*(xy_res**2)
		return [mean_x, mean_y, ice_area, num_ridges, ridge_area_all, ridge_areaL, mean_ridge_height_all, mean_ridge_heightL, mean_alt, mean_pitch, mean_roll, mean_vel, num_pts_section, stats_found]
	elif (stats_found==0):
		#a = ma.masked_all((0))
		#masked_val = mean(a)
		return [mean_x, mean_y, -999, 0,-999, -999, -999, -999, mean_alt, mean_pitch, mean_roll, mean_vel, num_pts_section, stats_found]

#-------------- ATM AND DMS PATHS------------------

datapath='./Data_output/'
rawdatapath = '../../../DATA/'
IBrawdatapath = rawdatapath+'/ICEBRIDGE/'
ATM_path = IBrawdatapath+'/ATM/ARCTIC/'
dms_path = IBrawdatapath+'/DMS/'
posAV_path =IBrawdatapath+'/POSAV/SEA_ICE/GR/'
figpath='./Figures/'

my_cmap=ro.perceptual_colormap("Linear_L", rawdatapath+'/OTHER/CMAPS/', reverse=1)


m=pyproj.Proj("+init=EPSG:3413")

calc_atm_stats=1

min_ridge_height = 0.2
pint=5
pwidth = 20
along_track_res=1000
pts_threshold=18000*(along_track_res/1000)
xy_res=2
num_points_req = 100/(xy_res**2)


year = 2011
days = 0
atm_file = 1
section=843


narrow=0
if (narrow==1):
	nstr='N'
	ATM_path = ATM_path+'/NARROW/'
else:
	nstr=''

ATM_year = ATM_path+str(year)+'/'
atm_files_year = glob(ATM_year+'*/')

atm_path_date = atm_files_year[days]
print 'ATM day:', atm_path_date

if (narrow==1):
	atm_files_in_day = glob(atm_path_date+'*.h5')
else:
	atm_files_in_day = ro.get_atm_files(atm_path_date, year)

#load POS file
posAV = loadtxt(posAV_path+str(year)+'_GR_NASA/sbet_'+str(atm_path_date[-9:-1])+'.out.txt', skiprows=1)
#GET POSITION OF PLANE AND 1km MARKERS FROM POSAV
xp, yp, dist, km_idxs, km_utc_times = ro.get_pos_sections(posAV, m, along_track_res)

print 'ATM file:', atm_files_in_day[atm_file], str(atm_file)+'/'+str(size(atm_files_in_day))

if (narrow==1):
	lonT, latT, elevationT, utc_timeT= ro.get_atmmerged(atm_files_in_day[atm_file], year, 1)
else:
	lonT, latT, elevationT, utc_timeT= ro.get_atmqih5(atm_files_in_day[atm_file], year, 1)


#IF SIZE OF DATA IS LESS THAN SOME THRESHOLD THEN DONT BOTHER ANALYZING
xT, yT = m(lonT, latT)
#GET POSAV INDICES COINCIDING WITH START AND END OF ATM FILE. ADD PLUS/MINUS 1 FOR SOME LEEWAY.
start_i = np.abs(km_utc_times - utc_timeT[0]).argmin()
end_i = np.abs(km_utc_times - utc_timeT[-1]).argmin()

for section in xrange(section, section+1):
	print section
	mean_x, mean_y, mean_alt, mean_pitch, mean_roll, mean_vel = ro.posav_section_info(m, posAV[km_idxs[section]:km_idxs[section+1]]	)
	print 'Mean altitude:', mean_alt
	print 'Mean pitch:', mean_pitch
	print 'Mean roll:', mean_roll

	if (abs(mean_alt-500)<200) & (abs(mean_pitch)<5) & (abs(mean_roll)<5):
		poly_path, vertices, sides = ro.get_pos_poly(xp, yp, km_idxs[section], km_idxs[section+1])

		xatm_km, yatm_km, elevation_km = ro.get_atm_poly(xT, yT, elevationT, km_utc_times, utc_timeT, poly_path, section)
		if (size(xatm_km)>0):
			lonDMS, latDMS = m(xatm_km[0], yatm_km[0], inverse=True)
			lonDMS_str = '%.2f' %lonDMS
			latDMS_str = '%.2f' %latDMS
			print lonDMS_str, latDMS_str
			xatm_km = xatm_km - np.amin(xatm_km)
			yatm_km = yatm_km - np.amin(yatm_km)
			spacing, minS, meanS, nearmaxS, maxS = ro.shot_spacing(xatm_km, yatm_km)

			mean_altSTR=' %.0f' % mean_alt
			mean_pitchSTR=' %.0f' % mean_pitch
			mean_rollSTR=' %.0f' % mean_roll
			mean_velSTR=' %.0f' % mean_vel
			meanSSTR = ' %.2f' % meanS
			nearmaxSSTR = ' %.2f' % nearmaxS
			maxSSTR = ' %.2f' % maxS

			plot_spacing()

