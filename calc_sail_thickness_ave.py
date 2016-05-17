############################################################## 
# Date: 20/01/16
# Name: calc_sail_thickness_ave.py
# Author: Alek Petty
# Description: Script to calculate sail height and IB ice thickness (10 km ave)
# Input requirements: indy results and IB thickness data

import matplotlib
matplotlib.use("AGG")
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import mpl_toolkits.basemap.pyproj as pyproj
from pylab import *
import IB_functions as ro
import numpy.ma as ma
from matplotlib import rc
from glob import glob
import os
from scipy import stats
from scipy.spatial import cKDTree as KDTree
from scipy.interpolate import griddata
from scipy import stats
import scipy.optimize as optimization

rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9
rcParams['legend.fontsize']=9
rcParams['font.size']=9
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})


def kdtree_clean1D(xx1d, yy1d, xS, yS, elevation1d, grid_res):
	#REMOVE DODGY ADDED DATA FROM THE REGRIDDING BASED ON KDTREE. 
	# dist is how far away the nearest neighbours are. 
	# need to decide on this threshold.
	# ONLY DO THIS FOR POINTS THAT HAVE ALREADY BEEN CLASSIFIED AS RIDGES
	grid_points = np.c_[xx1d, yy1d]
	tree = KDTree(np.c_[xS, yS])
	dist, _ = tree.query(grid_points, k=1)
	#dist = dist.reshape(xx2d.shape)
	elevation1d[dist > grid_res]=np.nan
	#elevation2d_KD=ma.masked_where(dist > grid_res, elevation2d)
	return elevation1d

def calc_mean_height_dist(year):
	mean_xptsB=[]
	mean_yptsB=[]
	mean_heightsR=[]
	sect_dists=[]

	#READ IN YEAR DATA BULK/INDY/COVARS
	xptsBT, yptsBT, swath_areaBT, num_labelsBT, sail_areaBT, sail_heightBT, sect_numBT, numptsBT = ro.get_bulk_ridge_stats(mib, mplot, year, datapath, section=1)
	xptsRT, yptsRT, lonRT, latRT, max_heightRT, sizeRT, sect_numRT = ro.get_indy_mean_max_height(mib, mplot, year, datapath, lonlat_section=1)

	#region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsRT, yptsRT), method='nearest')
	
	#if region<2:
	#		mask = where((region_maskR<11)&(max_heightRT<10)&(max_heightRT>0.2)&(lonRT>region_lonlat[0]) & (lonRT<region_lonlat[1]) & (latRT>region_lonlat[2]) & (latRT<region_lonlat[3]))
	#else:
	#	mask = where((region_maskR<11)&(max_heightRT<10)&(max_heightRT>0.2))
	
	#xptsRT=xptsRT[mask]
	#yptsRT=yptsRT[mask]
	#max_heightRT=max_heightRT[mask]
	#sizeRT=sizeRT[mask]
	#sect_numRT=sect_numRT[mask]
	sections=[]
	gd_sections=[]
	gd_distance=[]
	gd_ridges=[]
	for i in xrange(np.amin(sect_numBT), np.amax(sect_numBT)-10, num_kms):
		sections.append(1)
		sect1 = i
		sect2 = sect1+num_kms
		#BULK
		sect_idxsB = where((sect_numBT>=sect1) & (sect_numBT<sect2))[0]
		if (size(sect_idxsB)>=num_gd_sections):
			gd_sections.append(1)
			sect_dist = sqrt((xptsBT[sect_idxsB[-1]]-xptsBT[sect_idxsB[0]])**2 + (yptsBT[sect_idxsB[-1]]-yptsBT[sect_idxsB[0]])**2)
			sect_dists.append(sect_dist)
			
			if (sect_dist<num_kms*1500):
				gd_distance.append(1)
				#sails_areaB = sum(sail_areaBT[sect_idxsB])
				ice_areaB = sum(swath_areaBT[sect_idxsB])
				#INDY
				sect_idxsI = where((sect_numRT>=sect1) & (sect_numRT<sect2))[0]
				#DECIDE WHAT TO DO ABOUT NUMBER OF RIDGES NEEDED
				#print size(sect_idxsI)
				if (size(sect_idxsI)>=2):
					gd_ridges.append(1)
					mean_xptsB.append(mean(xptsBT[sect_idxsB]))
					mean_yptsB.append(mean(yptsBT[sect_idxsB]))
					mean_heightR = mean(max_heightRT[sect_idxsI])
					mean_heightsR.append(mean_heightR)
				#else:
				#	mean_heightsR.append(np.nan)
		
					
	return array(mean_xptsB), array(mean_yptsB), array(mean_heightsR), gd_sections, gd_distance, gd_ridges, sections,

def running_mean(xpts, ypts, var, window):
	mean_var=[]
	mean_xpts=[]
	mean_ypts=[]
	i=0
	j=0
	while j<size(xpts):
		#print i
		dist=0
		while dist<window:
			j+=1
			if j==size(xpts):
				break
			dist = sqrt((xpts[j]-xpts[i])**2 + (ypts[j]-ypts[i])**2)
		mean_var.append(np.mean(var[i:j]))
		mean_xpts.append(np.mean(xpts[i:j]))
		mean_ypts.append(np.mean(ypts[i:j]))
		i=j+0

	return array(mean_xpts), array(mean_ypts), array(mean_var)

def func(x, b):
    return b*(x**0.5)

def get_sail_thickness_ave(year, addsnow):
	mean_heightALL=[]
	thicknessALL=[]
	xpts_matchALL=[]
	ypts_matchALL=[]

	
	print year
	#mean_xpts, mean_ypts, mean_height = running_mean(xptsT, yptsT, max_heightT, window)
	mean_xpts, mean_ypts, mean_height,gd_sections, gd_distance, gd_ridges, sections = calc_mean_height_dist(year)

	mean_height.dump(savepath+'/Sailheightave_10km_'+str(year)+region_str+'nomatch.txt')
	mean_xpts.dump(savepath+'/xpts_10km_'+str(year)+region_str+'nomatch.txt')
	mean_ypts.dump(savepath+'/ypts_10km_'+str(year)+region_str+'nomatch.txt') 


	xptsIB, yptsIB, latsIB,lonsIB,thicknessIB, snowIB= ro.read_icebridgeALL(mplot, rawdatapath, year, mask_hi=1, mask_nonarctic=1)

	mean_xptsIB, mean_yptsIB, mean_thicknessIB = running_mean(xptsIB, yptsIB, thicknessIB, window)
	
	mean_thicknessIBR = griddata((mean_xptsIB, mean_yptsIB),mean_thicknessIB, (mean_xpts, mean_ypts), method='nearest', rescale=True)
	print 'size IB ave:', size(mean_thicknessIB)
	print 'size IB ave:', size(mean_thicknessIBR)
	mean_thicknessIBR = kdtree_clean1D(mean_xpts, mean_ypts, mean_xptsIB, mean_yptsIB, mean_thicknessIBR, kdtree_radius)
	#REMOVE LOW gridded thickness VALS
	mean_thicknessIBR[where(mean_thicknessIBR<0)] = np.nan 
	mean_xpts[where(mean_thicknessIBR<0)] = np.nan 
	mean_ypts[where(mean_thicknessIBR<0)] = np.nan 
	
	#FIND WHERE BOTH HAVE DATA
	match_data = where(~(np.isnan(mean_thicknessIBR+mean_height)))
	xpts_match = mean_xpts[match_data]
	ypts_match = mean_ypts[match_data]
	mean_thicknessIBR_match=mean_thicknessIBR[match_data]
	
	mean_height_match=mean_height[match_data]

	xpts_match.dump(savepath+'/xpts_10km_'+str(year)+region_str+'.txt') 
	ypts_match.dump(savepath+'/ypts_10km_'+str(year)+region_str+'.txt') 

	mean_height_match.dump(savepath+'/Sailheightave_10km_'+str(year)+region_str+'.txt') 

	mean_thicknessIBR_match.dump(savepath+'/IBthickness_10km_'+str(year)+region_str+'.txt') 


	return mean_thicknessIB, mean_thicknessIBR, mean_thicknessIBR_match, mean_height, mean_height_match,gd_sections, gd_distance, gd_ridges, sections


mplot = Basemap(projection='npstere',boundinglat=66,lon_0=0, resolution='l'  )
#-------------- IB Projection ------------------
mib=pyproj.Proj("+init=EPSG:3413")

thresh=20
fadd=''
ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd

figpath = './Figures/'
datapath = './Data_output/'+ftype+'/'
rawdatapath='../../../DATA/'
savepath=datapath+'/CORRELATIONS/INDY/'
if not os.path.exists(savepath):
	os.makedirs(savepath)

num_kms=10
num_gd_sections=3
window=num_kms*1000
kdtree_radius=5000
#--------------------------------------------------

region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)

region =2
if (region==0): 
	region_lonlat = [-150, 10, 81, 90]
	region_str='CA'
if (region==1): 
	region_lonlat = [-170, -120, 69, 79]
	region_str='BC'
if (region==2): 
	region_str=''

year=2015

mean_thicknessIB, mean_thicknessIBR, mean_thicknessIBR_match, mean_height, mean_height_match,gd_sections, gd_distance, gd_ridges, sections = get_sail_thickness_ave(year, addsnow=0)

print size(sections), size(gd_sections),size(gd_distance),size(gd_ridges)








