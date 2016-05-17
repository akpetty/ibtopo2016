############################################################## 
# Date: 20/01/16
# Name: calc_6distributionsBULK.py
# Author: Alek Petty
# Description: Script to calculate distributions of bulk topography stats (e.g. mean feature volume) across the MY/FY and CA/BC regions.
# Input requirements: Topography stats across all years
# Output: Distributions of bulk statistics (e.g. mean feature volunme)

import matplotlib
matplotlib.use("AGG")
# basemap import
from mpl_toolkits.basemap import Basemap, shiftgrid
# Numpy import
import numpy as np
import mpl_toolkits.basemap.pyproj as pyproj
from pylab import *
import IB_functions as ro
import numpy.ma as ma
from scipy.interpolate import griddata
import os

mplot = Basemap(projection='npstere',boundinglat=68,lon_0=0, resolution='l'  )
#rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = ['Computer Modern Roman']
#rcParams['text.usetex'] = True


def get_hist_year(region, type, year):
	hist=[]
	bins=[]

	if (region==0): 
		region_lonlat = [-150, 10, 81, 90]
		region_str='CA'
	if (region==1): 
		region_lonlat = [-170, -120, 69, 79]
		region_str='BC'

	#areavollonlat_big=1 or areavollonlat=1 to include/not include small features

	xptsT, yptsT, lonT, latT, sail_area_fracT, ridge_heightT, ice_volumeT= ro.get_bulk_ridge_stats(mib, mplot, year, datapath, areavollonlat=1)
	if (bulk_type==0):
		varT=np.copy(sail_area_fracT)
	elif (bulk_type==1):
		varT=np.copy(ice_volumeT)
	elif (bulk_type==2):
		varT=np.copy(ridge_heightT)

	region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)
	region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsT, yptsT), method='nearest')

	#ice_type, xptsA, yptsA = ro.get_ice_type_year(mplot, year-2009, res=1)
	ice_type, xptsA, yptsA = ro.get_mean_ice_type(mplot, rawdatapath, year, res=1)

	ice_typeR = griddata((xptsA.flatten(), yptsA.flatten()),ice_type.flatten(), (xptsT, yptsT), method='nearest')

	if (type==0):
		mask = where((ice_typeR<1.1) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
	if (type==1):
		mask = where((ice_typeR<0.6) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
	if (type==2):
		mask = where((ice_typeR<1.1) & (ice_typeR>0.9)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
	
	varT=varT[mask]

	
	histT, binsT = np.histogram(varT, bins=bin_vals)
	meanH = mean(varT)
	medianH = median(varT)
	stdH = std(varT)
	modeH = binsT[argmax(histT)] + bin_width/2.
	#hist.append(histT)
	#bins.append(binsT)

	return binsT, histT, meanH, medianH, modeH, stdH

def get_hist_allyears(region, type):
	hist=[]
	bins=[]

	if (region==0): 
		region_lonlat = [-150, 10, 81, 90]
		region_str='CA'
	if (region==1): 
		region_lonlat = [-170, -120, 69, 79]
		region_str='BC'

	varALL=[]
	for year in xrange(start_year, end_year+1):
		print year

		xptsT, yptsT, lonT, latT, sail_area_fracT, ridge_heightT, ice_volumeT= ro.get_bulk_ridge_stats(mib, mplot, year, datapath, areavollonlat=1)
		if (bulk_type==0):
			varT=np.copy(sail_area_fracT)
		elif (bulk_type==1):
			varT=np.copy(ice_volumeT)
		elif (bulk_type==2):
			varT=np.copy(ridge_heightT)

		region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)
		region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsT, yptsT), method='nearest')

		ice_type, xptsA, yptsA = ro.get_mean_ice_type(mplot, rawdatapath,year, res=1)
		ice_typeR = griddata((xptsA.flatten(), yptsA.flatten()),ice_type.flatten(), (xptsT, yptsT), method='nearest')

		if (type==0):
			mask = where((ice_typeR<1.1) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
		if (type==1):
			mask = where((ice_typeR<0.6) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
		if (type==2):
			mask = where((ice_typeR<1.1) & (ice_typeR>0.9)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
	
		varT=varT[mask]
		varALL.extend(varT)

	histT, binsT = np.histogram(varALL, bins=bin_vals)
	meanH = mean(varALL)
	medianH = median(varALL)
	stdH = std(varALL)
	modeH = binsT[argmax(histT)] + bin_width/2.
	#hist.append(histT)
	#bins.append(binsT)

	return binsT, histT, meanH, medianH, modeH, stdH

#--------------------------------------------------

#-------------- GET DMS Projection ------------------
mib=pyproj.Proj("+init=EPSG:3413")

bulk_type=1

if (bulk_type==0):
	bulk_str='area'
elif (bulk_type==1):
	bulk_str='vol'
elif (bulk_type==2):
	bulk_str='height'

thresh=20
fadd=''

ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd

datapath='./Data_output/'+ftype+'/'
figpath = './Figures/'
outpath= datapath+'DISTS/'
rawdatapath='../../DATA/'

if not os.path.exists(outpath):
	os.makedirs(outpath)


bin_width = 0.01
start_h=0.
print start_h
end_h=0.5
bin_vals=np.arange(start_h,end_h, bin_width)

start_year=2009
end_year=2014
num_years = end_year - start_year + 1

histALL=[]
binsALL=[]

statsT = np.zeros((num_years+1, 4))
years = np.arange(start_year, end_year+1)
years.astype('str')

statsALL = np.zeros(((num_years+1)*3, 9))

for t in xrange(3):
	for r in xrange(2):
		for year in xrange(start_year, end_year+1):
			print t, r, year
			binsT, histT, meanHT, medianHT, modeHT, stdHT = get_hist_year(r, t, year)
			binsT.dump(outpath+'binsT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(year)+bulk_str+'.txt') 
			histT.dump(outpath+'histT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(year)+bulk_str+'.txt') 
			statsT[year - start_year, 0] = meanHT
			statsT[year - start_year, 1] = stdHT
			#statsT[year - start_year, 1] = medianHT
			statsT[year - start_year, 2] = modeHT
			statsT[year - start_year, 3] = sum(histT)/1e4
			
		binsT, histT, meanHT, medianHT, modeHT, stdHT = get_hist_allyears(r, t)
		binsT.dump(outpath+'binsT_r'+str(r)+'_t'+str(t)+'_'+ftype+bulk_str+'ALLYEARS.txt') 
		histT.dump(outpath+'histT_r'+str(r)+'_t'+str(t)+'_'+ftype+bulk_str+'ALLYEARS.txt') 
		statsT[year - start_year+1,0] = meanHT
		statsT[year - start_year+1,1] = stdHT
		statsT[year - start_year+1,2] = modeHT
		statsT[year - start_year+1,3] = sum(histT)/1e4
		savetxt(outpath+'statsALL_r'+str(r)+'_t'+str(t)+'_'+ftype+bulk_str+'.txt', statsT, fmt='%.3f', header='Mean, SD, Mode, Number', delimiter='&')
		statsALL[t*(num_years+1):(t*(num_years+1))+num_years+1, 0]=np.arange(2009, 2016)
		statsALL[t*(num_years+1):(t*(num_years+1))+num_years+1,(r*4)+1:(r*4)+4+1 ]=statsT	
		#binsALL.append(binsT)
		#histALL.append(histT)

savetxt(outpath+'statsALL_CABCMYFY_'+ftype+bulk_str+'.txt', statsALL, fmt='& %.0f & %.2f (%.2f) & %.2f & %.2f & %.2f (%.2f) & %.2f & %.2f \\', header='Year, Mean (SD), Mode, Number, Mean (SD), Mode, Number')








