############################################################## 
# Date: 20/01/16
# Name: calc_6distributions.py
# Author: Alek Petty
# Description: Script to calculate distributions of individual topography stats (e.g. feature height) across the MY/FY and CA/BC regions.
# Input requirements: Topography stats across all years
# Output: Distributions of individual statistics (e.g. feature height)
# info: need to run with 0.1 and 0.3 bin width for normal and log plots

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



	xptsT, yptsT, lonT, latT, max_heightT, sizeT, sectionT = ro.get_indy_mean_max_height(mib, mplot, year, datapath, lonlat_section=1)
	
	region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)
	region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsT, yptsT), method='nearest')

	ice_type, xptsA, yptsA = ro.get_mean_ice_type(mplot, rawdatapath, year, res=1)
	#ice_type, xptsA, yptsA = ro.get_ice_type_year(mplot, year-2009, res=1)
	ice_typeR = griddata((xptsA.flatten(), yptsA.flatten()),ice_type.flatten(), (xptsT, yptsT), method='nearest')

	if (type==0):
		mask = where((ice_typeR<1.1) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8)&(max_heightT<10))
	if (type==1):
		mask = where((ice_typeR<0.6) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8)&(max_heightT<10))
	if (type==2):
		mask = where((ice_typeR<1.1) & (ice_typeR>0.9)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8)&(max_heightT<10))
	
	max_heightT=max_heightT[mask]

	
	histT, binsT = np.histogram(max_heightT, bins=bin_vals)
	meanH = mean(max_heightT)
	medianH = median(max_heightT)
	stdH = std(max_heightT)
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

	max_heightALL=[]
	for year in xrange(start_year, end_year+1):
		print year

		xptsT, yptsT, lonT, latT, max_heightT, sizeT, sectionT = ro.get_indy_mean_max_height(mib, mplot, year, datapath, lonlat_section=1)
		
		region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)
		region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsT, yptsT), method='nearest')

		#ice_type, xptsA, yptsA = ro.get_ice_type_year(mplot, year-2009, res=1)
		
		ice_type, xptsA, yptsA = ro.get_mean_ice_type(mplot, rawdatapath, year, res=1)

		ice_typeR = griddata((xptsA.flatten(), yptsA.flatten()),ice_type.flatten(), (xptsT, yptsT), method='nearest')

		if (type==0):
			mask = where((ice_typeR<1.1) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8)&(max_heightT<10))
		if (type==1):
			mask = where((ice_typeR<0.6) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8)&(max_heightT<10))
		if (type==2):
			mask = where((ice_typeR<1.1) & (ice_typeR>0.9)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8)&(max_heightT<10))
	
		max_heightT=max_heightT[mask]
		max_heightALL.extend(max_heightT)

	histT, binsT = np.histogram(max_heightALL, bins=bin_vals)
	meanH = mean(max_heightALL)
	medianH = median(max_heightALL)
	stdH = std(max_heightALL)
	modeH = binsT[argmax(histT)] + bin_width/2.
	#hist.append(histT)
	#bins.append(binsT)

	return binsT, histT, meanH, medianH, modeH, stdH

#--------------------------------------------------

#-------------- GET DMS Projection ------------------
mib=pyproj.Proj("+init=EPSG:3413")

thresh=20
fadd=''

ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd

datapath='./Data_output/'+ftype+'/'
figpath = './Figures/'
rawdatapath='../../../DATA/'
outpath= datapath+'DISTS/'

if not os.path.exists(outpath):
	os.makedirs(outpath)

bin_width = 0.3
start_h=0.
print start_h
end_h=10.0
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
			binsT.dump(outpath+'binsT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(year)+str(bin_width)+'.txt') 
			histT.dump(outpath+'histT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(year)+str(bin_width)+'.txt') 
			statsT[year - start_year, 0] = meanHT
			statsT[year - start_year, 1] = stdHT
			#statsT[year - start_year, 1] = medianHT
			statsT[year - start_year, 2] = modeHT
			statsT[year - start_year, 3] = sum(histT)/1e5
			
		binsT, histT, meanHT, medianHT, modeHT, stdHT = get_hist_allyears(r, t)
		binsT.dump(outpath+'binsT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(bin_width)+'ALLYEARS.txt') 
		histT.dump(outpath+'histT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(bin_width)+'ALLYEARS.txt') 
		statsT[year - start_year+1,0] = meanHT
		statsT[year - start_year+1,1] = stdHT
		statsT[year - start_year+1,2] = modeHT
		statsT[year - start_year+1,3] = sum(histT)/1e5
		savetxt(outpath+'statsALL_r'+str(r)+'_t'+str(t)+'_'+ftype+str(bin_width)+'.txt', statsT, fmt='%.3f', header='Mean, SD, Mode, Number', delimiter='&')
		statsALL[t*(num_years+1):(t*(num_years+1))+num_years+1, 0]=np.arange(2009, 2016)
		statsALL[t*(num_years+1):(t*(num_years+1))+num_years+1,(r*4)+1:(r*4)+4+1 ]=statsT	
		#binsALL.append(binsT)
		#histALL.append(histT)

savetxt(outpath+'statsALL_CABCMYFY_'+ftype+str(bin_width)+'.txt', statsALL, fmt='& %.0f & %.2f (%.2f) & %.2f & %.2f & %.2f (%.2f) & %.2f & %.2f \\', header='Year, Mean (SD), Mode, Number, Mean (SD), Mode, Number')













