############################################################## 
# Date: 20/01/16
# Name: calc_coast_dist_bulk.py
# Author: Alek Petty
# Description: Script to calculate distributions of bulk topography stats (e.g. mean feature volume) as a function of coastline proximity
# Input requirements: Topography stats across all years
# Output: Distributions of bulk statistics (e.g. mean feature volunme)
# Notes - need to repeat for two regions

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
from operator import truediv, mul
from netCDF4 import Dataset
import h5py

def get_ridge_data_sections(region):
	ice_volumeALL=[]
	xpts=[]
	ypts=[]

	if (region==0): 
		region_lonlat = [-150, 10, 81, 90]
		region_str='CA'
	if (region==1): 
		region_lonlat = [-170, -120, 69, 79]
		region_str='BC'

	for year in xrange(start_year, end_year+1):
		print year
		xptsT, yptsT, lonT, latT, sail_area_fracT, heightT, ice_volumeT= ro.get_bulk_ridge_stats(mib, mplot, year, datapath, areavollonlat=1)

		region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)
		region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsT, yptsT), method='nearest')

		mask = where((lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
		#mask = mask[0].tolist()
		lonT=lonT[mask]
		latT=latT[mask]
		xptsT=xptsT[mask]
		yptsT=yptsT[mask]
		ice_volumeT=ice_volumeT[mask]

		print year, size(mask)
		xpts.append(xptsT)
		ypts.append(yptsT)
		ice_volumeALL.append(ice_volumeT)

	return xpts, ypts, ice_volumeALL

#-------------- GET DMS Projection ------------------
mplot=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)
mib=pyproj.Proj("+init=EPSG:3413")

region=0
if (region==0):
	region_str='CA'
elif (region==1):
	region_str='BC'
	
thresh=20
fadd=''
ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd

#region=int(sys.argv[2])
print 'Ftype:', ftype
print 'Region:',region

rawdatapath = '../../DATA/'
datapath='./Data_output/'+ftype+'/'
outpath= datapath+'COAST/BULK/'
if not os.path.exists(outpath):
	os.makedirs(outpath)

start_year=2009
end_year=2014
num_years = end_year - start_year + 1

coast_file = h5py.File('../../DATA/OTHER/distance-to-coast_2m.nc','r')
lat_c = coast_file['lat'][:]
lon_c = coast_file['lon'][:]
z_c = coast_file['z'][:]

#REMOVE SOME OF THE COAST DISTANCE DATA THAT ISN'T NEEDED
lat_c=lat_c[4250:-1]
lon_c=lon_c[::20]
z_c=z_c[4250:-1, ::20]

xpts_c, ypts_c = mplot(*np.meshgrid(lon_c, lat_c))
xpts, ypts, sail_area_frac = get_ridge_data_sections(region)

for x in xrange(num_years):
	
	z_c_intT = griddata((xpts_c.flatten(), ypts_c.flatten()), z_c.flatten(),(xpts[x],ypts[x]),method='linear')
	mask_z=where(z_c_intT>0.)
	z_c_int_maT=z_c_intT[mask_z]
	z_c_int_maT.dump(outpath+'z_c_int_ma'+region_str+str(x+start_year)+'.txt') 
	sail_areaT=sail_area_frac[x][mask_z]
	sail_areaT.dump(outpath+'volume_ma'+region_str+str(x+start_year)+'.txt') 
	











