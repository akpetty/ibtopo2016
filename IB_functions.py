############################################################## 
# Date: 10/05/16
# Name: Alek Petty
# Description: Functions called by the various ibtopo Python scripts
# Extra: May not need all these libraries, so check before using if having issues.

from osgeo import osr, gdal
import numpy as np
from pyproj import Proj
from glob import glob
from pylab import *
from scipy import ndimage
from matplotlib import rc
from scipy.interpolate import griddata as griddatascipy
#from matplotlib.mlab import griddata
import time
import scipy.interpolate
import h5py
from scipy.spatial import cKDTree as KDTree
import os
import mpl_toolkits.basemap.pyproj as pyproj
from matplotlib.path import Path
from scipy import stats
import matplotlib.patches as patches
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from skimage.feature import peak_local_max
from skimage.morphology import watershed
from natsort import natsorted

rc("ytick",labelsize=10)
rc("xtick",labelsize=10)
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def get_bulk_ridge_stats(mib, mplot, year, ridge_stats_path, out_extra=0, lonlat=0, statsdata_compress=0, num_labels=0, areavol=0, section=0, section_vol_cor=0, areavollonlat=0, areavollonlat_big=0, getlevpercent=0):
	
	files = glob(ridge_stats_path+str(year)+'/bulk_ridge_stats_*.txt')
	# MAY NOT NEED THIS LINE NOW
	files=natsorted(files)
	print year, size(files)
	if (size(files)==0):
		print 'No data this year'

	data_bulk=np.array([]).reshape(0,18)
	for i in xrange(size(files)):
		dataT = load(files[i])
		data_bulk=vstack([data_bulk, dataT])
	#COMPRESS TO DATA WHERE THE PLANE WAS GOOD AND THERE WERE ENOUGH ATM POINTS IN THE SECTION
	good_mask = where((data_bulk[:, 16]==1)&(data_bulk[:, 7]<10)&(data_bulk[:, 6]<10))[0]
	data_bulk_good = data_bulk[good_mask]

	xptsT =data_bulk_good[:, 0]
	yptsT = data_bulk_good[:, 1]
	ice_areaT = data_bulk_good[:, 2]
	num_labelsT=data_bulk_good[:, 3]
	ridge_area_allT=data_bulk_good[:, 4]
	ridge_area_bigT=data_bulk_good[:, 5]
	ridge_height_allT=data_bulk_good[:, 6]
	ridge_height_bigT=data_bulk_good[:, 7]
	mean_altT=data_bulk_good[:, 8]
	mean_pitchT=data_bulk_good[:, 9]
	mean_rollT=data_bulk_good[:, 10]
	mean_velT=data_bulk_good[:, 11]
	num_ptsT=data_bulk_good[:, 12]
	levpercentT=data_bulk_good[:, 13]
	sectionnumT=data_bulk_good[:, 14]
	ridgefoundT=data_bulk_good[:, 15]
	pointsgoodT=data_bulk_good[:, 16]
	planegoodT=data_bulk_good[:, 17]

	lons, lats = mib(xptsT, yptsT, inverse=True)
	xptsP, yptsP = mplot(lons, lats)

	if (areavol==1):
		area_fracT = ridge_area_allT/ice_areaT
		volumeT = area_fracT*ridge_height_allT
		return array(xptsP), array(yptsP), array(area_fracT), array(volumeT)

	if (areavollonlat==1):
		area_fracT = ridge_area_allT/ice_areaT
		volumeT = area_fracT*ridge_height_allT
		return array(xptsP), array(yptsP), array(lons), array(lats), array(area_fracT), array(ridge_height_allT), array(volumeT)

	if (areavollonlat_big==1):
		area_fracT = ridge_area_bigT/ice_areaT
		volumeT = area_fracT*ridge_height_bigT
		return array(xptsP), array(yptsP), array(lons), array(lats), array(area_fracT), array(ridge_height_bigT), array(volumeT)

	if (getlevpercent==1):
		return array(xptsP), array(yptsP), array(levpercentT)

	if (num_labels==1):
		return array(xptsP), array(yptsP), array(num_labelsT)

	if (out_extra==1):
		return array(xptsP), array(yptsP), array(ice_areaT), array(num_labelsT), \
		array(ridge_area_allT), array(ridge_area_bigT), array(ridge_height_allT), \
		array(ridge_height_bigT), array(mean_altT), array(mean_pitchT), array(mean_rollT), \
		array(mean_velT), array(num_ptsT), array(ridgefoundT), array(pointsgoodT), array(planegoodT)
	
	if (section==1):
		return array(xptsP), array(yptsP), array(ice_areaT), array(num_labelsT), \
		array(ridge_area_bigT), array(ridge_height_bigT), array(sectionnumT).astype(int), array(num_ptsT)
	
	if (section_vol_cor==1):
		area_fracT = ridge_area_allT/ice_areaT
		volumeT = area_fracT*ridge_height_allT
		return array(xptsP), array(yptsP), array(volumeT), array(sectionnumT).astype(int)
	
	elif (lonlat==1):
		return array(xptsP), array(yptsP), array(lons), array(lats), array(ice_areaT), \
		array(num_labelsT), array(ridge_area_allT), array(ridge_area_bigT), array(ridge_height_allT), \
		array(ridge_height_bigT), array(ridgefoundT), array(pointsgoodT), array(planegoodT)
	
	else:
		return array(xptsP), array(yptsP), array(ice_areaT), array(num_labelsT), array(ridge_area_allT), \
		array(ridge_area_bigT), array(ridge_height_allT), array(ridge_height_bigT), array(ridgefoundT), \
		array(pointsgoodT), array(planegoodT)

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def return_xpts_ypts_fowler(path, m):
    lon_lat = loadtxt(path+'north_x_y_lat_lon.txt')
    lons = np.zeros((361, 361))
    lats = np.zeros((361, 361))
    lons = np.reshape(lon_lat[:, 3], (361, 361))
    lats = np.reshape(lon_lat[:, 2], (361, 361))
    xpts, ypts = m(lons, lats)
    return xpts, ypts, lons, lats

def get_era_winds(m, rawdatapath, start_year, end_year, start_month,  end_month, xyres):

	wind_data = Dataset(rawdatapath+'/WINDS/ERA/DATA/ERAI_WINDS_MONTHLY_1979-2014.nc', 'r')
	lats = wind_data.variables['latitude'][::-1]
	lons = wind_data.variables['longitude'][:]
	time = wind_data.variables['time'][:]/(24.*30.)
	time = time-time[0]
	time=time.reshape(time.shape[0]/12,12)
	u10 = wind_data.variables['u10'][:, ::-1, :]
	v10 = wind_data.variables['v10'][:, ::-1, :]
	u10=u10.reshape(u10.shape[0]/12, 12,u10.shape[1],u10.shape[2])
	v10=v10.reshape(v10.shape[0]/12, 12, v10.shape[1],v10.shape[2])

	u10_winter_mean= np.mean(u10[start_year-1979:end_year-1979+1,start_month:end_month+1], axis=tuple(range(0, 2)))
	v10_winter_mean= np.mean(v10[start_year-1979:end_year-1979+1,start_month:end_month+1], axis=tuple(range(0, 2)))


	u10_winter_meanS, lonsS = shiftgrid(180.,u10_winter_mean,lons,start=False)
	v10_winter_meanS, lonsS = shiftgrid(180.,v10_winter_mean,lons,start=False)

	u10_winter_meanSC, lonsSC = addcyclic(u10_winter_meanS, lonsS)
	v10_winter_meanSC, lonsSC = addcyclic(v10_winter_meanS, lonsS)

	xyres=100

	xvel,yvel,xptsW,yptsW = m.transform_vector(u10_winter_meanSC,v10_winter_meanSC,lonsSC,lats,xyres,xyres,returnxy=True,masked=True)
	wind_speed = sqrt((xvel**2) + (yvel**2))


	wind_speed = sqrt((xvel**2) + (yvel**2))
	return xptsW, yptsW, xvel, yvel, wind_speed

def get_box_xy(m, lonlat):
    lats = np.zeros((40))
    lats[0:10] = np.linspace(lonlat[3], lonlat[2], 10)
    lats[10:20] = np.linspace(lonlat[2], lonlat[2], 10)
    lats[20:30] = np.linspace(lonlat[2], lonlat[3], 10)
    lats[30:40] = np.linspace(lonlat[3], lonlat[3], 10)
    lons = np.zeros((40))
    lons[0:10] = np.linspace(lonlat[1], lonlat[1], 10)
    lons[10:20] = np.linspace(lonlat[1], lonlat[0], 10)
    lons[20:30] = np.linspace(lonlat[0], lonlat[0], 10)
    lons[30:40] = np.linspace(lonlat[0], lonlat[1], 10)
    xpts, ypts = m(lons, lats)

    return xpts, ypts

def perceptual_colormap(name, datapath, reverse=0):
    #cubeYF
    #cube1
    LinL = np.loadtxt(datapath+'{name}_0-1.csv'.format(name=name),delimiter=",")
    if (reverse==1):
        LinL=np.flipud(LinL)
    #blue 
    b3=LinL[:,2] # value of blue at sample n
    b2=LinL[:,2] # value of blue at sample n
    b1=np.linspace(0,1,len(b2)) # position of sample n - ranges from 0 to 1
    #green
    g3=LinL[:,1]
    g2=LinL[:,1]
    g1=np.linspace(0,1,len(g2))
    #red
    r3=LinL[:,0]
    r2=LinL[:,0]
    r1=np.linspace(0,1,len(r2))
    # creating list
    R=zip(r1,r2,r3)
    G=zip(g1,g2,g3)
    B=zip(b1,b2,b3)
    # transposing list
    RGB=zip(R,G,B)
    rgb=zip(*RGB)
    # creating dictionary
    k=['red', 'green', 'blue']
    data=dict(zip(k,rgb)) # makes a dictionary from 2 lists
    my_cmap = mpl.colors.LinearSegmentedColormap(name,data)
    return my_cmap

def get_mean_ice_type(mplot, rawdatapath, year, res=1):
	ice_type_path = rawdatapath+'/ICETYPE/OSISAF/GR/'+str(year)+'/'
	files = glob(ice_type_path+'*.nc')
	f = Dataset(files[0], 'r')
	lats = f.variables['lat'][::res, ::res]
	lons = f.variables['lon'][::res, ::res]
	xpts_type, ypts_type = mplot(lons, lats)
	#ice_type = np.zeros((6, lats.shape[0], lats.shape[1]))
	ice_type_days=[]
	for day in xrange(size(files)):
		f = Dataset(files[day], 'r')
		ice_typeT = f.variables['ice_type'][0,::res, ::res]

		#open water
		ice_typeT = np.where(ice_typeT==1, 0., ice_typeT)
		#first year ice
		ice_typeT = np.where(ice_typeT==2, 0.5, ice_typeT)
		#multi-year ice
		ice_typeT = np.where(ice_typeT==3, 1, ice_typeT)
		#ambiguous
		ice_typeT = np.where(ice_typeT==4, 0.75, ice_typeT)
		#set everything else to FY ice
		#ice_typeT = np.where((ice_typeT==-1), 0.5, ice_typeT)

		#EXTRAPOLATE FIRST YEAR (BC) AND MULTI_YEAR (CA) ALONG THE COAST
		#assume CA is MY and BC is FY
		#CA
		region_lonlatCA = [-150, 10, 81, 85]
		#BC
		region_lonlatBC = [-170, -120, 69, 79]
		ice_typeT = where((ice_typeT==-1)& (lons>region_lonlatBC[0]) & (lons<region_lonlatBC[1]) & (lats>region_lonlatBC[2]) & (lats<region_lonlatBC[3]), 0.5, ice_typeT )
		ice_typeT = where((ice_typeT==-1)& (lons>region_lonlatCA[0]) & (lons<region_lonlatCA[1]) & (lats>region_lonlatCA[2]) & (lats<region_lonlatCA[3]), 1., ice_typeT )
		
		ice_type_days.append(ice_typeT)
	ice_type = mean(ice_type_days, axis=0)
	return ice_type, xpts_type, ypts_type

def get_region_mask(datapath, mplot):
	header = 300
	datatype='uint8'
	file_mask = datapath+'/OTHER/region_n.msk'

	#8 - Arctic Ocean
	#9 - Canadian Archipelago
	#10 - Gulf of St Lawrence
	#11 - Land

	fd = open(file_mask, 'rb')
	region_mask = fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask[header:], [448, 304])

	mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
	mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
	lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	xpts, ypts = mplot(lons_mask, lats_mask)

	return region_mask, xpts, ypts

def calc_icebridge_flights(year, rawdatapath, loc):
    init_path = rawdatapath+'/ICEBRIDGE/POSAV/SEA_ICE/'+loc+'/'+str(year)+'_'+loc+'_NASA/'

    files = glob(init_path+'*.txt')
    lons=[]
    lats=[]

    for x in xrange(size(files)):
        data = loadtxt(files[x], skiprows=1)
        lats_t = data[::10, 1]
        lons_t = data[::10, 2]
        for x in xrange(size(lats_t)):
            lats.append(lats_t[x])
            lons.append(lons_t[x])

    return lons, lats

def calc_icebridge_flights_years(m,start_year,end_year, loc):
	xpts_all=[]
	ypts_all=[]
	for year in xrange(start_year, end_year+1, 1):
		print year
		lonsT, latsT = calc_icebridge_flights(year, 'GR')
		xptsT, yptsT=m(lonsT, latsT)
    	xpts_all.append(xptsT)
    	ypts_all.append(yptsT)
	return xpts_all, ypts_all

def get_atm_stats(year, atm_stats_path):
	
	#files = glob()
	files = glob(atm_stats_path+str(year)+'/atm_stats*.txt')
	print year, size(files)
	if (size(files)==0):
		print 'No data this year'

	min_spacingT=[]
	mean_spacingT=[]
	max_spacingT=[]
	nearmax_spacingT=[]

	for i in xrange(size(files)):
		#print i
		dataT = loadtxt(files[i])
		if (dataT.shape[0]>=5):
	    # [mean(xS), mean(yS), ice_area, size(label_nums_ma)-1, ridge_area_all, ridge_area_big, mean_ridge_height_all, mean_ridge_height_big]
			dataT=ma.masked_where(dataT==-999, dataT)


			min_spacingT.extend(dataT[:, 0])
			mean_spacingT.extend(dataT[:, 1])
			nearmax_spacingT.extend(dataT[:, 2])
			max_spacingT.extend(dataT[:, 3])

	return min_spacingT, mean_spacingT, nearmax_spacingT, max_spacingT


def read_icebridgeALL(mplot, rawdatapath, year, mask_hi=1, mask_nonarctic=1):
	lats_total=[] 
	lons_total=[]
	thickness_total=[]
	snow_thickness_total=[]
	if (year>2013):
		files = glob(rawdatapath+'/ICEBRIDGE/ICEBRIDGE_HI/QLOOK/'+str(year)+'*/*.txt')
	else:
		files = glob(rawdatapath+'/ICEBRIDGE/ICEBRIDGE_HI/IDCSI4/'+str(year)+'.*/*.txt')
	
	for x in xrange(size(files)):
		data = genfromtxt(files[x], delimiter=',', skip_header=1, dtype=str)
		# data is a table-like structure (a numpy recarray) in which you can access columns and rows easily
		lats = data[:, 0].astype(float)
		lons = data[:, 1].astype(float)
		thickness = data[:, 2].astype(float)
		snow_thickness = data[:, 7].astype(float)
		lats_total.extend(lats)
		lons_total.extend(lons)
		thickness_total.extend(thickness)
		snow_thickness_total.extend(snow_thickness)

	thickness_total=array(thickness_total)
	snow_thickness_total=array(snow_thickness_total)
	lats_total=array(lats_total)
	lons_total=array(lons_total)

	if (mask_hi==1):
		good_data=where((thickness_total>=0.)&(thickness_total<=20.))
		thickness_total = thickness_total[good_data]
		snow_thickness_total=snow_thickness_total[good_data]
		lats_total = lats_total[good_data]
		lons_total = lons_total[good_data]
	if (mask_nonarctic==1):
		region_mask, xptsM, yptsM = get_region_mask(rawdatapath, mplot)
		xptsIB, yptsIB = mplot(lons_total, lats_total)
		region_maskR = griddatascipy((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsIB, yptsIB), method='nearest')
		good_data = where((region_maskR==8))
		lats_total = lats_total[good_data]
		lons_total=lons_total[good_data]
		thickness_total=thickness_total[good_data]
		snow_thickness_total=snow_thickness_total[good_data]

	xpts,ypts = mplot(lons_total, lats_total)

	return xpts,ypts, lats_total, lons_total, thickness_total, snow_thickness_total

def get_indy_mean_max_height(mib, mplot, year, ridge_stats_path, lonlat=0, lonlat_section=0):

	files = glob(ridge_stats_path+str(year)+'/ridge_stats_*.txt')
	files=natsorted(files)
	if (size(files)==0):
		print 'No data this year'
	xpts=[]
	ypts=[]
	height=[]
	nearmax_height=[]
	max_height=[]
	sizeR=[]
	section_num=[]
	print 'Number files', size(files)
	for i in xrange(size(files)):
		#print i, size(files)
		#check file exists
		#if (os.stat(files[i]).st_size>0):
			#data = pd.read_csv(files[i], sep=' ', header=None)
		data = load(files[i])
		
		#check there is something in the data
		if (data.shape[0]>=1):
			#print size(data.shape)
			#check if just one row - extend doesn't like this for some reason so need to append
			# if (data.shape[0]==1):
			# 	print 'append:', i
			# 	xpts.append(data[0])
			# 	ypts.append(data[1])
			# 	height.append(data[4])
			# 	nearmax_height.append(data[5])
			# 	max_height.append(data[6])
			# else:
			xpts.extend(data[:,0])
			ypts.extend(data[:,1])
			height.extend(data[:,4])
			nearmax_height.extend(data[:,5])
			max_height.extend(data[:,6])
			sizeR.extend(data[:,7])
			section_num.extend(data[:,8])
			#else:
				#print 'size=0',i, data.shape[0]

	print 'num ridges', size(xpts)

	lons, lats = mib(xpts, ypts, inverse=True)

	xpts, ypts = mplot(lons, lats)
	if (lonlat==1):
		return array(xpts), array(ypts), array(lons), array(lats), array(height), array(nearmax_height), array(max_height), array(sizeR), array(section_num).astype(int)
	elif (lonlat_section==1):
		return array(xpts), array(ypts), array(lons), array(lats),array(max_height), array(sizeR), array(section_num).astype(int)
	else:
		return array(xpts), array(ypts), array(height), array(nearmax_height),array(max_height), array(sizeR)
		
def get_indy_covar(year, ridge_stats_path):

	files = glob(ridge_stats_path+str(year)+'/cov_matrix_*.txt')
	files=natsorted(files)
	if (size(files)==0):
		print 'No data this year'
	c1=[]
	c2=[]
	c3=[]
	c4=[]
	section_num=[]
	print 'Number files', size(files)
	for i in xrange(size(files)):
		data = load(files[i])

		#check there is something in the data
		if (data.shape[0]>=1):

			c1.extend(data[:,0])
			c2.extend(data[:,1])
			c3.extend(data[:,2])
			c4.extend(data[:,3])
			section_num.extend(data[:,4])

	return array(c1), array(c2), array(c3), array(c4), array(section_num).astype(int)

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
    return b*x**0.5

def get_indy_orientation(mib, mplot, year, ridge_stats_path, lonlat=0):
	
	files = glob(ridge_stats_path+str(year)+'/ridge_stats_*.txt')
	files=natsorted(files)
	if (size(files)==0):
		print 'No data this year'
	xpts=[]
	ypts=[]
	orientation=[]
	for i in xrange(size(files)):
		#print i, size(files)
		data = load(files[i])
		if (data.shape[0]>=1):
			#print size(data.shape)
			if (size(data.shape)==1):
				xpts.append(data[0])
				ypts.append(data[1])
				orientation.append(data[6])
			else:
				xpts.extend(data[:, 0])
				ypts.extend(data[:, 1])
				orientation.extend(data[:, 6])
	lons, lats = mib(xpts, ypts, inverse=True)

	xpts, ypts = mplot(lons, lats)
	if (lonlat==0):
		return array(xpts), array(ypts), array(orientation)
	else:
		return array(xpts), array(ypts), array(lons), array(lats), array(orientation)

def get_ibcao():
	bathy_file = Dataset('/Users/apetty/NOAA/IBCAO_V3_30arcsec_RR.grd','r')
	#bathy_file.variables()

	bathy = bathy_file.variables['z'][::10, ::20]
	lon_m = bathy_file.variables['x'][::20]
	lat_m = bathy_file.variables['y'][::10]
	xpts_m, ypts_m = m(*np.meshgrid(lon_m, lat_m))
	return xpts_m, ypts_m, lat_m, bathy

def get_ice_typeALL(res, mplot):
	ice_type_path = '/Users/apetty/NOAA/DATA/ICE_TYPE/OSISAF/GR/'
	files = glob(ice_type_path+'*.nc')
	f = Dataset(files[0], 'r')
	lats = f.variables['lat'][::res, ::res]
	lons = f.variables['lon'][::res, ::res]
	xpts_type, ypts_type = mplot(lons, lats)
	ice_type = np.zeros((6, lats.shape[0], lats.shape[1]))
	for x in xrange(6):
	    f = Dataset(files[x], 'r')
	    ice_type[x] = f.variables['ice_type'][0,::res, ::res]

	ice_type = np.where(ice_type<1.5, 0, ice_type)
	ice_type = np.where(ice_type==2, 0.25, ice_type)
	ice_type = np.where(ice_type==3, 0.5, ice_type)
	ice_type = np.where(ice_type==4, 0.375, ice_type)
	ice_type = ma.masked_where(ice_type<0.2, ice_type)
	return ice_type, xpts_type, ypts_type

def shot_spacing(x, y, out_min_only=0):
	min_spacing=[]

	for i in xrange(size(x)):
		xnew=np.delete(x, i)
		ynew=np.delete(y, i)
		dist = hypot(x[i]-xnew, y[i]-ynew)
		closest = np.amin(dist)
		#print closest
		min_spacing.append(closest)
	if (out_min_only==1):
		return min_spacing
	else:
		return min_spacing, amin(min_spacing), mean(min_spacing), np.percentile(min_spacing, 99), amax(min_spacing)

def get_atmmerged(atm_file, year, res, utc_time=1):

		atm = h5py.File(atm_file, 'r')
		elevation = atm['elevation'][::res]
		lon = atm['longitude'][:]
		lat = atm['latitude'][:]
		#pitch = atm['pitch'][:]
		#roll = atm['roll'][:]
		#azi =atm['azimuth'][:]
		#multiply by 1000 to put milliseconds in front of the decimal.
		gps_time=atm['time'][:]*1000
		atm.close()
		if (utc_time==1):
			utc_time = gpshhmmss_to_utc_seconds(gps_time, year, nooffset=1)
			return lon, lat, elevation, utc_time
		else:
			return lon, lat, elevation

def get_atmqih5(atm_file, year, res, utc_time=1):

	if (year>=2013):

		atm = h5py.File(atm_file, 'r')
		elevation = atm['elevation'][::res]
		lon = atm['longitude'][:]
		lat = atm['latitude'][:]
		pitch = atm['instrument_parameters']['pitch'][:]
		roll = atm['instrument_parameters']['roll'][:]
		azi =atm['instrument_parameters']['azimuth'][:]
		#multiply by 1000 to put milliseconds in front of the decimal.
		gps_time=atm['instrument_parameters']['time_hhmmss'][:]*1000
		atm.close()
		if (utc_time==1):
			utc_time = gpshhmmss_to_utc_seconds(gps_time, year)
			return lon, lat, elevation, utc_time
		else:
			return lon, lat, elevation

	else:
		fid = open(atm_file, 'rb')
		if (year<=2010):
			#BIG ENDIAN
			dtypeT='>i4'
		else:
			dtypeT='<i4'

		#header = np.fromfile(fid, dtype='>u4', count=3888)
		numofwords = np.fromfile(fid, dtype=dtypeT, count=1)/4
		blankarray = np.fromfile(fid, dtype=dtypeT, count=numofwords[0]-1)
		initialword = np.fromfile(fid, dtype=dtypeT, count=1)
		skipBytes = np.fromfile(fid, dtype=dtypeT, count=1)
		print skipBytes[0]
		if (skipBytes[0]>20000.):
			if (year==2009):
				skipBytes=[2832]
			elif (year==2010):
				skipBytes=[2928]
			elif (year==2011):
				skipBytes=[3888]
			elif (year==2012):
				skipBytes=[4176]

		fid.seek(0)
		fid.seek(skipBytes[0], os.SEEK_SET)

		data = np.fromfile(fid, dtype=dtypeT)
		data = data.reshape(-1, 12)
		atm=np.zeros((data.shape))
		atm[:, 0] = data[:, 0]/1000.
		atm[:, 1] = data[:, 1]/1000000.
		atm[:, 2] = data[:, 2]/1000000.
		atm[:, 3] = data[:, 3]/1000.
		atm[:, 4] = data[:, 4]
		atm[:, 5] = data[:, 5]
		atm[:, 6] = data[:, 6]/1000.
		atm[:, 7] = data[:, 7]/1000.
		atm[:, 8] = data[:, 8]/1000.
		atm[:, 9] = data[:, 9]/10.
		atm[:, 10] = data[:, 10]
		atm[:, 11] = data[:, 11]

		lat = atm[:, 1]
		lon = atm[:, 2]
		elevation = atm[:, 3]
		pitch = atm[:, 7]
		roll = atm[:, 8]
		gps_time = atm[:, 11]
		azi = atm[:, 6]
		#pulse_s = data[:, 4]
		#ref_s = data[:, 5]
		#azi = data[:, 6]/1000.
		#pdop = data[:, 9]/10.
		#p_width = data[:, 10]

		if (utc_time==1):
			utc_time = gpshhmmss_to_utc_seconds(gps_time, year)
			return lon, lat, elevation, utc_time
		else:
			return lon, lat, elevation

def get_atm_files(atm_path_date, year):
	if (year<=2012):
		return glob(atm_path_date+'*.qi')
	else:
		return glob(atm_path_date+'*.h5')

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

def get_dms_proj():
	#-------------- GET DMS Projection ------------------
	dms_path = '/Users/apetty/NOAA/ICEBRIDGE/DMS/2011/'
	dms_date = '20110323'
	dms_time = '17440152'
	image_path = glob(dms_path+'*'+dms_date+'_'+dms_time+'.tif')
	xT, yT, dms, geo = get_dms(image_path[0])
	spatialRef= osr.SpatialReference()
	spatialRef.ImportFromWkt(geo.GetProjectionRef())
	ProjDMS = Proj(spatialRef.ExportToProj4())
	return Proj(spatialRef.ExportToProj4())

def gpshhmmss_to_utc_seconds(gps_timeT, year, nooffset=0):
	
	if (year>2014):
		gps_utc_time=17
	elif (year==2014)|(year==2013)|(year==2012):
		gps_utc_time=16
	elif (year==2011) | (year==2010):
		gps_utc_time=15
	elif (year==2009):
		gps_utc_time=15
		#2009 is out by 15 seconds for some reason? Perhaps they changed the time twice??
		#ADDED AN OFFSET TO THE POSAV DATA INSTEAD.
	if (nooffset==1):
		gps_utc_time=0
	utc_timeT=[]
	print gps_utc_time
	for i in xrange(size(gps_timeT)):
		gps_time_str = "%09d" % int(gps_timeT[i])
		#gps_time_str = str(gps_timeT[i])
		utc_timeT.append((int(gps_time_str[0:2])*60*60)+(int(gps_time_str[2:4])*60)+(int(gps_time_str[4:6]))+(int(gps_time_str[6:9]))/1000.-gps_utc_time)

	return array(utc_timeT)

def grid_atm(xS, yS, xy_res):
	xxS = np.arange(np.amin(xS),np.amax(xS), xy_res)
	yyS = np.arange(np.amin(yS),np.amax(yS), xy_res)
	
	xx2d, yy2d = meshgrid(xxS, yyS)

	return xx2d, yy2d

def calc_level_ice(elevationS, pint, pwidth, min_ridge_height, lev_bounds=0):

	difs = [np.percentile(elevationS, i+ pwidth)-np.percentile(elevationS, i) for i in range(0, int(100)-pwidth, pint)]
	difs2= diff(array(difs))
	difnegs = where(difs2<0)[0]
	if (size(difnegs)>0):
		min_index = difnegs[-1]+1
	else:
		min_index =where(difs==min(difs))[0][0]
	level_elev = np.percentile(elevationS, (pint*min_index)+(pwidth/2))
	level_elevl = np.percentile(elevationS, (pint*min_index))
	level_elevu = np.percentile(elevationS, (pint*min_index)+pwidth)
	print 'Level ice elevation:', level_elev
	thresh = level_elev+min_ridge_height
	if (lev_bounds==1):
		return level_elev, level_elevl, level_elevu, min_index, thresh
	else:
		return level_elev, thresh, (pint*min_index)+(pwidth/2)

def calc_mean_spacing(xT, yT):
	#p = sqrt((xT-xT[0])**2 + (yT-yT[0])**2)
	pt_diffs = hypot(diff(xT), diff(yT))
	#pt_diffs = abs(np.diff(p))

	return mean(pt_diffs)

def grid_elevation(xS, yS,elevationS, xx2d, yy2d, thresh, kdtree=0):
	#INTERPOLATE ATM DATA ONTO THE 1m GRID>
	#MASK WHERE LESS THAN THE RIDGE ICE THRESHOLD

	#GIVE A BIT OF LEEWAY FOR THE INTERPOLATION
	#elevationS_Rnan =np.copy(elevationS)
	#elevationS_Rnan[elevationS<thresh-leeway]=np.nan

	elevation2d = griddatascipy((xS, yS), elevationS, (xx2d, yy2d), method='linear')
	if (kdtree==1):
		elevation2d = kdtree_clean(xx2d, yy2d, xS, yS, elevation2d)
	
	elevation2d_ma=ma.masked_where(isnan(elevation2d), elevation2d)

	elevation2d_ridge_ma=ma.masked_where(elevation2d_ma<thresh, elevation2d_ma)
	ridge_area = ma.count(elevation2d_ridge_ma)
	return elevation2d_ma, elevation2d_ridge_ma, ridge_area

def get_pos_sections(posAV, m, along_track_res):
	latp = posAV[:, 1]
	lonp = posAV[:, 2]
	xp, yp = m(lonp, latp)
	dist = hypot(diff(xp), diff(yp))
	cum_dist = cumsum(dist)
	km_indices=[]
	for seg_length in xrange(along_track_res, int(cum_dist[-1]), along_track_res):
		km_indices.append(np.abs(cum_dist - seg_length).argmin())
	km_utc_times = posAV[km_indices, 0]

	return xp, yp, cum_dist, km_indices, km_utc_times

def kdtree_clean(xx2d, yy2d, xS, yS, elevation2d):
	#REMOVE DODGY ADDED DATA FROM THE REGRIDDING BASED ON KDTREE. 
	# dist is how far away the nearest neighbours are. 
	# need to decide on this threshold.
	# ONLY DO THIS FOR POINTS THAT HAVE ALREADY BEEN CLASSIFIED AS RIDGES
	grid_points = np.c_[xx2d.ravel(), yy2d.ravel()]
	tree = KDTree(np.c_[xS, yS])
	dist, _ = tree.query(grid_points, k=1)
	dist = dist.reshape(xx2d.shape)
	elevation2d_KD=ma.masked_where(dist > 4, elevation2d)
	return elevation2d_KD

def get_labelsBIGNOWATER(elevation2d_ridge_ma, num_points_req):
	#GET BOOLEEN ARRAY WHERE TRUE IS VALID RIDGE DATA.

	elevation2d_ridge_mask = ~ma.getmask(elevation2d_ridge_ma)
	#SQURE OR DIAGONAL STRUCTIRE FOR LABELLING
	struct = ndimage.generate_binary_structure(2, 2)
	#scan across array and label connected componenets
	# returns a labelled array (integers for different labels) and the number of labels
	label_im, nb_labels = ndimage.label(elevation2d_ridge_mask, structure=struct)
	#find the unique label numbers and how many of each label are in the labelled array
	label_nums, label_num_counts = unique(label_im, return_counts=True)
	#make a note of the label numbers where there are less than a certain number of points in that label
	#labels_small = label_nums[where(label_num_counts<num_points_req)]

	label_imBIG = np.copy(label_im)
	mask_size = label_num_counts < num_points_req
	remove_pixel = mask_size[label_imBIG]

	label_imBIG[remove_pixel] = 0
	labels = np.unique(label_imBIG)
	label_imBIG = np.searchsorted(labels, label_imBIG)

	#return the new compressed, masked label array and the label numbers
	return label_imBIG

def get_labels(elevation2d_ridge_ma, xy_res, min_ridge_size, min_ridge_height):
	#GET BOOLEEN ARRAY WHERE TRUE IS VALID RIDGE DATA.
	num_points_req=min_ridge_size/(xy_res**2)
	elevation2d_ridge_mask = ~ma.getmask(elevation2d_ridge_ma)
	#SQURE OR DIAGONAL STRUCTIRE FOR LABELLING
	struct = ndimage.generate_binary_structure(2, 2)
	#scan across array and label connected componenets
	# returns a labelled array (integers for different labels) and the number of labels
	label_im, nb_labels = ndimage.label(elevation2d_ridge_mask, structure=struct)
	#find the unique label numbers and how many of each label are in the labelled array
	label_nums, label_num_counts = unique(label_im, return_counts=True)
	#make a note of the label numbers where there are less than a certain number of points in that label
	#labels_small = label_nums[where(label_num_counts<num_points_req)]

	label_imBIG = np.copy(label_im)
	mask_size = label_num_counts < num_points_req
	remove_pixel = mask_size[label_imBIG]

	label_imBIG[remove_pixel] = 0
	labels = np.unique(label_imBIG)
	label_imBIG = np.searchsorted(labels, label_imBIG)

	new_labels = np.zeros((label_im.shape))
	label_num_new=0
	for label_num in xrange(1, np.amax(label_imBIG)+1):
		elevation2d_ridge_maT=ma.masked_where(label_imBIG!=label_num, elevation2d_ridge_ma)
		imageT = elevation2d_ridge_maT.filled(0) - min_ridge_height
		local_maxi = peak_local_max(imageT, threshold_rel=0.25, min_distance=25/xy_res,exclude_border=False,indices=False)
		markers = ndimage.label(local_maxi)[0]

		if (np.amax(markers>0)):
			labelsT = watershed(-elevation2d_ridge_maT, markers)
			labelsT = ma.masked_where(ma.getmask(elevation2d_ridge_maT), labelsT)
			new_labels[where(labelsT!=0)]=labelsT[where(labelsT!=0)]+int(np.amax(new_labels))
		else:
			new_labels[where(label_imBIG==label_num)]=int(np.amax(new_labels)+1)

	new_labels=new_labels.astype('int')

	#return the new compressed, masked label array and the label numbers
	return new_labels

def calc_ridge_stats(elevation2d_ridge_ma, num_ridges, label_im, xx2d, yy2d, level_elev, section_num, calc_orientation=0):
	#GET RIDGE STATS
	#first create an empty array of appropriate size (number of labels by number of stats)
	#x, y, xst, ystd, mean height, max height, orientation
	
	ridge_stats = ma.masked_all((num_ridges, 9))
	cov_matrix = ma.masked_all((num_ridges, 5))

	#make an empty gridded array to be filled with the ridge heights relative to level ice surface
	ridge_height_mesh = ma.masked_all((xx2d.shape))
	for i in xrange(1, num_ridges+1):
		#print i
		#get aray indices of each valid (big) label
		index_t = where(label_im==i)
		#height wrt to the lowest somehting percentile elevation
		ridge_height_mesh[index_t] = np.mean(elevation2d_ridge_ma[index_t]) - level_elev
		#mean x position  of ridge
		ridge_stats[i-1, 0] = mean(xx2d[index_t])
		#mean y position  of ridge
		ridge_stats[i-1, 1] = mean(yy2d[index_t])
		#mean x std of ridge points
		ridge_stats[i-1, 2] = std(xx2d[index_t])
		#mean y std of ridge points
		ridge_stats[i-1, 3] = std(yy2d[index_t])
		#mean height of ridge relative to level ice surface
		ridge_stats[i-1, 4] = np.mean(elevation2d_ridge_ma[index_t]) - level_elev
		#max (95th percentile) height of ridge relative to level ice surface
		ridge_stats[i-1, 5] =  np.percentile(elevation2d_ridge_ma[index_t], 95) - level_elev

		ridge_stats[i-1, 6] = np.amax(elevation2d_ridge_ma[index_t])-level_elev
		#only want one coordinate size for number of points!
		ridge_stats[i-1, 7] = np.size(index_t[0])
		#section number ridge belongs to
		ridge_stats[i-1, 8] = section_num
		#CALCULATE ORIENTATION OF EACH RIDGE.
		if (calc_orientation==1):
			cov_matrix[i-1, 0:4]=np.cov(xx2d[index_t], yy2d[index_t]).flatten()
			cov_matrix[i-1, 4] = section_num

	return ridge_stats, ridge_height_mesh, cov_matrix, index_t
	
def calc_bulk_stats(stats_found, num_pts_section):

	if (stats_found==1):	
		ice_area = ma.count(elevation2d)*(xy_res**2)
		ridge_area_all = ma.count(elevation2d_ridge_ma)*(xy_res**2)
		mean_ridge_height_all = np.mean(elevation2d_ridge_ma) - level_elev
		mean_ridge_heightL = np.mean(ridge_height_mesh)
		ridge_areaL = ma.count(ridge_height_mesh)*(xy_res**2)
		return [mean_x, mean_y, ice_area, size(label_numsL)-1, ridge_area_all, ridge_areaL, mean_ridge_height_all, mean_ridge_heightL, mean_alt, mean_pitch, mean_roll, num_pts_section, stats_found]
	elif (stats_found==0):
		#a = ma.masked_all((0))
		#masked_val = mean(a)
		return [mean_x, mean_y, -999, -999,-999, -999, -999, -999, mean_alt, mean_pitch, mean_roll, num_pts_section, stats_found]

def posav_section_info(m, posAV_section):
	mean_lat=mean(posAV_section[:, 1])
	mean_lon=mean(posAV_section[:, 2])
	mean_x, mean_y = m(mean_lon, mean_lat)
	mean_alt=mean(posAV_section[:, 3])
	mean_vel=mean(posAV_section[:, 4])
	mean_pitch=mean(posAV_section[:, 5])
	mean_roll=mean(posAV_section[:, 6]) 

	return mean_x, mean_y, mean_alt, mean_pitch, mean_roll, mean_vel

def get_atm_poly(xT, yT, elevationT, km_utc_times, utc_timeT, poly_path, i):
	pos_time1 = km_utc_times[i]-4
	pos_time2 = km_utc_times[i+1]+4
	xpos_rough = xT[where((utc_timeT>pos_time1)&(utc_timeT<pos_time2))]
	ypos_rough = yT[where((utc_timeT>pos_time1)&(utc_timeT<pos_time2))]
	elevation_rough = elevationT[where((utc_timeT>pos_time1)&(utc_timeT<pos_time2))]
	coords_rough = np.dstack((xpos_rough,ypos_rough))
	inpoly = poly_path.contains_points(coords_rough[0])

	xpos_km = xpos_rough[where(inpoly==True)]
	ypos_km = ypos_rough[where(inpoly==True)]
	elevation_km = elevation_rough[where(inpoly==True)]

	return xpos_km, ypos_km, elevation_km

def get_pos_poly(xp, yp, idx1, idx2):

	gradx = xp[idx1+1]-xp[idx1-1]
	grady = yp[idx1+1]-yp[idx1-1]
	dx1 = 300*cos(arctan(grady/gradx)+(pi/2))
	dy1 = 300*sin(arctan(grady/gradx)+(pi/2))

	xp1 = xp[idx1]-dx1
	xp2 = xp[idx1]+dx1
	yp1 = yp[idx1]-dy1
	yp2 = yp[idx1]+dy1

	x1y1=[xp1, yp1]
	x2y2=[xp2, yp2]

	gradx = xp[idx2+1]-xp[idx2-1]
	grady = yp[idx2+1]-yp[idx2-1]
	dx2 = 300*cos(arctan(grady/gradx)+(pi/2))
	dy2 = 300*sin(arctan(grady/gradx)+(pi/2))

	xp3 = xp[idx2]+dx2
	xp4 = xp[idx2]-dx2
	yp3 = yp[idx2]+dy2
	yp4 = yp[idx2]-dy2

	x3y3=[xp3, yp3]
	x4y4=[xp4, yp4]

	hypot1 = hypot(xp2-xp1, yp2-yp1)
	hypot2 = hypot(xp3-xp2, yp3-yp2)
	hypot3 = hypot(xp4-xp3, yp4-yp3)
	hypot4 = hypot(xp1-xp4, yp1-yp4)

	return Path([x1y1, x2y2, x3y3, x4y4, x1y1], closed=True), [x1y1, x2y2, x3y3, x4y4, x1y1], [hypot1, hypot2, hypot3, hypot4]

