############################################################## 
# Date: 20/01/16
# Name: plot_thickness_sail_scatter.py
# Author: Alek Petty
# Description: Script to plot scatter of sail height and IB ice thickness (10 km ave)
# Input requirements: indy results and IB thickness data
# Output: scatter plot

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


def read_icebridge_uncALL(mplot, rawdatapath, year, mask_hi=1, mask_nonarctic=1):
	lats_total=[] 
	lons_total=[]
	thickness_total=[]
	thickness_unc_total=[]
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
		thickness_unc = data[:, 3].astype(float)
		snow_thickness = data[:, 7].astype(float)
		lats_total.extend(lats)
		lons_total.extend(lons)
		thickness_total.extend(thickness)
		thickness_unc_total.extend(thickness_unc)
		snow_thickness_total.extend(snow_thickness)

	thickness_total=array(thickness_total)
	thickness_unc_total=array(thickness_unc_total)
	snow_thickness_total=array(snow_thickness_total)
	lats_total=array(lats_total)
	lons_total=array(lons_total)

	if (mask_hi==1):
		good_data=where((thickness_total>=0.)&(thickness_total<=20.))
		thickness_total = thickness_total[good_data]
		thickness_unc_total=thickness_unc_total[good_data]
		snow_thickness_total=snow_thickness_total[good_data]
		lats_total = lats_total[good_data]
		lons_total = lons_total[good_data]
	if (mask_nonarctic==1):
		xptsIB, yptsIB = mplot(lons_total, lats_total)
		region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsIB, yptsIB), method='nearest')
		good_data = where((region_maskR==8))
		lats_total = lats_total[good_data]
		lons_total=lons_total[good_data]
		thickness_total=thickness_total[good_data]
		thickness_unc_total=thickness_unc_total[good_data]
		snow_thickness_total=snow_thickness_total[good_data]

	xpts,ypts = mplot(lons_total, lats_total)

	return xpts,ypts, lats_total, lons_total, thickness_total, thickness_unc_total,snow_thickness_total

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

	region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsRT, yptsRT), method='nearest')
	
	if region<2:
		mask = where((region_maskR<11)&(max_heightRT<10)&(max_heightRT>0.2)&(lonRT>region_lonlat[0]) & (lonRT<region_lonlat[1]) & (latRT>region_lonlat[2]) & (latRT<region_lonlat[3]))
	else:
		mask = where((region_maskR<11)&(max_heightRT<10)&(max_heightRT>0.2))
	
	xptsRT=xptsRT[mask]
	yptsRT=yptsRT[mask]
	max_heightRT=max_heightRT[mask]
	sizeRT=sizeRT[mask]
	sect_numRT=sect_numRT[mask]


	for i in xrange(np.amin(sect_numBT), np.amax(sect_numBT)-10, num_kms): #sect_numBT[0], sect_numBT[-1]-num_kms, num_kms):
		sect1 = i
		sect2 = sect1+num_kms
		#BULK
		sect_idxsB = where((sect_numBT>=sect1) & (sect_numBT<sect2))[0]
		if (size(sect_idxsB)>=num_gd_sections):
			sect_dist = sqrt((xptsBT[sect_idxsB[-1]]-xptsBT[sect_idxsB[0]])**2 + (yptsBT[sect_idxsB[-1]]-yptsBT[sect_idxsB[0]])**2)
			sect_dists.append(sect_dist)
			if (sect_dist<num_kms*1200):
				mean_xptsB.append(mean(xptsBT[sect_idxsB]))
				mean_yptsB.append(mean(yptsBT[sect_idxsB]))
				#sails_areaB = sum(sail_areaBT[sect_idxsB])
				ice_areaB = sum(swath_areaBT[sect_idxsB])
				#INDY
				sect_idxsI = where((sect_numRT>=sect1) & (sect_numRT<sect2))[0]
				#DECIDE WHAT TO DO ABOUT NUMBER OF RIDGES NEEDED
				if (size(sect_idxsI)>=2):
					mean_heightR = mean(max_heightRT[sect_idxsI])
					mean_heightsR.append(mean_heightR)
				else:
					mean_heightsR.append(np.nan)
		
					
	return array(mean_xptsB), array(mean_yptsB), array(mean_heightsR)

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

def running_mean_unc(xpts, ypts, var, window):
	mean_var_unc=[]
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
		unc=sqrt(np.sum(var[i:j]**2))/size(var[i:j])
		mean_var_unc.append(unc)
		mean_xpts.append(np.mean(xpts[i:j]))
		mean_ypts.append(np.mean(ypts[i:j]))
		i=j+0

	return array(mean_xpts), array(mean_ypts), array(mean_var_unc)

def func(x, b):
    return b*(x**0.5)

def correlate_sails_IB(addsnow):
	mean_heightALL=[]
	thicknessALL=[]
	thickness_uncIBALL=[]
	xpts_matchALL=[]
	ypts_matchALL=[]

	trend_years=[]
	intercept_years=[]
	r_years=[]
	prob_years=[]
	stderr_years=[]
	p_coef_years=[]

	thickness_years=[]
	height_years=[]
	for year in xrange(start_year, end_year+1):
		print year
		#mean_xpts, mean_ypts, mean_height = running_mean(xptsT, yptsT, max_heightT, window)
		mean_xpts, mean_ypts, mean_height = calc_mean_height_dist(year)

		xptsIB, yptsIB, latsIB,lonsIB,thicknessIB, thickness_uncIB,snowIB= read_icebridge_uncALL(mplot, rawdatapath, year, mask_hi=1, mask_nonarctic=1)

		mean_xptsIB, mean_yptsIB, mean_thicknessIB = running_mean(xptsIB, yptsIB, thicknessIB, window)
		mean_xptsIB, mean_yptsIB, mean_thickness_uncIB = running_mean(xptsIB, yptsIB, thickness_uncIB, window)
	
		mean_thicknessIBR = griddata((mean_xptsIB, mean_yptsIB),mean_thicknessIB, (mean_xpts, mean_ypts), method='nearest', rescale=True)
		mean_thicknessIBR = griddata((mean_xptsIB, mean_yptsIB),mean_thicknessIB, (mean_xpts, mean_ypts), method='nearest', rescale=True)
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

		if (addsnow==1):
			mean_xptsIBSNOW, mean_yptsIBSNOW, mean_snowIB = ro.running_mean(xptsIB, yptsIB, snowIB, window)
			mean_snowIBR = griddata((mean_xptsIB, mean_yptsIB),mean_snowIB, (mean_xpts, mean_ypts), method='nearest', rescale=True)
			mean_snowIBR = kdtree_clean1D(mean_xpts, mean_ypts, mean_xptsIB, mean_yptsIB, mean_snowIBR, kdtree_radius)
			#mean_snowIBR[where(np.isnan(mean_thicknessIBR))] = np.nan
			mean_snowIBR_match=mean_snowIBR[match_data]
			mean_height_match = mean_height_match+mean_snowIBR_match
			mean_snowIBR_match.dump(savepath+'/IBsnow_10km_'+str(year)+region_str+'.txt') 

		p= optimization.curve_fit(func, mean_thicknessIBR_match, mean_height_match, [1.])
		trendT, interceptT, rT, probT, stderrT = stats.linregress(p[0][0]*sqrt(mean_thicknessIBR_match), mean_height_match) 
		p_coef_years.append(p[0][0])
		trend_years.append(trendT)
		intercept_years.append(interceptT)
		r_years.append(rT)
		prob_years.append(probT)
		stderr_years.append(stderrT)

		thickness_years.append(mean_thicknessIBR_match)
		height_years.append(mean_height_match)

		xpts_matchALL.extend(xpts_match)
		ypts_matchALL.extend(ypts_match)
		mean_heightALL.extend(mean_height_match)
		thicknessALL.extend(mean_thicknessIBR_match)
		thickness_uncIBALL.extend(mean_thickness_uncIB)

	thicknesses = [np.mean(thickness_years[x]) for x in xrange(num_years)]
	heights = [np.mean(height_years[x]) for x in xrange(num_years)]
	pALL= optimization.curve_fit(func, thicknessALL, mean_heightALL, [1.])
	trend2, intercept2, r2, prob2, stderr2 = stats.linregress(pALL[0][0]*sqrt(thicknessALL), mean_heightALL) 

	return mean_heightALL, thicknessALL, thickness_years, height_years, r_years, r2, pALL[0][0], p_coef_years[0:6], thickness_uncIBALL


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

my_cmap=ro.perceptual_colormap("Linear_L", rawdatapath+'/OTHER/CMAPS/', reverse=1)


num_kms=10
num_gd_sections=0.3*num_kms
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

start_year=2009
end_year=2014
num_years = end_year - start_year +1
thickness_line= np.arange(0.,8., 0.1)


mean_heightALL_NOSNOW, thicknessALL_NOSNOW, thickness_yearsNOSNOW, height_yearsNOSNOW, ryears_NOSNOW, r2_NOSNOW, pALL_NOSNOW, p_coef_yearsNOSNOW, thickness_uncIBALL = correlate_sails_IB(addsnow=0)
#mean_heightALL_SNOW, thicknessALL_SNOW, thickness_yearsSNOW, height_yearsSNOW, ryears_SNOW, r2_SNOW, pALL_SNOW, p_coef_yearsSNOW = correlate_sails_IB(addsnow=1)
ryears_NOSNOW.append(r2_NOSNOW)
p_coef_yearsNOSNOW.append(pALL_NOSNOW)

thickness_uncIBALL_sq=[thickness_uncIBALL[x]**2 for x in xrange(size(thickness_uncIBALL))]

thickness_unc_tot=sqrt(np.sum(thickness_uncIBALL_sq))/size(thickness_uncIBALL)


std_err=sqrt(sum(((mean_heightALL_NOSNOW/pALL_NOSNOW)**2-thicknessALL_NOSNOW)**2)/size(thicknessALL_NOSNOW))

years = [x for x in xrange(start_year, end_year+1)]
years.append(0000)
data_NOSNOW= zip(years, ryears_NOSNOW, p_coef_yearsNOSNOW)

savetxt(savepath+'/Correlation_data_nosnow'+region_str+'.txt', data_NOSNOW, fmt= '%04d &%.2f &%.2f \\')
			
colors=['m', 'r', 'g', 'b', 'c', 'y']
fig = figure(figsize=(3.5,2.2))
ax1=subplot(1, 1, 1)
plot(thickness_line, pALL_NOSNOW*sqrt(thickness_line), 'k')

for x in xrange(end_year+1-start_year):
	scatter(thickness_yearsNOSNOW[x], height_yearsNOSNOW[x],color=colors[x],s=5, marker='x', alpha=0.5)
	plot(thickness_line, p_coef_yearsNOSNOW[x]*sqrt(thickness_line), colors[x])
	r_str = '%.2f' % ryears_NOSNOW[x]
	p_str = '%.2f' % p_coef_yearsNOSNOW[x]
	ax1.annotate(str(x+2009)+'\n(r:'+r_str+' b:'+p_str+')', xy=(1.015, 0.9-(0.15*x)), xycoords='axes fraction', horizontalalignment='left',color=colors[x])

std_err_str = '%.2f' % std_err
ax1.annotate(r'$\sigma_r$:'+std_err_str+' m', xy=(0.7, 0.07), xycoords='axes fraction', horizontalalignment='left',color='k')

r_strALL_NOSNOW = '%.2f' % r2_NOSNOW
p_strALL_NOSNOW = '%.2f' % pALL_NOSNOW
ax1.annotate('All'+'\n(r:'+r_strALL_NOSNOW+' b:'+p_strALL_NOSNOW+')', xy=(1.015, 0.00), xycoords='axes fraction', horizontalalignment='left',color='k')

ax1.set_ylabel('Feature height (m)', labelpad=2)
ax1.set_xlabel('Ice thickness (m)')
ax1.set_ylim(0., 2.5)
ax1.set_xlim(0, 8)
#ax1.annotate('(a) no snow', xy=(0.05, 0.93), xycoords='axes fraction', horizontalalignment='left',color='k')
subplots_adjust(left = 0.115, right = 0.75, bottom = 0.165, top = 0.96, wspace=0.22)
savefig(figpath+'figure12.png', dpi=300)















