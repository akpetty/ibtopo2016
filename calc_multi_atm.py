############################################################## 
# Date: 20/01/16
# Name: calc_multi_atm.py
# Author: Alek Petty
# Description: Main script to calculate sea ice topography from IB ATM data
# Input requirements: ATM data, PosAV data (for geolocation)
# Output: topography datasets

import matplotlib
matplotlib.use("AGG")
import IB_functions as ro
import mpl_toolkits.basemap.pyproj as pyproj
from osgeo import osr, gdal
from pyproj import Proj
from glob import glob
from pylab import *
from scipy import ndimage
from matplotlib import rc
#from scipy.interpolate import griddata 
from matplotlib.mlab import griddata
import time
import scipy.interpolate
import h5py
from scipy.spatial import cKDTree as KDTree
import os

def calc_bulk_stats():
		
	ice_area=-999
	ridge_area_all=-999
	mean_ridge_height_all=-999
	mean_ridge_heightL=-999
	ridge_areaL=-999
	num_ridges_out=-999
	levpercent_out=-999
	num_pts_section=-999
	# IF SECTION GOOD THEN GET ICE SWATH AREA
	if (points_good==1):
		ice_area = ma.count(elevation2d)*(xy_res**2)
		levpercent_out=levpercent
	# IF SECTION GOOD AND HAVE SOME RIDGING THEN ASSIGN TOTAL RIDGE AREA AND ELEVATION
	if ((points_good==1)&(found_ridges==1)):
		ridge_area_all = ma.count(elevation2d_ridge_ma)*(xy_res**2)
		mean_ridge_height_all = np.mean(elevation2d_ridge_ma) - level_elev

	# IF SECTION GOOD AND WE HAVE NO RIDGING (AREA OF RIDGING = 0) THEN ASSIGN ZERO RIDGE AREA HEIGHT
	if ((points_good==1)&(found_ridges==0)):
		ridge_area_all = 0.
		mean_ridge_height_all = 0.
	
	#IF GOOD SECTION BUT NO BIG RIDGES THEN SET THESE VALUES TO ZERO
	if ((points_good==1)&(found_big_ridge==0)):
		mean_ridge_heightL=0.
		ridge_areaL=0.
		num_ridges_out=0
	# IF WE FOUND SOME BIG RIDGES THENA SSIGN BIG RIDGE AREA HEIGHT AND NUMBER
	if ((points_good==1)&(found_big_ridge==1)):	
		mean_ridge_heightL = np.mean(ridge_height_mesh)
		ridge_areaL = ma.count(ridge_height_mesh)*(xy_res**2)
		num_ridges_out = num_ridges


	return [mean_x, mean_y, ice_area, num_ridges_out, ridge_area_all, ridge_areaL, mean_ridge_height_all, mean_ridge_heightL, mean_alt, mean_pitch, mean_roll, mean_vel, num_pts_section,levpercent_out, section_num, found_ridges, points_good, plane_good]
	
#-------------- ATM AND DMS PATHS------------------

datapath='./Data_output/'
rawdatapath = '../../../DATA/ICEBRIDGE/'
ATM_path = rawdatapath+'/ATM/ARCTIC/'
posAV_path =rawdatapath+'/POSAV/SEA_ICE/GR/'
#posAV_path ='/Volumes/TBOLT_HD_PETTY/POSAV/'
m=pyproj.Proj("+init=EPSG:3413")
#FREE PARAMETERS

min_ridge_height = 0.2
along_track_res=1000
pwidth=20
pint=5
xy_res=2
start_year=2009
end_year=2009
min_ridge_size=100

sh=0
if (sh==1):
	print 'Ridge threshold:', sys.argv[1]
	print 'Along track res:',sys.argv[2]
	print 'xy res:',sys.argv[3]
	print 'Start year:',sys.argv[4]
	print 'End year:',sys.argv[5]
	min_ridge_height = float(sys.argv[1])
	along_track_res = int(sys.argv[2])
	xy_res = int(sys.argv[3])
	start_year=int(sys.argv[4])
	end_year=int(sys.argv[5])

pts_threshold=15000
num_points_req = min_ridge_size/(xy_res**2)
section_num=0

print 'Num points req', num_points_req

ftype = str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm'
outpath = datapath+ftype+'/'


for year in xrange(start_year, end_year+1):
	ATM_year = ATM_path+str(year)+'/'
	atm_files_year = glob(ATM_year+'/*/')
	#for days in xrange():
	for days in xrange(size(atm_files_year)):
		atm_path_date = atm_files_year[days]
		print 'ATM day:', atm_path_date
		atm_files_in_day = ro.get_atm_files(atm_path_date, year)
		#load POS file
		posAV = loadtxt(posAV_path+str(year)+'_GR_NASA/sbet_'+str(atm_path_date[-9:-1])+'.out.txt', skiprows=1)
		#GET POSITION OF PLANE AND 1km MARKERS FROM POSAV
		xp, yp, dist, km_idxs, km_utc_times = ro.get_pos_sections(posAV, m, along_track_res)

		for atm_file in xrange(size(atm_files_in_day)):
			atm_statsALL=np.array([]).reshape(0,3)
			ridge_statsALL=np.array([]).reshape(0,9)
			covarALL=np.array([]).reshape(0,5)
			bulk_statsALL=np.array([]).reshape(0,18)
			print 'ATM file:', atm_files_in_day[atm_file], str(atm_file)+'/'+str(size(atm_files_in_day))
			lonT, latT, elevationT, utc_timeT= ro.get_atmqih5(atm_files_in_day[atm_file], year, 1)
			#IF SIZE OF DATA IS LESS THAN SOME THRESHOLD THEN DONT BOTHER ANALYZING
			if (size(utc_timeT)<100):
				break
			xT, yT = m(lonT, latT)
			#GET POSAV INDICES COINCIDING WITH START AND END OF ATM FILE. ADD PLUS/MINUS 1 FOR SOME LEEWAY.
			start_i = np.abs(km_utc_times - utc_timeT[0]).argmin()

			end_i = np.abs(km_utc_times - utc_timeT[-1]).argmin()
			print 'START/END:', start_i, end_i

			for i in xrange(start_i -1, end_i + 1):
				section_num+=1
				found_ridges=0
				found_big_ridge=0
				plane_good=0
				points_good=0
				ridge_statsT = np.array([]).reshape(0,9)
				cov_matrix = np.array([]).reshape(0,5)
				#label_numsL=np.array(0)
				mean_x, mean_y, mean_alt, mean_pitch, mean_roll, mean_vel = ro.posav_section_info(m, posAV[km_idxs[i]:km_idxs[i+1]]	)
				print '    '
				print str(i)+'/'+str(end_i + 1)
				print 'Mean altitude:', mean_alt
				print 'Mean pitch:', mean_pitch
				print 'Mean roll:', mean_roll
				print 'Mean vel:', mean_vel
				
				if (abs(mean_alt-500)<200) & (abs(mean_pitch)<5) & (abs(mean_roll)<5):
					plane_good=1
					
					poly_path, vertices, sides = ro.get_pos_poly(xp, yp, km_idxs[i], km_idxs[i+1])
					
					xatm_km, yatm_km, elevation_km = ro.get_atm_poly(xT, yT, elevationT, km_utc_times, utc_timeT, poly_path, i)
					num_pts_section = size(xatm_km)
					print 'Num pts in section:', size(xatm_km)
					#if there are more than 15000 pts in the 1km grid (average of around 20000) then proceed
					
					if  (num_pts_section>pts_threshold):	
						points_good=1
						#ro.plot_atm_poly(m, xatm_km, yatm_km, elevation_km, poly_path, i, out_path, year)
						
						#GET ATM GRID
						xx2d, yy2d = ro.grid_atm(xatm_km, yatm_km, xy_res)
						print 'Grid:', size(xx2d[0]), size(xx2d[1])
						# CALCULATE THE LEVEL ICE SURFACE USING THE CUMULATIVE DISTRIBUTION
						#THRESH IS THE LEVEL ICE PLUS RIDGED ICE ELEVATION
						level_elev, thresh, levpercent = ro.calc_level_ice(elevation_km, pint, pwidth, min_ridge_height)

						#level_elev, thresh, levpercent = ro.calc_level_ice(elevation_km, pwidth, min_ridge_height)

						elevation2d, elevation2d_ridge_ma, ridge_area = ro.grid_elevation(xatm_km, yatm_km,elevation_km, xx2d, yy2d, thresh, kdtree=1)
						elevation2d_ridge_maL =elevation2d_ridge_ma-level_elev
						#IF THERE IS EVEN A LITTLE BIT OF RIDGING (might not actually be enough for a big areal ridge from the labelling) then proceed to clean up data.
						if (ridge_area>0):
							found_ridges=1
							#CLEAN UP DATA WITH KDTREE AROUND RIDGE POINTS
							#REMOVE FOR PRELIMINARY STUDIES AS COMPUTATIONALLY EXPENSIVE!
							#elevation2d_ridge_ma = kdtree_clean()
							#GET RIDGE LABELS - MAKE SURE RIDGES ARE ABOVE CERTAIN SIZE, DECIDED BY NUM_PTS_REQ
							label_im  = ro.get_labels(elevation2d_ridge_maL, xy_res, min_ridge_size, min_ridge_height)
							# DECIDE IF WE WANT TO CALCULATE RIDGE ORIENTATION OR NOT.
							if (np.amax(label_im)>=1):
								found_big_ridge=1
								print 'Found Ridge!'
								print 'Number of labels:', np.amax(label_im)
								num_ridges = np.amax(label_im)
								#GET RIDGE STATS IF WE DID FIND A RIDGE
								ridge_statsT, ridge_height_mesh, cov_matrix, indexT = ro.calc_ridge_stats(elevation2d_ridge_ma, num_ridges, label_im, xx2d, yy2d, level_elev, section_num, calc_orientation=1)
								
						#CALCULATE BULK STATISTICS AS WE HAVE VALID NUMBER OF POINTS WITHIN THE SECTION
						
					else:
						print 'No data - WHY?! --------------'
						print 'Num pts in section:', size(xatm_km)
				#ASSIGN BULK STATISTICS AS WE HAVE NOT CARRIED OUT RIDGE CALCULATION AS PLANE IS DOING FUNNY THINGS
				bulk_statsT = calc_bulk_stats()
				ridge_statsALL = vstack([ridge_statsALL, ridge_statsT])
				covarALL = vstack([covarALL, cov_matrix])
				bulk_statsALL = vstack([bulk_statsALL, bulk_statsT])
			
			if not os.path.exists(outpath+str(year)):
				os.makedirs(outpath+str(year))

			ridge_statsALL.dump(outpath+str(year)+'/ridge_stats_'+str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm_poly'+str(atm_path_date[-9:-1])+'_f'+str(atm_file).zfill(3)+'.txt')
			covarALL.dump(outpath+str(year)+'/cov_matrix_'+str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm_poly'+str(atm_path_date[-9:-1])+'_f'+str(atm_file).zfill(3)+'.txt')
			bulk_statsALL.dump(outpath+str(year)+'/bulk_ridge_stats_'+str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm_poly'+str(atm_path_date[-9:-1])+'_f'+str(atm_file).zfill(3)+'.txt')
 			
 			#CAN OUTPUT AS TEXT FILES INSTEAD - BIGGER BUT CAN OPEN RAW
			#savetxt(outpath+str(year)+'/ridge_stats_'+str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm_poly'+str(atm_path_date[-9:-1])+'_f'+str(atm_file)+'.txt', ridge_statsALL)
			#savetxt(outpath+str(year)+'/cov_matrix_'+str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm_poly'+str(atm_path_date[-9:-1])+'_f'+str(atm_file)+'.txt', covarALL)
			#savetxt(outpath+str(year)+'/bulk_ridge_stats_'+str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm_poly'+str(atm_path_date[-9:-1])+'_f'+str(atm_file)+'.txt', bulk_statsALL)


