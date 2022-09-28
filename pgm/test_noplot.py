#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import os, os.path
import sys

import numpy as np
import amv as amv

from pylab import *

undef = -9999.

# ./test_noplot.py /home/daves/intercomparison2021/BRZ/4th_AMVIC_INPE_Test_2b_final.txt  BRZexp22  'BRZ Exp22' CQI 60 69

# Pick up parameters
file_amv=sys.argv[1]
source=sys.argv[2]
qc_type=sys.argv[4]
qc_min=int(sys.argv[5])
print("QI min {0}".format(qc_min))
qc_max=int(sys.argv[6])
print("QI max {0}".format(qc_max))

if qc_type == 'QINF' :
   qc_pos=16  # 0 is first position
   qc_str='QINF:'+str(qc_min) + '-' + str(qc_max)
   qc_txt='.qinf_'+str(qc_min) + '-' + str(qc_max)
   
if qc_type == 'QIWF' :
   qc_pos=17
   qc_str='QIWF:'+str(qc_min) + '-' + str(qc_max)
   qc_txt='.qiwf_'+str(qc_min) + '-' + str(qc_max)
   
if qc_type == 'CQI' :
   qc_pos=18
   qc_str='CQI:'+str(qc_min) + '-' + str(qc_max)
   qc_txt='.cqi_'+str(qc_min) + '-' + str(qc_max)

fig_title=sys.argv[3] + qc_str


#files
omb_u_fig='omb_u.'+source+qc_txt+'.png'
omb_v_fig='omb_v.'+source+qc_txt+'.png'
omb_s_fig='omb_s.'+source+qc_txt+'.png'
omb_vd_fig='omb_vd.'+source+qc_txt+'.png'
main_fig='bfit.'+source+qc_txt+'.png'
amv_fig='amv_location.'+source+qc_txt+'.png'
bfit_fig='bfit_location.'+source+qc_txt+'.png'


print("Reading amv data file: "  + file_amv)
#delimiter=';'
delimiter=','
#amv data is mcidas west positive -180 to 180, data is flipped and rotated to 0-360 longitude
#amv_data = amv.read_txt(file_amv,qc_pos,qc_min,qc_max,delimiter=delimiter,lonflip=True,lon0=180.,shift=-1)
#amv_data = amv.read_txt(file_amv,qc_pos,qc_min,qc_max,delimiter=delimiter,lonflip=True,lon0=180.)
amv_data = amv.read_txt(file_amv,qc_pos,qc_min,qc_max,delimiter=delimiter,lonflip=None,lon0=180.)
#lay_loc = [(amv_data[4,:] >=700.)]
#print lay_loc

print("Reading forecast file")
#file_fcst = '/home/daves/BestFit/DecodedForecast_20120917070613Z_20120917120000Z_12_O_MPFS03'
#fcst_data = amv.read_DecodedForecast_MSG(file_fcst)
#file_fcst = '/home/snebuda/icomp/data/ECM_EI_AN_20160721_PL.grb'
#file_fcst = '/data/rdworak/Intercomp/model/ECM_EI_AN_20160721_PL.grb'
# file_fcst = '/home/daves/intercomparison2021/ERA5/ERA5_UV_prs_20191020_hourly.grib'
#file_fcst = '/Users/sreiner/Documents/stratus/datafiles/ERA5_UV_prs_20191020_hourly.grib'
file_fcst = '/data/sreiner/datafiles/ERA5_UV_prs_20191020_hourly.grib'
#file_fcst = '/data/rdworak/Intercomp/model/ECM_EI_AN_20160721_PL.grb'
#datetime=2016072112
datetime=2019102012
#forecast data is 0-360 Longitude which works for H8 AMV data
#no handling of lon mismatch in amv.locate
fcst_data, dt_out = amv.read_MSG_grib(file_fcst,datetime=datetime)

print("Finding grid location")
grid_i,grid_j = amv.locate(amv_data,fcst_data)
#print("value of i".format(grid_i[1]))
#print("value of j".format(grid_j[1]))
print("Finding best fit")

amv_num = amv_data.shape[1]
amv_u=np.empty(amv_num)
amv_v=np.empty(amv_num)
bfit_u=np.empty(amv_num)
bfit_v=np.empty(amv_num)
bfit_prs=np.empty(amv_num)
bfit_flag=np.empty(amv_num)
bfit_spd=np.empty(amv_num)
bfit_dir=np.empty(amv_num)
bg_u=np.empty(amv_num)
bg_v=np.empty(amv_num)
bg_spd=np.empty(amv_num)
bg_dir=np.empty(amv_num)

print("Number of AMV {0}".format(amv_num))

amv_spd = amv_data[0,:]
amv_dir = amv_data[1,:]

n = 0
while (n < amv_num):
    amv_single = amv_data[:,n]
    fcst_profile = fcst_data[grid_i[n],grid_j[n],:,:]

    verbose = False
    bfit_u[n],bfit_v[n],bfit_prs[n],bfit_flag[n] = amv.bestfit(amv_single,fcst_profile,verbose)
    if bfit_prs[n] != undef:
        bfit_spd[n],bfit_dir[n] = amv.spddir(bfit_u[n],bfit_v[n])
    else:
        bfit_spd[n] = undef
        bfit_dir[n] = undef

    bg_u[n],bg_v[n] = amv.bg(amv_single,fcst_profile,verbose)
    bg_spd[n],bg_dir[n] = amv.spddir(bg_u[n],bg_v[n])
    amv_u[n],amv_v[n] = amv.uvcomp(amv_spd[n],amv_dir[n])

    n += 1

amv_prs = amv_data[4,:]
amv_lat = amv_data[2,:]
amv_lon = amv_data[3,:]
amv_qc = amv_data[5,:]

# Compute change in pressure for AMV best fit (sometimes there are bestfit heights, but not original match of AMV to background
# when the AMV height is higher than the lowest background pressure)
#bfit_loc = [bfit_prs !=undef]
#bfit_loc = [(bfit_prs !=undef) & (bg_spd !=undef)]
bfit_loc = np.logical_and(bfit_prs != undef, bg_spd != undef)
#bfit_loc = [(bfit_prs !=undef) & (bg_spd !=undef) & (amv_prs < 700.) & (amv_prs > 400.)]
tmp = bfit_prs[bfit_loc]
bfit_num = tmp.shape[0]
dp = bfit_prs[bfit_loc] - amv_prs[bfit_loc] 
dp_x = amv_lon[bfit_loc]
dp_y = amv_lat[bfit_loc]
dp_x_none = amv_lon[bfit_prs == undef]
dp_y_none = amv_lat[bfit_prs == undef]
dp_x_down = amv_lon[bfit_prs > amv_prs]
dp_y_down = amv_lat[bfit_prs > amv_prs]
dp_x_up = amv_lon[(bfit_prs != undef) & (bfit_prs < amv_prs)]
dp_y_up = amv_lat[(bfit_prs != undef) & (bfit_prs < amv_prs)]

#set new background fields as if amv was moved to best fit pressure
bg_u_new = bg_u.copy()
bg_u_new[bfit_loc] = bfit_u[bfit_loc]
bg_v_new = bg_v.copy()
bg_v_new[bfit_loc] = bfit_v[bfit_loc]
bg_spd_new = bg_spd.copy()
bg_spd_new[bfit_loc] = bfit_spd[bfit_loc]

#Compute OMB for U, V, Spd, Vector at AMV prs and best fit prs
u_omb_bf = amv_u[bfit_loc] - bfit_u[bfit_loc] 
v_omb_bf = amv_v[bfit_loc] - bfit_v[bfit_loc] 
spd_omb_bf = amv_spd[bfit_loc] - bfit_spd[bfit_loc] 
vd_omb_bf = np.sqrt( (amv_u[bfit_loc]-bfit_u[bfit_loc])**2 + (amv_v[bfit_loc]-bfit_v[bfit_loc])**2 )

u_omb_bf_orig = amv_u[bfit_loc] - bg_u[bfit_loc] 
v_omb_bf_orig = amv_v[bfit_loc] - bg_v[bfit_loc] 
spd_omb_bf_orig = amv_spd[bfit_loc] - bg_spd[bfit_loc] 
vd_omb_bf_orig = np.sqrt( (amv_u[bfit_loc]-bg_u[bfit_loc])**2 + (amv_v[bfit_loc]-bg_v[bfit_loc])**2 )

u_omb_bg = amv_u[bg_u !=undef] - bg_u[bg_u !=undef] 
v_omb_bg = amv_v[bg_v !=undef] - bg_v[bg_v !=undef] 
spd_omb_bg = amv_spd[bg_spd !=undef] - bg_spd[bg_spd !=undef] 
vd_omb_bg = np.sqrt( (amv_u[bg_spd !=undef]-bg_u[bg_spd !=undef])**2 + (amv_v[bg_spd !=undef]-bg_v[bg_spd !=undef])**2 )

u_omb_new = amv_u[bg_u_new !=undef] - bg_u_new[bg_u_new !=undef] 
v_omb_new = amv_v[bg_v_new !=undef] - bg_v_new[bg_v_new !=undef] 
spd_omb_new = amv_spd[bg_spd_new !=undef] - bg_spd_new[bg_spd_new !=undef] 
vd_omb_new = np.sqrt( (amv_u[bg_spd_new !=undef]-bg_u_new[bg_spd_new !=undef])**2 + (amv_v[bg_spd_new !=undef]-bg_v_new[bg_spd_new !=undef])**2 )

mean_vd_bg = np.mean(vd_omb_bg)
std_vd_bg = np.std(vd_omb_bg)
rms_vd_bg = np.sqrt( mean_vd_bg**2 + std_vd_bg**2 )

mean_vd_new = np.mean(vd_omb_new)
std_vd_new = np.std(vd_omb_new)
rms_vd_new = np.sqrt( mean_vd_new**2 + std_vd_new**2 )

#Write out stats to a file
#fo = open("amv_vd_stats.txt","a")
fo = open('stats_'+source+qc_txt+'.txt',"a")
statsvd = "{0} {1} Total Number, Best Fit Number, VD OMB Mean,RMSE and VD OMB after fit Mean, RMSE: {2} {3} {4:.2f} {5:.2f} {6:.2f} {7:.2f}\n" \
.format(file_amv,qc_str,amv_num,bfit_num,mean_vd_bg,rms_vd_bg,mean_vd_new,rms_vd_new)
fo.write( statsvd);
fo.close()


frac = np.empty(4)
flag_title=['Found','Not Constrained','No sufficient minimum','No forecast pressure match']
for f in range (4):
    frac[f] = float((bfit_flag == f).sum()) / float(amv_num)

