#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import os, os.path
import sys

import numpy as np
import amv as amv

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from matplotlib.colors import LogNorm
from pylab import *

undef = -9999.

# Pick up parameters
file_amv=sys.argv[1]
source=sys.argv[2]
qc_type=sys.argv[4]
#qc_min=int(sys.argv[5])
qc_min=80
#print qc_min
#qc_max=int(sys.argv[6])
qc_max=100
#print qc_max

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
delimiter=','
#delimiter=None
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
file_fcst = '/Users/sreiner/Documents/stratus/datafiles/ERA5_UV_prs_20191020_hourly.grib'
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
amv_qc  = amv_data[5,:]

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
fo = open("amv_vd_stats.txt","a")
statsvd = "{0} {1} Total Number, Best Fit Number, VD OMB Mean,RMSE and VD OMB after fit Mean, RMSE: {2} {3} {4:.2f} {5:.2f} {6:.2f} {7:.2f}\n" \
.format(file_amv,qc_str,amv_num,bfit_num,mean_vd_bg,rms_vd_bg,mean_vd_new,rms_vd_new)
fo.write( statsvd);
fo.close()


frac = np.empty(4)
flag_title=['Found','Not Constrained','No sufficient minimum','No forecast pressure match']
for f in range (4):
    frac[f] = float((bfit_flag == f).sum()) / float(amv_num)

# Page 1


fig = plt.figure(figsize=(12,12))

xmin=-40.
xmax= 40.
statx=0.05
staty=0.9
bins=np.arange(81) - 40.
label = 'U Component AMV-Forecast [m/s]'

ax=plt.subplot(2,2,1)
title = 'Before Fit'
values = u_omb_bg
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,2)
title = 'After Fit'
values = u_omb_new
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,3)
title = 'Subset of BFIT points Before Fit'
values = u_omb_bf_orig
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,4)
title = 'Subset of BFIT points After Fit'
values = u_omb_bf
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

plt.figtext(0.5,0.96,fig_title,ha='center',color='black',weight='bold',size='large')

#plt.show()
plt.gcf().set_size_inches(13, 13)
plt.savefig(omb_u_fig)
plt.clf()

print("wrote figure file: "  + omb_u_fig)

# Page 2


fig = plt.figure(figsize=(12,12))

xmin=-20.
xmax= 20.
statx=0.05
staty=0.9
bins=np.arange(81)*0.5 - 20.
label = 'V Component AMV-Forecast [m/s]'

ax=plt.subplot(2,2,1)
title = 'Before Fit'
values = v_omb_bg
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,2)
title = 'After Fit'
values = v_omb_new
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,3)
title = 'Subset of BFIT points Before Fit'
values = v_omb_bf_orig
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,4)
title = 'Subset of BFIT points After Fit'
values = v_omb_bf
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

plt.figtext(0.5,0.96,fig_title,ha='center',color='black',weight='bold',size='large')

#plt.show()
plt.gcf().set_size_inches(13, 13)
plt.savefig(omb_v_fig)
plt.clf()

print("wrote figure file: "  + omb_v_fig)

# Page 3


fig = plt.figure(figsize=(12,12))

xmin=-40.
xmax= 40.
statx=0.05
staty=0.9
bins=np.arange(81) - 40.
label = 'Speed AMV-Forecast [m/s]'

ax=plt.subplot(2,2,1)
title = 'Before Fit'
values = spd_omb_bg
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,2)
title = 'After Fit'
values = spd_omb_new
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,3)
title = 'Subset of BFIT points Before Fit'
values = spd_omb_bf_orig
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,4)
title = 'Subset of BFIT points After Fit'
values = spd_omb_bf
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='left',va='center',transform=ax.transAxes)

plt.figtext(0.5,0.96,fig_title,ha='center',color='black',weight='bold',size='large')

#plt.show()
plt.gcf().set_size_inches(13, 13)
plt.savefig(omb_s_fig)
plt.clf()

print("wrote figure file: "  + omb_s_fig)

# Page 4


fig = plt.figure(figsize=(12,12))

xmin=0.
xmax= 40.
statx=0.9
staty=0.9
bins=np.arange(81)*0.5
label = 'Vector Difference AMV-Forecast [m/s]'

ax=plt.subplot(2,2,1)
title = 'Before Fit'
values = vd_omb_bg
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='right',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,2)
title = 'After Fit'
values = vd_omb_new
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='right',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,3)
title = 'Subset of BFIT points Before Fit'
values = vd_omb_bf_orig
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='right',va='center',transform=ax.transAxes)

ax=plt.subplot(2,2,4)
title = 'Subset of BFIT points After Fit'
values = vd_omb_bf
n,bins,patches = plt.hist(values,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
plt.xlabel(label)
plt.ylabel('Number')
plt.title(title)
stats="Mean {0:.2f} STD {1:.2f}".format(np.mean(values),np.std(values))
text(statx,staty,stats,ha='right',va='center',transform=ax.transAxes)

plt.figtext(0.5,0.96,fig_title,ha='center',color='black',weight='bold',size='large')

#plt.show()
plt.gcf().set_size_inches(13, 13)
plt.savefig(omb_vd_fig)
plt.clf()

print("wrote figure file: "  + omb_vd_fig)


# Page 5
x=amv_data[3,:]
y=amv_data[2,:]


fig = plt.figure(figsize=(12,12))

plt.subplot(3,2,1)
num_bins=50
n,bins,patches = plt.hist(x,num_bins,facecolor='green',alpha=0.5)
n,bins,patches = plt.hist(dp_x,num_bins,facecolor='orange',alpha=0.5)
plt.xlim(210.,360.)
#plt.xlim(-60.,60.)
#plt.ylim(0.,600.)
plt.xlabel('Longitudes')
plt.ylabel('Number')
plt.title(r'Green all AMV, Yellow found best fit')

plt.subplot(3,2,3)
num_bins=50
n,bins,patches = plt.hist(y,num_bins,facecolor='green',alpha=0.5)
n,bins,patches = plt.hist(dp_y,num_bins,facecolor='orange',alpha=0.5)
plt.xlim(-60.,60.)
#plt.ylim(0.,600.)
plt.xlabel('Latitudes')
plt.ylabel('Number')
#plt.title(r'Green all AMV, Yellow found best fit')

plt.subplot(4,2,7)
#not enough flag=3 to plot for these test cases
labels=flag_title[0:3]
sizes=frac[0:3]
colors=['lightblue','orange','lightgreen','red']
plt.pie(sizes,labels=labels,colors=colors,autopct='%.0f%%')
plt.axis('equal')

plt.subplot(2,2,2)
num_bins=50
n,bins,patches = plt.hist(dp,num_bins,facecolor='blue',alpha=0.5)
plt.xlim(-300.,300.)
#plt.ylim(0.,100.)
plt.xlabel('dp')
plt.ylabel('Number')
plt.title(r'Histogram of AMV Best Fit - Original Pressure')

plt.subplot(2,2,4)
plt.scatter(dp_x_none,dp_y_none,s=1,color='0.8')
plt.scatter(dp_x_up,dp_y_up,s=1,color='r')
plt.scatter(dp_x_down,dp_y_down,s=1,color='b')
#m = Basemap(projection='cyl',llcrnrlat=-80.,urcrnrlat=80.,llcrnrlon=-130.,urcrnrlon=-20.,resolution='c')
m = Basemap(projection='cyl',llcrnrlat=-80.,urcrnrlat=80.,llcrnrlon=230.,urcrnrlon=340.,resolution='c')
m.drawcoastlines()
plt.title(r'Grey no fit, Red BFIT up, Blue BFIT down ')

plt.figtext(0.5,0.96,fig_title,ha='center',color='black',weight='bold',size='large')

#plt.show()
plt.gcf().set_size_inches(13, 13)
plt.savefig(main_fig)
plt.clf()

print("wrote figure file: "  + main_fig)


# Page 6

fig = plt.figure(figsize=(12,12))

ax=plt.subplot(2,1,1)
lat = amv_lat
lon = amv_lon
prs = amv_prs
low_lat = lat[prs >=700.]
low_lon = lon[prs >=700.]
mid_lat = lat[(prs < 700.) & (prs >400.)]
mid_lon = lon[(prs < 700.) & (prs >400.)]
hig_lat = lat[prs <=400.]
hig_lon = lon[prs <=400.]
plt.scatter(low_lon,low_lat,s=1,color='darkblue')
plt.scatter(mid_lon,mid_lat,s=1,color='g')
plt.scatter(hig_lon,hig_lat,s=1,color='orange')
#m = Basemap(projection='cyl',llcrnrlat=-80.,urcrnrlat=80.,llcrnrlon=-130.,urcrnrlon=-20.,resolution='c')
#m = Basemap(projection='cyl',llcrnrlat=-80.,urcrnrlat=80.,llcrnrlon=60.,urcrnrlon=220.,resolution='c')
m = Basemap(projection='cyl',llcrnrlat=-80.,urcrnrlat=80.,llcrnrlon=230.,urcrnrlon=340.,resolution='c')
m.drawcoastlines()
plt.xlabel(r'Blue 700 hPa, Green 700-400 hPa, Yellow Above 400 hPa')
plt.title(r'AMV Location')


ax=plt.subplot(2,2,3)
plt.gca().invert_yaxis()
xmin=-60.
xmax=60.
ymin=1000.
ymax=100.
plt.xlim(xmin,xmax)
#plt.ylim(ymin,ymax)
#plt.scatter(lat,prs,s=1,color='0.8')
colors=['purple','blue','cyan','green','yellow','orange','red','magenta','black']
spd_min=0.
spd_max=10.
count=0
while (spd_max < 100):
    lat_bin=lat[(amv_spd>=spd_min) & (amv_spd<spd_max)]
    prs_bin=prs[(amv_spd>=spd_min) & (amv_spd<spd_max)]
    plt.scatter(lat_bin,prs_bin,s=1,color=colors[count])
    count+=1
    spd_min=spd_max
    spd_max=spd_max+10.

plt.xlabel('Latitude [deg]')
plt.ylabel('Pressure [hPa]')
plt.title(r'Color indicates speed - purple slow, red fast')


ax=plt.subplot(2,2,4)
xmin=0.
xmax=100.
bins=np.arange(101)
n,bins,patches = plt.hist(amv_qc,bins=bins,facecolor='green',alpha=0.5)
plt.xlim(xmin,xmax)
#plt.ylim(0.,600.)
plt.xlabel('QI')
plt.ylabel('Number')

plt.figtext(0.5,0.96,fig_title,ha='center',color='black',weight='bold',size='large')

plt.gcf().set_size_inches(13, 13)
plt.savefig(amv_fig)
plt.clf()
print("wrote figure file: "+amv_fig)

# Page 7

plt.subplot(2,2,1)
bfp = bfit_prs[(bfit_prs !=undef)]
amvp = amv_prs[(bfit_prs !=undef)]

dp_low = bfp[(amvp>=700.)] - amvp[(amvp>=700.)] 
dp_mid = bfp[(amvp<700.) & (amvp>400.)] - amvp[(amvp<700.) & (amvp>400.)] 
dp_hig = bfp[(amvp<=400.)] - amvp[(amvp<=400.)] 

#num_bins=50
#n,bins,patches = plt.hist(dp,num_bins,facecolor='blue',alpha=0.5)
bins=arange(61)*10. - 300.
hist1,binedges = np.histogram(dp_low,bins=bins)
x1=0.5*(binedges[1:]+binedges[:-1])
hist2,binedges = np.histogram(dp_mid,bins=bins)
x2=0.5*(binedges[1:]+binedges[:-1])
hist3,binedges = np.histogram(dp_hig,bins=bins)
x3=0.5*(binedges[1:]+binedges[:-1])
plt.plot(x1,hist1,'darkblue',x2,hist2,'g',x3,hist3,'orange')
plt.xlim(-300.,300.)
#plt.ylim(0.,100.)
plt.xlabel('BFIT pressure - AMV pressure [hPa]')
plt.ylabel('Number')
plt.title(r'Blue - low, Green - mid, Yellow - high')

plt.subplot(2,2,2)
lat = amv_lat[bfit_prs != undef]
lon = amv_lon[bfit_prs != undef]
prs = amv_prs[bfit_prs != undef]
low_lat = lat[prs >=700.]
low_lon = lon[prs >=700.]
mid_lat = lat[(prs < 700.) & (prs >400.)]
mid_lon = lon[(prs < 700.) & (prs >400.)]
hig_lat = lat[prs <=400.]
hig_lon = lon[prs <=400.]

plt.scatter(low_lon,low_lat,s=1,color='darkblue')
plt.scatter(mid_lon,mid_lat,s=1,color='g')
plt.scatter(hig_lon,hig_lat,s=1,color='orange')
#m = Basemap(projection='cyl',llcrnrlat=-80.,urcrnrlat=80.,llcrnrlon=-130.,urcrnrlon=-20.,resolution='c')
#m = Basemap(projection='cyl',llcrnrlat=-80.,urcrnrlat=80.,llcrnrlon=60.,urcrnrlon=220.,resolution='c')
m = Basemap(projection='cyl',llcrnrlat=-80.,urcrnrlat=80.,llcrnrlon=230.,urcrnrlon=340.,resolution='c')
m.drawcoastlines()
plt.xlabel(r'Blue Below 700 hPa, Green 700-400 hPa, Yellow Above 400 hPa')
plt.title(r'BFIT AMV Location')

ax=plt.subplot(2,2,3)
plt.gca().invert_yaxis()
xmin=-60.
xmax=60.
ymin=1000.
ymax=100.
lat_none = amv_lat[bfit_prs == undef]
prs_none = amv_prs[bfit_prs == undef]
lat_down = amv_lat[bfit_prs > amv_prs]
prs_down = amv_prs[bfit_prs > amv_prs]
lat_up = amv_lat[(bfit_prs != undef) & (bfit_prs < amv_prs)]
prs_up = amv_prs[(bfit_prs != undef) & (bfit_prs < amv_prs)]
plt.xlim(xmin,xmax)
#plt.ylim(ymin,ymax)
plt.scatter(lat_none,prs_none,s=1,color='0.8')
plt.scatter(lat_up,prs_up,s=1,color='r')
plt.scatter(lat_down,prs_down,s=1,color='b')
plt.xlabel('Latitude [deg]')
plt.ylabel('Pressure [hPa]')
plt.title(r'Red BFIT up, Blue BFIT down, Grey no BFIT')

plt.subplot(2,2,4)
plt.gca().invert_yaxis()
xmin = 0.
xmax = 70.
ymin = 1000.
ymax = 100.
spd_none = amv_spd[bfit_prs == undef]
prs_none = amv_prs[bfit_prs == undef]
spd_down = amv_spd[bfit_prs > amv_prs]
prs_down = amv_prs[bfit_prs > amv_prs]
spd_up = amv_spd[(bfit_prs != undef) & (bfit_prs < amv_prs)]
prs_up = amv_prs[(bfit_prs != undef) & (bfit_prs < amv_prs)]

plt.xlim(xmin,xmax)
#plt.ylim(ymin,ymax)
plt.scatter(spd_none,prs_none,s=1,color='0.8')
plt.scatter(spd_down,prs_down,s=1,color='b')
plt.scatter(spd_up,prs_up,s=1,color='r')
plt.title(r'Red BFIT up, Blue BFIT down, Grey no BFIT')
plt.xlabel(r'Speed [m/s]')
plt.ylabel(r'Pressure [hPa]')

plt.figtext(0.5,0.96,fig_title,ha='center',color='black',weight='bold',size='large')


plt.gcf().set_size_inches(13, 13)
plt.savefig(bfit_fig)
plt.clf()

print("wrote figure file: "+bfit_fig)

