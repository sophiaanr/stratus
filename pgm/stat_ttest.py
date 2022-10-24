#!/usr/bin/env python3

import numpy as np
from scipy import stats
import amv as amv
import sys
undef = -9999

file_amv = sys.argv[1]
qc_min = [60, 70, 80, 90]
qc_max = [69, 79, 89, 100]
qc_pos = 18
source = sys.argv[2]
print(type(source))

# place to hold the ttest data
ttest_data = []

file_fcst = '/data/sreiner/datafiles/ERA5_UV_prs_20191020_hourly.grib'
datetime = 2019102012

fcst_data, dt_out = amv.read_MSG_grib(file_fcst, datetime=datetime)

delimiter = ','

for i in range(len(qc_min)):
    amv_data = amv.read_txt(file_amv, qc_pos, qc_min[i], qc_max[i], delimiter=delimiter, lonflip=None, lon0=180.)

    print("Finding grid location")
    grid_i, grid_j = amv.locate(amv_data, fcst_data)

    print("Finding best fit")
    amv_num = amv_data.shape[1]
    amv_u = np.empty(amv_num)
    amv_v = np.empty(amv_num)
    bfit_u = np.empty(amv_num)
    bfit_v = np.empty(amv_num)
    bfit_prs = np.empty(amv_num)
    bfit_flag = np.empty(amv_num)
    bfit_spd = np.empty(amv_num)
    bfit_dir = np.empty(amv_num)
    bg_u = np.empty(amv_num)
    bg_v = np.empty(amv_num)
    bg_spd = np.empty(amv_num)
    bg_dir = np.empty(amv_num)

    print("Number of AMV {0}".format(amv_num))

    amv_spd = amv_data[0, :]
    amv_dir = amv_data[1, :]

    n = 0
    while n < amv_num:
        amv_single = amv_data[:, n]
        fcst_profile = fcst_data[grid_i[n], grid_j[n], :, :]

        verbose = False
        bfit_u[n], bfit_v[n], bfit_prs[n], bfit_flag[n] = amv.bestfit(amv_single, fcst_profile, verbose)
        if bfit_prs[n] != undef:
            bfit_spd[n], bfit_dir[n] = amv.spddir(bfit_u[n], bfit_v[n])
        else:
            bfit_spd[n] = undef
            bfit_dir[n] = undef

        bg_u[n], bg_v[n] = amv.bg(amv_single, fcst_profile, verbose)
        bg_spd[n], bg_dir[n] = amv.spddir(bg_u[n], bg_v[n])
        amv_u[n], amv_v[n] = amv.uvcomp(amv_spd[n], amv_dir[n])

        n += 1

    amv_prs = amv_data[4, :]
    amv_lat = amv_data[2, :]
    amv_lon = amv_data[3, :]
    amv_qc = amv_data[5, :]

    # Compute change in pressure for AMV best fit (sometimes there are bestfit heights, but not original match of AMV to background
    # when the AMV height is higher than the lowest background pressure)
    # bfit_loc = [bfit_prs !=undef]
    # bfit_loc = [(bfit_prs !=undef) & (bg_spd !=undef)]
    bfit_loc = np.logical_and(bfit_prs != undef, bg_spd != undef)
    # bfit_loc = [(bfit_prs !=undef) & (bg_spd !=undef) & (amv_prs < 700.) & (amv_prs > 400.)]
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

    # set new background fields as if amv was moved to best fit pressure
    bg_u_new = bg_u.copy()
    bg_u_new[bfit_loc] = bfit_u[bfit_loc]
    bg_v_new = bg_v.copy()
    bg_v_new[bfit_loc] = bfit_v[bfit_loc]
    bg_spd_new = bg_spd.copy()
    bg_spd_new[bfit_loc] = bfit_spd[bfit_loc]

    # Compute OMB for U, V, Spd, Vector at AMV prs and best fit prs
    u_omb_bf = amv_u[bfit_loc] - bfit_u[bfit_loc]
    v_omb_bf = amv_v[bfit_loc] - bfit_v[bfit_loc]
    spd_omb_bf = amv_spd[bfit_loc] - bfit_spd[bfit_loc]
    vd_omb_bf = np.sqrt((amv_u[bfit_loc] - bfit_u[bfit_loc]) ** 2 + (amv_v[bfit_loc] - bfit_v[bfit_loc]) ** 2)

    u_omb_bf_orig = amv_u[bfit_loc] - bg_u[bfit_loc]
    v_omb_bf_orig = amv_v[bfit_loc] - bg_v[bfit_loc]
    spd_omb_bf_orig = amv_spd[bfit_loc] - bg_spd[bfit_loc]
    vd_omb_bf_orig = np.sqrt((amv_u[bfit_loc] - bg_u[bfit_loc]) ** 2 + (amv_v[bfit_loc] - bg_v[bfit_loc]) ** 2)

    u_omb_bg = amv_u[bg_u != undef] - bg_u[bg_u != undef]
    v_omb_bg = amv_v[bg_v != undef] - bg_v[bg_v != undef]
    spd_omb_bg = amv_spd[bg_spd != undef] - bg_spd[bg_spd != undef]
    vd_omb_bg = np.sqrt(
        (amv_u[bg_spd != undef] - bg_u[bg_spd != undef]) ** 2 + (amv_v[bg_spd != undef] - bg_v[bg_spd != undef]) ** 2)

    u_omb_new = amv_u[bg_u_new != undef] - bg_u_new[bg_u_new != undef]
    v_omb_new = amv_v[bg_v_new != undef] - bg_v_new[bg_v_new != undef]
    spd_omb_new = amv_spd[bg_spd_new != undef] - bg_spd_new[bg_spd_new != undef]
    vd_omb_new = np.sqrt((amv_u[bg_spd_new != undef] - bg_u_new[bg_spd_new != undef]) ** 2 + (
                amv_v[bg_spd_new != undef] - bg_v_new[bg_spd_new != undef]) ** 2)

    mean_vd_bg = np.mean(vd_omb_bg)
    std_vd_bg = np.std(vd_omb_bg)
    rms_vd_bg = np.sqrt(mean_vd_bg ** 2 + std_vd_bg ** 2)

    mean_vd_new = np.mean(vd_omb_new)
    std_vd_new = np.std(vd_omb_new)
    rms_vd_new = np.sqrt(mean_vd_new ** 2 + std_vd_new ** 2)

    ttest_data.append(vd_omb_bg)

    # Write out stats to a file
    # fo = open("amv_vd_stats.txt","a")

    if len(ttest_data) > 1:
        t, p = stats.ttest_ind(ttest_data[i-1], ttest_data[i], equal_var=True)
    else:
        t, p = undef, undef

    with open(f'stats_{source}_CQI.txt', 'a') as f:
        statsvd = f'{file_amv} CQI{qc_min[i]}-{qc_max[i]}: Total Number, Best Fit Number, VD OMB Mean, VD OMB STD, RMSE, t-test (VD OMB), p (VD OMB): ' \
                f'{amv_num} {bfit_num} {mean_vd_bg:.2f} {std_vd_bg:.2f} {rms_vd_bg:.2f} {t:.2f} {p:.4f}\n'
        f.write(statsvd)

