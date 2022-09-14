
#!/usr/bin/env python

import os, os.path
import sys
import numpy as np
import math
from math import radians, sin, cos, acos
import pygrib


def read_txt(file,qc_pos,qc_min,qc_max,delimiter=None,lonflip=None,lon0=None):

#lonflip not equal to None will flip the sign on longitude
#lon0 two options, 0= -180 to 180 lon
#                  180 = 0 to 360 lon
#                  any other value, lon will not be changed

    undef = -9999.0


# count the valid AMV
    amv_num = 0
    with open(file,'r') as f:
        for line in f:
            try:
                amv_list=parse_amv(line,qc_pos,qc_min,qc_max,delimiter=delimiter)
                if (amv_list[0] != undef): 
                    amv_num +=1
            except:
                continue
    f.closed 

# allocate np arrays
    print('number of amv in read_txt = {0}'.format(amv_num))
    amv_spd = np.empty(amv_num)
    amv_dir = np.empty(amv_num)
    amv_prs = np.empty(amv_num)
    amv_lat = np.empty(amv_num)
    amv_lon = np.empty(amv_num)
    amv_qc  = np.empty(amv_num)
    
#TargetID;Latitude;Longitude;TSize;SSize;Speed;Direction;Height;LLC;ModelSpeed;ModelDir;Albedo;MaxCorr;TM;HeightError;HAM;QI;QIF;CQI

# read again to fill np arrays
    count = 0
    with open(file,'r') as f:
        for line in f:
            try:
                amv_list=parse_amv(line,qc_pos,qc_min,qc_max,lonflip=lonflip,lon0=lon0,delimiter=delimiter)
                if (amv_list[0] != undef): 
                    amv_spd[count] = amv_list[0]
                    amv_dir[count] = amv_list[1]
                    amv_lat[count] = amv_list[2]
                    amv_lon[count] = amv_list[3]
                    amv_prs[count] = amv_list[4]
                    amv_qc[count] = amv_list[5]
                    count +=1
            except:
                continue
    f.closed 

# place in one variable for convenience
    amv_data = np.vstack((amv_spd,amv_dir,amv_lat,amv_lon,amv_prs,amv_qc))

    return amv_data



def parse_amv(line,qc_pos,qc_min,qc_max,lonflip=None,lon0=None,delimiter=None):

# contain format of text file and data checking in one function

    undef = -9999.0
    #print('in parse delimiter={0}'.format(delimiter))

#TargetID;Latitude;Longitude;TSize;SSize;Speed;Direction;Height;LLC;ModelSpeed;ModelDir;Albedo;MaxCorr;TM;HeightError;HAM;QI;QIF;CQI
#2018 intercomparison

#last intercomparison
    tlat = 1 ; vlat=[-61.,61.]
    #tlon = 2 ; vlon=[-61.,61.]
    tlon = 2 ; vlon=[-180.,360.]
    tspd = 5 ; vspd=[0.,150.]
    tdir = 6 ; vdir=[0.,361.]
    tprs = 7 ; vprs=[10.,1020.]
    tqc  = qc_pos ; vqc =[qc_min,qc_max]

    if delimiter==None:
        token=line.split()
    else:
        token=line.split(delimiter)
 
    lat=undef
    lon=undef
    spd=undef
    dir=undef
    prs=undef
    qc =undef
    #print('lat {0} lon {1} spd {2} dir {3} prs {4} qc {5}'.format(token[tlat],token[tlon],token[tspd],token[tdir],token[tprs],token[qc_pos]))

# basic valid data check, could add QI or other flag check here
    valid=True
    if (float(token[tlat]) < vlat[0]) or (float(token[tlat]) > vlat[1]):
        valid=False
    if (float(token[tlon]) < vlon[0]) or (float(token[tlon]) > vlon[1]):
        valid=False
    if (float(token[tspd]) < vspd[0]) or (float(token[tspd]) > vspd[1]):
        valid=False
    if (float(token[tdir]) < vdir[0]) or (float(token[tdir]) > vdir[1]):
        valid=False
    if (float(token[tprs]) < vprs[0]) or (float(token[tprs]) > vprs[1]):
        valid=False
    if (float(token[tqc])  < vqc[0])  or (float(token[tqc])  > vqc[1]):
        valid=False


    if (valid):
        lat = float(token[tlat]) 
        lon = float(token[tlon])
        spd = float(token[tspd])
        dir = float(token[tdir])
        prs = float(token[tprs])
        qc  = float(token[tqc])
    if lonflip!=None:
        lon = -lon
    if lon0==0. and lon>180.:
        lon=lon-360.
    if lon0==180. and lon<0.:
        lon=lon+360.

    amv_list = [spd,dir,lat,lon,prs,qc]

    return amv_list



def write_txt(file,amv_data):

    amv_num  = np.size(amv_data,axis=1)

    with open(file,'w') as f:
        out_string=file+" spd,dir,lat,lon,prs"
        f.write(out_string)
        f.write("\n")
        i = 0
        while (i<amv_num):
            out_string="{0} {1} {2} {3} {4} {5}".format(amv_data[0,i],amv_data[1,i],amv_data[2,i],amv_data[3,i],amv_data[4,i],amv_data[5,i])
            f.write(out_string)
            f.write("\n")

            i +=1
    f.closed 

    return 

def read_DecodedForecast_MSG(file):

    """Usage: DecFcst = read_DecodedForecast_MSG(file)
    where, e.g.:
    file = '[...]/DecodedForecast_20120503070606Z_20120503120000Z_12_V_MPFS07'
    """

#   	Decoded forecast header (4,776 bytes)
#	Decoded_Forecast_Header = BYTARR(4776)

#   	Decoded forecast data point (40 bytes)
#	Decoded_Forecast_Point = {Latitude      : FLOAT(0), $
#                             Longitude     : FLOAT(0), $
#                             Pressure      : FLOAT(0), $
#                             Geopotential  : FLOAT(0), $
#                             Temperature   : FLOAT(0), $
#                             WVMixingRatio : FLOAT(0), $
#                             O3MixingRatio : FLOAT(0), $
#                             WindSpeed     : FLOAT(0), $
#                             WindDirection : FLOAT(0), $
#                             DewPointTemp  : FLOAT(0)}
#   Decoded forecast data array (29,160,000 bytes)
#   Decoded_Forecast_Array = REPLICATE(Decoded_Forecast_Point, 40L * 135L * 135L)

# This could input defined
    num_x = 135   
    num_y = 135
    num_var = 10
    num_lev = 40
    num_pt = num_x * num_y
    arrout = np.zeros((num_x,num_y,10,40))  

#   Open file and skip header
    fo = open(file)
    fo.seek(4776)

    print('Program is reading decoded forecast data:{0}'.format(file))

#   '>f4' is float which is big-endian
    forecast=np.fromfile(fo,dtype=('>f4'))
    fo.close()
    fore=np.reshape(forecast,(num_pt,num_lev,num_var))

    latitude      = fore[:,:,0]
    longitude     = fore[:,:,1]
    pressure      = fore[:,:,2]
    geopotential  = fore[:,:,3]
    temperature   = fore[:,:,4]
    wvmixingratio = fore[:,:,5]
    o3mixingratio = fore[:,:,6]
    windspeed     = fore[:,:,7]
    winddirection = fore[:,:,8]
    dewpointtemp  = fore[:,:,9]
	
    for j in range(num_pt):

#hardwired for this input file
        x = np.floor(longitude[j,0] + 67)
        y = np.floor(latitude[j,0] + 67)

        arrout[x,y,0,:] = pressure[j,:]
        arrout[x,y,1,:] = geopotential[j,:]
        arrout[x,y,2,:] = temperature[j,:]
        arrout[x,y,3,:] = wvmixingratio[j,:]
        arrout[x,y,4,:] = o3mixingratio[j,:]
        arrout[x,y,5,:] = windspeed[j,:]
        arrout[x,y,6,:] = winddirection[j,:]
        arrout[x,y,7,:] = dewpointtemp[j,:]
        arrout[x,y,8,:] = latitude[j,:]
        arrout[x,y,9,:] = longitude[j,:]

    return arrout


def locate(amv_data,fcst_data):

    amv_num  = np.size(amv_data,axis=1)
    grid_i = np.zeros(amv_num,dtype=int)
    grid_j = np.zeros(amv_num,dtype=int)
    fcst_lat  = fcst_data[:,:,8,0]
    fcst_lon  = fcst_data[:,:,9,0]
    #toplon=361.
    #max_lon = max(fcst_lon.all,toplon)
    #print('Longitudes'.format(max_lon))
    numx = np.size(fcst_data,axis=0)
    numy = np.size(fcst_data,axis=1)

    n = 0
    while (n<amv_num):
        #print('n={0}'.format(n))
        amv_lat  = amv_data[2,n]
        amv_lon  = amv_data[3,n] 
        grid_diff = (amv_lat-fcst_lat)**2 + (amv_lon-fcst_lon)**2
        indx= np.argmin(grid_diff)
        indx_2d=np.unravel_index(indx,(numx,numy))
        grid_i[n] = indx_2d[0]
        grid_j[n] = indx_2d[1]
        n +=1
    
    return grid_i,grid_j 



def bestfit(amv_data,fcst_data,verbose):

    """Finds the background model best fit pressure associated with the AMV.
       The model best-fit pressure is the height (in pressure units) where the
       vector difference between the observed AMV and model background is a
       minimum.  This calculation may only work approximately 1/3 of the time.

       Reference:
       Salonen et al (2012), "Characterising AMV height assignment error by 
       comparing best-fit pressure statistics from the Met Office and ECMWF 
       System."  Proceedings of the 11th International Winds Workshop, 
       Auckland, New Zealand, 20-24 February 2012.

       Input contained in amv_data and fcst_data:
       amv_spd -  AMV speed m/s 
       amv_dir -  AMV direction deg
       amv_prs -  AMV pressure hPa
       fcst_spd - (level) forecast speed m/s
       fcst_dir - (level) forecast direction (deg)
       fcst_prs - (level) forecast pressure (hPa)

       Output contained in bf_data:
       SatwindBestFitU - AMV best fit U component m, unconstrained value is undef
       SatwindBestFitV - AMV best fit V component m, unconstrained value is undef
       bfit_prs - AMV best fit pressure m/s, unconstrained value is undef
       flag - 0 found, 1 not contrained, 2 vec diff minimum not met, 3 failed to find suitable fcst pressure match

       History:
       10/2012 - Steve Wanzong - Created in Fortran
       10/2013 - Sharon Nebuda - rewritten for python
    """
    undef = -9999.0

    amv_spd = amv_data[0]
    amv_dir = amv_data[1]
    amv_prs = amv_data[4]
    amv_lat = amv_data[2]
    amv_lon = amv_data[3]

    fcst_spd = fcst_data[5,:]
    fcst_dir = fcst_data[6,:]
    fcst_prs = fcst_data[0,:]
    fcst_lat = fcst_data[8,0]
    fcst_lon = fcst_data[9,0]

    fcst_num_levels = fcst_spd.shape[0]

#   verbose = True
#   verbose = False

    SatwindBestFitPress = undef
    SatwindBestFitU = undef
    SatwindBestFitV = undef

    PressDiff = 150.                      # pressure above and below AMV to look for fit
    TopPress = 50.                        # highest level to allow search

    flag = 3
    bf_data = np.vstack((undef,undef,undef,flag))

    if (amv_prs<TopPress):
        if (verbose):
            print('AMV location lat,lon,prs ({0},{1},{2}) is higher than pressure {3}'.format(amv_lat,amv_lon,amv_prs,TopPress))
        return bf_data

#Calculate the pressure +/- 150 hPa from the AMV pressure.
    PressMax = amv_prs + PressDiff
    PressMin = max((amv_prs-PressDiff),TopPress)

#1d array of indicies to consider for best fit location
    kk = np.where((fcst_prs<PressMax) & (fcst_prs>PressMin))
    if (len(kk[0]) ==0):
        if (verbose):
            print('AMV location lat,lon,prs ({0},{1},{2}) failed to find fcst prs around AMV'.format(amv_lat,amv_lon,amv_prs))
        return bf_data

#Diagnostic field: Find the model minimum speed and maximum speed within PressDiff of the AMV.
    if (verbose):
        SatwindMinSpeed = min(fcst_spd[kk])
        SatwindMaxSpeed = max(fcst_spd[kk])

#Compute U anv V for both AMVs and forecast
    amv_uwind = -amv_spd * np.sin(math.radians(amv_dir))
    amv_vwind = -amv_spd * np.cos(math.radians(amv_dir))
#   fcst_uwind = -fcst_spd[:] * np.sin(math.radians(fcst_dir[:]))
#   fcst_vwind = -fcst_spd[:] * np.cos(math.radians(fcst_dir[:]))
    dr=0.017453
    fcst_uwind = -fcst_spd * np.sin(dr*fcst_dir)
    fcst_vwind = -fcst_spd * np.cos(dr*fcst_dir)

#Calculate the vector difference between the AMV and model background at all levels. 
    VecDiff = np.sqrt((amv_uwind - fcst_uwind) ** 2 + (amv_vwind - fcst_vwind) ** 2)

#Find the model level of best-fit pressure, from the minimum vector difference.
    MinVecDiff = min(VecDiff[kk])
    imin=-1
    for i, item in enumerate(VecDiff):
        if MinVecDiff == VecDiff[i]:
            if i in kk[0]:
                imin = i

    if (imin ==-1 ):
        if (verbose):
            print('AMV location lat,lon,prs ({0},{1},{2}) failed to find min vector difference in layers around AMV'.format(amv_lat,amv_lon,amv_prs))

        return bf_data


#Use a parabolic fit to find the best-fit pressure.
#p2 - Minimized model pressure at level imin (hPa)
#v2 - Minimized vector difference at level imin (m/s)
#p1 - 1 pressure level lower in atmosphere than p2
#p3 - 1 pressure level above in atmosphere than p2
#v1 - Vector difference 1 pressure level lower than p2
#v3 - Vector difference 1 pressure level above than p2

    p2 = fcst_prs[imin]
    v2 = VecDiff[imin]

# assumes fcst data level 0 at surface and (fcst_num_levels-1) at model top
#if bottom model level
    if imin == 0:
        SatwindBestFitPress = p2
    else:
        p3 = fcst_prs[imin+1]
        p1 = fcst_prs[imin-1]
        v3 = VecDiff[imin+1]
        v1 = VecDiff[imin-1]

#if top of allowed region
        if p3 < TopPress:
            SatwindBestFitPress = p2

#check not collinear
        elif (v1 != v2 and v2 != v3):
            SatwindBestFitPress = p2 - (0.5 * 
            ((((p2 - p1) * (p2 - p1) * (v2 - v3)) - ((p2 - p3) * (p2 - p3) * (v2 - v1))) / 
            (((p2 - p1) * (v2 - v3)) - ((p2 - p3) * (v2 - v1)))))
            if (SatwindBestFitPress < p3) or (SatwindBestFitPress > p1):
                if (verbose):
                    print('Best Fit not found between two pressure layers')
                    print('SatwindBestFitPress {0} p1 {1} p2 {2} p3 {3} imin {4}'.format(SatwindBestFitPress,p1,p2,p3,imin))
                SatwindBestFitPress = p2
        else:
            SatwindBestFitPress = p2

#Find best fit U and V by linear interpolation.
    if p2 == SatwindBestFitPress:
        SatwindBestFitU = fcst_uwind[imin]
        SatwindBestFitV = fcst_vwind[imin]
    else:
        if p2 < SatwindBestFitPress:
            LevBelow = imin - 1
            LevAbove = imin
            Prop = (SatwindBestFitPress - p1) / (p2 - p1)
        else:
            LevBelow = imin
            LevAbove = imin + 1      
            Prop = (SatwindBestFitPress - p2) / (p3 - p2)

        SatwindBestFitU = fcst_uwind[LevBelow] * (1.0 - Prop) + fcst_uwind[LevAbove] * Prop
        SatwindBestFitV = fcst_vwind[LevBelow] * (1.0 - Prop) + fcst_vwind[LevAbove] * Prop

    

# Check to see if the best fit pressure is constrained.

    SatwindGoodConstraint = 0
    flag = 2

    if MinVecDiff <= 4.0:
      
        SatwindGoodConstraint = 1
        flag = 1
        
        for ilev in range(fcst_num_levels):
            if fcst_prs[ilev] >= TopPress: 
          
                if ((fcst_prs[ilev] < (SatwindBestFitPress - 100.)) or  \
                   (fcst_prs[ilev] > (SatwindBestFitPress + 100.))) and  \
                   (VecDiff[ilev] <= (MinVecDiff + 2.0)):
                   SatwindGoodConstraint = 0
           
    if SatwindGoodConstraint == 1:
        bfit_prs = SatwindBestFitPress
        bfit_u = SatwindBestFitU
        bfit_v = SatwindBestFitV
        flag = 0
    else:
        bfit_prs = undef
        bfit_u = undef
        bfit_v = undef


    if (verbose):
        print('*** AMV best-fit ***')
        print('AMV -> p/minspd/maxspd: {0} {1} {2}'.format(amv_prs,SatwindMinSpeed,SatwindMaxSpeed))
        print('Bestfit -> p1,p2,p3,v1,v2,v3: {0} {1} {2} {3} {4} {5}'.format(p1,p2,p3,v1,v2,v3))
        print('Bestfit -> pbest,bfu,bfv,amvu,amvv,bgu,bgv: {0} {1} {2} {3} {4} {5} {6}'.format(
        SatwindBestFitPress,SatwindBestFitU,SatwindBestFitV,amv_uwind,amv_vwind,fcst_uwind[imin],fcst_vwind[imin]))
        print('Good Constraint: {0}'.format(SatwindGoodConstraint))
        print('Minimum Vector Difference: {0}'.format(VecDiff[imin]))
        print('Vector Difference Profile: ')
        print(VecDiff)
        print('Pressure Profile: ')
        print(fcst_prs)

        if (abs(SatwindBestFitU - amv_uwind) > 4.0) or  (abs(SatwindBestFitV - amv_vwind) > 4.0):
            print('U Diff: {0}'.format(abs(SatwindBestFitU - amv_uwind)))
            print('V Diff: {0}'.format(abs(SatwindBestFitV - amv_vwind)))

    bf_data = np.vstack((bfit_u,bfit_v,bfit_prs,flag))

    return bf_data




def bg(amv_data,fcst_data,verbose):

    """Finds the background model U and V components at the AMV pressure level
       5/2014 - Sharon Nebuda 
    """
    undef = -9999.0

    amv_prs = amv_data[4]

    fcst_spd = fcst_data[5,:]
    fcst_dir = fcst_data[6,:]
    fcst_prs = fcst_data[0,:]

    fcst_num_levels = fcst_spd.shape[0]

    bg_u = undef
    bg_v = undef

    bg_data = np.vstack((bg_u,bg_v))

#Compute U anv V for forecast
#   fcst_uwind = -fcst_spd[:] * np.sin(math.radians(fcst_dir[:]))
#   fcst_vwind = -fcst_spd[:] * np.cos(math.radians(fcst_dir[:]))
    dr=0.017453
    fcst_uwind = -fcst_spd * np.sin(dr*fcst_dir)
    fcst_vwind = -fcst_spd * np.cos(dr*fcst_dir)

    if (verbose):
        print(fcst_uwind)

# assumes fcst data level 0 at surface and (fcst_num_levels-1) at model top
# data not well behaved (level 0 = level 1, level at top = 0 or -9999)
# hard wired this search 
#   if (fcst_prs[1] < fcst_prs[0]):
    k = 0
    while (k < (fcst_num_levels-1)):

        if ((fcst_prs[k] >= amv_prs) and (amv_prs >= fcst_prs[k+1])):

            LevBelow = k
            LevAbove = k + 1      
            Prop = (amv_prs - fcst_prs[LevAbove]) / (fcst_prs[LevBelow] - fcst_prs[LevAbove])
            
            bg_u = fcst_uwind[LevBelow] * Prop + fcst_uwind[LevAbove] * (1.-Prop)
            bg_v = fcst_vwind[LevBelow] * Prop + fcst_vwind[LevAbove] * (1.-Prop)
            k=fcst_num_levels
            if (verbose):
                print('{0} {1} {2} {3} {4}'.format(k,fcst_uwind[LevBelow],fcst_uwind[LevAbove],Prop,bg_u))

        k+=1


# assumes fcst data level (fcst_num_levels-1) at surface and 0 at model top
#   else: 
#       print('k=0 prs {0} k=1 prs {1}'.format(fcst_prs[0],fcst_prs[1])
#       k = fcst_num_levels-1
#       while (k > 0 and fcst_prs[k] > -999.):
#   
#           if ((fcst_prs[k] >= amv_prs) and (amv_prs >= fcst_prs[k-1])):
#   
#               LevBelow = k 
#               LevAbove = k - 1
#               Prop = (amv_prs - fcst_prs[LevAbove]) / (fcst_prs[LevAbove] - fcst_prs[LevBelow])
#   
#               bg_u = fcst_uwind[LevBelow] * Prop + fcst_uwind[LevAbove] * (1.-Prop)
#               bg_v = fcst_vwind[LevBelow] * Prop + fcst_vwind[LevAbove] * (1.-Prop)
#               print('shouldnt be here {0} {1} {2} {3} {4}'.format(k,fcst_prs[k],amv_prs,bg_u,bg_v))
#   
#           k-=1


#   if (bg_u == undef):
#       print('Did not find U {0} {1} {2}'.format(amv_prs, fcst_prs[0], fcst_prs[fcst_num_levels-1])
    bg_data = np.vstack((bg_u,bg_v))

    return bg_data


def spddir(ucomp,vcomp):

    """Computes speed and direction (wind barb convention) from U and V wind components
    """

    undef = -9999.

    if (ucomp != undef and vcomp !=undef):
        speed = np.sqrt( (ucomp*ucomp) + (vcomp*vcomp) )
        r2d = 180./np.pi
        direction = 90. - np.arctan(vcomp/ucomp) * r2d
    else: 
        speed = undef
        direction = undef

    return speed, direction


def uvcomp(speed,direction):

    """Computes U and V components from speed and direction
    """
    undef = -9999.

    if (speed != undef and direction !=undef):
        uwind = -speed * np.sin(math.radians(direction))
        vwind = -speed * np.cos(math.radians(direction))
    else:
        uwind = undef 
        vwind = undef

    return uwind, vwind


def read_MSG_grib(file,datetime=None):

    """Usage: Fcst, timeout = read_MSG_grib(file,time=time)  
       Returns a pressure, lat, lon, speed, direction in a 
       padded array for backwards compatibility
       time in YYYYMMDDHH format, long integer
       number of levels and number of time samples is hardwired
       if time is not specified, first time is returned
       if time is does not match file times, first time is returned
    """

    grbs=pygrib.open(file)
    uall=grbs.select(name='U component of wind')
    vall=grbs.select(name='V component of wind')
    #or grb=grbs.select(name='Temperature')[0]  # for one 2d grib message
    #without [0] get all levels and all times in a list of numpy arrays of complex data type
    #print(grb[0].keys())
    #print(grb[0].level)
    #print(grb[0].year)
    #print(grb[0])
    
    nlev = 37
    ntim = 24
    
    dt_file = np.zeros(ntim,dtype=int)
    for n in range(ntim):
        nm = n*nlev  # grib message number
        date = uall[nm].dataDate
        time = uall[nm].hour
        #string_datetime = '{0}{1:2d}'.format(date1,time1)
        dt_file[n] = date*100 + time

    nm=1
    dt_out = dt_file[0]
    if datetime!=None:
        for n in range(ntim):
            if datetime == dt_file[n]:
                nm = n
                dt_out = dt_file[n]
        
    nstart = nlev * nm

    prs = np.zeros(nlev)
    for n in range(nlev):
        prs[n] = uall[n].level
    
    lats, lons = uall[0].latlons()
    dims = lats.shape
    
    nlon = dims[0]
    nlat = dims[1]
    
    fcst_data=np.zeros((nlon,nlat,10,nlev))
    uwnd=np.zeros((nlon,nlat))
    vwnd=np.zeros((nlon,nlat))
    wspd=np.zeros((nlon,nlat))
    wdir=np.zeros((nlon,nlat))
    for n in range(nlev):
        ns = nstart + n
        nl = nlev - n - 1 # flip so data is surface to model top
        uwnd[:,:]=uall[ns].values  # same as uall[ns]['values']
        vwnd[:,:]=vall[ns].values  
        wspd = np.square(uwnd) + np.square(vwnd)
        wspd = np.sqrt(wspd) 
        wdir = 270.-np.arctan2(vwnd,uwnd)*180./np.pi
        fix = np.where(wdir>= 360)
        wdir[fix] = wdir[fix] - 360.

        fcst_data[:,:,0,nl]=prs[n]
        fcst_data[:,:,5,nl]=wspd[:,:]
        fcst_data[:,:,6,nl]=wdir[:,:]
        fcst_data[:,:,8,nl]=lats[:,:]
        fcst_data[:,:,9,nl]=lons[:,:]
    
    grbs.close()

    return fcst_data, dt_out
    
# old decoder
#       arrout[x,y,0,:] = pressure[j,:]
#       arrout[x,y,1,:] = geopotential[j,:]
#       arrout[x,y,2,:] = temperature[j,:]
#       arrout[x,y,3,:] = wvmixingratio[j,:]
#       arrout[x,y,4,:] = o3mixingratio[j,:]
#       arrout[x,y,5,:] = windspeed[j,:]
#       arrout[x,y,6,:] = winddirection[j,:]
#       arrout[x,y,7,:] = dewpointtemp[j,:]
#       arrout[x,y,8,:] = latitude[j,:]
#       arrout[x,y,9,:] = longitude[j,:]



#uwnd, vwnd, prs, lons, lats, timeout = amv.read_MSG_grib_uv(file,datetime=datetime)
def read_MSG_grib_uv(file,datetime=None):

    """Usage: uwnd,vwnd,prs,lons,lats, timeout = read_MSG_grib_uv(file,time=time)  
       Returns a pressure, lat, lon, speed, direction in a 
       padded array for backwards compatibility
       time in YYYYMMDDHH format, long integer
       number of levels and number of time samples is hardwired
       if time is not specified, first time is returned
       if time is does not match file times, first time is returned
    """

    grbs=pygrib.open(file)
    uall=grbs.select(name='U component of wind')
    vall=grbs.select(name='V component of wind')
    #or grb=grbs.select(name='Temperature')[0]  # for one 2d grib message
    #without [0] get all levels and all times in a list of numpy arrays of complex data type
    #print(grb[0].keys())
    #print(grb[0].level)
    #print(grb[0].year)
    #print(grb[0])
    
    nlev = 37
    ntim = 24
    
    dt_file = np.zeros(ntim,dtype=int)
    for n in range(ntim):
        nm = n*nlev  # grib message number
        date = uall[nm].dataDate
        time = uall[nm].hour
        #string_datetime = '{0}{1:2d}'.format(date1,time1)
        dt_file[n] = date*100 + time

    nm=1
    dt_out = dt_file[0]
    if datetime!=None:
        for n in range(ntim):
            if datetime == dt_file[n]:
                nm = n
                dt_out = dt_file[n]
        
    nstart = nlev * nm

    prs = np.zeros(nlev)
    for n in range(nlev):
        prs[n] = uall[n].level
    
    lats, lons = uall[0].latlons()
    dims = lats.shape
    
    nlon = dims[0]
    nlat = dims[1]
    
    uwnd=np.zeros((nlon,nlat,nlev))
    vwnd=np.zeros((nlon,nlat,nlev))
    for n in range(nlev):
        ns = nstart + n
        nl = nlev - n - 1 # flip so data is surface to model top
        uwnd[:,:,nl]=uall[ns].values  # same as uall[ns]['values']
        vwnd[:,:,nl]=vall[ns].values  

    grbs.close()

    return uwnd,vwnd,prs,lons,lats, dt_out
    
# old decoder
#       arrout[x,y,0,:] = pressure[j,:]
#       arrout[x,y,1,:] = geopotential[j,:]
#       arrout[x,y,2,:] = temperature[j,:]
#       arrout[x,y,3,:] = wvmixingratio[j,:]
#       arrout[x,y,4,:] = o3mixingratio[j,:]
#       arrout[x,y,5,:] = windspeed[j,:]
#       arrout[x,y,6,:] = winddirection[j,:]
#       arrout[x,y,7,:] = dewpointtemp[j,:]
#       arrout[x,y,8,:] = latitude[j,:]
#       arrout[x,y,9,:] = longitude[j,:]


# https://medium.com/@petehouston/calculate-distance-of-two-locations-on-earth-using-python-1501b1944d97
def great_circle(lat1, lon1, lat2, lon2):
    """Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)"""

    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    return 6371 * (
        acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
    )

