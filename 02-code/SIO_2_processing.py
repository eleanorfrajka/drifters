"""
SCRIPT 2
Edit raw data and compute velocity components.

Steps:
- discard data where GPS has failed (flag based on lat threshold) and outside 
the North Atlantic box 
- apply LOWESS (Locally Weighted Scatterplot Smoother) method to lat/lon 
(Elipot et al, 2016)
- regrid on an hourly grid:
    > 'linear' for continuous variables
    > 'nearest neighbour' for discrete ones
- compute velocity using lowess filtered lat/lon and apply a threshold
- save as netcdf file

Last modified: 5 Sep 2021
"""

####################---------   LOCAL FUNCTIONS    ---------####################

def num_ibad(var, ibad):
    """
    var  : xarray dataset that contains 'time' and 'flag' variables
    ibad : int or float; value of the (bad) flag

    Returns the total number of flagged points. It also print the total number
    of points and the number of flagged ones.

    <!> This function breaks if the flag and time variable names change.
    """
    all_npts = var.time.size
    ibad_npts = var.flag.where(var.flag==ibad, drop=True).size 

    print("Total # of points: %s" % str(all_npts))
    print("Total # of flagged points: %s" % str(ibad_npts))
    
    return int(ibad_npts)

def latlon_extremes(var, igood):
    """
    var   : xarray dataset that containg the 'flag' and lat/lon variables.
    igood : int or float; value of the good flag

    Prints the min/max of lat and lon.

    <!> This function breaks if the lat/lon variable names change.
    """
    # Select lat/lon where flag is good
    bool_cond = var.flag.values==igood

    print('lat min: ', var.GPS_Latitude_deg.values[bool_cond].min())
    print('lat max: ', var.GPS_Latitude_deg.values[bool_cond].max())
    print('lon min: ', var.GPS_Longitude_deg.values[bool_cond].min())
    print('lon max: ', var.GPS_Longitude_deg.values[bool_cond].max())


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import modules

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import datetime
import gsw
import os
import glob
import sys
import re

# Local import
from SIO_wrap import dir_tree, fnames
from SIO_wrap.lowess import LatLonLocalWess
from SIO_wrap import jlab_python as jlab

################################################################################
####################-----------   USER EDITS    ------------####################

# Path for the output data
data_dir = dir_tree.dir_out

#--------------
# Flag value - good=1, bad=4
ibad = 4 
igood = 1

# FLAG 1: Latitude values less than this threshold are flagged
latbad_threshold = 0

# FLAG 2: North Atlantic box; values outside this box are flagged
# Defined from -180 to 180 (lon) and -90 to 90 (lat)
na_lonlim = [-80, 30]
na_latlim = [40, 80]

#--------------
# Speed threshold; values outside [-3, 3] m/s are flagged
threshold_ms = 3 # in m/s

#--------------
# LOWESS params
poly_order = 1
bandwidth = 2

#--------------
# Time formats
tstamp_strftime = '%Y%m%d'  # Filename timestamp
timcol_strftime = '%Y-%m-%d %H:%M:%S'  # Convert text to datetime format 

# Reference date for computing time in seconds
# Can use an earlier time reference if data start before 2000
ref_time = datetime.datetime(2000, 1, 1)

#--------------
# List of variable names split between float/int types based on whether the 
# variables are continuous or discrete, respectively.
# <!> If the names of variables change, update the lists by printing a list of 
# all the names from the raw datafile: list(xarrayDataset.keys())

integ_vars = ['Drogue_cnts', 'GPS_HDOP', 'GPS_FixDelay', 'GPS_TTFF', 
            'GPS_NumSat', 'SBD_Transmit_Delay', 'SBD_Retries']

float_vars = ['GPS_Latitude_deg', 'GPS_Longitude_deg', 'SST_degC',
            'SLP_mB', 'Battery_volts', 'llon', 'llat']

################################################################################
####################----------   load raw data    ----------####################

# Extract a list with the names of existing raw data files.
existing_files = glob.glob(os.path.join(data_dir, fnames.fname_rawdata + '*'))

# ~ ~ print update ~ ~ 
if len(existing_files) > 0:
    print("Existing raw data files: \n%s\n" % existing_files)
else:
    sys.exit("No previous raw data files.\n")

# ~ ~ filenaming convention ~ ~
# If there are multiple files with raw data (i.e. non-updated datasets), select 
# the latest one updated.
# The file names are distinguished by the timestamp appended to the filename 
# and has <tstamp_strftime> format (see 'user edits' section).
# The data are cropped such that the last day is fully sampled (spans 0h-23h).
# The timestamp in the filename is the latest downloaded fully sampled day.

# Extract the timestamp part of the filename(s) in a list
tstamp = [date for file in existing_files 
            for date in re.findall("(\d{8})", file)]

# Convert to datetime and pick the most recent timestamp
tstamp_date = pd.to_datetime(tstamp, format=tstamp_strftime)
fname_timestamp = tstamp[tstamp_date.argmax()]

# Load the raw file with the latest timestamp
ds_fname = f"{fnames.fname_rawdata}{fname_timestamp}.nc"
ds_fpath = os.path.join(data_dir, ds_fname)

print("Opening file: %s\n" % ds_fpath)
ds_raw = xr.open_dataset(ds_fpath)

# Create field FLAG and start by labelling all data 'good' (flag=igood)
ds_raw["flag"] = ('n', igood * np.ones(ds_raw.time.shape, dtype=np.int8))

# Total number of points 
total_points = int(ds_raw.time.size)

################################################################################
####################-------------   FLAG 1    --------------####################
# Flag samples where lat < latbad_threshold. 
# >> There are some points where lat=-90 probably where the GPS failed

ds_raw.flag.values[ds_raw.GPS_Latitude_deg.values<latbad_threshold] = ibad

# Print update on flagged data
print("\n> > Flags - stage 1 < <\n")
num_flags1 = num_ibad(ds_raw, ibad)
num_percent1 = 100 * num_flags1 / total_points
print("Flagged data as percentage: %s\n" % str(num_percent1))

latlon_extremes(ds_raw, igood)

################################################################################
####################-------------   FLAG 2    --------------####################
# Flag data in San Diego/outside of North Atlantic

lat_logical = np.logical_or(ds_raw.GPS_Latitude_deg<min(na_latlim),
                            ds_raw.GPS_Latitude_deg>max(na_latlim))
ds_raw.flag.values[lat_logical] = ibad

lon_logical = np.logical_or(ds_raw.GPS_Longitude_deg<min(na_lonlim),
                            ds_raw.GPS_Longitude_deg>max(na_lonlim))
ds_raw.flag.values[lon_logical] = ibad

# Print update on flagged data
print("\n> > Flags - stage 1 + 2 < <\n")
num_flags2 = num_ibad(ds_raw, ibad)
num_percent2 = 100 * num_flags2 / total_points
print("Total flagged data as percentage: %s\n" % str(num_percent2))

num_flags2_diff = num_flags2 - num_flags1
print("Flagged data at this step only: %s\n" % str(num_flags2_diff))

latlon_extremes(ds_raw, igood)

# Velocity flags here to remove GPS data. - note from Eleanor


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Remove flagged data, otherwise it might affect the lowess filtering

print("Removing flagged data. \n")
#ds = ds_raw.where(ds_raw.flag==1, drop=True) # this is too slow

pd_raw = ds_raw.to_pandas()
pd2 = pd_raw[pd_raw.flag==1]

# convert back to xarray
ds = pd2.to_xarray()

print("Drop the <flag> field for now. \n")
ds1 = ds.drop_vars('flag')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print("Convert time to number of seconds relative to the reference time: %s\n"
    % str(ref_time))

dtime_sec = (pd.to_datetime(ds1.time.values) - ref_time).total_seconds()
ds1['time_seconds'] = ("n", dtime_sec)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find the unique list of Platform IDs (i.e. drifters)
PID = list(set(ds1.Platform_ID.values))

# Create entries for lowess lat/lon
ds1['llat'] = ('n', np.ones(len(ds1.time)))
ds1['llon'] = ('n', np.ones(len(ds1.time)))

# Empty lists for storing the hourly data and new column names
hourly_arr = []
arr_names = []

# Iterate through all the drifters 
for i in range(len(PID)):
    print("PID: %s" % str(PID[i]))

    drift_i = ds1.where(ds1.Platform_ID==PID[i], drop=True)

    # Change the coordinate to time instead of index n
    drift_ii = drift_i.assign_coords(time=drift_i.time)
    drift_ii = drift_ii.drop_vars('n')
    drift_ii = drift_ii.drop_vars('Platform_ID')

    # make sure time axis is ascending
    drift_ii = drift_ii.sortby('time', ascending=True)
    
    # <!> If lat/lon variable names changes this breaks
    lon = drift_ii.GPS_Longitude_deg.values
    lat = drift_ii.GPS_Latitude_deg.values
    time_sec = drift_ii.time_seconds.values

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # LOWESS: Locally Weighted Scatterplot Smoother (Elipot et al 2016)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print("Applying the lowess method to lat and lon...")
    llat, llon = LatLonLocalWess(time_sec, lon, lat, poly_order, bandwidth)

    drift_ii['llon'] = ('n', llon)
    drift_ii['llat'] = ('n', llat)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # hourly interpolation
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # List of variable names
    varnames_list = list(drift_ii.keys())

    # Sampled time in datetime format
    time = pd.to_datetime(drift_ii.time.values)

    # Round up the time to the nearest hour
    time_rounded = time.round('60min')
    # > start time
    time_start = time_rounded[0]
    # > for end time check the rounded value is not higher than sampled one
    if time_rounded[-1] == time[-1]:
        time_end = time_rounded[-1]
    else:
        time_end = time_rounded[-2]

    # Hourly time grid
    tgrid = pd.date_range(time_start, time_end, freq='1H')

    # Grid time in seconds (with same time ref as time_seconds field)
    tgrid_sec = (tgrid - ref_time).total_seconds()

    # Create temporary array to store hourly interp data
    r_len = len(integ_vars) + len(float_vars)
    c_len = len(tgrid)
    arr = np.ones((r_len, c_len))
    
    # Iterate through the variable names and set the interpolation method 
    # based on the dtype of the data
    k = 0
    for varname in varnames_list:
        if (varname in integ_vars) or (varname in float_vars):
            if varname in integ_vars:
                interp_method='nearest'
            elif varname in float_vars:
                interp_method = 'linear'

            var = drift_ii[varname].values
            fc_interp = interp1d(time_sec, var, interp_method, fill_value='extrapolate')
            arr[k, :] = fc_interp(tgrid_sec)
            
            # Store names of variables in the order they are saved in the array
            # These are the same for every iteration so only do it once
            if i == 0:
                arr_names.append(varname)

            k += 1
        else:
            print("Variables not interpolated hourly: < %s >\n" % varname)
 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # velocity calculation
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # find index of llat and llon
    llat_h = arr[arr_names.index("llat"), :] 
    llon_h = arr[arr_names.index("llon"), :]  

    u, v = jlab.latlon2uv(tgrid_sec, llat_h, llon_h)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Append fields to the hourly array (e.g. velocity and fields that don't 
    # need interpolating)
    # Array structure will look like (top-down): 
    # < Time, Deployment time, PlatID, arr_names[..], u, v >

    # 1. Add the Platform ID field
    plat_id = np.ones(len(tgrid)) * PID[i]
    arr2 = np.vstack((plat_id, arr))

    # 2. Create field 'deployment date'
    deployment = np.ones(len(tgrid)) * tgrid_sec[0]
    arr2 = np.vstack((deployment, arr2))

    # 3. Add the time field
    arr2 = np.vstack((tgrid_sec, arr2))
    

    # 4.1 Add the velocity component fields
    arr2 = np.vstack((arr2, u))
    arr2 = np.vstack((arr2, v))

    # 4.2 Add row of ones for velocity flags
    flag_vel = np.ones(len(tgrid))
    flag_vel[abs(u) > threshold_ms] = ibad 
    flag_vel[abs(v) > threshold_ms] = ibad 

    arr2 = np.vstack((arr2, flag_vel))
    
    if i == 0:
        hourly_arr = arr2

        arr2_names = ['Platform_ID'] + arr_names
        arr2_names = ['Deployment_sec'] + arr2_names
        arr2_names = ['time_seconds'] + arr2_names

        arr2_names = arr2_names + ['u', 'v', 'flag']
    else:
        hourly_arr = np.hstack((hourly_arr, arr2))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# convert to a dataset and add attributes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# xarray ds
hds_pd = pd.DataFrame(hourly_arr.T, columns=arr2_names)


hds = hds_pd.to_xarray()
# convert times to datetime

hds_time = [(datetime.timedelta(seconds=hds.time_seconds.values[j]) 
            + ref_time) for j in range(len(hds.index))]
hds_depl = [(datetime.timedelta(seconds=hds.Deployment_sec.values[j]) 
            + ref_time) for j in range(len(hds.index))]

hds['time'] = ('index', hds_time)
hds['Deployment_date'] = ('index', hds_depl)

# attributes for Dataset
hds.time_seconds.attrs["units"] = "seconds since %s" % str(ref_time.strftime(timcol_strftime))
hds.llon.attrs["long_name"] = "longitude_lowess"
hds.llat.attrs["long_name"] = "latitude_lowess"
hds.u.attrs["units"] = "m/s"
hds.v.attrs["units"] = "m/s"
hds.u.attrs["long_name"] = "zonal_velocity"
hds.v.attrs["long_name"] = "meridional_velocity"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# save dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Filename and path of dataset
data_fname = f"{fnames.fname_data}{fname_timestamp}.nc"
data_fpath = os.path.join(data_dir, data_fname)

# Save dataset to netcdf file
hds.to_netcdf(data_fpath)

print("File saved: %s" % data_fpath)
print("Script 2 finished. \n")
