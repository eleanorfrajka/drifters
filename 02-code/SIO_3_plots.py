"""
SCRIPT 3
Plot drifter trajectories
- compare lat/lon with lowess lat/lon
- velocity for each drifter

Last modified: 6 Sep 2021
"""

# Import modules

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

import datetime
import os
import glob
import sys
import re

# Local import
from SIO_wrap import dir_tree, fnames

####################---------    FUNCTIONS    ---------####################

def plot_lowess(pid):
    """
    Plot lat/lon and the lowess filtered values
    """
    print("PID: %s" % str(pid))

    lon = hds.GPS_Longitude_deg[hds.Platform_ID==pid].values
    lat = hds.GPS_Latitude_deg[hds.Platform_ID==pid].values

    llon = hds.llon[hds.Platform_ID==pid].values
    llat = hds.llat[hds.Platform_ID==pid].values

    fig, ax = plt.subplots()
    ax.scatter(lon, lat, c='k')
    ax.plot(llon, llat, c='orange')

def plot_vel(pid, igood):
    """
    pid     : float, Platform ID
    igood   : int or float, good flag

    Plot velocity components.
    """

    print("PID: %s" % str(pid))

    llon = hds.llon[hds.Platform_ID==pid].values
    llat = hds.llat[hds.Platform_ID==pid].values

    u = hds.u[hds.Platform_ID==pid].values
    v = hds.v[hds.Platform_ID==pid].values
    vel_flag = hds.flag[hds.Platform_ID==pid].values

    flag_cond = vel_flag==igood

    for vel, vel_name in zip((u, v), ('u', 'v')):
        fig, ax = plt.subplots()
        cs_vel = ax.scatter(llon[flag_cond], llat[flag_cond],
                            c=vel[flag_cond], cmap='plasma', s=1, 
                            marker='s', vmin=-.2, vmax=.2)
        fig.colorbar(cs_vel, ax=ax, extend='both')
        ax.set_title('%s (m/s)' % vel_name) 

################################################################################
####################-----------   USER EDITS    ------------####################
# Path for the output data
data_dir = dir_tree.dir_out

igood = 1

#--------------
# Time formats
tstamp_strftime = '%Y%m%d'
timcol_strftime = '%Y-%m-%d %H:%M:%S'

# Reference date for computing time in seconds
ref_time = datetime.datetime(2000, 1, 1)

################################################################################
####################----------   load data    ----------####################

# Extract a list with the names of existing raw data files.
existing_files = glob.glob(os.path.join(data_dir, fnames.fname_data + '*'))

# ~ ~ print update ~ ~ 
if len(existing_files) > 0:
    print("Existing data files: \n%s\n" % existing_files)
else:
    sys.exit("No previous data files.\n")

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
ds_fname = f"{fnames.fname_data}{fname_timestamp}.nc"
ds_fpath = os.path.join(data_dir, ds_fname)

print("Opening file: %s\n" % ds_fpath)
hds = xr.open_dataset(ds_fpath)

# Extract unique values of the platform ID (to select separate drifters)
PID = list(set(hds.Platform_ID.values))

print("Index, Platform ID, Deployment Time")
for i in range(len(PID)):
    depl_time = hds.Deployment_date[hds.Platform_ID==PID[i]][0]
    depl_time = depl_time.dt.strftime(timcol_strftime).values
    print(i, int(PID[i]), depl_time)

indx = int(input("Enter the index of a drifter to plot trajectory: "))

plt.ion()
plot_lowess(PID[indx])
plot_vel(PID[indx], igood)