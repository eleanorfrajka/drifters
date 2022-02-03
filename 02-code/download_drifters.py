import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import requests
import os
import glob
import yaml
#import re
#import sys
import datetime
from datetime import timedelta
from io import StringIO
#from bs4 import BeautifulSoup

# Local import 
# > Make sure SIO_wrap dir is on the same path as this script.
from SIO_wrap import dir_tree, fnames

from setdir import *

#=============================================================================
#. Local functions
#=============================================================================
def load_one_drifter(pidstr, dnld_config, download_start_date): 
    
    from SIO_wrap import fnames
    
    #--------------------------------
    # Configuration details
    #--------------------------------
    base_url = dnld_config['base_url']
    username = dnld_config['username']
    password = dnld_config['password']
    
    tstr = f'start_date={download_start_date}'
    pidwrapstr = '&platform_id='

    full_url = base_url+tstr+pidwrapstr+pidstr

    #--------------------------------
    # Request file
    #--------------------------------
    resp = requests.get(full_url, auth=(username, password))
    # Print the response code (200 is good. If you get something else, may be
    # a password problem)

    # To print content: resp.content
    aa = resp.content.decode("utf-8")
    data_df = pd.read_csv(StringIO(aa))
    
    #--------------------------------
    # Clean up column names
    #--------------------------------
    tmp = data_df.columns.str.strip()
    tmp = tmp.str.replace(" ", "", regex=True)
    tmp = tmp.str.replace("-", "_", regex=True)
    tmp = tmp.str.replace("(", "_", regex=True)
    tmp = tmp.str.replace(')', '', regex=True)
    data_df.columns = tmp
    
    #--------------------------------
    # Remove the </br> column (IF IT EXISTS)
    #--------------------------------
    if '</br>' in data_df.columns:
        data_df = data_df.drop(columns='</br>')
        # data_df.dtypes
    
    #--------------------------------
    # Convert the time column to a timestamp
    #--------------------------------
    time_colname = 'Timestamp_UTC'
    data_df[time_colname] = pd.to_datetime(data_df[time_colname],
                                       format=timcol_strftime) 
    
    # Prep to convert xarray
    data_df2 = data_df
    data_df2["time"] = data_df2["Timestamp_UTC"].values
    data_df2 = data_df2.set_index("time")
    data_df2 = data_df2.drop(columns='Timestamp_UTC')
    #--------------------------------
    # Convert to xarray
    #--------------------------------
    ds = data_df2.to_xarray()

    # Sort by time ascending
    ds = ds.sortby('time', ascending=True)

    return ds


#=============================================================================
#. Get user defined parmaeters from a config file
#=============================================================================
with open('download_config.yml') as f:
    dnld_config = yaml.safe_load(f)

# Download start date must have format yyyy-mm-dd. Default is set to the
# beginning of the TERIFIC project, i.e 2019-12-04.
download_start_date = dnld_config['download_start_date']
print("\nDefault download start date: %s\n" % download_start_date)

# How many data points need to be available in order to update the data file?
# For hourly data, 24 data points is a full day.
num_update = dnld_config['num_update']

if 0:
    keyslist = list(dnld_config.keys())
    for keyname in keyslist:
        print(dnld_config[keyname])
        tmp = dnld_config[keyname]
        print(type(tmp))
        exec(keyname+'='+tmp)


        
#=============================================================================
#. Should not need to change these
#=============================================================================

# Path where data are saved. Can be changed in file SIO_wrap/dir_tree.py
data_dir = dir_tree.dir_out

# String formatting for time for:
#   - the download url, 
#   - appending to the filename
# 	- the data time column, respectively.
url_strf = '%Y-%m-%d'
tstamp_strf = '%Y%m%d'
mtime_strf = '%Y-%m-%dT%H:%M:%S'
timcol_strftime = '%Y-%m-%d %H:%M:%S'

# Get the list of Platform IDs
PID = pd.read_csv(cat_proc_path('PID_list.txt'), header='infer', index_col=0)


# Compare when the file was last updated to now
todaystr = datetime.datetime.now().strftime(url_strf)


#=============================================================================
#. Run the update doesn't work on 300234068342280 and 300234066513050
#=============================================================================


counter = 0
#4 and 57
#------------------------------------------------------------------------
# Loop through all the PID in TERIFIC
#------------------------------------------------------------------------
for i in range(0,len(PID)):
    # Get a single platform ID from the full list
    PID1 = (PID["PID"].values)[i]
    pidstr = PID1.astype('str')
    print(str(counter)+'. '+pidstr)

    # Use the same file name with just pid and then the PID
    fname = 'pid'+str(PID1)+'*.nc'
    existing_files = glob.glob(cat_raw_path(fname))

    #------------------------------------------------------
    # If there are no existing data files, then download
    # all data and save to netcdf
    #------------------------------------------------------
    if (len(existing_files)==0) | (PID1==300234068342280):
        print('       - no existing files')
        ds_new = load_one_drifter(pidstr, dnld_config,
                                  download_start_date)
        #------------------------------------------------------
        # Create the attributes
        maxtime = ds_new.time.max().values
        maxtimestr = pd.to_datetime(maxtime).strftime(mtime_strf)
        
        PID1_from_data = ds_new.Platform_ID.values[-1]
        pidstr_from_data = PID1_from_data.astype('str')
        if not PID1_from_data==PID1:
            print('========== PID mismatch: '+pidstr+'~='+pidstr_from_data)
        attr_dict = {"Platform_ID": PID1_from_data,
                    "End Time": maxtimestr,
                    "Project": dnld_config['project_name'],
                    "Originator": dnld_config['operator_name'],
                    "Institution": dnld_config['institution_name'],
                    "Date created": todaystr,
                    }

        ds_new = ds_new.assign_attrs(attr_dict)

        #------------------------------------------------------
        # Save file to raw 
        fname = 'pid'+pidstr_from_data+'.nc'
        outfile_with_path = cat_raw_path(fname)
        ds_new.to_netcdf(cat_raw_path(fname), mode='w')
        print(str(counter)+'. pid('+pidstr_from_data+
              ') - New end:'+maxtimestr)
        print('---- PID: '+str(PID1_from_data)+' and filename: '+fname)
        counter += 1
        
    #------------------------------------------------------
    # If there are existing data files, check whether it
    # needs to be updated - download and append
    #------------------------------------------------------
    if len(existing_files) > 0:

        #------------------------------------------------
        # Load the newest file in the list
        latest_file = max(existing_files, key=os.path.getctime)
        ds_previous = xr.open_dataset(latest_file)
        PID1_from_data = ds_previous.attrs['Platform_ID']
        if not PID1_from_data==PID1:
            print('--------- MISMATCH: PID and ds_previous')

        #------------------------------------------------
        # Check whether the file needs to be updated
        maxtimestr_prev = ds_previous.attrs['End Time']
        maxtime_prev = datetime.datetime.strptime(maxtimestr_prev, mtime_strf)
        # End date (without time)
        download_update = maxtime_prev.strftime(url_strf)
        print('       - existing files, ending '+download_update)
        

        #------------------------------------------------
        # If the last date is the same as today, then do not proceed
        if not todaystr==download_update:
            #------------------------------------------------
            # Download drifter data since the previous date
            ds_update = load_one_drifter(pidstr, dnld_config, download_update)

            # If there were at least 24 data points since the last data file,
            # then append the PID to the update list
            if len(ds_update["time"]) > num_update:
                print('       - new data '+str(len(ds_update["time"])))

                # Remove platform ID as a long DataArray - 
                # it's now an attribute
#                ds_update = ds_update.drop('Platform_ID')
                PID1_from_data = ds_update.Platform_ID.values[-1]
                pidstr_from_data = PID1_from_data.astype('str')
                if not PID1_from_data==PID1:
                    print('========== PID mismatch: '
                          +pidstr+'~='+pidstr_from_data)


                # Merge the old and new data
                ds_new = xr.merge([ds_update, ds_previous], 
                                  combine_attrs="drop")

                #------------------------------------------------------
                # Create the attributes
                maxtime = ds_update.time.max().values
                maxtimestr = pd.to_datetime(maxtime).strftime(mtime_strf)

                attr_dict = {"Platform_ID": PID1,
                            "End Time": maxtimestr,
                            "Project": dnld_config['project_name'],
                            "Originator": dnld_config['operator_name'],
                            "Institution": dnld_config['institution_name'],
                            "Date created": todaystr,
                            }

                ds_new = ds_new.assign_attrs(attr_dict)

                #------------------------------------------------------
                # Save file to raw 
                fname = 'pid'+pidstr+'.nc'

                outfile_with_path = cat_raw_path(fname)
                ds_new.to_netcdf(cat_raw_path(fname), mode='w')
                print(str(counter)+'. pid('+pidstr+') - Ended:'\
                      +download_update
                      +', New end:'+maxtimestr)
                counter += 1
            else:
                print('       - no new data')




print('No more files to update')