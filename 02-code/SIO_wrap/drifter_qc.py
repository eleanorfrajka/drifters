import numpy as np
import pandas as pd

##################---------   STATUS FUNCTIONS    ---------##################

def latlon_extremes(var, igood):
    """
    var   : xarray dataset that containg the 'flag' and lat/lon variables.
    igood : int or float; value of the good flag

    Prints the min/max of lat and lon.

    <!> This function breaks if the lat/lon variable names change.
    """
    # Select lat/lon where flag is good
    bool_cond = var.flag.values==igood

    if np.count_nonzero(bool_cond):
        lonmin = var.GPS_Longitude_deg.values[bool_cond].min().astype('str')
        latmin = var.GPS_Latitude_deg.values[bool_cond].min().astype('str')
        lonmax = var.GPS_Longitude_deg.values[bool_cond].max().astype('str')
        latmax = var.GPS_Latitude_deg.values[bool_cond].max().astype('str')
    
        print('Lon in ('+lonmin+', '+lonmax+'), Lat in ('+latmin+', '+latmax+')')
        
    else:
        print('No good values remaining')
#    print('lat min: ', var.GPS_Latitude_deg.values[bool_cond].min())
#    print('lat max: ', var.GPS_Latitude_deg.values[bool_cond].max())
#    print('lon min: ', var.GPS_Longitude_deg.values[bool_cond].min())
#    print('lon max: ', var.GPS_Longitude_deg.values[bool_cond].max())
    

##################---------   QC FUNCTIONS    ---------##################

def flag_gps_dropout(ds_raw, varname, mode_val, badflag):
    # Find where the variable has the value given by mode_val, and flag them
    ds_raw.flag.values[ds_raw[varname]==mode_val] = badflag
    
    numflags = num_ibad(ds_raw, badflag)
    return ds_raw, numflags

def flag_gps_region(ds_raw, varname, varlim, badflag):
    # Find where the variable is outside the region given by
    # varlim
    lat_logical = np.logical_or(ds_raw[varname]<min(varlim),
                                ds_raw[varname]>max(varlim))
    ds_raw.flag.values[lat_logical] = badflag
    numflags = num_ibad(ds_raw, badflag)
    
    return ds_raw, numflags


def flag_vel(ds_raw, varu, varv, varmax, badflag):
    
    uvel = ds_raw[varu]
    vvel = ds_raw[varv]
    velmag = np.sqrt(np.square(uvel) + np.square(vvel))

    ds_raw.flag.values[velmag > varmax] = badflag
    
    numflags = num_ibad(ds_raw, badflag)
    
    return ds_raw, numflags

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

#    print("Total # of points: %s" % str(all_npts))
#    print("Total # of flagged points: %s" % str(ibad_npts))
    
    return int(ibad_npts)


def create_hourly_grid(time_da, ref_time):

    # Compute the hourly time grid
    # Create the regular time grid
    time_as_series = time_da.to_series()
    dt = time_as_series.diff().mean()

    # Sampled time in datetime format
    time1 = pd.to_datetime(time_da.values)

    # Round up the time to the nearest hour
    time_rounded = time1.round('60min')
    # > start time
    time_start = time_rounded[0]
    # > for end time check the rounded value is not higher than sampled one
    if time_rounded[-1] == time1[-1]:
        time_end = time_rounded[-1]
    else:
        time_end = time_rounded[-2]

    # Hourly time grid
    tgrid = pd.date_range(time_start, time_end, freq='1H')
    tgrid_sec = (tgrid - ref_time).total_seconds()

    return tgrid, tgrid_sec

