import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats
from SIO_wrap import jlab_python as jlab
from SIO_wrap.lowess import LatLonLocalWess
from scipy.interpolate import interp1d





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

def flag_gps_dropout(ds_qc, varname, mode_val, badflag):
    # Find where the variable has the value given by mode_val, and flag them
    ds_qc.flag.values[ds_qc[varname]==mode_val] = badflag
    
    numflags = num_ibad(ds_qc, badflag)
    return ds_qc, numflags

def flag_gps_region(ds_qc, varname, varlim, badflag):
    # Find where the variable is outside the region given by
    # varlim
    lat_logical = np.logical_or(ds_qc[varname]<min(varlim),
                                ds_qc[varname]>max(varlim))
    ds_qc.flag.values[lat_logical] = badflag
    numflags = num_ibad(ds_qc, badflag)
    
    return ds_qc, numflags


def flag_vel(ds_qc, varu, varv, varmax, badflag):
    
    uvel = ds_qc[varu]
    vvel = ds_qc[varv]
    velmag = np.sqrt(np.square(uvel) + np.square(vvel))

    ds_qc.flag.values[velmag > varmax] = badflag
    
    numflags = num_ibad(ds_qc, badflag)
    
    return ds_qc, numflags

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

##################---------   Processing    ---------##################

def drifter_flagbad(ds_raw, fields_to_remove, lonname, latname,
                    good_flagval, ref_time,
                    bad_gps_flagval, na_latlim, na_lonlim,
                   uvelname, vvelname, val_threshold, bad_vel_flagvel):
    #--------------------------------------
    # Total number of points 
    total_points = int(ds_raw.time.size)

    # Compute a time vector in seconds
    dtime_sec = (pd.to_datetime(ds_raw.time.values) - 
                 ref_time).total_seconds()
    ds_raw['time_seconds'] = ("time", dtime_sec)
    

    # Find the modes for lat and lon
    mode1 = stats.mode(ds_raw[latname])
    lat_mode = mode1.mode
    mode_count = mode1.count
    mode1 = stats.mode(ds_raw[lonname])
    lon_mode = mode1.mode
    
    # Create a dictionary of attributes
    qc_attr_dict = {"total_points_orig": total_points}

    ########################################################################
    ######## Step 2. Create ds_qc the quality-controlled GPS fixes
    """

    ds_qc: Flag + remove bad GPS data
    - Where the GPS has dropped out (= mode),
    - Where it is outside the North Atlantic box (na_lonlim, na_latlim)
    - Where the velocity is unrealistic (larger than 3 m/s)
    - Where the time series of longitude and latitude (compared to a one-dimensional five-point median filter)
        is more than five standard deviations from the five-point median)

    Compute the u and v velocities on the QC GPS fixes
    """

    ds_qc = ds_raw


    ds_qc = ds_qc.drop(fields_to_remove, errors='ignore')

    # ----------------------------------------------------
    # Create the FLAG field, initialised to 'good'
    ds_qc["flag"] = ('time', good_flagval * np.ones(ds_qc.time.shape, 
                                                     dtype=np.int8))


    print(type(ds_qc))
    # ------------------------------------------------------
    # Flag bad gps (-90 lat, -180 lon)
    if mode_count > 100:
        # Kind of arbitrary but only do this if there are a lot of mode value
        ds_qc, numflags = flag_gps_dropout(ds_qc, latname, 
                                           lat_mode, bad_gps_flagval)

        ds_qc, numflags = flag_gps_dropout(ds_qc, lonname, 
                                           lon_mode, bad_gps_flagval)

        latlon_extremes(ds_qc, good_flagval)

    # Count data discarded + save to dictionary
    numflags1 = num_ibad(ds_qc, bad_gps_flagval)
    qc_attr_dict["GPS_dropouts"] = numflags1

    # ---------------------------
    # Flag gps locations outside the region (at Scripps)
    ds_qc, numflags = flag_gps_region(ds_qc, latname,
                                      na_latlim, bad_gps_flagval)
    ds_qc, numflags = flag_gps_region(ds_qc, lonname,
                                      na_lonlim, bad_gps_flagval)

    latlon_extremes(ds_qc, good_flagval)


    # Count data discarded + save to dictionary
    numflags2 = num_ibad(ds_qc, bad_gps_flagval) - numflags1
    qc_attr_dict["GPS_out_of_region"] = numflags2

    if (total_points - numflags1 - numflags2)==0:
        # Check whether there are any good data left
        print('No processed data')
        
        doneit = 0
        return doneit, ds_raw, qc_attr_dict
    else:

        # ----------------------------------------------------
        # Remove flagged data, otherwise it might affect the lowess filtering
        ds_qc = ds_qc.where(ds_qc.flag==good_flagval, drop=True) 

        # ----------------------------------------------------
        # ----------------------------------------------------
        # Calculate velocity - using jlab.latlon2uv
        GPSlat = ds_qc[latname].to_numpy()
        GPSlon = ds_qc[lonname].to_numpy()
        print('doing forward')
        u_orig, v_orig = jlab.latlon2uv_forward_mine(ds_qc['time_seconds'], GPSlat, GPSlon)

        # Add field to the ds_qc dataset
        ds_qc[uvelname] = ('time', u_orig)
        ds_qc[vvelname] = ('time', v_orig)

        # ---------------------------
        # Flag the too-high velocities (> 3 m/s)
        ds_qc, numflags3 = flag_vel(ds_qc, uvelname, 
                                    vvelname, val_threshold,
                                    bad_vel_flagvel)
        # Save to attribute dictionary
        qc_attr_dict["velocity_exceed_threshold"] = numflags3

        # Removing flagged data
        ds_qc = ds_qc.where(ds_qc.flag==good_flagval, drop=True) 

        print('Flagged '+str(numflags1)+' for GPS dropouts, '+str(numflags2)+' for region, '+str(numflags3)+' for velo')
    
        doneit = 1
        return doneit, ds_qc, qc_attr_dict
        
def drifter_filter(ds_qc, lonname, latname,
                   poly_order, bandwidth):

    ######## Step 3. Run Lowess filter
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Create entries for lowess lat/lon
    time1 = ds_qc.time.values
    ds_lowess = xr.Dataset(coords=dict(time=(["time"], time1)))
    ds_lowess = ds_lowess.assign_attrs(ds_qc.attrs)

    ds_lowess[latname] = ('time', np.ones(len(ds_qc.time)))
    ds_lowess[lonname] = ('time', np.ones(len(ds_qc.time)))


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # LOWESS: Locally Weighted Scatterplot Smoother (Elipot et al 2016)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    time_sec = ds_qc['time_seconds'].values
    lon = ds_qc[lonname].values
    lat = ds_qc[latname].values
    llat, llon = LatLonLocalWess(time_sec, lon, lat, poly_order, 
                                 bandwidth)

    ds_lowess[latname] = ('time', llat)
    ds_lowess[lonname] = ('time', llon)

    # Cannot use the same code to Lowess filter the velocities since it
    # treats latitude and longitude differently.
    return ds_lowess
    
    
def drifter_hourly(ds_lowess, ds_qc, ds_raw, tgrid_hourly, tgrid_sec,
                  integ_vars, float_vars, lonname, latname, uvelname, vvelname,
                  qc_attr_dict):

    ######## Step 4. Interpolate hourly
     # Revised from Oana's method (using arrays) to use dataset in xarray.


    # Empty lists for storing the hourly data and new column names
    hourly_arr = []
    arr_names = []


    total_points_qc = int(ds_raw.time.size)
    qc_attr_dict["qc_points"] = total_points_qc


    # Create the new hourly xarray dataset
    ds_hourly = xr.Dataset(coords=dict(time=(["time"], tgrid_hourly)))
    ds_hourly = ds_hourly.assign_attrs(ds_lowess.attrs)
    ds_hourly
    time_sec = ds_qc.time_seconds
    time_sec_raw = ds_raw.time_seconds

    # Get the variable names
    varnames_list = list(ds_raw.keys())
    varnames_qc_list = list(ds_qc.keys())
    varnames_lowess_list = list(ds_lowess.keys())

    # Create temporary array to store hourly interp data
    c_len = len(tgrid_hourly)
    arr = np.ones((c_len,))

    # Iterate through the variable names and set the interpolation method 
    # based on the dtype of the data
    k = 0
    for varname in set(varnames_list + varnames_qc_list):
        if (varname in varnames_lowess_list):
            # Lowess filtered variables (GPS lat and lon)
            if varname in integ_vars:
                interp_method = 'nearest'
            elif varname in float_vars:
                interp_method = 'linear'

            var = ds_lowess[varname].values
            fc_interp = interp1d(time_sec, var, interp_method,
                                bounds_error=False)
            # Need to fill the lat/lon with NaN
            arr = fc_interp(tgrid_sec)

            ds_hourly[varname] = ('time', arr)
            #            print('Using lowess for '+varname)

        elif (varname in varnames_qc_list):
            # QC drifter variables (velocity)
            if varname in integ_vars:
                interp_method = 'nearest'
            elif varname in float_vars:
                interp_method = 'linear'

            var = ds_qc[varname].values
            fc_interp = interp1d(time_sec, var, interp_method,
                                 bounds_error=False)
            arr = fc_interp(tgrid_sec)

            ds_hourly[varname] = ('time', arr)
             #           print('Using qc for '+varname)


        elif (varname in integ_vars) or (varname in float_vars):
            # Raw drifter variables
            if varname in integ_vars:
                interp_method = 'nearest'
            elif varname in float_vars:
                interp_method = 'linear'

            var = ds_raw[varname].values
            fc_interp = interp1d(time_sec_raw, var, interp_method,
                                 bounds_error=False)
#                                 fill_value='extrapolate')
            arr = fc_interp(tgrid_sec)

            ds_hourly[varname] = ('time', arr)
        else:
            print("Variables not interpolated hourly: < %s >\n" % varname)

    # Chop to first non-NaN latitude and longitude values
    lat1 = ds_hourly[latname]
    ibad = np.isnan(lat1)
    ds_hourly = ds_hourly.drop_isel(time=ibad)
    lon1 = ds_hourly[lonname]
    ibad = np.isnan(lon1)
    ds_hourly = ds_hourly.drop_isel(time=ibad)


    # Use attributes to store units and names for DataArrays
    ds_hourly[lonname].attrs["long_name"] = "longitude_lowess"
    ds_hourly[latname].attrs["long_name"] = "latitude_lowess"
    ds_hourly[uvelname].attrs["units"] = "m/s"
    ds_hourly[vvelname].attrs["units"] = "m/s"
    ds_hourly[uvelname].attrs["long_name"] = "zonal_velocity"
    ds_hourly[vvelname].attrs["long_name"] = "meridional_velocity"
    
    return ds_hourly

