"""
SCRIPT 1
Download raw data from website and store a netcdf file.

<!> dependency <!> 
Dir SIO_wrap/ contains files:
	> dir_tree.py 	: setting up directory tree
	> fnames.py 	: file name convention

TERIFIC Data are available since 2019-12-04.

Scrape data from web and stitch them together with existing data.
If there is existing data, the script looks for the last time sample from the 
previous download and updates dataset from there.

Filename timestamp: latest day that is fully sampled [from 00h to 23h] 

Last modified: 5 Sep 2021
"""

####################---------   LOCAL FUNCTIONS    ---------####################
def mysplit(s, delim=None):
	"""
	Parameters
	----------
	s 		: string 
	delim 	: delimiter is a string, default separator is None 

	The built-in split() method splits a string into a list but does not 
	ignore empty strings and when applying it for the TERIFIC GDP data 
	it was creating an additional empty column.  This functions removes the
	empty strings. [last checked: Aug 2021]

	Returns
	-------
	List of strings where the elements are the substrings separated 
	by the specifed delimiter.
	"""
	return [x for x in s.split(delim) if x]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import modules

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import requests
import os
import glob
import re
import sys
from datetime import timedelta

# Local import 
# > Make sure SIO_wrap dir is on the same path as this script.
from SIO_wrap import dir_tree, fnames

# Contents of namespaces
#dir(dir_tree)
#dir(fnames)

####################-----------   USER EDITS    ------------####################

# Path where data are saved. Can be changed in file SIO_wrap/dir_tree.py
data_dir = dir_tree.dir_out

# SIO username and password
username = "uk-noc"
password = "noc-drifter"

# Download URL (main body without the start date)
download_url = ("https://gdp.ucsd.edu/cgi-bin/projects/uk-noc/"
				"drifter.py?start_date=") 

# Download start date must have format yyyy-mm-dd. Default is set to the
# beginning of the TERIFIC project, i.e 2019-12-04.
download_start_date = "2019-12-04"
print("\nDefault download start date: %s\n" % download_start_date)

# String formatting for time for:
#   - the download url, 
#   - appending to the filename
# 	- the data time column, respectively.
url_strftime = '%Y-%m-%d'
tstamp_strftime = '%Y%m%d'
timcol_strftime = '%Y-%m-%d %H:%M:%S'


################################################################################
################################################################################

# Extract a list with the names of existing raw data files.
existing_files = glob.glob(os.path.join(data_dir, fnames.fname_rawdata + '*'))

# ~ ~ print update ~ ~ 
if len(existing_files) > 0:
	print("Existing raw data files: \n%s" % existing_files)
else:
	print("No previous raw data files.\n")

# ~ ~ filenaming convention ~ ~
# If there are multiple files with raw data (i.e. non-updated datasets), select 
# the latest one updated.
# The file names are distinguished by the timestamp appended to the filename 
# and has <tstamp_strftime> format (see 'user edits' section).
# The data are cropped such that the last day is fully sampled (spans 0h-23h).
# The timestamp in the filename is the latest downloaded fully sampled day.

if len(existing_files) > 0:

	# Extract the timestamp part of the filename(s) in a list
	tstamp = [date for file in existing_files 
				for date in re.findall("(\d{8})", file)]

	# Convert to datetime and pick the most recent timestamp
	tstamp_date = pd.to_datetime(tstamp, format=tstamp_strftime)
	prev_fname_timestamp = tstamp[tstamp_date.argmax()]

	# Load the previously updated file
	prev_fname = f"{fnames.fname_rawdata}{prev_fname_timestamp}.nc"
	prev_fpath = os.path.join(data_dir, prev_fname)
	prev_ds = xr.open_dataset(prev_fpath)

	# Make sure the time variable is sorted in ascending order
	#prev_time = prev_ds.time.sortby(prev_ds.time)

	# Set download start date to +1 day from the last day of previous dataset
	latest_date = tstamp_date.max()
	download_start_date = (latest_date 
							+ timedelta(days=1)).strftime(url_strftime)

	print("Download start date changed to: %s\n" % download_start_date)

# Combine body of download link with start date
data_url = download_url + download_start_date

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Scrape data from website 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

print("Scraping data from the website starting from: %s \n.......\n" 
		% download_start_date)

data_html = requests.get(data_url, auth=(username, password))
# If print(data_html) outputs <Response [200]> then code worked
# To print content: data_html.content


# Extract text from the html page
print("Parsing web data ...")

data_soup = BeautifulSoup(data_html.text, "html.parser") 
data_text = data_soup.text


# Split the text after every newline character '\n' into separate rows
print("Splitting text into separate rows ...")

data_rows = data_text.splitlines()

# Further split the data into columns, but first:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Edit the header
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#
# The first row contains the header.  This is repeated for every drifter but
# because we can ID the drifters by the Platform ID, we remove the header lines.
# We use the header to name the columns of data, but first we remove unwanted
# characters (spaces/parantheses/dashes). 

header_raw = data_rows[0]

print("\nHeader format before processing:\n%s\n" % header_raw)

# Remove header lines
data_rows_clean = [x for x in data_rows if header_raw not in x]

# Remove unwanted characters from the header. 
# [!!!] These might change if there are new columns added/names change. 
header = header_raw.replace(" ", "")	
header = header.replace("-",  "_")
header = header.replace("(", "_")
header = header.replace(")", "")

# Split the header into columns
col_names = mysplit(header, ',')

print("\nHeader after removing unwanted spaces and characters:\n%s\n" % col_names)

# check header matches the lists of integer/float names
#if all(item in col_names for item in integ_vars)==False:
#	print("List of integer var names does not match the column names")
#	print("Check variable %s" % str(integ_vars))
#if all(item in col_names for item in float_vars)==False:
#	print("List of float var names does not match the column names")
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Store data in a pandas dataframe 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# > Split each row of data into columns; delimiter: comma and space (', ')

print("Create a pandas dataframe....")
df = []		# create a list

for i in range(len(data_rows_clean)):
	df.append(mysplit(data_rows_clean[i], ', ')) 

# > Assign each column a name using the edited header
data_df = pd.DataFrame(df, columns=col_names)

#print("Data stored in a pandas dataframe. Data fields: \n%s" % data_df.keys())

# > Change the formatting of the time column from text to datetime[64]
# extract name of column that contains time
#time_colname = data_df.filter(like=('Time' or 'time')).columns 
time_colname = 'Timestamp_UTC'
data_df[time_colname] = pd.to_datetime(data_df[time_colname],
											format=timcol_strftime)
#print(data_df)

# > Sort rows by time
print("Sorting rows by time ..\n")
data_df = data_df.sort_values(by=time_colname)

#data_df = data_df.iloc[:1200000]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# convert pandas dataset to xarray dataset 
# (easier to save as netcdf file)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define a dictionary and first populate it with the float variables.
# Treat the time variable separately because it has type datetime64[ns].
#
# Dictionary uses an ordinate (n) as a coordinate; decided not to use the time 
# because it does not have unique values although this can be changed.
print("Converting DataFrame to xarray Dataset..\n")
dd = {}		# create an empty dict

for coln in col_names:
	if coln != time_colname:
		var = (pd.to_numeric(data_df[coln]))
		#print(var)
		dd[coln] = ("n", data_df[coln].astype(var.dtype).values)


# changed the time variable name to 'time'
dd["time"] = ("n", data_df[time_colname].values)

# xarray dataset
ds = xr.Dataset(dd,
	coords={"n" : np.arange(len(data_df[time_colname].values))})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# crop data so that the last day is fully sampled and there are no overlaps
# when the data are updated; basically discard the last day if it's incomplete
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

end_datetime = pd.to_datetime(ds.time.values[-1])
end_datestr = end_datetime.strftime(url_strftime)

penultimate_datetime = end_datetime - timedelta(days=1)
penultimate_datestr = penultimate_datetime.strftime(url_strftime)

if download_start_date == end_datestr:
	sys.exit("No updated data. Last full day available is %s" 
		% penultimate_datestr)


cutoff_date = pd.to_datetime(end_datestr +" 00:00:00", format=timcol_strftime)

ds_crop = ds.where(ds.time<cutoff_date, drop=True)

# timestamp for filename
fname_timestamp = penultimate_datetime.strftime(tstamp_strftime)

# stitch together the files
if len(existing_files) > 0:
	print("Stitch updated dataset with the previous one. \n")
	# use previously opened dataset (prev_ds)
	# put both datasets in a list
	d = []
	d.append(prev_ds)
	d.append(ds_crop)

	# merge list into a dataset
	new_ds = xr.concat(d, dim='n')

else:
	new_ds = ds_crop


# Filename and path of (updated) dataset
update_fname = f"{fnames.fname_rawdata}{fname_timestamp}.nc"
update_fpath = os.path.join(data_dir, update_fname)

# Save dataset to netcdf file
new_ds.to_netcdf(update_fpath)

print("File updated/saved: %s" % update_fpath)
print("Script 1 finished. \n")
