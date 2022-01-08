"""
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 DIRECTORY TREE and PATHS
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 terific_nocs/
		- drifter_aux_data/			: auxiliary datasets
		- drifter_out/				: scripts output saved here
		- drifter_scripts/			: main scripts
			-- SIO_wrap/			: helper scripts
 
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Last edited: 6 Sep 2021
Edited os.makedir to os.makedirs (maybe associated with python change 3.9.4 to 3.9.5?)
"""

# Import modules
import os

####################-----------   USER EDITS    ------------####################

# Specify the path where to create the directory structure
path = "/Users/eddifying/Python/"

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# name of main directory (joined with the path)
maindir = os.path.join(path, 'drifters/')

# subdirectories (with path)
dir_scr = os.path.join(maindir, '02-code/')
dir_sio = os.path.join(dir_scr, 'SIO_wrap/')
dir_data = os.path.join(maindir,'01-data/')
dir_out = os.path.join(dir_data, '02-intermediate/')
dir_aux = os.path.join(dir_data, '04-aux/')

dir_list = (dir_scr, dir_sio, dir_aux, dir_out)


# - - - check/create DIRECTORY TREE - - - - 

def create_dir(path_dir):
    """
    Check if a directory 'path_dir' exists
    and if not, create it.
    """
    if os.path.isdir(path_dir):
    	print("%s already exists" % path_dir)
    else:
    	os.makedirs(path_dir)
    	print("%s created" % path_dir)


if os.path.isdir(maindir):
	print("%s already exists" % maindir)
	for dir_var in dir_list:
		create_dir(dir_var)
else:
	os.makedirs(maindir)
	print("%s created" % maindir)

	for dir_var in dir_list:
		os.makedirs(dir_var)
		print("%s created" % dir_var)
