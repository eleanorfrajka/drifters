import os
import datetime
import shutil
import requests

## Left off here
# https://stackoverflow.com/questions/56950987/download-file-from-url-and-save-it-in-a-folder-python/56951038
def download_file(url, folder_name):
    local_filename = url.split('/')[-1]
    path = os.path.join("/{}/{}".format(folder_name, local_filename))
    with requests.get(url, stream=True) as r:
        with open(path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

    return local_filename
##

def create_figdir():
    #figdir = '../03-results/figures/yyyy-mm/';
    yyyymm = datetime.datetime.now().strftime('%Y-%m')
    figdir = os.path.join('..','03-results','figures',yyyymm)
    if not os.path.isdir(figdir):
        os.makedirs(figdir, exist_ok=True)
        #os.mkdir(figdir)
    return figdir
        
def save_figure(fig,figname):
    figdir = create_figdir()
    figname_with_ext = 'fig_'+figname+'.png'
    fig.savefig(os.path.join(figdir,figname_with_ext))

def cat_data_path(data_folder,filename):
    input_file = os.path.join('..',data_folder,filename)
    return input_file

def cat_raw_path(fname):
    output_file_with_path = os.path.join('..','01-data','01-raw',fname)
    return output_file_with_path

def cat_interim_path(fname):
    output_file_with_path = os.path.join('..','01-data','02-intermediate',fname)
    return output_file_with_path

def cat_proc_path(fname):
    output_file_with_path = os.path.join('..','01-data','03-processed',fname)
    return output_file_with_path