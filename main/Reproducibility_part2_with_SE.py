"""
Created on Mon Nov 27 09:22:13 2023

@author: jaebeom
"""

import numpy as np
from astropy.io import fits
import os

def make_param_file(parameter_set):
    fm = open('%s.params' %(cname), 'w')
    for i in range(len(parameter_set)):
        fm.write('%s\n' %(parameter_set[i]))
    fm.close()
    return 0
    
def run_se(data, header, config_file, param_file):
    header['NAXIS'] == 2   
    os.system('mkdir catalogue')
    for i in range(len(data)):
        fits.writeto('tmp_%s.fits' %(i), data[i], header, overwrite=True)
        os.system('sex tmp_%s.fits -c %s -CATALOG_NAME cat_%s.dat -PARAMETERS_NAME %s -FILTER %s \
                  -CLEAN %s -PHOT_APERTURES %s -BACK_SIZE %s -BACK_FILTERSIZE %s\
                  -BACKPHOTO_TYPE %s -CHECKIMAGE_TYPE %s -CHECKIMAGE_NAME obj_only_%s.fits' \
                  %(i, config_file, i, param_file, FILTER, CLEAN, PHOT_APERTURES, BACK_SIZE, BACK_FILTERSIZE, \
                    BACKPHOTO_TYPE, CHECKIMAGE_TYPE,i))
        os.system('rm -rf tmp*.fits')
        os.system('mv cat_*.dat catalogue')
    return 0
    
def make_obj_cube(original_data, original_header):
    data = original_data
    header = original_header
    source_cube = np.zeros(data.shape[0]*data.shape[1]*data.shape[2]).\
        reshape(data.shape[0],data.shape[1],data.shape[2])
    for i in range(len(data)):
        f1 = fits.open('obj_only_%s.fits' %(i))
        source_cube[i] = np.array(f1[0].data, dtype='float32')
        print ('obj_%s complete...' %(i))
    source_cube = np.array(source_cube, dtype='float32')
    fits.writeto('source_only_cube_BS%s_BF%s.fits' %(BACK_SIZE, BACK_FILTERSIZE), source_cube, header, overwrite=True)
    os.system('rm -rf obj_*.fits')
    return source_cube

def obj_sub_cube(original_cube, source_only_cube):
    new_cube = original_cube - source_only_cube
    fits.writeto('source_sub_cube_BS%s_BF%s.fits' %(BACK_SIZE, BACK_FILTERSIZE), new_cube, header)
    return new_cube
  
    
#================= Path setting
datapath = '/home/jaebeom/desktop/sdc/data/ZW3.msw_image_pbcor_900crop.fits'
workpath = '/home/jaebeom/desktop/sdc/se_share'
#==============================================================================

#================= Source extractor parameter setting part
cname = 'default'       # Configuration and parameter file name
FILTER = 'N'
CLEAN = 'N'
PHOT_APERTURES = 20
BACK_SIZE = 20
BACK_FILTERSIZE = 1
BACKPHOTO_TYPE = 'LOCAL' #set GLOBAL or LOCAL
CHECKIMAGE_TYPE = 'OBJECTS'
#==============================================================================
  

f = fits.open(datapath)
data = f[0].data
header = f[0].header
os.chdir(workpath)
os.system('sex -d > %s.conf' %(cname))
param_set = ['NUMBER', 'MAG_APER', 'MAGERR_APER', 'X_IMAGE', 'Y_IMAGE', ]

config_file = '%s.conf' %(cname)
param_file = '%s.params' %(cname)

make_param_file(param_set)
run_se(data, header, config_file, param_file)
source_cube = make_obj_cube(data, header)
obj_sub_cube(data, source_cube)
