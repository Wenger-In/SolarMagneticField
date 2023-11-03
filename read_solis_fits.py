from astropy.io import fits
import numpy as np
import os

solis_dir = 'E:/Research/Data/SOLIS/'
target_name = '_000_int-mas_dim-900.fits'

# 遍历文件夹
for root, dirs, files in os.walk(solis_dir):
    for file in files:
        if file.endswith(target_name):
            file_path = os.path.join(root, file)
            fits_file = fits.open(file_path)
            fits_file.info()
            primary_hdu = fits_file[0] 
            data = primary_hdu.data
            fits_file.close()
            save_name = 'cr' + str(file[17:21]) + '.fits'
            fits.writeto(solis_dir + save_name, data, overwrite=True)
