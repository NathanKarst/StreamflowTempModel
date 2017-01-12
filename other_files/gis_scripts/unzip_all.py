# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 13:38:05 2016

Script to unzip all PRISM precip files within 'local' directory
@author: daviddralle
"""

from zonal_stats import zonal_stats as zs
import zipfile
import gdal, rasterio, fiona
import os
home = os.environ['HOME']
local = '/Dropbox/research/streamflow_temp/data/prism_monthly'
local_data = home + local + '/ppt'

#years = os.listdir(local_data)
years = ['2015','2016']
for year in years:
        try: 
            test = float(year)
        except ValueError:
            continue
        
        zips =os.listdir(local_data + '/' + year)
        for zip_file in zips:
                if 'stable' in zip_file:
                    zip_ref = zipfile.ZipFile(local_data + '/' + year + '/' + zip_file, 'r')
                    #os.mkdir(zip_file[:-4])
                    zip_ref.extractall(local_data + '/' + year)
                    zip_ref.close()
                    os.remove(local_data + '/' + year + '/' + zip_file)
                else:
                    os.remove(local_data + '/' + year + '/' + zip_file)
        


