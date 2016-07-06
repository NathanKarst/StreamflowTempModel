# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 13:38:05 2016

Script to unzip any PRISM files within 'raw_data' directory. 
All data should be 'stable', daily PRISM data. Do not use
provisional data. 

@author: daviddralle
"""

import zipfile
import os
from os.path import dirname


def main():
    parent_dir = dirname(dirname(os.getcwd()))
    raw_data_dir = os.path.join(parent_dir,'raw_data')
    prism_vars = ['ppt','tmean','tmin','tmax']
    
    for prism_var in prism_vars:
        years = os.listdir(os.path.join(raw_data_dir,prism_var))
        if len(years)==1:
            print 'Warning: Could not find any folders containing ' + prism_var + ' data.'
            continue
            
        for year in years[1:]:
            print 'unzipping ' + year + ' ' + prism_var + ' data.'
            zips =os.listdir(os.path.join(raw_data_dir,prism_var,year))
            for zip_file in zips:
                    if 'stable' in zip_file:
                        zip_ref = zipfile.ZipFile(os.path.join(raw_data_dir,prism_var,year,zip_file), 'r')
                        #os.mkdir(zip_file[:-4])
                        zip_ref.extractall(os.path.join(raw_data_dir,prism_var,year))
                        zip_ref.close()
                        os.remove(os.path.join(raw_data_dir,prism_var,year,zip_file))
                    else:
                        os.remove(os.path.join(raw_data_dir,prism_var,year,zip_file))
        

if __name__ == '__main__': 
    main()

