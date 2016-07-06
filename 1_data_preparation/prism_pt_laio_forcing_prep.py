# -*- coding: utf-8 -*-
'''
Script to extract REW hydrologic forcing data from PRISM precip and mean temp datasets 
for use in the Priestley-Taylor ET model. Data is generated for use in the
vadose zone model developed in Laio et al. 2009.

Output is an hdf5 database of forcing data
stored in the model_data folder of the model parent directory. Data is 
indexed by rew_x, where x is the REW id. 

The lat/long of each REW centroid is used to fetch the value at each 
gridpoint of the climate forcing surface. 

'''


import pandas as pd
import gdal
import os
import numpy as np
import fiona
import shapely
from shapely import geometry
from os.path import dirname
import glob
import sys



def main(): 
    #extract these prism variables to all REWs for all times
    prism_vars = ['ppt','tmean']
    parent_dir = dirname(dirname(os.getcwd()))
    sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','lib'))
    import zonal_stats as zs
    
    try:
        basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    except RuntimeError:
        print 'Cannot find basins shapefile. Please make sure basins shapefile is located in \n the model directory under /raw_data/basins_poly'

    #Get list of lat/long centroid tuples for each REW
    #index in pos_list corresponds to REW ID number
    fc = fiona.open(basins)
    pos_list = []
    shapefile_record = fc.next()
    
    #dictionary for converting IDs to labels in topology
    translate_labels ={}
    for shapefile_record in fc:
        geo = shapefile_record['geometry']
        shape = shapely.geometry.asShape(shapefile_record['geometry'])
        long_point = shape.centroid.coords.xy[0][0]
        lat_point = shape.centroid.coords.xy[1][0]
        pos = (long_point, lat_point)
        pos_list.append(pos)
        translate_labels[shapefile_record['id']] = shapefile_record['properties']['value']
    
    #initialize hdf5 forcing database in model_data folder
    forcing_store = pd.HDFStore(os.path.join(parent_dir,'model_data','forcing.h5'))
    num_REW = len(pos_list)
    
    #dictionary that will hold list of rew forcing lists for each forcing variable
    forcing_dict = {}
    dates_list = []

    for prism_var in prism_vars:
        years = os.listdir(os.path.join(parent_dir,'raw_data',prism_var))
        if len(years)==1:
            'Error: No data exists for ' + prism_var
            return
            
        vals_list = [[] for i in range(num_REW)]
        date_list=[]
        #for all years of forcing variable, load each day of raster data and extract to all REWs
        for year in years[1:]:
            print '...processing ' + str(year) + ' ' + prism_var + ' data...'
            raster_list = glob.glob(os.path.join(parent_dir,'raw_data',prism_var,year,'*.bil'))
            for rast in raster_list:
                date_list.append(rast[-16:][:8])
                raster_file = os.path.join(parent_dir,'raw_data',prism_var,year,rast)
                gdata = gdal.Open(raster_file)
                gt = gdata.GetGeoTransform()
                data = gdata.ReadAsArray().astype(np.float)
                gdata = None
                for rew_id in range(num_REW): 
                    pos = pos_list[rew_id]
                    x = int((pos[0] - gt[0])/gt[1])
                    y = int((pos[1] - gt[3])/gt[5])
                    vals_list[rew_id].append(data[y, x])
                    
        forcing_dict[prism_var]=vals_list
        dates_list.append(date_list)
  
    
    #check that for each forcing variable, the dates are the same
    unique_dates = [list(i) for i in set(tuple(i) for i in dates_list)]
    if len(unique_dates)>1:
        print 'Forcing data dates do not all match, please check raw data'
        return
  
    #initialize a list of pandas forcing dataframes for each REW
    rng = pd.date_range(unique_dates[0][0],unique_dates[0][-1])
    rew_dfs = [ pd.DataFrame(data=None,index=rng) for i in range(num_REW) ]
    
    #get array of radiative forcing data for each REW
    #assumes rn folder is filled with 366 daily integrated net rad 
    #forcing data generated from DEM. For each day/REW, the forcing is
    #computed as the mean forcing within the REW polygon. 
    rad_files = os.listdir(os.path.join(parent_dir,'raw_data','rn'))
    rad_files = [x for x in rad_files if '.' not in x]
    #store radiation data in numpy array of size num_REW by len(rng)
    rad_array = np.empty([num_REW,366])
    day=1
    for rad_file in rad_files:
        stats = zs.zonal_stats(basins, os.path.join(parent_dir,'raw_data','rn',rad_file))
        rad_array[:,day-1]=np.array([x['mean'] for x in stats])
        day+=1
    
    #for each forcing variable, add the forcing data to each REW dataframe
    for prism_var in prism_vars:
        for rew_id in range(num_REW):
            rew_dfs[rew_id][prism_var] = forcing_dict[prism_var][rew_id]
    
    #for each REW, plug each year of radiative data into dataframe
    for rew_id in range(num_REW):
        df[rew_id]['rn']=np.nan
        for year in years:
            df[rew_id].rn[year]=rad_array[rew_id,0:(len(df[year])-1)]
            
            
    for rew_id in range(num_REW):
        forcing_store['rew_' + str(rew_id)] = rew_dfs[rew_id]
        

            
            
if __name__ == '__main__': 
    main()




