'''
Script to extract REW hydrologic forcing data from PRISM precip, mean temp, and tdmean data
for use in the Priestley-Taylor ET model. Data is generated for use vadose zone models 
of the type developed in Laio et al. 2009.

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
import pickle
from functools import partial

parent_dir = dirname(dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','lib'))
import zonal_stats as zs
import meteolib as meteo
import evaplib as evap



def main():
    """ Extract hydrologic forcing data. 
    
    Script to extract REW hydrologic forcing data from PRISM precip, mean temp, and tdmean data
    for use in the Priestley-Taylor ET model. Data is generated for vadose zone models
    of the type developed in Laio et al. 2009.

    The lat/long of each REW centroid is used to fetch the value at each 
    generatedridpoint of the climate forcing surface. 
    
    Dataframe columns include:
        - 
    
    Args: 
        - None
        
    Returns:
        - None
    """
    prism_vars = ['ppt','tmean', 'tdmean']
    
    #Get basins shapefile. 
    try:
        basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    except RuntimeError:
        print 'Cannot find basins shapefile. Please make sure basins shapefile is located in \n the model directory under /raw_data/basins_poly'


    #Get list of lat/long centroid tuples for each REW
    fc = fiona.open(basins)
    shapefile_record = fc.next()
    pos_dict={}
    for shapefile_record in fc:
        shape = shapely.geometry.asShape(shapefile_record['geometry'])
        long_point = shape.centroid.coords.xy[0][0]
        lat_point = shape.centroid.coords.xy[1][0]
        pos = (long_point, lat_point)
        pos_dict[int(shapefile_record['properties']['cat'])]=pos
        
    
    #initialize hdf5 forcing database in model_data folder
    #forcing_store = pd.HDFStore(os.path.join(parent_dir,'model_data','forcing.h5'))

    
    #dictionary that will hold list of rew forcing lists for each forcing variable
    forcing_dict = {}
    dates_dict = {}

    for prism_var in prism_vars:
        years = os.listdir(os.path.join(parent_dir,'raw_data',prism_var))
        
        #key of vals_dict is rew id, corresponding list 
        #is forcing data in dateorder of date_list
        vals_dict = {k:[] for k in pos_dict.keys()}
        date_list=[]
        #for all years of forcing variable, load each day of raster data and extract to all REWs
        for year in years:
            try: 
                int(year)
            except ValueError:
                continue

            print '...processing ' + str(year) + ' ' + prism_var + ' data...'

            #Data should be in .bil format from PRISM
            #Data for each REW is extracted to REW centroid lat/long, 
            #Not zonally averaged over the REW area. Zonal averaging can be used, 
            #but it is computationally intensive 
            raster_list = glob.glob(os.path.join(parent_dir,'raw_data',prism_var,year,'*.bil'))
            for rast in raster_list:
                date_list.append(rast[-16:][:8])
                raster_file = os.path.join(parent_dir,'raw_data',prism_var,year,rast)
                gdata = gdal.Open(raster_file)
                gt = gdata.GetGeoTransform()
                data = gdata.ReadAsArray().astype(np.float)
                gdata = None
                for rew_id in pos_dict.keys(): 
                    pos = pos_dict[rew_id]
                    x = int((pos[0] - gt[0])/gt[1])
                    y = int((pos[1] - gt[3])/gt[5])
                    vals_dict[rew_id].append(data[y, x])
        
        #store forcing variable data for each rew 
        #and corresponding list of daily dates          
        forcing_dict[prism_var]=vals_dict
        dates_dict[prism_var]=date_list
  
    
    #check that for each forcing variable, the dates are the same
    unique_dates = [list(i) for i in set(tuple(i) for i in dates_dict.values())]
    if len(unique_dates)>1:
        print 'Forcing data dates do not all match, please check raw data'
        return
  
    #initialize a list of pandas forcing dataframes for each REW
    rng = pd.date_range(unique_dates[0][0],unique_dates[0][-1])
    rew_dfs_dict = {k:pd.DataFrame(data=None,index=rng) for k in pos_dict.keys()}
    
    #for each date in rew_dfs_dict, fill radiative forcing data
    #currently, this accessses an rn data frame generated by ipython notebook 
    rew_dfs_dict = get_rad(rew_dfs_dict, parent_dir)

    #for each forcing variable, add the forcing data to each REW dataframe
    for prism_var in prism_vars:
        #need to convert ppt from mm/day to cm/day for PRISM data
        if prism_var=='ppt':
            div = 10.0
        else:
            div = 1.0
        for rew_id in pos_dict.keys():
            rew_dfs_dict[rew_id][prism_var] = np.array(forcing_dict[prism_var][rew_id])/div

    #for each rew, compute ETmax
    #...to do this, first compute pressure...
    rew_elevations = get_rew_elevations(basins, parent_dir)
    rew_pressures = get_rew_pressures(rew_elevations)
    #...then compute relative humidity and use priestley taylor for PET
    for rew_id in pos_dict.keys():
        rh = get_rh(rew_dfs_dict[rew_id]['tmean'], rew_dfs_dict[rew_id]['tdmean'])

        #priestley taylor equation for PET
        #needs mean temp, relative humidity, atmospheric pressure, net radiation, 
        #and soil heat flux (assumed zero at daily timescale)
        pet = evap.Ept(
            rew_dfs_dict[rew_id]['tmean'],
            rh, 
            rew_pressures[rew_id]*np.ones(len(rh)), 
            rew_dfs_dict[rew_id]['rn'],
            np.zeros(len(rh))
            )
        #convert to cm/day
        rew_dfs_dict[rew_id]['pet'] = pet/10.0

            
    pickle.dump( rew_dfs_dict, open( os.path.join(parent_dir,'model_data','rew_forcing_dict.p'), "wb" ) )
        

def get_rh(tmean, tdmean):
    a = 17.625
    b = 243.04
    return 100*np.exp(a*tdmean/(b+tdmean))/np.exp(a*tmean/(b+tmean))

def get_rad(rew_dfs_dict, parent_dir):
    #for now, uses data frame in /raw_data/rn folder
    #assumes that rn already has AT LEAST date range of PRISM data
    df = pickle.load( open( os.path.join(parent_dir, 'raw_data','rn', 'rn.p'), "rb" ) )
    rew_id = rew_dfs_dict.keys()[0]
    start = rew_dfs_dict[rew_id].index[0]
    stop = rew_dfs_dict[rew_id].index[-1]
    for rew_id in rew_dfs_dict.keys():
        rew_dfs_dict[rew_id]['rn'] = df[start:stop]

    return rew_dfs_dict
          
def get_rew_elevations(basins, parent_dir):
    #get mean rew elevation
    dem_file = os.listdir(os.path.join(parent_dir,'raw_data','dem'))
    dem_file = [x for x in dem_file if 'dem' in x]
    elev_stats = zs.zonal_stats(basins, os.path.join(parent_dir,'raw_data','dem',dem_file[0]))
    
    rew_elevations = {}
    translated = translate_to_rew_id(parent_dir)

    for stat in elev_stats: 
        rew_elevations[translated[stat['fid']]] = stat['mean']

    return rew_elevations

def get_rew_pressures(rew_elevations):
    #elevations in meters, output pressure in Pa
    #using simple elevation/pressure relationship
    rew_pressures = {}
    for key in rew_elevations.keys():
        rew_pressures[key] = 1000*101.325*((293-0.0065*rew_elevations[key])/293)**5.26

    return rew_pressures


def translate_to_fid(parent_dir):
    #function that translates feature id's to rew id's
    #returns a dictionary where the key is fid and value is rew id
    basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    translated = {}

    fc = fiona.open(basins)
    shapefile_record = fc.next()
    for shapefile_record in fc:
        translated[int(shapefile_record['properties']['cat'])] = int(shapefile_record['id'])
    return translated

def translate_to_rew_id(parent_dir):
    #function that translates rew id's to feature id's
    basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    translated = {}

    fc = fiona.open(basins)
    shapefile_record = fc.next()
    for shapefile_record in fc:
        translated[int(shapefile_record['id'])] = int(shapefile_record['properties']['cat'])
    return translated

if __name__ == "__main__":
    main()

    # BELOW CAN BE USED FOR QUERYING RN RASTER USING REWS
    # rad_files = os.listdir(os.path.join(parent_dir,'raw_data','rn'))
    #     rad_files = [x for x in rad_files if '.' not in x]
    #     #store radiation data in numpy array of size num_REW by len(rng)
    #     rad_dict = {k: [] for k in pos_dict.keys()}
    #     counter=1
    #     for rad_file in rad_files:
    #         print 'processing radiation data, DOY: ', counter
    #         stats = zs.zonal_stats(basins, os.path.join(parent_dir,'raw_data','rn',rad_file))
    #         for stat in stats:
    #             rad_dict[stat['fid']].append(stat['mean'])

    #         counter=counter+1
            
            





