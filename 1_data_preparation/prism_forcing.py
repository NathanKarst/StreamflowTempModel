
# coding: utf-8

# # Forcing Extraction Tutorial
# 
# In this notebook, we'll present an example of converting PRISM precipitation, mean temperature, and mean dew point data into a format that can be used easily in the types of vadose and groundwater zone models we cover in the [custom zone models tutorial](custom_zone_models.ipynb). 

# In[1]:

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


# We first need to specify the geographic areas of interest. We'll assume that these regions are specified via shapefiles located in the `raw_data` folder. 

# In[2]:

def main():

    try:
        basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    except RuntimeError:
        print 'Cannot find basins shapefile. Please make sure basins shapefile is located in \n the model directory under /raw_data/basins_poly'


    # For simplicity, we'll assume that the forcing along the spatial extent of each basin is constant and equal to the forcing observed at its centroid. We could also consider implementing zonal averaging or some other more faithful approximation, but this would be much more computationally intensive, and quite likely unnecessary for lower resolution climate datasets. 

    # In[3]:

    fc = fiona.open(basins)
    shapefile_record = fc.next()
    pos_dict={}
    for shapefile_record in fc:
        shape = shapely.geometry.asShape(shapefile_record['geometry'])
        long_point = shape.centroid.coords.xy[0][0]
        lat_point = shape.centroid.coords.xy[1][0]
        pos = (long_point, lat_point)
        pos_dict[int(shapefile_record['properties']['cat'])]=pos


    # With the list of locations in hand, we can begin unpacking the forcing data. Here we'll focus on three time series: daily precipitation (`ppt`), mean daily temperature (`tmean`), and mean daily dew point (`tdmean`). 
    # 
    # Our main data structure here will be a dictionary whose keys are forcing data types and whose values are themselves dictionaries. These inner dictionaries have REW IDs as keys, and forcing data at particular instances of time as values.
    # 
    # We track time in a separate dictionary, again keyed based on forcing data. The values here are lists of dates at which each type of forcing data was observed. 
    # 
    # Note: we assume that each forcing time series is stored exactly as would be downloaded from the PRISM download page. A simple bash script for downloading daily PRISM data can be found [HERE](https://github.com/daviddralle/downloadPrism). For the purposes of the model, the data should be stored in the `raw_data` folder in a foldername equal to the variable name, for instance precipitation (`ppt`) is in `raw_data/ppt`.

    # In[4]:

    forcing_dict = {}
    dates_dict = {}

    prism_vars = ['ppt','tmean', 'tdmean']

    for prism_var in prism_vars:
        #print "Extracting " + str(prism_var) + " data..."
        years = os.listdir(os.path.join(parent_dir,'raw_data',prism_var))

        # Initialize empty data structures -- a dict keyed on REW ID with empty lists as values, 
        # and an empty list that will hold dates at which prism_var was observed
        vals_dict = {k:[] for k in pos_dict.keys()}
        date_list=[]
    
        #for all years of forcing variable, load each day of raster data and extract to all REWs
        for year in years:
            try: 
                int(year)
            except ValueError:
                continue



            # Assuming that data is in .bil format from PRISM.
            # Data for each REW is extracted to REW centroid lat/long. 
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
                
            #print str(year) + ' ' + prism_var + ' data processed.'                

        forcing_dict[prism_var] = vals_dict
        dates_dict[prism_var] = date_list


    # Before using the data for modeling, we need to ensure that the observation dates for each forcing time series align. We leverage Python's sets to make sure that all the datelists are in fact the same.

    # In[5]:

    unique_dates = [list(i) for i in set(tuple(i) for i in dates_dict.values())]
    if len(unique_dates)>1:
        print 'Forcing data dates do not all match, please check raw data!'
#     else:
#         print 'Dates match!'


    # In general, we will house all timeseries data in Pandas dataframes in order to leverage some very nice resampling functionality. For now, we'll just initialize one empty dataframe for each REW. Each of these dataframes will be indexed by date and will have one column for each forcing timeseries. 

    # In[6]:

    rng = pd.date_range(unique_dates[0][0],unique_dates[0][-1])
    rew_dfs_dict = {k:pd.DataFrame(data=None,index=rng) for k in pos_dict.keys()}

    #print('REW IDs: ' + str(pos_dict.keys()))

    # We can first set the average daily net radiation, `rn`, in units of W/m$^2$. We'll assume that this has already been somehow computed (since net radiation is not a PRISM variable) and placed as a [pickled](https://docs.python.org/2/library/pickle.html) time series, `raw_data/rn.p`. Note that the net radiation time series must include at least the date range of the other forcing data, but could also include more. For this tutorial, the `rn` dataset is not spatially distributed and is considered the same for each REW.

    # In[7]:

    df = pickle.load( open( os.path.join(parent_dir, 'raw_data','rn', 'rn.p'), "rb" ) )
    rew_id = rew_dfs_dict.keys()[0]

    #get start stop dates of extracted PRISM data. Pull rn data from these dates. 
    start = rew_dfs_dict[rew_id].index[0]
    stop = rew_dfs_dict[rew_id].index[-1]

    for rew_id in rew_dfs_dict.keys():
        rew_dfs_dict[rew_id]['rn'] = df[start:stop]


    # We can now simply cycle over the PRISM variables we enumerated earlier. Here, we might have to convert units to be in line with the main modeling effort's assumption that lengths are in units of cm and time is in units of days. 

    # In[8]:

    for prism_var in prism_vars:
        #need to convert ppt from mm/day to cm/day for PRISM data
        if prism_var=='ppt': div = 10.0
        else: div = 1.0
        for rew_id in pos_dict.keys():
            rew_dfs_dict[rew_id][prism_var] = np.array(forcing_dict[prism_var][rew_id])/div


    # The last remaining forcing timeseries required by the model presented in this tutorial is daily average potential evapotranspiration (in units of cm/day). Here, we'll use the Priestly-Taylor equation to compute this forcing.   Before we get started, we need to build some helper functions. 

    # In[9]:



    # We finally have all the pieces to compute PET. We'll use the Priestley-Taylor implementation from the `evap` module. We assume that the soil heat flux is zero at the daily timescale. Once all forcing timeseries have been computed, we save to `rew_forcing.p` in the `model_data` directory. 

    # In[10]:

    # For now, just set elevation to 250m. This is only to get REW pressure
    # for PET calculation. In the future, we'll have pressure data 
    # and this step will be unecessary. 
    rew_elevations = {i:250 for i in pos_dict.keys()}
    rew_pressures = get_rew_pressures(rew_elevations)
    for rew_id in pos_dict.keys():
        rh = get_rh(rew_dfs_dict[rew_id]['tmean'], rew_dfs_dict[rew_id]['tdmean'])

        pet = evap.Ept(
            rew_dfs_dict[rew_id]['tmean'],
            rh, 
            rew_pressures[rew_id]*np.ones(len(rh)), 
            rew_dfs_dict[rew_id]['rn'],
            np.zeros(len(rh))
            )
    
        rew_dfs_dict[rew_id]['pet'] = pet/10.0 # convert to cm/day


    # # This is new!!!! Need to document

    # In[11]:

    # the forcing dictionary must have climate_group id's as keys. 
    # for the time being, we will merge REWs within each climate group
    # using the mean of the forcings for all REWs in the climate group. 
    rew_config = pickle.load( open( os.path.join(parent_dir,'model_data','rew_config.p'), "rb" ) )
    
    #print('REW IDs identified by rew_config: ' + str(rew_config.keys()))
    climate_group_forcing = {}
    for climate_group in set([rew_config[i]['climate_group'] for i in rew_config.keys()]):
        rew_ids_in_climate_group = [i for i in rew_config.keys() if rew_config[i]['climate_group']==climate_group]
        rew_dfs  = [rew_dfs_dict[rew_id] for rew_id in rew_ids_in_climate_group]
        df = pd.concat(rew_dfs)
        df = df.groupby(df.index).mean()
        climate_group_forcing[climate_group] = df

    
    
    pickle.dump( climate_group_forcing, open( os.path.join(parent_dir,'model_data','climate_group_forcing.p'), "wb" ) )




# compute relative humidity using temp and dewpoint
def get_rh(tmean, tdmean):
    a = 17.625
    b = 243.04
    return 100*np.exp(a*tdmean/(b+tdmean))/np.exp(a*tmean/(b+tmean))


def get_rew_pressures(rew_elevations):
    #elevations in meters, output pressure in Pa
    #using simple elevation/pressure relationship
    rew_pressures = {}
    for key in rew_elevations.keys():
        rew_pressures[key] = 1000*101.325*((293-0.0065*rew_elevations[key])/293)**5.26

    return rew_pressures

# translate feature IDs to REW IDs
def translate_to_fid(parent_dir):
    basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    translated = {}

    fc = fiona.open(basins)
    shapefile_record = fc.next()
    for shapefile_record in fc:
        translated[int(shapefile_record['properties']['cat'])] = int(shapefile_record['id'])
    return translated

# translate REW IDs to feature IDs
def translate_to_rew_id(parent_dir):
    basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    translated = {}

    fc = fiona.open(basins)
    shapefile_record = fc.next()
    for shapefile_record in fc:
        translated[int(shapefile_record['id'])] = int(shapefile_record['properties']['cat'])
    return translated
    
    
if __name__ == '__main__': main()
