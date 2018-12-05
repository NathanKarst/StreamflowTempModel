import numpy as np
import os
import pickle
from datetime import date
import time
import pandas as pd
import gdal
import fiona
import shapely
from shapely import geometry,ops
from os.path import dirname
import glob
import sys
from functools import partial
import pyproj
import geopandas as gp

parent_dir = dirname(dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','2_hillslope_discharge'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','3_channel_routing'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','4_temperature'))
from vadoseZone import LaioVadoseZone, PorporatoVadoseZone, PorporatoPreferentialVadoseZone, SimpleRockMoistureZone, PreferentialRockMoistureZone
from groundwaterZone import GroundwaterZone, Melange, NonlinearReservoir, NonlinearReservoir, TwoLinearReservoir, TwoParallelLinearReservoir, LinearToNonlinearReservoir
from temperature import LagrangianSimpleTemperatureChengHeatedGW, ImplicitEulerWesthoff
from channel import TrapezoidalChannel
  
def rew_params():
    """ Write REW parameter groups and channel parameter files. 
    
    This function writes three files to the model_data directory:
        - parameter_group_params.p: parameters for each REW parameter group, as specified in model_data/rew_config.p under the 'parameter_group' key. All REWs in a given parameter group will have identical parameters. 
        - parameter_ranges.p: ranges for some subset of parameters specified for each parameter group. This file can be used fo model calibration. Parameters whose ranges are not specified are assumed constant.
        - channel_params.p: parameters for each REW's channel model. 
    Args:
        - None
        
    Returns: 
        - None. 
    """

    parent_dir = dirname(dirname(os.getcwd()))
    rew_config = pickle.load( open( os.path.join(parent_dir,'model_data','rew_config.p'), "rb" ) )
    rews = rew_config.keys()
    parameter_groups = set([rew_config[i]['parameter_group'] for i in rews])

    # Model SF @ MIRANDA - Calibrated Melange and Coastal Belt using Dry and Elder only
    # Parameter group 2 = Melange
    # Parameter group 1 = Coastal Belt
    parameter_group_params = {
    2:{'smax':19.18, 'b':2.41, 'vz':PorporatoVadoseZone, 'capacity':2.69, 'gz':Melange, 'storageVZ':10, 'storageGZ':1, 'eta':0.75, 'a':0.10},
    1:{'smaxS':2.579, 'smaxR':30.3, 'eta':1.0, 'alpha':0.296, 'res2': 1.0, 'res1': 1.0, 'gz': LinearToNonlinearReservoir , 'b': 2.09,'storageS': 1.0, 'a': 0.00216, 'k12': 0.366, 'storageR': 1.0, 'f': 0.136, 'k1': 0.308, 'vz': PreferentialRockMoistureZone },
                           }  
    parameter_ranges = {i:{'k12':(0,0.4), 'k1':(0.5,4.0), 'eta':(0.2, 1.0)} for i in parameter_groups}
    channel_params = {i:{'mannings_n':0.05, 'e':0.0108, 'f':0.3759, 'g':2.0, 'h':-0.1029, 'volume':1.0, 'model':TrapezoidalChannel} for i in rews}
    channel_params_ranges = {i:{} for i in rews}
    temperature_params = {i:{'model':ImplicitEulerWesthoff, 'alphaw':0.06 ,'cp':4186.0, 'kh':6.0,'rho':1000.0,'sigma':5.67e-8, 'temperature':11.0, 'kf':1.84, 'tau0':0.89, 'ktau':1.8, 'Tgw_offset':12.0} for i in rews}
    temperature_params_ranges = {i:{'tau0':(0.1,2), 'ktau':(0.01,10.0),'kf':(.1,20.0),'kh':(1.0,20.0), 'alphaw':(0.05, 0.3)} for i in rews}


    # Save config files
    pickle.dump( parameter_group_params, open( os.path.join(parent_dir,'model_data','parameter_group_params.p'), "wb" ) )
    pickle.dump( channel_params, open( os.path.join(parent_dir,'model_data','channel_params.p'), "wb" ) )
    pickle.dump( temperature_params, open( os.path.join(parent_dir,'model_data','temperature_params.p'), "wb" ) )
    pickle.dump( parameter_ranges, open( os.path.join(parent_dir,'model_data','parameter_ranges.p'), "wb" ) )
    pickle.dump( temperature_params_ranges, open( os.path.join(parent_dir,'model_data','temperature_params_ranges.p'), "wb" ) )
    pickle.dump( channel_params_ranges, open( os.path.join(parent_dir,'model_data','channel_params_ranges.p'), "wb" ) )

def model_config(outputFilename='model_config.p'): 
    """ Write model configuration file.
    
    This function writes a dictionary containing model configuration information to model_config.p in the model_data directory. 
    
    Keys of the dictionary include: 
        - start_date (datetime instance): Start date for running model 
        - stop_date (datetime instance): Stop date for running model
        - spinup_date (datetime instance): Date after which model has spun up. Only calibrate model after this date.
        - Tmax (float): Number of days of simulation
        - dt_hillslope (float): Hillslope simulation timestep in days
        - dt_channel (float): Channel simulation timestep in days
        - resample_freq_hillslope (float): Frequency at which to resample forcing data depending on timestep to solve the hillslope model
        - resample_freq_channel (float): Frequency at which to resample forcing data depending on timestep to solve the channel model
        - timestamps_hillslope (datetimes): Times at which to simulate hillslope dynamics
        - timestamps_channel (datetimes): Times at which to simulate channel dynamics

    Args: 
        - None
        
    Returns: 
        - None
    
    """
    #start/stop dates for running model

    # start_date = date(1981, 1, 1)             
    # spinup_date = date(1984, 10, 1)
    # stop_date = date(2017, 9, 30)
    start_date = date(2010, 1, 1)             
    spinup_date = date(2012, 10, 1)
    stop_date = date(2015, 9, 30)

    Tmax = 1.0*(stop_date - start_date).days

    #hillslope timestep information
    dt_hillslope = 1/8.
    resample_freq_hillslope = str(int(dt_hillslope*24*60)) + 'T'

    #channel timestep information
    dt_channel = 4./1440.
    resample_freq_channel = str(int(dt_channel*24*60)) + 'T'

    #temperature timestep information
    dt_temperature = 64./1440.
    resample_freq_temperature = str(int(dt_temperature*24*60)) + 'T'
    
    parent_dir = dirname(dirname(os.getcwd()))
    model_config = {'spinup_date':spinup_date, 
                    'start_date':start_date, 
                    'stop_date':stop_date, 
                    'Tmax':Tmax,
                    'dt_hillslope':dt_hillslope, 
                    'resample_freq_hillslope':resample_freq_hillslope, 
                    'dt_channel':dt_channel, 
                    'resample_freq_channel':resample_freq_channel, 
                    'dt_temperature':dt_temperature,
                    'resample_freq_temperature':resample_freq_temperature,
                    }

    pickle.dump( model_config, open( os.path.join(parent_dir,'model_data',outputFilename), "wb" ) )


def rew_config():
    """ Write REW configuration file. 
    
    This function converts GIS data related topology and physical characteristics into a dictionary, where each key is an REW id. 
    This dictionary is written to rew_config.p in the model_data directory.
    
    Each value in the dictionary is another dictionary with the attributes of each REW. 
    
    Keys (attributes) include:
        - next_stream (int): index of REW containing child stream; set to -1 if no such stream exists.
        - prev_str01 (int): index of REW containing parent stream; set to 0 if no such stream exists.
        - prev_str02 (int): index of REW containing parent stream; set to 0 if no such stream exists.
        - strahler (int): Horton-Strahler number of REW's stream
        - shreve (int): Shreve number of REW's stream
        - length (float): channel length
        - flow_accum (float): total upstream accumulated area
        - out_dist (float): distance to outlet 
        - elev_drop (float): eleveation drop from REW to outlet
        - elevation (float): elevation of REW at centroid
        - gradient (float): gradient along REW's channel
        - parameter_group: REW parameter group
        - climate_group: REW climate group  
        - group: REW group (which is just an ordered pair of the climate and parameter group numbers)
    
    Args: 
        - None
        
    Returns:
        - None
    """
    parent_dir = dirname(dirname(os.getcwd()))
    top_file = os.path.join(parent_dir,'raw_data','topology','topology.csv')
    df = pd.read_csv(top_file)
    basin_file = os.path.join(parent_dir,'raw_data','topology','basin.csv')
    df_basins = pd.read_csv(basin_file)
    
    #remove duplicates, cleanup dataframe
    df = df.drop_duplicates()
    df = df[df.cat != 0]
    df['next_stream']=df['next_stream'].astype(int)
    df['prev_str01']=df['prev_str01'].astype(int)
    df['prev_str02']=df['prev_str02'].astype(int)
    df['strahler']=df['strahler'].astype(int)
    df['shreve']=df['shreve'].astype(int)

    
    #get area of each REW into table 
    df['rew']=df['cat']

    #make sure gradient is nonzero for all channels
    for i,row in df.iterrows():
        if row.gradient==0:
            prev_2 = row['prev_str02']
            prev_1 = row['prev_str01']
            if (prev_2==0)|(prev_1==0):
                prev_grad2 = 0.001
                prev_grad1 = 0.001
            else:
                prev_grad2 = float(df.gradient.loc[df.rew==prev_2])
                prev_grad1 = float(df.gradient.loc[df.rew==prev_1])
                
            prev_grad = np.min([prev_grad2, prev_grad1])
            print('Zero gradient in REW %s! Set to minimum gradient of adjacent, upstream REWs.\n\n'%(int(row.rew)))
            df.ix[i,'gradient'] = prev_grad

    rew_config = df[['rew','next_stream','prev_str01','prev_str02','strahler','shreve','length','flow_accum','out_dist','elev_drop','gradient']].set_index('rew')
    del df_basins['label']
    df_basins = df_basins.drop_duplicates('cat')
    df_basins = df_basins.set_index('cat')
    rew_config = pd.concat([rew_config, df_basins], axis=1)

    #Assign parameter groups and climate groups
    #for the time being, assume each REW is in its own climate group
    #Each REW "Group" is the combination of its climate AND parameter group
    rew_config = _get_parameter_groups(rew_config, parent_dir)
    rew_config = _get_elevations(rew_config, parent_dir)
    rew_config = _get_climate_groups(rew_config, parent_dir) 
    rew_config = _get_interception_factor(rew_config)
    rew_config['group'] = [item for item in zip(rew_config['parameter_group'], rew_config['climate_group'])]

    #Get basins from shapefile to check that the REW ids and basins ids match    
    try:
        basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    except RuntimeError:
        print('Cannot find basins shapefile. Please make sure basins shapefile is located in \n the model directory under /raw_data/basins_poly')
    #often there will be a couple extra dangling microbasins that need to be dropped
    basins_gdf = gp.read_file(basins)
    cols = basins_gdf.columns
    basins_gdf['shapelyarea'] = basins_gdf.geometry.area
    basins_gdf = basins_gdf.sort_values('shapelyarea').drop_duplicates(subset=['cat'], keep='last')
    basins_gdf = basins_gdf[cols]
    basins_gdf.to_file(basins)
    basins_list = list(basins_gdf['cat'].values)
    basins_list.sort()
    
    x = list(df['rew'])
    x.sort()
    if x!=basins_list:
        print('REW IDs do not match basins IDs. REW config file cannot be written. \n Please clear out model data folders and re-run extract stream basins script.')        

    #get areas in cm^2, lengths in cm
    rew_config['area_sqcm'] = rew_config['area_sqkm']*10**10
    rew_config['upstream_area'] = 0
    for rew_id in rew_config.index:
        rew_config.loc[rew_config.index==rew_id,'upstream_area'] = _get_upstream_contributing(rew_config, rew_id)
    rew_config['length']=rew_config['length']*100
    rew_config['out_dist']=rew_config['out_dist']*100
    rew_config['elev_drop']=rew_config['elev_drop']*100

    #write centroids to file
    _write_coords(parent_dir=parent_dir)

    #for the sake of consistency with other config files, write as a dictionary
    rew_config = rew_config.to_dict('index')

    #save config dataframe into model_data folder
    pickle.dump( rew_config, open( os.path.join(parent_dir,'model_data','rew_config.p'), "wb" ) )
    
def _get_parameter_groups(rew_config, parent_dir):
    """ Fetch parameter groups from parameter groups raster. 
    
    If the parameter groups .tif is not located in the folder, 
    it is assumed that all REWs belong to the same parameter group
    Args:
        - rew_config 
        - parent_dir
        
    Returns: 
        - updated rew_config
    """

    raster_file = os.path.join(parent_dir,'raw_data','parameter_groups','parameter_groups.tif')
    try:
        gdata = gdal.Open(raster_file)
        gt = gdata.GetGeoTransform()
    except:
        for rew_id in rew_config.index: 
            #rew_config.set_value(rew_id, 'parameter_group', 1)
            rew_config.at[rew_id, 'parameter_group'] = 1
        return rew_config

    data = gdata.ReadAsArray().astype(np.float)
    gdata = None
    pos_dict = _get_coords(parent_dir)
    rew_config['parameter_group'] = 0

    for rew_id in rew_config.index: 
        pos = pos_dict[rew_id]
        x = int((pos[0] - gt[0])/gt[1])
        y = int((pos[1] - gt[3])/gt[5])
        try: 
            rew_config.at[rew_id, 'parameter_group'] = data[y,x]
            #rew_config.set_value(rew_id, 'parameter_group', data[y, x])
            
        except:  #default set parameter group to melange
            rew_config.at[rew_id, 'parameter_group'] = 2
            #rew_config.set_value(rew_id, 'parameter_group', 2)

    return rew_config


def _get_elevations(rew_config, parent_dir):
    """ Fetch REW elevations. 
    
    Args:
        - rew_config 
        - parent_dir
        
    Returns: 
        - updated rew_config
    """

    raster_file = os.path.join(parent_dir,'raw_data','dem','dem.tif')
    try:
        gdata = gdal.Open(raster_file)
    except:
        return [0]*len(rew_config)

    gt = gdata.GetGeoTransform()
    data = gdata.ReadAsArray().astype(np.float)
    gdata = None
    pos_dict = _get_coords(parent_dir)
    rew_config['elevation'] = 0

    for rew_id in rew_config.index: 
        pos = pos_dict[rew_id]
        x = int((pos[0] - gt[0])/gt[1])
        y = int((pos[1] - gt[3])/gt[5])
        rew_config.at[rew_id, 'elevation'] = data[y,x]
        #rew_config.set_value(rew_id, 'elevation', data[y, x])
        
    return rew_config

def _get_climate_groups(rew_config, parent_dir):
    """ Fetch climate groups from climate groups raster. 
    
    If the climate groups .tif is not located in the folder, 
    it is assumed that all REWs belong to the same parameter group
    Args:
        - rew_config 
        - parent_dir
        
    Returns: 
        - updated rew_config
    """
    raster_file = os.path.join(parent_dir,'raw_data','climate_groups','climate_groups.tif')
    try:
        # if climate groups are specified, load the climate groups raster
        gdata = gdal.Open(raster_file)
        gt = gdata.GetGeoTransform()
    except:
        # if there is no climate group raster, each REW gets its own unique climate group
        for rew_id in rew_config.index: 
            rew_config.at[rew_id, 'climate_group'] = int(rew_id)
        return rew_config

    # if climate group raster exists, extract the climate group for each REW
    data = gdata.ReadAsArray().astype(np.float)
    gdata = None
    pos_dict = _get_coords(parent_dir)
    rew_config['climate_group'] = 0
    for rew_id in rew_config.index: 
        pos = pos_dict[rew_id]
        x = int((pos[0] - gt[0])/gt[1])
        y = int((pos[1] - gt[3])/gt[5])
        try: 
            #rew_config.set_value(rew_id, 'climate_group', data[y, x])
            rew_config.at[rew_id, 'climate_group'] = data[y, x]

        except:  #default set parameter group to melange
            rew_config.at[rew_id, 'climate_group'] = 2
            #rew_config.set_value(rew_id, 'climate_group', 2)

        
    return rew_config

def _get_interception_factor(rew_config):
    # TODO: Make this a raster operation; interception factor can be specified as a raster map
    rew_config['interception_factor'] = 0.0
    for rew_id in rew_config.index: 
        if rew_config['parameter_group'].loc[rew_id]==2:
            #rew_config.set_value(rew_id, 'interception_factor', 0.1)
            rew_config.at[rew_id, 'interception_factor'] = 0.1
        elif rew_config['parameter_group'].loc[rew_id]==1:
            rew_config.at[rew_id, 'interception_factor'] = 0.4
    return rew_config


def _get_upstream_contributing(rew_config, rew_id):
    if (rew_config.prev_str02.loc[rew_id]==0) and (rew_config.prev_str01.loc[rew_id]==0):
        return rew_config.area_sqcm.loc[rew_id]
    else:
        return rew_config.area_sqcm.loc[rew_id] + _get_upstream_contributing(rew_config, rew_config.prev_str02.loc[rew_id]) + _get_upstream_contributing(rew_config, rew_config.prev_str01.loc[rew_id])

def _write_coords(parent_dir):
    basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    basins_shape = gp.GeoDataFrame.from_file(basins)
    # basins_shape['coords'] = basins_shape['geometry'].apply(lambda x: x.representative_point().coords[:])
    coords = basins_shape['geometry'].apply(lambda x: x.representative_point().coords[:])
    x = [c[0][0] for c in coords]
    y = [c[0][1] for c in coords]
    basins_shape['x'] = x
    basins_shape['y'] = y

    basins_shape_latlon = basins_shape.to_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    latlons = basins_shape_latlon['geometry'].apply(lambda x: x.representative_point().coords[:])
    lon = [c[0][0] for c in latlons]
    lat = [c[0][1] for c in latlons]
    basins_shape['lon'] = lon
    basins_shape['lat'] = lat

    basins_shape.set_index('cat').drop(['label','area_sqkm','geometry'],axis=1).to_csv(os.path.join(parent_dir,'raw_data','basins_centroids','points.csv'))

def _get_coords(parent_dir):
    basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    fc = fiona.open(basins)
    shapefile_record = fc.next()
    pos_dict={}
    for shapefile_record in fc:
        shape = shapely.geometry.asShape(shapefile_record['geometry'])
        long_point = shape.centroid.coords.xy[0][0]
        lat_point = shape.centroid.coords.xy[1][0]
        pos = (long_point, lat_point)
        pos_dict[int(shapefile_record['properties']['cat'])]=pos
    return pos_dict


def main():
    model_config()
    rew_config()
    rew_params()
    
if __name__ == '__main__': main()
    
    