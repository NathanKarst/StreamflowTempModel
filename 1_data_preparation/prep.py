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
from vadoseZone import LaioVadoseZone, PorporatoVadoseZone, SimpleRockMoistureZone
from groundwaterZone import GroundwaterZone, NonlinearReservoir, NonlinearReservoir, TwoLinearReservoir, TwoParallelLinearReservoir
from temperature import SimpleTemperature
from channel import SimpleChannel
  
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
    spinup_date = date(2013, 07, 01)             
    start_date = date(2012, 7, 01)
    stop_date = date(2014, 12, 30)
    
    Tmax = 1.0*(stop_date - start_date).days

    #hillslope timestep information
    dt_hillslope = 1/4.
    t_hillslope = np.linspace(0,Tmax,np.ceil(Tmax/dt_hillslope)+1)
    resample_freq_hillslope = str(int(dt_hillslope*24*60)) + 'T'
    timestamps_hillslope = pd.date_range(start_date, stop_date, freq=resample_freq_hillslope)

    #channel timestep information
    dt_channel = 1/1440.
    t_channel = np.linspace(0, Tmax, np.ceil(Tmax/dt_channel)+1)
    resample_freq_channel = str(int(dt_channel*24*60)) + 'T'
    timestamps_channel = pd.date_range(start_date, stop_date, freq=resample_freq_channel)

    #temperature timestep information
    dt_temperature = 1/1440.
    t_temperature = np.linspace(0, Tmax, np.ceil(Tmax/dt_temperature)+1)
    resample_freq_temperature = str(int(dt_temperature*24*60)) + 'T'
    timestamps_temperature = pd.date_range(start_date, stop_date, freq=resample_freq_temperature)
    
    parent_dir = dirname(dirname(os.getcwd()))

    model_config = {'spinup_date':spinup_date, 
                    'start_date':start_date, 
                    'stop_date':stop_date, 
                    'Tmax':Tmax,
                    'dt_hillslope':dt_hillslope, 
                    't_hillslope':t_hillslope, 
                    'resample_freq_hillslope':resample_freq_hillslope, 
                    'dt_channel':dt_channel, 
                    't_channel':t_channel,
                    'resample_freq_channel':resample_freq_channel, 
                    'timestamps_channel':timestamps_channel, 
                    'timestamps_hillslope':timestamps_hillslope,
                    'dt_temperature':dt_temperature,
                    't_temperature':t_temperature,
                    'resample_freq_temperature':resample_freq_temperature,
                    'timestamps_temperature':timestamps_temperature}

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
        - gradient (float): gradient along REW's channel
        - parameter_group: REW parameter group
        - climate_group: REW climate group  
        - group: REW group (which is just an ordered pair of the climate and parameter group numbers)
    
    Args: 
        - None
        
    Returns:
        - None
    """
    
    #get topology file generated by grass script
    #print '\n'
    #print 'Fetching REW configuration data...'
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

    #Check to make sure there are no ridiculous slope values; if anything is zero, set to upstream neighbor
    #df['gradient'].loc[df.gradient==0] = df['gradient'].loc[df.rew==df['prev_str02'].loc[df.gradient==0]]
    
    
    #get area of each REW into table 
    df['rew']=df['cat']
    rew_config = df[['rew','next_stream','prev_str01','prev_str02','strahler','shreve','length','flow_accum','out_dist','elev_drop','gradient']].set_index('rew')
    del df_basins['label']
    df_basins.set_index('cat', inplace=True)
    rew_config = pd.concat([rew_config, df_basins], axis=1)

 
    #Assign parameter groups and climate groups
    #for the time being, assume each REW is in its own climate group
    rew_config['parameter_group']=get_parameter_groups(rew_config)
    rew_config['climate_group']=range(len(rew_config))
    rew_config['group'] = zip(rew_config['parameter_group'], rew_config['climate_group'])


    #Get basins from shapefile to check that the REW ids and basins ids match
    try:
        basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    except RuntimeError:
        print 'Cannot find basins shapefile. Please make sure basins shapefile is located in \n the model directory under /raw_data/basins_poly'
    fc = fiona.open(basins)
    shapefile_record = fc.next()
    basins_list = []
    for shapefile_record in fc:
        rew_idx = int(shapefile_record['properties']['cat'])
        basins_list.append(rew_idx)
        if rew_idx == -1: continue
        shape = shapely.geometry.asShape(shapefile_record['geometry'])


    #get areas in cm^2, lengths in cm
    rew_config['area_sqcm'] = rew_config['area_sqkm']*10**10
    rew_config['upstream_area'] = 0
    for rew_id in rew_config.index:
        rew_config.loc[rew_config.index==rew_id,'upstream_area'] = _get_upstream_contributing(rew_config, rew_id)
    rew_config['length']=rew_config['length']*100
    rew_config['out_dist']=rew_config['out_dist']*100
    rew_config['elev_drop']=rew_config['elev_drop']*100


    #Print results of REW setup
    x = list(df['rew'])
    x.sort()
    basins_list.sort()
    #print('REW IDs used: ' + str(x))
    #print('Corresponding basin IDs: ' + str(basins_list))
    if x!=basins_list:
        print 'REW IDs do not match basins IDs. REW config file cannot be written. \n Please clear out model data folders and re-run extract stream basins script.'        
#     print('Total number of REW IDs used: %d'%len(x))
#     print 'Total number of unique REW parameter group(s): ' + str(len(set(rew_config['parameter_group'])))
#     print 'Total number of unique REW climate group(s): ' + str(len(set(rew_config['climate_group'])))
#     print 'Total number of unique REW group(s): ' + str(len(set(rew_config['group']))) + '\n'

    #write centroids to file
    _write_coords(parent_dir=parent_dir)

    #for the sake of consistency with other config files, write as a dictionary
    rew_config = rew_config.to_dict('index')

    #save config dataframe into model_data folder
    pickle.dump( rew_config, open( os.path.join(parent_dir,'model_data','rew_config.p'), "wb" ) )
    
def get_parameter_groups(rew_config):
    #to be used later for grouping REWs to save VadoseZone compute time
    #for now, just return a single group, 0
    return [0]*len(rew_config)


def _get_upstream_contributing(rew_config, rew_id):
    if (rew_config.prev_str02.loc[rew_id]==0) and (rew_config.prev_str01.loc[rew_id]==0):
        return rew_config.area_sqcm.loc[rew_id]
    else:
        return rew_config.area_sqcm.loc[rew_id] + _get_upstream_contributing(rew_config, rew_config.prev_str02.loc[rew_id]) + _get_upstream_contributing(rew_config, rew_config.prev_str01.loc[rew_id])

def _write_coords(parent_dir):
    basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    basins_shape = gp.GeoDataFrame.from_file(basins)
    basins_shape['coords'] = basins_shape['geometry'].apply(lambda x: x.representative_point().coords[:])
    basins_shape['coords'] = [coords[0] for coords in basins_shape['coords']]

    basins_shape.set_index('cat').drop(['label','area_sqkm','geometry'],axis=1).to_csv(os.path.join(parent_dir,'raw_data','basins_centroids','points.csv'))
    
    
def rew_params():
    """ Write REW parameter groups and channel parameter files. 
    
    This function writes three files to the model_data directory:
        - parameter_group_params.p: parameters for each REW parameter group, as specified in model_data/rew_config.p under the 'parameter_group' key. All REWs in a given parameter group will have identical parameters. 
        - parameter_ranges.p: ranges for some subset of parameters specified for each parameter group. This file can be used for model calibration. Parameters whose ranges are not specified are assumed constant.
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

    #NONLINEAR GROUNDWATER RESERVOIR, Simple rock moisture vadose zone
    parameter_group_params = {i:{'ET':0, 'emax':0.5, 'leakage':0, 'nR':0.05, 'nS':.5, 's0R':.2, 's0S':.3,'stR':.6,'stS':.5 ,'sfc':0.51, 'zrR':1000, 'zrS': 50, 'f':.7, 'storageR':0,'storageS':0, 'storageVZ':0,'storageGZ':0,'discharge':0,'a':0.0001064, 'b':3, 'vz':SimpleRockMoistureZone, 'gz':NonlinearReservoir} for i in parameter_groups}
    parameter_ranges = {i:{'zrR':(500,3000), 'zrS':(20,100), 'a':(.0001, .001), 'nR':(.01, .1)} for i in parameter_groups}
    channel_params = {i:{'mannings_n':0.03, 'e':0.01, 'f':0.39, 'volume':0, 'model':SimpleChannel} for i in rews}
    temperature_params = {i:{'cp':4186.0, 'eps':1.0, 'Tgw':14.0, 'alphaw':0.15, 'rho':1000.0, 'kh':20.0,'sigma':5.67e-8, 'temperature':15.0, 'model':SimpleTemperature} for i in rews}


    pickle.dump( parameter_group_params, open( os.path.join(parent_dir,'model_data','parameter_group_params.p'), "wb" ) )
    pickle.dump( channel_params, open( os.path.join(parent_dir,'model_data','channel_params.p'), "wb" ) )
    pickle.dump( temperature_params, open( os.path.join(parent_dir,'model_data','temperature_params.p'), "wb" ) )
    pickle.dump( parameter_ranges, open( os.path.join(parent_dir,'model_data','parameter_ranges.p'), "wb" ) )
    
    
def main():
    model_config()
    rew_config()
    rew_params()
    
if __name__ == '__main__': main()
    
    