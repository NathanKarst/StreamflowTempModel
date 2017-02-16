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
from vadoseZone import LaioVadoseZone, PorporatoVadoseZone, SimpleRockMoistureZone, PreferentialRockMoistureZone
from groundwaterZone import GroundwaterZone, LinearToNonlinearMelange, Melange, NonlinearReservoir, NonlinearReservoir, TwoLinearReservoir, TwoParallelLinearReservoir, LinearToNonlinearReservoir
from temperature import SimpleTemperature, LagrangianSimpleTemperature, EulerianWesthoff, LaxWendroffWesthoff, LagrangianSimpleTemperatureTriangular
from channel import SimpleChannel, NoChannel
  
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

    spinup_date = date(2015, 12, 05)             
    start_date = date(2013, 07, 01)
    stop_date = date(2016, 7, 30)
    
    Tmax = 1.0*(stop_date - start_date).days

    #hillslope timestep information
    dt_hillslope = 1/4.
    # t_hillslope = np.linspace(0,Tmax,np.ceil(Tmax/dt_hillslope)+1)
    resample_freq_hillslope = str(int(dt_hillslope*24*60)) + 'T'
    # timestamps_hillslope = pd.date_range(start_date, stop_date, freq=resample_freq_hillslope)

    #channel timestep information
    dt_channel = 2/1440.
    # t_channel = np.linspace(0, Tmax, np.ceil(Tmax/dt_channel)+1)
    resample_freq_channel = str(int(dt_channel*24*60)) + 'T'
    # timestamps_channel = pd.date_range(start_date, stop_date, freq=resample_freq_channel)

    #temperature timestep information
    dt_temperature = 16./1440.
    # t_temperature = np.linspace(0, Tmax, np.ceil(Tmax/dt_temperature)+1)
    resample_freq_temperature = str(int(dt_temperature*24*60)) + 'T'
    # timestamps_temperature = pd.date_range(start_date, stop_date, freq=resample_freq_temperature)
    
    parent_dir = dirname(dirname(os.getcwd()))

    model_config = {'spinup_date':spinup_date, 
                    'start_date':start_date, 
                    'stop_date':stop_date, 
                    'Tmax':Tmax,
                    'dt_hillslope':dt_hillslope, 
                    #'t_hillslope':t_hillslope, 
                    'resample_freq_hillslope':resample_freq_hillslope, 
                    'dt_channel':dt_channel, 
                    # 't_channel':t_channel,
                    'resample_freq_channel':resample_freq_channel, 
                    # 'timestamps_channel':timestamps_channel, 
                    # 'timestamps_hillslope':timestamps_hillslope,
                    'dt_temperature':dt_temperature,
                    # 't_temperature':t_temperature,
                    'resample_freq_temperature':resample_freq_temperature,
                    # 'timestamps_temperature':timestamps_temperature
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
    #Each REW "Group" is the combination of its climate AND parameter group
    rew_config = _get_parameter_groups(rew_config, parent_dir)
    rew_config = _get_elevations(rew_config, parent_dir)
    rew_config = _get_climate_groups(rew_config, parent_dir) 
    rew_config['group'] = zip(rew_config['parameter_group'], rew_config['climate_group'])
    rew_config['interception_factor'] = _get_interception_factor(rew_config)


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

    x = list(df['rew'])
    x.sort()
    basins_list.sort()
    if x!=basins_list:
        print 'REW IDs do not match basins IDs. REW config file cannot be written. \n Please clear out model data folders and re-run extract stream basins script.'        

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
    except:
        return [0]*len(rew_config)

    gt = gdata.GetGeoTransform()
    data = gdata.ReadAsArray().astype(np.float)
    gdata = None
    pos_dict = _get_coords(parent_dir)
    rew_config['parameter_group'] = 0

    for rew_id in rew_config.index: 
        pos = pos_dict[rew_id]
        x = int((pos[0] - gt[0])/gt[1])
        y = int((pos[1] - gt[3])/gt[5])
        try: 
            rew_config.set_value(rew_id, 'parameter_group', data[y, x])
        except:  #default set parameter group to melange
            rew_config.set_value(rew_id, 'parameter_group', 2)

        
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
        rew_config.set_value(rew_id, 'elevation', data[y, x])
        
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
        gdata = gdal.Open(raster_file)
    except:
        return [0]*len(rew_config)

    gt = gdata.GetGeoTransform()
    data = gdata.ReadAsArray().astype(np.float)
    gdata = None
    pos_dict = _get_coords(parent_dir)
    rew_config['climate_group'] = 0

    for rew_id in rew_config.index: 
        pos = pos_dict[rew_id]
        x = int((pos[0] - gt[0])/gt[1])
        y = int((pos[1] - gt[3])/gt[5])

        try:
            rew_config.set_value(rew_id, 'climate_group', data[y, x])
        except:
            rew_config.set_value(rew_id, 'climate_group', 2)
        
    return rew_config

def _get_interception_factor(rew_config):
    return [0.2]*len(rew_config)

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


    # daymet_angeloPPT, calibrated with logNSE on Savio
    # # FOR ELDER CREEK: Linear into nonlinear reservoir, preferential rock moisture vadose zone
    parameter_group_params = {i:{'zrS': 75., 'zrR': 869.7, 'alpha':0.119, 'res2': 1.0, 'res1': 1.0, 'gz': LinearToNonlinearReservoir , 'nR': 0.073, 'b': 2.063, 'stS': 0.6, 'storageS': 1.0, 'nS': 0.4, 'a': 0.0037, 'k12': 0.486, 'storageR': 100.0, 'f': 0.798, 's0R': 0.343, 's0S': 0.19, 'k1': 0.2477, 'stR': 0.698, 'vz': PreferentialRockMoistureZone } for i in parameter_groups}          
    parameter_ranges = {i:{ 'zrR':(500.,1200.),'k1':(0.2,0.4),'k12':(0.3,0.5),'nR':(0.01,0.4),'f':(.1,.9),'s0R':(0,.4),'stR':(0.1,0.9), 'b':(1.8,2.5), 'alpha':(.05,.95),'a':(.0005,.01)} for i in parameter_groups}
    channel_params = {i:{'mannings_n':0.1, 'e':0.01, 'f':0.39, 'volume':1.0, 'model':SimpleChannel} for i in rews}
    channel_params_ranges = {i:{'mannings_n':(.03,.15)} for i in rews}
    temperature_params = {i:{'windspeed':1.0,'c1':1.0, 'c2':1.0, 'cp':4186.0, 'eps':0.95, 'Tgw':11.0, 'alphaw':0.05, 'rho':1000.0, 'kh':5.5969,'sigma':5.67e-8, 'temperature':11.0, 'model':EulerianWesthoff} for i in rews}
    temperature_params_ranges = {i:{'kh':(0.1,20.0), 'c1':(0.1,3.0), 'c2':(0.1,3.0)} for i in rews}

    # David's "no routing" temp model test
    parameter_group_params = {i:{'zrS': 75., 'zrR': 869.7, 'alpha':0.119, 'res2': 1.0, 'res1': 1.0, 'gz': LinearToNonlinearReservoir , 'nR': 0.073, 'b': 2.063, 'stS': 0.6, 'storageS': 1.0, 'nS': 0.4, 'a': 0.0037, 'k12': 0.486, 'storageR': 100.0, 'f': 0.798, 's0R': 0.343, 's0S': 0.19, 'k1': 0.2477, 'stR': 0.698, 'vz': PreferentialRockMoistureZone } for i in parameter_groups}          
    parameter_ranges = {i:{ 'zrR':(500.,1200.),'k1':(0.2,0.4),'k12':(0.3,0.5),'nR':(0.01,0.4),'f':(.1,.9),'s0R':(0,.4),'stR':(0.1,0.9), 'b':(1.8,2.5), 'alpha':(.05,.95),'a':(.0005,.01)} for i in parameter_groups}
    channel_params = {i:{'volume':1.0, 'model':NoChannel} for i in rews}
    channel_params_ranges = {i:{ } for i in rews}
    temperature_params = {i:{'mannings_n':0.1, 'windspeed':1.0,'thetahalf':10600000000.0, 'thetamax':50.0*3.14/180, 'cp':4186.0, 'eps':0.95, 'Tgw':11.0, 'alphaw':0.05, 'rho':1000.0, 'kh':5.5969,'sigma':5.67e-8, 'temperature':11.0, 'model':LagrangianSimpleTemperatureTriangular} for i in rews}
    temperature_params_ranges = {i:{'kh':(0.1,20.0), 'c1':(0.1,3.0), 'c2':(0.1,3.0)} for i in rews}

    # Giving things a shot on Dry Creek, Melange site. 
    parameter_group_params = {i:{'gz':LinearToNonlinearMelange, 'vz': PorporatoVadoseZone, 'zr':46.58, 'sw':0.145, 'sfc':0.352, 'n':0.272, 'a':0.086, 'b':2.5, 'capacity':7.7, 'res1':1.0, 'res2':1.0,'k12':0.5, 'k1':3.6, 'storageVZ':1.0} for i in parameter_groups}          
    parameter_ranges = {i:{ 'zr':(5.0, 50.0), 'sw':(0.1,0.3), 'sfc':(0.3,0.65), 'n':(0.1,0.5), 'a':(0.1, 0.001), 'k12':(0.05,1.0), 'k1':(0.5,5.0), 'b':(1.5, 4.0), 'capacity':(1.0, 10.0)} for i in parameter_groups}
    channel_params = {i:{'volume':1.0, 'model':NoChannel} for i in rews}
    channel_params_ranges = {i:{ } for i in rews}
    temperature_params = {i:{'mannings_n':0.1, 'windspeed':1.0,'thetahalf':10600000000.0, 'thetamax':50.0*3.14/180, 'cp':4186.0, 'eps':0.95, 'Tgw':11.0, 'alphaw':0.05, 'rho':1000.0, 'kh':5.5969,'sigma':5.67e-8, 'temperature':11.0, 'model':LagrangianSimpleTemperatureTriangular} for i in rews}
    temperature_params_ranges = {i:{'kh':(0.1,20.0), 'c1':(0.1,3.0), 'c2':(0.1,3.0)} for i in rews}


    pickle.dump( parameter_group_params, open( os.path.join(parent_dir,'model_data','parameter_group_params.p'), "wb" ) )
    pickle.dump( channel_params, open( os.path.join(parent_dir,'model_data','channel_params.p'), "wb" ) )
    pickle.dump( temperature_params, open( os.path.join(parent_dir,'model_data','temperature_params.p'), "wb" ) )
    pickle.dump( parameter_ranges, open( os.path.join(parent_dir,'model_data','parameter_ranges.p'), "wb" ) )
    pickle.dump( temperature_params_ranges, open( os.path.join(parent_dir,'model_data','temperature_params_ranges.p'), "wb" ) )
    pickle.dump( channel_params_ranges, open( os.path.join(parent_dir,'model_data','channel_params_ranges.p'), "wb" ) )

    

def KML_to_params():
    try:
        basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    except RuntimeError:
        print 'Cannot find basins shapefile. Please make sure basins shapefile is located in \n the model directory under /raw_data/basins_poly'

    # basins are stored in one projection, but lithology is in another
    proj = pyproj.Proj('+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs ')

    # get REW centroids first
    fc = fiona.open(basins)
    shapefile_record = fc.next()
    pos_dict={}
    for shapefile_record in fc:
        shape = shapely.geometry.asShape(shapefile_record['geometry'])
        long_point = shape.centroid.coords.xy[0][0]
        lat_point = shape.centroid.coords.xy[1][0]
        pos = tuple(proj(long_point,lat_point,inverse=True))
        pos_dict[int(shapefile_record['properties']['cat'])]=pos
        
    print(pos_dict)
    
    # now get soil types
    careList = ['Yg','CB','C','UM','QTw']
    fc = fiona.open(os.path.join(parent_dir,'raw_data','geo','cali_geo.shp'))
    geo_dict = dict([(key,None) for key in pos_dict.keys()])
    for shapefile_record in fc:
        shape = shapely.geometry.shape(shapefile_record['geometry'])
        soilType = str(shapefile_record['properties']['Name'])
        if soilType not in careList: continue
        
        # if a centroid is in a given soil type region, associate the centroid with soil type
        done = []
        for key,val in pos_dict.items():
            if key in done: continue
            if shape.contains(shapely.geometry.Point(val)): 
                geo_dict[key] = soilType
                done.append(key)
                
    print(geo_dict)
        
    


def main():
    model_config()
    rew_config()
    rew_params()
    
if __name__ == '__main__': main()
    
    