import os
import sys
from os.path import dirname
parent_dir = dirname(dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','2_hillslope_discharge'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','3_channel_routing'))
import random
random.seed()
from vadoseZone import *
import glob
from groundwaterZone import *
from REW import REW
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import pickle
from datetime import date
import pandas as pd
import numpy as np
import geopandas as gp
import time
from tqdm import tnrange, tqdm_notebook
import sys
import shapely
import fiona
import folium
from multiprocessing import Process, Queue, current_process, freeze_support
from ast import literal_eval as make_tuple

parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir, 'StreamflowTempModel', '1_data_preparation'))
rew_config = pickle.load( open( os.path.join(parent_dir,'model_data','rew_config.p'), "rb" ) )
climate_group_forcing = pickle.load( open( os.path.join(parent_dir,'model_data','climate_group_forcing.p'), "rb" ) )
parameter_group_params = pickle.load( open( os.path.join(parent_dir,'model_data','parameter_group_params.p'), "rb" ))
model_config = pickle.load( open( os.path.join(parent_dir, 'model_data', 'model_config.p'), 'rb'))
parameter_ranges = pickle.load( open( os.path.join(parent_dir, 'model_data', 'parameter_ranges.p'), 'rb'))
start_date = model_config['start_date']
stop_date = model_config['stop_date']
spinup_date = model_config['spinup_date']
Tmax = model_config['Tmax']
dt = model_config['dt_hillslope']
t = model_config['t_hillslope']
resample_freq_hillslope = model_config['resample_freq_hillslope']
timestamps_hillslope = model_config['timestamps_hillslope']


def get_groups_to_calibrate(subwatershed_name):
    shapefile_path = os.path.join(parent_dir, 'raw_data','watershed_poly', subwatershed_name)

    # get coordinate tuples corresponding to each REW
    # these are stored in the x, y columns of points
    points = pd.read_csv( os.path.join(parent_dir, 'raw_data', 'basins_centroids', 'points.csv')).set_index('cat')

    # check to see which REWs fall within sub-watershed
    ids_in_subwatershed = []
    with fiona.open(shapefile_path) as fiona_collection:
        for shapefile_record in fiona_collection:
            # note: the shapefile record must be of type polygon, not multi-polygon
            # i.e. the sub-watershed must be a single polygon
            shape = shapely.geometry.Polygon( shapefile_record['geometry']['coordinates'][0] )

            for index, row in points.iterrows(): 
                point =  shapely.geometry.Point(row.x, row.y)
                if shape.contains(point):
                    ids_in_subwatershed.append(index)

    ids_in_subwatershed = list(set(ids_in_subwatershed))

    # if no REWs found inside sub-watershed, 
    # assume the sub-watershed is contained within a single REW. 
    # Here, find the id of that REW
    if len(ids_in_subwatershed)==0:
        subwatershed_shape = gp.GeoDataFrame.from_file(shapefile_path)
        basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
        with fiona.open(basins) as fiona_collection:
            for shapefile_record in fiona_collection:
                shape    = shapely.geometry.Polygon(shapefile_record['geometry']['coordinates'][0])
        if shape.contains(subwatershed_shape['geometry'].loc[0].centroid):
            ids_in_subwatershed.append(shapefile_record['properties']['cat'])


    groups_to_calibrate = []
    for rew_id in ids_in_subwatershed:
        groups_to_calibrate.append(rew_config[rew_id]['group'])

    print('REWs %s are located within the calibration sub-watershed' % str(ids_in_subwatershed))
    print('The groups %s will be run for calibration purposes' % str(groups_to_calibrate))
    return groups_to_calibrate, ids_in_subwatershed

def generate_parameter_sets(N, parameter_group_params, parameter_ranges):
    parameter_realz = []
    for i in range(N):
        parameter_group_params_current = {}
        for w in parameter_group_params.keys():
            parameter_group_params_current[w] = parameter_group_params[w].copy()

        for j, parameter_group in enumerate(parameter_ranges.keys()):
                for k, parameter in enumerate(parameter_ranges[parameter_group].keys()):
                    new_value = random.random()*(parameter_ranges[parameter_group][parameter][1] - parameter_ranges[parameter_group][parameter][0]) + parameter_ranges[parameter_group][parameter][0]
                    parameter_group_params_current[parameter_group][parameter] = new_value
        parameter_realz.append(parameter_group_params_current)

    return parameter_realz
    
    
    
def calibrate(arguments):
    calibration_data_filename, groups_to_calibrate, ids_in_subwatershed, parameter_realz, objective_function, minimize_objective_function, cpu = arguments
    
    # Load calibration data
    calibration_data = pickle.load( open(os.path.join(parent_dir,'calibration_data',calibration_data_filename)))
    calibration_data = calibration_data[spinup_date:stop_date]
    
    N = len(parameter_realz)

    # for each parameter realization
    objs = []
    best_fit = pd.DataFrame({'modeled':np.zeros(len(timestamps_hillslope))}, index=timestamps_hillslope).resample('D').mean()
    if minimize_objective_function: 
        objs_curr = np.inf
    else:
        objs_curr = -np.inf

    best_index = -1
    desc = "Core #%s"%(cpu)
    for i in tqdm_notebook(range(N), desc=desc):
        solved_groups = {}
        parameter_group_params = {}
        parameter_group_params = parameter_realz[i]

        solved_group_hillslopes_dict = {}
        for group_id in groups_to_calibrate:

            parameter_group_id = group_id[0]
            climate_group_id = group_id[1]

            vz = parameter_group_params[parameter_group_id]['vz'](**parameter_group_params[parameter_group_id])
            gz = parameter_group_params[parameter_group_id]['gz'](**parameter_group_params[parameter_group_id])    

            rew = REW(vz, gz,  **{'pet':climate_group_forcing[climate_group_id].pet, 'ppt':climate_group_forcing[climate_group_id].ppt, 'aspect':90})

            # storageVZ    = np.zeros(np.size(t))
            # storageGZ     = np.zeros(np.size(t))
            discharge       = np.zeros(np.size(t))
            leakage         = np.zeros(np.size(t))
            # ET              = np.zeros(np.size(t))

            # Resample pet and ppt to integration timestep
            ppt = np.array(rew.ppt[start_date:stop_date].resample(resample_freq_hillslope).ffill())
            pet = np.array(rew.pet[start_date:stop_date].resample(resample_freq_hillslope).ffill())

            # Solve group hillslope
            for l in range(len(t)):
                rew.vz.update(dt,**{'ppt':ppt[l],'pet':pet[l]})
                # storageVZ[l] = rew.vz.storageVZ
                leakage[l]      = rew.vz.leakage
                # ET[l]           = rew.vz.ET   
                rew.gz.update(dt,**{'leakage':leakage[l]})
                # storageGZ[l] = rew.gz.storageGZ
                discharge[l] = rew.gz.discharge

            # resample as daily data
            solved_groups[group_id] = pd.DataFrame({'discharge':discharge}, index=timestamps_hillslope).resample('D').mean()

        total_area = 0
        for rew_id in ids_in_subwatershed:
            total_area += rew_config[rew_id]['area_sqkm']

        name = str(i) + 'discharge'
        solved_subwatershed = pd.DataFrame({name:np.zeros(len(timestamps_hillslope))}, index=timestamps_hillslope).resample('D').mean()

        solved_subwatershed_array = np.zeros(int(len(solved_subwatershed)))
        for rew_id in ids_in_subwatershed:
            solved_subwatershed_array += rew_config[rew_id]['area_sqkm']/total_area*solved_groups[rew_config[rew_id]['group']]['discharge']

        solved_subwatershed[name] = solved_subwatershed_array
        objs.append(objective_function(solved_subwatershed[name][spinup_date:stop_date],calibration_data['runoff'][spinup_date:stop_date]))
        
        if minimize_objective_function:
            if objs[i]<objs_curr:
                objs_curr = objs[i]
                best_index = i
                best_fit = solved_subwatershed[name].copy()
                print('Min objective function value so far is: ' + str(objs_curr))
        else:
            if objs[i]>objs_curr:
                objs_curr = objs[i]
                best_index = i
                best_fit = solved_subwatershed[name].copy()
    
    return (best_fit, objs_curr, best_index)

