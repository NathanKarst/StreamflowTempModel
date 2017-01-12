import os
import sys
from os.path import dirname
parent_dir = dirname(dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','2_hillslope_discharge'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','3_channel_routing'))
import random
from vadoseZone import *
import glob
from groundwaterZone import *
from REW import REW
import numpy as np
import pickle
from datetime import date
from datetime import datetime
import pandas as pd
import numpy as np
import geopandas as gp
import time
import sys
import copy
import shapely
import fiona
from ast import literal_eval as make_tuple
import multiprocessing as mp
parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir, 'StreamflowTempModel', '1_data_preparation'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','lib'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','4_temperature'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','3_channel_routing'))

rew_config = pickle.load( open( os.path.join(parent_dir,'model_data','rew_config.p'), "rb" ) )
climate_group_forcing = pickle.load( open( os.path.join(parent_dir,'model_data','climate_group_forcing.p'), "rb" ) )
model_config = pickle.load( open( os.path.join(parent_dir, 'model_data', 'model_config.p'), 'rb'))
temperature_params = pickle.load( open( os.path.join(parent_dir, 'model_data', 'temperature_params.p'), 'rb'))
temperature_params_ranges = pickle.load( open( os.path.join(parent_dir, 'model_data', 'temperature_params_ranges.p'), 'rb'))
hill_groups = pickle.load( open( os.path.join(parent_dir,'model_data','solved_hillslope_discharge.p'), "rb" ) )
solved_channel_routing = pickle.load( open( os.path.join(parent_dir,'model_data','solved_channel_routing.p'), "rb" ) )
channel_params = pickle.load( open( os.path.join(parent_dir,'model_data','channel_params.p'), "rb" ))
radiation = pickle.load( open(os.path.join(parent_dir, 'raw_data', 'radiation', 'radiation.p'),'rb') )
ta_ea = pickle.load( open(os.path.join(parent_dir, 'raw_data', 'ta_ea', 'ta_ea.p'),'rb') )
# mean_rew_lpis = pickle.load( open(os.path.join(parent_dir, 'raw_data', 'mean_rew_lpis', 'mean_rew_lpis.p'),'rb') )
mean_rew_lpis = {1:.1,2:.1,3:.1}

#start/stop dates for running model  
#spinup date is the date after start_date for which we assume model is finished spinning up         
start_date = model_config['start_date']
stop_date = model_config['stop_date']
spinup_date = model_config['spinup_date']
Tmax = model_config['Tmax']
dt = model_config['dt_temperature']
resample_freq_channel = model_config['resample_freq_channel']
resample_freq_hillslope = model_config['resample_freq_hillslope']
resample_freq_temperature = model_config['resample_freq_temperature']
t = np.linspace(0,Tmax,np.ceil(Tmax/dt)+1)
timestamps_hillslope = pd.date_range(start_date, stop_date, freq=resample_freq_hillslope)
timestamps_channel = pd.date_range(start_date, stop_date, freq=resample_freq_channel)
timestamps_temperature = pd.date_range(start_date, stop_date, freq=resample_freq_temperature)

def main(argv):
    if argv[0]=='False':
        minimize_objective_function = False
    else: 
        minimize_objective_function = True

    N = int(argv[1])
    subwatershed_calibration_name = argv[2]
    subwatershed_name = argv[3]

    ids_in_subwatershed = get_rews_to_calibrate(subwatershed_name)

    cores = mp.cpu_count()
    print('There are %s cores on this machine, \n%s model runs will be performed on each core'%(str(cores), str(N)))
    calibration_data = pickle.load( open(os.path.join(parent_dir,'calibration_data',subwatershed_calibration_name)))
    calibration_data = calibration_data[spinup_date:stop_date]

    arguments = []
    for cpu in range(0,cores):
        arguments.append((subwatershed_calibration_name, ids_in_subwatershed, N, objective_function, minimize_objective_function, cpu))

    pool = mp.Pool()

    # results is a list of 3-tuples. Each 3-tuple includes, in order, 
    # the best fit model run time series
    # the objective function value of the best run, and the index of the best run 
    # parameter set
    results = pool.map(calibrate, arguments)
    cpu_objs = [results[i][1] for i in range(cores)]
    if minimize_objective_function:
        cpu_best = np.argmin(cpu_objs)
    else:
        cpu_best = np.argmax(cpu_objs)

    best_fit = results[cpu_best][0]
    best_objective = results[cpu_best][1]
    best_params = copy.deepcopy(results[cpu_best][2])

    print('With an objective function value of %0.2f, the best parameter set is:' % (best_objective))
    print(best_params)

# Nash sutcliffe efficiency. Should be maximized for best fit. 
def objective_function(modeled, observed):
    inds = ((modeled != 0) & (observed != 0))
    if np.sum(modeled)<0.01:
        return -9999.0
    elif np.isnan(np.sum(modeled)):
        return -9999.0
    else:
        return 1-np.sum((observed.loc[inds]-modeled.loc[inds])**2)/np.sum((observed.loc[inds]-np.mean(observed.loc[inds]))**2)


def get_rews_to_calibrate(subwatershed_name):
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
    # assume the sub-watershed is contained within a single REW
    # find the id of that REW
    if len(ids_in_subwatershed)==0:
        subwatershed_shape = gp.GeoDataFrame.from_file(shapefile_path)
        basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
        with fiona.open(basins) as fiona_collection:
            for shapefile_record in fiona_collection:
                shape    = shapely.geometry.Polygon(shapefile_record['geometry']['coordinates'][0])
        if shape.contains(subwatershed_shape['geometry'].loc[0].centroid):
            ids_in_subwatershed.append(shapefile_record['properties']['cat'])


    print('REWs %s are located within the calibration sub-watershed' % str(ids_in_subwatershed))
    return ids_in_subwatershed

def generate_parameter_set(parameters, parameter_ranges, outlet_id):
    parameters_current = copy.deepcopy(parameters)
    
    for k, parameter in enumerate(parameter_ranges[outlet_id].keys()):
                new_value = random.random()*(parameter_ranges[outlet_id][parameter][1] - parameter_ranges[outlet_id][parameter][0]) + parameter_ranges[outlet_id][parameter][0]
                for rew_id in parameters.keys():
                    parameters_current[rew_id][parameter] = new_value

    return parameters_current

  
def calibrate(arguments):
    calibration_data_filename, ids_in_subwatershed, N, objective_function, minimize_objective_function, cpu = arguments
    
    # Load calibration data
    calibration_data = pickle.load( open(os.path.join(parent_dir,'calibration_data',calibration_data_filename)))
    calibration_data = calibration_data[spinup_date:stop_date]
    
    best_fit = pd.DataFrame({'modeled':np.zeros(len(timestamps_hillslope))}, index=timestamps_hillslope).resample('D').mean()
    best_parameter_set = {}
    if minimize_objective_function: 
        objs_curr = np.inf
        best_obj = np.inf
    else:
        objs_curr = -np.inf
        best_obj = -np.inf
    best_index = -1

    #Find REW with largest total contributing area
    outlet_id = ids_in_subwatershed[0]
    area_max = rew_config[outlet_id]['upstream_area']
    for rew_id in ids_in_subwatershed:
        if rew_config[rew_id]['upstream_area']>area_max:
            outlet_id = rew_id

    desc = "Core #%s"%(cpu)
    for i in range(N):

        channel_network = {}
        for rew_id in ids_in_subwatershed: 
            args = rew_config[rew_id].copy()
            args.update(channel_params[rew_id])
            channel_network[rew_id] = args['model'](rew_id=rew_id, **args)


        shreves = [rew_config[rew_id]['shreve'] for rew_id in ids_in_subwatershed]
        rewQueue = [rew_id for (shreve,rew_id) in sorted(zip(shreves,ids_in_subwatershed))]
        network_temps = {}        
        parameters_current = generate_parameter_set(temperature_params, temperature_params_ranges, outlet_id)

        temperature_network = {}
        for rew_id in ids_in_subwatershed: 
            args = rew_config[rew_id].copy()
            args.update(parameters_current[rew_id])
            temperature_network[rew_id] = args['model'](rew_id=rew_id, **args)

        for rew_id in rewQueue:
            
            shreve  = rew_config[rew_id]['shreve']
            group_id = rew_config[rew_id]['group']
            climate_group_id = group_id[1]
            rew_df = climate_group_forcing[climate_group_id]
            
            Lin = np.array(radiation[rew_id]['Lin'][start_date:stop_date].resample(resample_freq_temperature).ffill())
            Sin = np.array(radiation[rew_id]['Sin'][start_date:stop_date].resample(resample_freq_temperature).ffill())
            
            # Get other forcing data, re-sample to current timestamps
            # everything converted to meters/seconds/kelvin/kg
            temp_ea = ta_ea[rew_id].resample(resample_freq_temperature).ffill()
            ppt = 1.15741e-7*np.array(climate_group_forcing[climate_group_id][start_date:stop_date].ppt.resample(resample_freq_temperature).ffill())
            Ta = np.array(temp_ea['ta'][start_date:stop_date])
            hillslope_discharge = pd.DataFrame({'discharge':hill_groups[group_id]['discharge']}, index=hill_groups[group_id].index)
            hillslope_overlandFlow = pd.DataFrame({'overlandFlow':hill_groups[group_id]['overlandFlow']}, index=hill_groups[group_id].index)
            hillslope_volumetric_overlandFlow = 1.15741e-11*np.array(hillslope_overlandFlow[start_date:stop_date].overlandFlow.resample(resample_freq_temperature).ffill())*rew_config[rew_id]['area_sqcm']
            hillslope_volumetric_discharge = 1.15741e-11*np.array(hillslope_discharge[start_date:stop_date].discharge.resample(resample_freq_temperature).ffill())*rew_config[rew_id]['area_sqcm']
            volume = 1e-6*np.array(solved_channel_routing[rew_id][start_date:stop_date].volumes.resample(resample_freq_temperature).ffill())
            volumetric_discharge = 1.15741e-11*np.array(solved_channel_routing[rew_id][start_date:stop_date].volumetric_discharge.resample(resample_freq_temperature).ffill())
            width = 0.01*channel_network[rew_id].width
            length = 0.01*channel_network[rew_id].length

            temp = np.zeros(np.shape(t))
            
             #get upstream discharges, upstream temperatures
            if shreve == 1:
                vol_1 = np.zeros(np.shape(t))
                vol_2 = np.zeros(np.shape(t))

                temp_1 = np.zeros(np.shape(t))
                temp_2 = np.zeros(np.shape(t))
            else:
                upstream_1 = rew_config[rew_id]['prev_str01']
                upstream_2 = rew_config[rew_id]['prev_str02']

                vol_1 = 1.15741e-11*np.array(solved_channel_routing[upstream_1][start_date:stop_date].volumetric_discharge.resample(resample_freq_temperature).ffill())
                vol_2 = 1.15741e-11*np.array(solved_channel_routing[upstream_2][start_date:stop_date].volumetric_discharge.resample(resample_freq_temperature).ffill())

                temp_1 = np.array(network_temps[upstream_1].reindex(index=timestamps_temperature,fill_value=np.nan).interpolate(method='linear'))
                temp_2 = np.array(network_temps[upstream_2].reindex(index=timestamps_temperature,fill_value=np.nan).interpolate(method='linear'))


            for l in range(len(t)):
                tempArgs = {'vol_1':vol_1[l],'temp_1':temp_1[l],'vol_2':vol_2[l],'temp_2':vol_2[l],
                            'hillslope_volumetric_discharge':hillslope_volumetric_discharge[l], 'hillslope_volumetric_overlandFlow':hillslope_volumetric_overlandFlow[l], 
                            'volumetric_discharge':volumetric_discharge[l], 'volume':volume[l], 'Ta':Ta[l], 'Lin':Lin[l], 'Sin':Sin[l], 
                            'ppt':ppt[l], 'width':width,'length':length}
            
                temperature_network[rew_id].update(dt, **tempArgs)
                temp[l]=temperature_network[rew_id].temperature

            network_temps[rew_id] = pd.DataFrame(temp,index=timestamps_temperature,columns=['temperature'])


        rng = calibration_data.index
        solved_outlet = network_temps[outlet_id].temperature.reindex(rng, method='nearest')

        objs_curr = objective_function(solved_outlet[spinup_date:stop_date],calibration_data['temperature'][spinup_date:stop_date])

        if minimize_objective_function:
            if objs_curr<best_obj:
                best_obj = objs_curr
                best_fit = solved_outlet.copy()
                best_parameter_set = copy.deepcopy(parameters_current)
        else:
            if objs_curr>best_obj:
                best_obj = objs_curr
                best_fit = solved_outlet.copy()
                best_parameter_set = copy.deepcopy(parameters_current)

    return (best_fit, best_obj, best_parameter_set)


if __name__ == '__main__':
    main(sys.argv[1:])




