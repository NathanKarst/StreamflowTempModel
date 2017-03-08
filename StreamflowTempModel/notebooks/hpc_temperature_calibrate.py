

import pandas as pd
import gdal
import os
import datetime
from datetime import timedelta
import numpy as np
import fiona
import shapely
from shapely import geometry
from os.path import dirname
import glob
import sys
import pickle
import copy
import random
from functools import partial
parent_dir = dirname(dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','lib'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','4_temperature'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','3_channel_routing'))
from temperature import SimpleTemperature, LagrangianSimpleTemperature
from channel import SimpleChannel
import zonal_stats as zs
import meteolib as meteo
import evaplib as evap
from ast import literal_eval as make_tuple
import multiprocessing as mp


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

    print(subwatershed_calibration_name)

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
    modeled = modeled.resample('D').mean()
    observed = observed.resample('D').mean()
    inds = ((modeled != 0) & (observed != 0))&((modeled.index.month>=5)&(modeled.index.month<=9))
    # if np.sum(modeled)<0.01:
    #     return -9999.0
    # elif np.isnan(np.sum(modeled)):
    #     return -9999.0
    # else:
    #     return 1-np.sum((observed.loc[inds]-modeled.loc[inds])**2)/np.sum((observed.loc[inds]-np.mean(observed.loc[inds]))**2)

    if np.sum(modeled)<0.01:
        return np.inf
    elif np.isnan(np.sum(modeled)):
        return np.inf
    else:
        return np.sqrt(1.0/len(observed.loc[inds])*np.sum((observed.loc[inds] - modeled.loc[inds])**2))


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
        sys.stdout.write('\rWorking on iteration %d out of %d \n' % (i+1,N))
        sys.stdout.flush()
        
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
            width = channel_network[rew_id].width
            length = channel_network[rew_id].length
            
            Lin = np.array(radiation[rew_id]['Lin'][start_date:stop_date].resample(resample_freq_temperature).ffill())
            Sin = np.array(radiation[rew_id]['Sin'][start_date:stop_date].resample(resample_freq_temperature).ffill())
            doy = np.array(radiation[rew_id]['doy'][start_date:stop_date].resample(resample_freq_temperature).ffill())

            temp_ea = ta_ea[rew_id].resample(resample_freq_temperature).interpolate()
            ppt_daily = climate_group_forcing[climate_group_id][start_date:stop_date].ppt
            ppt = np.array(climate_group_forcing[climate_group_id][start_date:stop_date].ppt.resample(resample_freq_temperature).ffill())
            ta = np.array(temp_ea['ta'][start_date:stop_date])
            ea = np.array(temp_ea['ea'][start_date:stop_date])
            hillslope_discharge = pd.DataFrame({'discharge':hill_groups[group_id]['discharge']}, index=hill_groups[group_id].index)
            hillslope_overlandFlow = pd.DataFrame({'overlandFlow':hill_groups[group_id]['overlandFlow']}, index=hill_groups[group_id].index)
            
            hillslope_volumetric_overlandFlow = np.array(hillslope_overlandFlow[start_date:stop_date].overlandFlow.resample(resample_freq_temperature).ffill())*rew_config[rew_id]['area_sqcm']
            hillslope_volumetric_discharge = np.array(hillslope_discharge[start_date:stop_date].discharge.resample(resample_freq_temperature).ffill())*rew_config[rew_id]['area_sqcm']
            hillslope_volumetric_discharge_daily = hillslope_discharge[start_date:stop_date].discharge*rew_config[rew_id]['area_sqcm']
            
            volumetric_discharge = np.array(solved_channel_routing[rew_id][start_date:stop_date].volumetric_discharge.resample(resample_freq_temperature).ffill())
            volumetric_discharge_daily = solved_channel_routing[rew_id][start_date:stop_date].volumetric_discharge
            temp = np.zeros(np.shape(t))
            
            #get upstream discharges, upstream temperatures
            vol_1_daily = 0
            vol_2_daily = 0
            if shreve == 1:
                vol_1 = np.zeros(np.shape(t))
                vol_2 = np.zeros(np.shape(t))

                temp_1 = np.zeros(np.shape(t))
                temp_2 = np.zeros(np.shape(t))
            else:
                upstream_1 = rew_config[rew_id]['prev_str01']
                upstream_2 = rew_config[rew_id]['prev_str02']

                vol_1 = np.array(solved_channel_routing[upstream_1][start_date:stop_date].volumetric_discharge.resample(resample_freq_temperature).ffill())
                vol_2 = np.array(solved_channel_routing[upstream_2][start_date:stop_date].volumetric_discharge.resample(resample_freq_temperature).ffill())

                vol_1_daily = solved_channel_routing[upstream_1][start_date:stop_date].volumetric_discharge[0]
                vol_2_daily = solved_channel_routing[upstream_2][start_date:stop_date].volumetric_discharge[0]
                
                temp_1 = np.array(network_temps[upstream_1].temperature)
                temp_2 = np.array(network_temps[upstream_2].temperature)

            # Now get volumes in channel link. 
            volume = np.array(solved_channel_routing[rew_id][start_date:stop_date].volumes.resample(resample_freq_temperature).interpolate())
            start_temp_model = int(1/dt*(len(pd.date_range(start_date,spinup_date))-365))
            w = 0
            for l in range(len(t)):
                if l<start_temp_model:
                    temp[l] = temperature_network[rew_id].temperature
                else:  
                    varyArgs = ['doy','vol_1','temp_1','vol_2','temp_2','hillslope_volumetric_discharge', 'hillslope_volumetric_overlandFlow', 'volumetric_discharge', 'volume', 'ta', 'Lin', 'Sin', 'ppt', 'ea']
                    constArgs = ['width','length']
                    tempArgs = {}
                    for arg in varyArgs: tempArgs[arg] = copy.copy(locals()[arg][l])
                    for arg in constArgs: tempArgs[arg] = copy.copy(locals()[arg])
                    
                    temp[l]=temperature_network[rew_id].temperature
                    temperature_network[rew_id].update(dt, **tempArgs)
                    

            network_temps[rew_id] = pd.DataFrame(temp,index=timestamps_temperature,columns=['temperature'])


        rng = calibration_data.index
        solved_outlet = network_temps[outlet_id].temperature.reindex(rng, method='nearest')
        objs_curr = objective_function(solved_outlet[spinup_date:stop_date],calibration_data['temperature'][spinup_date:stop_date])
        print(objs_curr)
        print('\n')
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
    print best_obj
    return (best_fit, best_obj, best_parameter_set)


if __name__ == '__main__':
    main(sys.argv[1:])




