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
import mpld3
import time
import sys
import copy
import shapely
import fiona
from ast import literal_eval as make_tuple
import hillslope_calibration
import multiprocessing as mp

parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir, 'StreamflowTempModel', '1_data_preparation'))
from prep import rew_params
rew_params()

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

subwatershed_name = 'elder'
groups_to_calibrate, ids_in_subwatershed = hillslope_calibration.get_groups_to_calibrate(subwatershed_name + str('.shp'))

# Nash sutcliffe efficiency. Should be maximized for best fit. 
def objective_function(modeled, observed):
    inds = ((modeled != 0) & (observed != 0))
    if np.sum(modeled)<0.01:
        return -9999.0
    elif np.isnan(np.sum(modeled)):
        return -9999.0
    else:
        return 1-np.sum((observed.loc[inds]-modeled.loc[inds])**2)/np.sum((observed.loc[inds]-np.mean(observed.loc[inds]))**2)
minimize_objective_function = False


#specify the number of parameter sets to generate
N = 250000
cores = mp.cpu_count()
print('There are %s cores on this machine, \n%s model runs will be performed on each core'%(str(cores), str(N)))



subwatershed_calibration_name = 'elder_runoff' + '.p'
calibration_data = pickle.load( open(os.path.join(parent_dir,'calibration_data',subwatershed_calibration_name)))
calibration_data = calibration_data[spinup_date:stop_date]

parameters_per_core = {}
arguments = []
for cpu in range(0,cores):
    parameters_per_core[cpu] = hillslope_calibration.generate_parameter_sets(N, parameter_group_params, parameter_ranges)
    arguments.append((subwatershed_calibration_name, groups_to_calibrate, ids_in_subwatershed, parameters_per_core[cpu], objective_function, minimize_objective_function, cpu))

pool = mp.Pool()

# results is a list of 3-tuples. Each 3-tuple includes, in order, 
# the best fit model run time series
# the objective function value of the best run, and the index of the best run 
# parameter set
results = pool.map(hillslope_calibration.calibrate, arguments)
cpu_objs = [results[i][1] for i in range(cores)]
if minimize_objective_function:
    cpu_best = np.argmax(cpu_objs)
else:
    cpu_best = np.argmin(cpu_objs)

best_index = results[cpu_best][2]
best_fit = results[cpu_best][0]
best_objective = results[cpu_best][1]

print('With an objective function value of %0.2f, the best parameter set is:' % (np.max(best_objective)))
print(parameters_per_core[cpu_best][best_index])




