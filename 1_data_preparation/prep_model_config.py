import numpy as np
import os
import pickle
from datetime import date
import time
parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
import pandas as pd
  
def main(): 
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
    dt_channel = 1/360.
    t_channel = np.linspace(0, Tmax, np.ceil(Tmax/dt_channel)+1)
    resample_freq_channel = str(int(dt_channel*24*60)) + 'T'
    timestamps_channel = pd.date_range(start_date, stop_date, freq=resample_freq_channel)


    model_config = {'spinup_date':spinup_date, 'start_date':start_date, 'stop_date':stop_date, 'Tmax':Tmax,
    				'dt_hillslope':dt_hillslope, 't_hillslope':t_hillslope, 'resample_freq_hillslope':resample_freq_hillslope, 'dt_channel':dt_channel, 't_channel':t_channel,
                    'resample_freq_channel':resample_freq_channel, 'timestamps_channel':timestamps_channel, 'timestamps_hillslope':timestamps_hillslope}

    pickle.dump( model_config, open( os.path.join(parent_dir,'model_data','model_config.p'), "wb" ) )



if __name__ == '__main__':
	main()