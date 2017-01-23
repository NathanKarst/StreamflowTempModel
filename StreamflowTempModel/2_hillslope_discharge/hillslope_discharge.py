
# coding: utf-8

# # Hillslope Discharge Tutorial
# 
# This notebook outlines how to access parameters and forcing data across a collection of REWs in order to compute hillslope discharge. We pay special attention to dealing with REW parameter and climate groups, which potentially simplify calibration procedures. In this example, we lay out the general organization to deal with multiple parameter groups but will in fact only consider one.  

# In[1]:

import os
from matplotlib import pylab
import sys
from os.path import dirname
parent_dir = dirname(dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','2_hillslope_discharge'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','3_channel_routing'))

from vadoseZone import *
from groundwaterZone import *
from REW import REW
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import pickle
from datetime import date
import pandas as pd
import numpy as np
import time
import sys


def main():

    # First, we'll need to set up some data structures to hold REW parameters and forcing data. We'll assume that these have already been computed and have been stored in the `model_data` subfolder of the parent folder. For information on the creation of these data, refer to the [model parameterization](./parameterize_and_configure.ipynb) or [data preparation](./prism_forcing.ipynb) notebooks. 
    # 
    # The `rew_config` dictionary has as its keys the list of all REWs to be modeled in the simulation. The value of each key contains data primarily related to physical topological attributes of each REW. The most important field in this example will be `group`, which specifies which collection of parameters we will be using to configure each REW. See the description of `group_params` below for more details.
    # 
    # The `rew_forcing` dictionary contains a number of forcing time series. The only one we will be using here is precipiation, `ppt`. 
    # 
    # The `group_params` dictionary has as its keys the list of all REW parameter groups. For each group, we can specific both scalars like field capacity or maximum evapotranspiration, and classes like specific groundwater and vadose zone mdoels. **The parameters enumerated for each group must include those needed to fully populate the chosen vadose zone and groundwater zone models, and the naming conventions between `group_params` and each class must match.** For instance, if the chosen vadose zone model has a field capacity attribute named `sfc`, then each value in the `group_params` dictionary must have a `sfc` key. Each REW will be instantiated using the parameters corresponding to the its `group` field in `rew_config`. This allows us to perform expensive calibration procedures only a few representative REW "types" (e.g., coastal belt covered in a Douglas fir), rather than on each REW individually. 
    # 
    # The `model_config` dictionary contains the high level parameters needed to spool up and execute the time domain simulation. Values include time stamps at which hillslope discharge will be computed, as well as the time step lengths and resampling frequencies necessary to map forcing data onto the desired simulation timescales. 

    # In[2]:

    # Load config files, forcing file, and paramters for each group
    parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))

    sys.path.append(os.path.join(parent_dir, 'StreamflowTempModel', '1_data_preparation'))

    # These dictionaries contain the all the data we'll need to instantiate 
    rew_config = pickle.load( open( os.path.join(parent_dir,'model_data','rew_config.p'), "rb" ) )
    climate_group_forcing = pickle.load( open( os.path.join(parent_dir,'model_data','climate_group_forcing.p'), "rb" ) )
    parameter_group_params = pickle.load( open( os.path.join(parent_dir,'model_data','parameter_group_params.p'), "rb" ))
    model_config = pickle.load( open( os.path.join(parent_dir, 'model_data', 'model_config.p'), 'rb'))


    # Before starting the simulation, it's convenient to load some of the model configuration data into the local workspace, rather than continually pulling from the `model_config` dictionary. 
    # 
    # The most important quantities here are the `start_date`, `stop_date` and `spinup_date`. The first two are self-explanatory; the last specifies the date at which we assume transients related, for instance, to starting with empty groundwater and vadose zone stocks have exited the system. 
    # 
    # Note that `model_config` typically contains datetime indices, time steps, and resampling frequencies related to both the hillslope dischage _and_ channel routing models. Since we'll only be dealing with the hillslope part of the process in this example, we'll simplify the notation. 

    # In[3]:

    start_date = model_config['start_date']
    stop_date = model_config['stop_date']
    spinup_date = model_config['spinup_date']
    Tmax = model_config['Tmax']
    dt = model_config['dt_hillslope']
    t = np.linspace(0,Tmax,np.ceil(Tmax/dt)+1)
    resample_freq_hillslope = model_config['resample_freq_hillslope']
    timestamps_hillslope = pd.date_range(start_date, stop_date, freq=resample_freq_hillslope)


    # With the forcing data, REW parameterizations, and overarching model description in hand, we're finally ready to run the model. We assume that in terms of hillslope discharge, no REW depends on any other, and so we can simply simulate discharge sequentially along the list of REWs. For each REW, we
    # 
    # * Instantiate a `REW` instance using the parameterization provided in `group_params`. For each REW, the appropriate group ID is located in `rew_config[REW ID]['group']`;
    # * Intialize empty vadose zone and groundwater zone stocks, as well as empty time series related to the major fluxes of each zone;
    # * Populate the attributes of both the vadose and groundwater zones. At REW instantiation, both the vadose zone and groundwater zone are empty. Therefore, the `group_params` dictionary **must** specify all quantities necessary for model execution, and these quantities must have the same variable name in both `group_params` and the corresponding model class; 
    # * Resample forcing data at the hillslope discharge simulation timescale specified in `model_config`;
    # * Simulate hillslope discharge at each time specified in the `t_hillsope` field of `model_config`;
    # * Store the flux and stock time series data in a Pandas data frame for each REW, and bundle these Pandas data frames into a dictionary of solved hillslopes indexed by its group ID;
    # * Save these solved hillslopes to the pickle file `solved_group_hillslopes` in the `model_data folder`.

    # In[4]:

    group_ids = [rew_config[i]['group'] for i in rew_config.keys()]
    solved_group_hillslopes_dict = {}

    for group_id in group_ids:

        parameter_group_id = group_id[0]
        climate_group_id = group_id[1]
        
        vz = parameter_group_params[parameter_group_id]['vz'](**parameter_group_params[parameter_group_id])
        gz = parameter_group_params[parameter_group_id]['gz'](**parameter_group_params[parameter_group_id])    
    
        rew = REW(vz, gz,  **{'pet':climate_group_forcing[climate_group_id].pet, 'ppt':climate_group_forcing[climate_group_id].ppt, 'aspect':90})

        storageVZ    = np.zeros(np.size(t))
        storageGZ     = np.zeros(np.size(t))
        discharge       = np.zeros(np.size(t))
        leakage         = np.zeros(np.size(t))
        ET              = np.zeros(np.size(t))
        overlandFlow = np.zeros(np.size(t))

        # Resample pet and ppt to integration timestep
        ppt = np.array(rew.ppt[start_date:stop_date].resample(resample_freq_hillslope).ffill())
        pet = np.array(rew.pet[start_date:stop_date].resample(resample_freq_hillslope).ffill())

        # Solve group hillslope
        for i in range(len(t)):
            rew.vz.update(dt,**{'ppt':ppt[i],'pet':pet[i]})
            storageVZ[i] = rew.vz.storageVZ
            leakage[i]      = rew.vz.leakage
            ET[i]           = rew.vz.ET   
            rew.gz.update(dt,**{'leakage':leakage[i]})
            storageGZ[i] = rew.gz.storageGZ
            discharge[i] = rew.gz.discharge
            overlandFlow[i] = rew.vz.overlandFlow
            
        totalIn = ppt.cumsum()*dt + storageVZ[0] + storageGZ[0] + discharge[0]*dt + ET[0]*dt
        totalPresent = storageVZ + storageGZ
        totalOut = discharge.cumsum()*dt + ET.cumsum()*dt
        
        balance = (totalOut + totalPresent)/totalIn 
        print(max(balance))
        print(min(balance))


        # Save all results as daily data. 
        solved_hillslope = pd.DataFrame({'storageVZ':storageVZ, 'leakage':leakage, 'ET':ET, 'storageGZ':storageGZ, 'discharge':discharge, 'overlandFlow':overlandFlow}, index=timestamps_hillslope)
        solved_group_hillslopes_dict[group_id] = solved_hillslope.resample('D').last()
        
    

    
    

    pickle.dump( solved_group_hillslopes_dict, open( os.path.join(parent_dir,'model_data','solved_hillslope_discharge.p'), "wb" ) )
    
if __name__ == '__main__': main()