import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os
import pickle
from datetime import date
import pandas as pd
import spotpy
import time
from channel import SimpleChannel

#must import all hillslope modules to set models for each group
import sys
from os.path import dirname
parent_dir = dirname(dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','2_hillslope_discharge'))
from vadoseZone import LaioVadoseZone, PorporatoVadoseZone
from REW import REW
from groundwaterZone import GroundwaterZone, NonlinearReservoir, NonlinearReservoir, TwoLinearReservoir, TwoParallelLinearReservoir

  
def main(): 

    hill_groups = pickle.load( open( os.path.join(parent_dir,'model_data','solved_group_hillslopes_dict.p'), "rb" ) )
    rew_config = pickle.load( open( os.path.join(parent_dir,'model_data','rew_config.p'), "rb" ) )
    group_params = pickle.load( open( os.path.join(parent_dir,'model_data','group_params.p'), "rb" ))
    model_config = pickle.load( open( os.path.join(parent_dir, 'model_data', 'model_config.p'), 'rb'))
    rew_forcing_dict = pickle.load( open( os.path.join(parent_dir,'model_data','rew_forcing_dict.p'), "rb" ) )
    channel_params = pickle.load( open( os.path.join(parent_dir,'model_data','channel_params.p'), "rb" ))

    #start/stop dates for running model  
    #spinup date is the date after start_date for which we assume model is finished spinning up         
    start_date = model_config['start_date']
    stop_date = model_config['stop_date']
    spinup_date = model_config['spinup_date']
    Tmax = model_config['Tmax']
    dt = model_config['dt_channel']
    t = model_config['t_channel']
    resample_freq_channel = model_config['resample_freq_channel']
    resample_freq_hillslope = model_config['resample_freq_hillslope']
    timestamps_hillslope = model_config['timestamps_hillslope']
    timestamps_channel = model_config['timestamps_channel']

    #create a dictionary of channels, one entry for each REW
    channel_network = {}
    for rew_id in rew_config.index:
        mannings_n = channel_params[rew_id]['mannings_n']
        e = channel_params[rew_id]['e']
        f = channel_params[rew_id]['f']
        volume = channel_params[rew_id]['volume']
        channel_network[rew_id] = SimpleChannel(rew_id=rew_id, rew_config=rew_config, volume=volume, mannings_n=mannings_n, e=e, f=f)


    network_volumetric_discharges = {}
    network_volumes = {}

    t0 = time.time()
    #depth first traversal of network to solve rew flow
    for order in range(1,rew_config.shreve.max()+1):
        for rew_id in rew_config.loc[rew_config.shreve==order].index:
            if order==1:
                #print 'Order = ' + str(order) + '; rew_id = ' + str(rew_id)
                group_id = rew_config.group.loc[rew_id]
                ppt = np.array(rew_forcing_dict[rew_id][start_date:stop_date].ppt.resample(resample_freq_channel).ffill())
                hillslope_discharge = pd.DataFrame({'discharge':hill_groups[group_id]['discharge']}, index=hill_groups[group_id].index)
                hillslope_volumetric_discharge = np.array(hillslope_discharge[start_date:stop_date].discharge.resample(resample_freq_channel).ffill())*rew_config.area_sqcm.loc[rew_id]
                volumetric_discharge = np.zeros(np.size(t))
                volumes = np.zeros(np.size(t))
                approx = 0
                for i in range(len(t)):
                    result = channel_network[rew_id].update(dt, upstream_volumetric_discharge=0, hillslope_volumetric_discharge=hillslope_volumetric_discharge[i] , ppt=ppt[i])
                    approx = (approx or result)
                    volumetric_discharge[i]=channel_network[rew_id].volumetric_discharge
                    volumes[i] = channel_network[rew_id].volume

                if approx==1: 
                    print '\nWarning: Numerical instability encountered. Consider decreasing timestep size. \nDischarge for REW ' + str(rew_id) + ' had to be approximated for some timesteps. \n'

                network_volumetric_discharges[rew_id]=pd.DataFrame({'volumetric_discharge':volumetric_discharge}, index=timestamps_channel).resample('D').mean()
                network_volumes[rew_id] = pd.DataFrame({'volumes':volumes}, index=timestamps_channel).resample('D').mean()
            else:
                #print 'Order = ' + str(order) + '; rew_id = ' + str(rew_id)
                group_id = rew_config.group.loc[rew_id]
                ppt = np.array(rew_forcing_dict[rew_id][start_date:stop_date].ppt.resample(resample_freq_channel).ffill())
                hillslope_discharge = pd.DataFrame({'discharge':hill_groups[group_id]['discharge']}, index=hill_groups[group_id].index)
                hillslope_volumetric_discharge = np.array(hillslope_discharge[start_date:stop_date].discharge.resample(resample_freq_channel).ffill())*rew_config.area_sqcm.loc[rew_id]
                volumetric_discharge = np.zeros(np.size(t))
                up = np.array(network_volumetric_discharges[rew_config.prev_str01.loc[rew_id]].volumetric_discharge.resample(resample_freq_channel).ffill() + network_volumetric_discharges[rew_config.prev_str02.loc[rew_id]].volumetric_discharge.resample(resample_freq_channel).ffill())
                volumes = np.zeros(np.size(t))
                approx = 0
                for i in range(len(t)):
                    result = channel_network[rew_id].update(dt, upstream_volumetric_discharge=up[i], hillslope_volumetric_discharge=hillslope_volumetric_discharge[i] , ppt=ppt[i])
                    approx = (approx or result)
                    volumetric_discharge[i]=channel_network[rew_id].volumetric_discharge
                    volumes[i] = channel_network[rew_id].volume

                if approx==1: 
                    print '\nWarning: Numerical instability encountered. Consider decreasing timestep size. \nDischarge for REW ' + str(rew_id) + ' had to be approximated for some timesteps. \n'

                network_volumetric_discharges[rew_id]=pd.DataFrame({'volumetric_discharge':volumetric_discharge}, index=timestamps_channel).resample('D').mean()
                network_volumes[rew_id] = pd.DataFrame({'volumes':volumes}, index=timestamps_channel).resample('D').mean()

    print 'Kinematic wave solution took ' + str((time.time()-t0)) + ' seconds'
    pickle.dump( network_volumetric_discharges, open( os.path.join(parent_dir,'model_data','network_volumetric_discharges.p'), "wb" ) )
    pickle.dump( network_volumes, open( os.path.join(parent_dir,'model_data','network_volumes.p'), "wb" ) )


    #plot data
    discharge = np.array(network_volumetric_discharges[2].volumetric_discharge/rew_config.upstream_area.loc[2])
    elder_df = pickle.load( open( os.path.join(parent_dir, 'calibration_data', 'elder_2010_2015.p'), 'rb'))
    elder_runoff_df = elder_df['runoff'][start_date:stop_date]
    elder_runoff = np.array(elder_df['runoff'][start_date:stop_date])
    outflows_df = pd.DataFrame({'data':elder_runoff, 'modeled':discharge}, index = elder_runoff_df.index)
    outflows_df.plot()
    plt.show()

            

        # else:
        #     for rew_id in rew_config.loc[rew_config.shreve==order].index:
        #         group_id = rew_config.group.loc[rew_id]
        #         ppt = np.array(rew_forcing_dict[rew_id][start_date:stop_date].ppt.resample(resample_freq).ffill())
        #         hillslope_volumetric_discharge = hill_groups[group_id]['discharge']*rew_config.area_sqcm.loc[rew_id]
        #         volumetric_discharge = np.zeros(np.size(t))
        #         up_discharge = network_volumetric_discharges[rew_config.prev_str01.loc[rew_id]]+network_volumetric_discharges[rew_config.prev_str02.loc[rew_id]]
        #         for i in range(len(t)):
        #             channel_network[rew_id].update(dt, upstream_volumetric_discharge=up_discharge[i], hillslope_volumetric_discharge=hillslope_volumetric_discharge[i] , ppt=ppt[i])
        #             volumetric_discharge[i]=channel_network[rew_id].volumetric_discharge

        #         network_volumetric_discharges[rew_id]=volumetric_discharge











if __name__ == '__main__': main()