from vadoseZone import LaioVadoseZone, PorporatoVadoseZone
from groundwaterZone import GroundwaterZone, NonlinearReservoir, NonlinearReservoir, TwoLinearReservoir, TwoParallelLinearReservoir
from REW import REW
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import os
import pickle
from datetime import date
import pandas as pd
import numpy as np
import spotpy
import time
import sys
calibrate = False
  
def main(): 

    #load config files, forcing file, and paramters for each group
    parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
    if (not calibrate):
        sys.path.append(os.path.join(parent_dir, 'StreamflowTempModel', '1_data_preparation'))
        from prep_rew_params import main as pparam 
        pparam()

    rew_config = pickle.load( open( os.path.join(parent_dir,'model_data','rew_config.p'), "rb" ) )
    rew_forcing_dict = pickle.load( open( os.path.join(parent_dir,'model_data','rew_forcing_dict.p'), "rb" ) )
    group_params = pickle.load( open( os.path.join(parent_dir,'model_data','group_params.p'), "rb" ))
    model_config = pickle.load( open( os.path.join(parent_dir, 'model_data', 'model_config.p'), 'rb'))


    #start/stop dates for running model  
    #spinup date is the date after start_date for which we assume model is finished spinning up         
    start_date = model_config['start_date']
    stop_date = model_config['stop_date']
    spinup_date = model_config['spinup_date']
    Tmax = model_config['Tmax']
    dt = model_config['dt_hillslope']
    t = model_config['t_hillslope']
    resample_freq_hillslope = model_config['resample_freq_hillslope']
    timestamps_hillslope = model_config['timestamps_hillslope']

    #for each group, create and solve REW with all corresponding output variables. 
    solved_group_hillslopes_dict = {}
    for group_id in group_params.keys():
        #get an REW id corresponding to the group so that we can use its forcing data
        #perhaps sometime down the line we take an average over the group members or something
        #or just ensure that each member of a group has same forcing as other members of same group
        #for now, just take the first one with that group id 
        rew_id = rew_config.loc[rew_config.group==group_id].index[0]
        rew = REW(group_params[group_id]['vz'], group_params[group_id]['gz'],  **{'pet':rew_forcing_dict[rew_id].pet, 'ppt':rew_forcing_dict[rew_id].ppt, 'aspect':90})
        
        storage    = np.zeros(np.size(t))
        groundwater     = np.zeros(np.size(t))
        discharge       = np.zeros(np.size(t))
        leakage         = np.zeros(np.size(t))
        ET              = np.zeros(np.size(t))

        #parameterize vz and gz
        for param in [x for x in rew.vz.__dict__.keys() if (x!='rew')]:
            setattr(rew.vz, param, group_params[group_id][param])
        for param in [x for x in rew.gz.__dict__.keys() if (x!='rew')]:
            setattr(rew.gz, param, group_params[group_id][param])

        #resample pet and ppt to integration timestep
        ppt = np.array(rew.ppt[start_date:stop_date].resample(resample_freq_hillslope).ffill())
        pet = np.array(rew.pet[start_date:stop_date].resample(resample_freq_hillslope).ffill())
        
        #Solve REW hillslope
        for i in range(len(t)):
            rew.vz.update(dt,**{'ppt':ppt[i],'pet':pet[i]})
            storage[i] = rew.vz.storage
            leakage[i]      = rew.vz.leakage
            ET[i]           = rew.vz.ET   
            rew.gz.update(dt,**{'leakage':leakage[i]})
            groundwater[i] = rew.gz.groundwater
            discharge[i] = rew.gz.discharge
    
        #Save all results as daily data. 
        solved_hillslope = pd.DataFrame({'storage':storage, 'leakage':leakage, 'ET':ET, 'groundwater':groundwater, 'discharge':discharge}, index=timestamps_hillslope)
        solved_group_hillslopes_dict[group_id] = solved_hillslope.resample('D').mean()

    pickle.dump( solved_group_hillslopes_dict, open( os.path.join(parent_dir,'model_data','solved_group_hillslopes_dict.p'), "wb" ) )



    #for calibration, only return discharge data after spinup period
    if calibrate:
        return discharge[-len(pd.date_range(spinup_date, stop_date, freq=resample_freq_hillslope)):]
    else:
        # YtdWater       = np.cumsum(ppt*dt)
        # YtdET          = np.cumsum(ET*dt)
        # YtdDischarge   = np.cumsum(discharge*dt)
        # plt.plot(t,storage/YtdWater,label='Soil Storage')
        # plt.plot(t,groundwater/YtdWater,label='Groundwater')    
        # plt.plot(t,YtdET/YtdWater,label='Cumulative ET')
        # plt.plot(t,YtdDischarge/YtdWater,label='Cumulative Discharge')
        # plt.plot(t,(storage+YtdET+groundwater+YtdDischarge)/YtdWater,label='Total')
        # plt.ylabel('Portion of Water Balance []')
        # plt.xlabel('Time [d]')
        # plt.legend()
        # plt.ylim(0,1.25)
        # plt.show()

        discharge = np.array(solved_group_hillslopes_dict[0].discharge)
        sns.set_style('dark')
        elder_df = pickle.load( open( os.path.join(parent_dir, 'calibration_data', 'elder_2010_2015.p'), 'rb'))
        outflows_df = pd.DataFrame({ 'data':np.array(elder_df['runoff'][start_date:stop_date]),'modeled':discharge}, index = elder_df[start_date:stop_date].index)
        ax = outflows_df.plot()
        fig = ax.get_figure()
        fig.axes[0].set_ylabel('q [cm/day]')
        ax.grid(True)
        plt.show()
        # ax = outflows_df.plot(secondary_y=['ppt'])
        # fig = ax.get_figure()
        # fig.axes[1].invert_yaxis()
        # fig.axes[0].set_ylabel('q [cm/day]')
        # fig.axes[1].set_ylabel('ppt [cm/day]')
        # ax.grid(True)
        # plt.show()
        # plt.tight_layout()
        # plt.savefig('~/Desktop/test.pdf')

        # plt.figure()
        # plt.plot(t,YtdET)
        # plt.ylabel('cumulative ET [cm]')
        # plt.xlabel('Time [d]')
        # nse = spotpy.objectivefunctions.nashsutcliff(elder_runoff[-len(pd.date_range(spinup_date, stop_date, freq=resample_freq)):],discharge[-len(pd.date_range(spinup_date, stop_date, freq=resample_freq)):])
        # print nse
        # plt.show()

    

def memoryAddress(input): return hex(id(input))

if __name__ == '__main__': main()



        



