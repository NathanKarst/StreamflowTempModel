import os
import numpy as np
from os.path import dirname
import glob
import pickle

#must import all hillslope modules to set models for each group
import sys
parent_dir = dirname(dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','2_hillslope_discharge'))
from vadoseZone import LaioVadoseZone, PorporatoVadoseZone
from groundwaterZone import GroundwaterZone, NonlinearReservoir, NonlinearReservoir, TwoLinearReservoir, TwoParallelLinearReservoir


def main():
    """ Write REW groups and channel parameter files. 
    
    This function writes three files to the model_data directory:
        - group_params.p: parameters for each REW group, as specified in model_data/rew_config.p under the 'group' key. All REWs in a given group will have identical parameters. 
        - param_ranges.p: ranges for some subset of parameters specified for each REW. This file can be used for model calibration. Parameters whose ranges are not specified are assumed to be constant.
        - channel_params.p: parameters for each REW's channel model. 
    Args:
        - None
        
    Returns: 
        - None. 
    """


    parent_dir = dirname(dirname(os.getcwd()))
    rew_config = pickle.load( open( os.path.join(parent_dir,'model_data','rew_config.p'), "rb" ) )
    

    #NONLINEAR GROUNDWATER RESERVOIR, PORPORATO VADOSE ZONE
    group_params = {i:{'ET':0, 'emax':0.5, 'leakage':0, 'n':0.29, 'sfc':0.51, 'storage':0, 'sw':0.30, 'zr':76.4, 
                    'discharge':0, 'groundwater':0, 'a':0.0001064, 'b':3, 'vz':PorporatoVadoseZone, 'gz':NonlinearReservoir} for i in set(rew_config['group'])}

    param_ranges = {i:{'zr':(50,300), 'sw':(.1,.4), 'a':(.0001,.01), 'sfc':(.4,.9),'n':(.2,.8)} for i in set(rew_config['group'])}


    # SINGLE LINEAR RESERVOIR MODEL, PORPORATO VADOSE ZONE
    # group_params = {i:{'ET':.1, 'emax':0.5, 'leakage':0.1, 'n':0.43, 'sfc':0.5, 'storage':0.1, 'sw':0.19, 'zr':46, 
    #                 'discharge':0.1, 'groundwater':0.1, 'k':0.5, 'groundwaterMax':200} for i in set(rew_config['group'])}
    # param_ranges = {i:{'zr':(0,200), 'sw':(.1,.4), 'k':(.1,2), 'sfc':(.4,.9),'n':(.2,.8)} for i in set(rew_config['group'])}


    #TWO RESERVOIR GROUNDWATER, PORPORATO VADOSE ZONE
    # group_params = {i:{'ET':0, 'emax':0.5, 'leakage':0.0, 'n':0.29, 'sfc':0.4, 'storage':0.0, 'sw':0.18, 'zr':95, 
    #                 'discharge':0.0, 'groundwater':0.0, 'k1':0.35, 'k2':0.038, 'k12':.99, 'res1':0.0, 'res2':0.0} for i in set(rew_config['group'])}
    # param_ranges = {i:{'zr':(50,300), 'sw':(.1,.4), 'k1':(.1,1), 'k2':(0.01,1), 'k12':(.1,1), 'sfc':(.4,.9),'n':(.2,.8)} for i in set(rew_config['group'])}


    #TWO PARALLEL LINEAR RESERVOIR, PORPORATOR VADOSE ZONE
    # group_params = {i:{'ET':0, 'emax':0.5, 'res1':0, 'res2':0, 'leakage':0, 'n':0.29, 'sfc':0.51, 'storage':0, 'sw':0.30, 'zr':76.4, 
    #                 'discharge':0, 'groundwater':0, 'k1':0.5, 'k2':0.05, 'f1':0.5, 'vz':PorporatoVadoseZone, 'gz':TwoParallelLinearReservoir} for i in set(rew_config['group'])}
    # param_ranges = {i:{'zr':(50,300), 'sw':(.1,.4), 'k1':(.1,1), 'sfc':(.4,.9),'n':(.2,.8), 'k2':(0.01,0.4), 'f1':(.01,.99)} for i in set(rew_config['group'])}


    channel_params = {i:{'mannings_n':0.03, 'e':0.01, 'f':0.39, 'volume':0} for i in set(rew_config.index)}


    pickle.dump( group_params, open( os.path.join(parent_dir,'model_data','group_params.p'), "wb" ) )
    pickle.dump( channel_params, open( os.path.join(parent_dir,'model_data','channel_params.p'), "wb" ) )
    pickle.dump( param_ranges, open( os.path.join(parent_dir,'model_data','param_ranges.p'), "wb" ) )

    
    
if __name__ == '__main__': 
    main()