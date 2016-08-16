
from vadoseZone import LaioVadoseZone, PorporatoVadoseZone
from groundwaterZone import GroundwaterZone, LinearReservoir, TwoLinearReservoir
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
import sandboxDralle
parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))


class spot_setup(object):
    def __init__(self,group_params, param_ranges):
        self.params = []
        for group in param_ranges.keys():
            for param in param_ranges[group].keys():
                self.params.append(
                    spotpy.parameter.Uniform(
                        param, param_ranges[group][param][0], param_ranges[group][param][1]
                    )
                )
        

    def parameters(self):
        #in the future, change this to deal with multiple groups of parameters...right now, group = 0 is only group
        param_realization = spotpy.parameter.generate(self.params)
        return param_realization
                

    def simulation(self, vector):
        parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
        group_params = pickle.load( open( os.path.join(parent_dir,'model_data','group_params.p'), "rb" ))
        for i, param_object in enumerate(self.params):
            param_name = param_object.name
            param_value = vector[i]
            group_params[0][param_name] = param_value

        pickle.dump( group_params, open( os.path.join(parent_dir,'model_data','group_params.p'), "wb" ) )
        simulations=sandboxDralle.main()
        return simulations   

     
    def evaluation(self):
        model_config = pickle.load( open( os.path.join(parent_dir, 'model_data', 'model_config.p'), 'rb'))
        elder_df = pickle.load( open( os.path.join(parent_dir, 'calibration_data', 'elder_2010_2015.p'), 'rb'))
        spinup_date = model_config['spinup_date']
        stop_date = model_config['stop_date']
        elder_runoff = np.array(elder_df['runoff'][spinup_date:stop_date].resample(model_config['resample_freq']).ffill())
        return elder_runoff  
        
    def objectivefunction(self,simulation,evaluation):
        objectivefunction = -spotpy.objectivefunctions.rmse((simulation),(evaluation))
        return objectivefunction