


from matplotlib import pyplot as plt
import numpy as np
import os
import pickle
from datetime import date
import pandas as pd
import numpy as np
import time
import sys

## add hillslope discharge folder to path
parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, os.path.join(parent_dir,'StreamflowTempModel','2_hillslope_discharge'))

from vadoseZone import LaioVadoseZone, PorporatoVadoseZone
from groundwaterZone import GroundwaterZone, LinearReservoir
from REW import REW

def main():


    parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
    rew_config = pickle.load( open( os.path.join(parent_dir,'model_data','rew_config.p'), "rb" ) )
    rew_forcing_dict = pickle.load( open( os.path.join(parent_dir,'model_data','rew_forcing_dict.p'), "rb" ) )
    group_params = pickle.load( open( os.path.join(parent_dir,'model_data','group_params.p'), "rb" ))
    model_config = pickle.load( open( os.path.join(parent_dir, 'model_data', 'model_config.p'), 'rb'))

    solved_hillslope_discharge = pickle.load(open( os.path.join(parent_dir,'model_data','solved_group_hillslopes_dict.p'), "rb"))
    
    outlets = rew_config[rew_config['next_stream'] == -1].index.tolist()
    adj = dict((int(index),[int(row['prev_str01']),int(row['prev_str02'])]) for index,row in rew_config.iterrows())
    ordered = dfs(outlets,[0],adj)[1:] ## remove initial 0 from the list
    ordered = ordered[::-1] ## go from strahler 1 to strahler n, instead of other way round
    
    
    
def dfs(queue,haveSeen,adj):
    if not queue: return haveSeen
    current = queue.pop()
    haveSeen.append(current)
    queue += [next for next in adj[current] if next not in haveSeen]
    
    return dfs(queue,haveSeen,adj)
        



if __name__ == '__main__':  main()