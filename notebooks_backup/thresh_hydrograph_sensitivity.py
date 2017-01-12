import pandas as pd
import gdal
import os
import numpy as np
import fiona
import shapely
from shapely import geometry
from os.path import dirname
import glob
import sys
import pickle
from functools import partial
import geopandas as gp
from matplotlib import pyplot as plt
import seaborn as sns

parent_dir = dirname(dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','lib'))
import zonal_stats as zs
import meteolib as meteo
import evaplib as evap

sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','1_data_preparation'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','2_hillslope_discharge'))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','3_channel_routing'))

import prism_forcing as forcing
import prep
import hillslope_discharge
import channel_routing


def main():

    for i in range(5):
        thresh = 160000/2**i

        print('Threshold: %d'%thresh)

        print('\tIdentifying basins...')
        os.system('./../1_data_preparation/extract_stream_basins_topology_standalone.sh ' + str(thresh) + ' ../..')
        
        print('\tPreparing basic data structures...')
        prep.main()
        rew_config = pickle.load(open(os.path.join(parent_dir,'model_data','rew_config.p'), "rb" ))   
        print('\t\tTotal REWs identified: %d'%len(rew_config.keys()))
        
        print('\tExtracting forcing datasets...')
        forcing.main()
        
        print('\tComputing hillslope discharge...')
        hillslope_discharge.main()
        
        print('\tRouting discharge through channel network...')
        channel_routing.main()        
        routed = pickle.load(open( os.path.join(parent_dir,'model_data','solved_channel_routing.p'), "rb" ))   
        
        outletREWId = [rew_id for rew_id,value in rew_config.items() if value['next_stream']==-1][0]
#         outletGroup = rew_config[outletREWId]['group']
        print('\t\tOutlet REW ID found: %d'%outletREWId)
        
        plt.plot(routed[outletREWId]['volumetric_discharge'],label=str(thresh))
    
    plt.legend()
    plt.show()






def plotBasins(thresh=''):
    basins = glob.glob(os.path.join(parent_dir,'raw_data','basins_poly','*.shp'))[0]
    basins_shape = gp.GeoDataFrame.from_file(basins)
    basins_shape['coords'] = basins_shape['geometry'].apply(lambda x: x.representative_point().coords[:])
    basins_shape['coords'] = [coords[0] for coords in basins_shape['coords']]

    ax1 = basins_shape.plot()
#    for idx, row in basins_shape.iterrows():
#        print_str = 'REW' + str(row['cat'])
#        plt.annotate(s=print_str, xy=row['coords'],
#                     horizontalalignment='center',fontsize=15)

    streams = glob.glob(os.path.join(parent_dir,'raw_data','streams_poly','*.shp'))[0]
    streams_shape = gp.GeoDataFrame.from_file(streams)
    streams_shape.plot(ax=ax1,color='blue')

    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])
    ax1.patch.set_facecolor('white')
    
    plt.savefig('basins_'+thresh+'.png')


if  __name__ == '__main__': main()