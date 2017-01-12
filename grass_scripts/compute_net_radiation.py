#!/usr/bin/python
"""
Created on Wed Jun 29 20:10:09 2016
Script to compute net radiation for all days of the year for given DEM
Must be run from within GRASS BASH terminal 

arguments: 
-timestep in hours (float)
-dem (string)
-start day (int 1-365)
-stop day (int 1-365)



@author: daviddralle
"""

import sys
from grass.pygrass.modules import Module
import os

def main(argv):
    print 'computing '
    step=float(argv[0])
    elevin=argv[1]
    start_day=int(argv[2])
    stop_day=int(argv[3])
    parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
    out_folder=os.path.join(parent_dir,"raw_data","rn")
    
    #compute irradiance for each day
    for day in range(start_day,stop_day):
        rsun = Module("r.sun", elevation=elevin,step=step,glob_rad='rad.global.'+str(day),overwrite=True,run_=False, day=day)
        rsun.run()
        #r.out.gdal geotiff does not include datum; use AAIGrid format for output raster
        rout = Module("r.out.gdal",overwrite=True,f=True,input='rad.global.'+str(day),output=os.path.join(out_folder,'rn'+str(day)),format="AAIGrid")
        gremove = Module("g.remove",f=True,type='raster',name='rad.global.'+str(day))

if __name__ == "__main__":
    main(sys.argv[1:])




