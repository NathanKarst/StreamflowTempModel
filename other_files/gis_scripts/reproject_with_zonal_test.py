# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import geopandas as gpd
import fiona 
import rasterio as rs
from zonal_stats import zonal_stats

driver = 'ESRI Shapefile'
raster_path = "/Users/daviddralle/Dropbox/research/streamflow_temp/data/elder_test/dem/imgn40w124_13.img"
shape_path = '/Users/daviddralle/Dropbox/research/streamflow_temp/data/elder_test/elder_polygon_prism/prism_elder_poly.shp'
elder_poly = gpd.read_file(shape_path)
src = rasterio.open(raster_path)


elder_poly = elder_poly.to_crs(src.crs)
elder_poly.to_file(shape_path,driver=driver)

stats = zonal_stats(shape_path,raster_path)






