#!/bin/bash
####################################################
# extract_stream_basins_topology.sh
# by David Dralle, 2016
# purpose: linux shell script to extract stream 
# network topology from DEM for ModCZ. This script 
# must be run with the corresponding stream network
# extraction Jupyter notebook in the /notebooks 
# subdirectory. 
####################################################

#clear out any existing data
rm $MODEL/raw_data/topology/basin.csv
rm $MODEL/raw_data/topology/topology.csv
rm $MODEL/raw_data/basins_poly/*
rm $MODEL/raw_data/streams_poly/*

M=dem
r.in.gdal input=$MODEL/raw_data/dem/dem.tif output=dem_unfilled
r.mapcalc "dem_unfilled_nulled = if(dem_unfilled ,dem_unfilled , null() , null() )"
g.region raster=dem_unfilled_nulled -p
g.region res=50 -ap
r.resamp.stats --overwrite input=dem_unfilled_nulled output=$M

g.region raster=dem
g.region -p > info.txt
nsres="$(grep "nsres" info.txt)"
nsres=${nsres:6}
ewres="$(grep "ewres" info.txt)"
ewres=${ewres:6}
rm info.txt

# convert m^2 threshold into #cells threshold
cell_area=$(echo "$nsres*$ewres" | bc)
THRESH=$(echo "$THRESHMETERS/$cell_area" | bc)
echo "Cell area is $cell_area"
echo "Cells threshold is $THRESH"

# get accumulation raster and extract stream network
ACCUMSTRING="accum_$THRESHMETERS"
DIRSTRING="dir_$THRESHMETERS"
STREAMSTRING="stream_$THRESHMETERS"
r.watershed -a --overwrite elevation=$M accumulation=$ACCUMSTRING
r.stream.extract --overwrite elevation=$M threshold=$THRESH stream_length=5 stream_raster=$STREAMSTRING stream_vector=stream_vector_temp direction=$DIRSTRING

#uncomment and install r extensions, which are not pre-installed with grass 7.0.x
#g.extension r.stream.order
#g.extension r.stream.basins

#get the maximum accumulation point and use as watershed outlet for rest of analysis
#r.watershed
# STATS=$(r.describe -r $ACCUMSTRING)
# searchstring="-"
# MAX=${STATS#*$searchstring}
# MAX=$(echo ${MAX%.*})
#create watershed outlet point, get east/north coords of outlet
# r.mapcalc --overwrite "newmap = if($ACCUMSTRING < $MAX, null(), $ACCUMSTRING)"
# r.out.xyz --overwrite input=newmap output="$MODEL/aux_data/watershed_outlet"
# EAST=$(tr '|' '\n' < $MODEL/aux_data/watershed_outlet | sed -n 1p)
# NORTH=$(tr '|' '\n' < $MODEL/aux_data/watershed_outlet | sed -n 2p)
#get watershed polygon corresponding to outlet
# WATERSTRING="watershed_$THRESH"
# r.water.outlet --overwrite input=$DIRSTRING output=$WATERSTRING coordinates=$EAST,$NORTH
#if necessary
# r.mask --overwrite raster=$WATERSTRING

# get topologic properties of stream network
STREAMVECT="stream_vect_$THRESHMETERS"
r.stream.order --overwrite accumulation=$ACCUMSTRING elevation=$M stream_rast=$STREAMSTRING direction=$DIRSTRING stream_vect=$STREAMVECT
v.out.ogr --overwrite -c input=$STREAMVECT type=line output="$MODEL/raw_data/streams_poly"

# streamvector to points, export
v.to.rast --overwrite input=$STREAMVECT output=streamrastertemp type=line use=val value=1
POINTSVECT="points_$THRESHMETERS"
r.to.vect --overwrite input=streamrastertemp output=$POINTSVECT type=point
v.out.ogr --overwrite -c input=$POINTSVECT type=point output="$MODEL/raw_data/streams_points"



# export table with stream topology info
db.out.ogr --overwrite input=$STREAMVECT output="$MODEL/raw_data/topology/topology.csv"

#get raster of basins corresponding to stream network
BASINSTRING="basins_$THRESHMETERS"
BASINVECT="basins_vect_$THRESHMETERS"
BASINVECTCLEAN="basins_vect_clean_$THRESHMETERS"
r.stream.basins --overwrite direction=$DIRSTRING stream_rast=$STREAMSTRING basins=$BASINSTRING

#convert basins raster into basins polygons, rebuild topology of vector layer
r.to.vect --overwrite -v input=$BASINSTRING output=$BASINVECT type=area
v.build --overwrite map=$BASINVECT option=build

#NOTE! NEEDS TO BE RE-WRITTEN SO THAT THRESH=MINIMUM MAPPING UNIT + EPSILON. 
#The point of v.clean is to get rid of dangly pieces of raster around the edge of the map
v.clean --overwrite input=$BASINVECT output=$BASINVECTCLEAN type=area tool=rmarea thresh=$cell_area
v.db.addcolumn map=$BASINVECTCLEAN columns="area_sqkm DOUBLE PRECISION"
v.to.db map=$BASINVECTCLEAN option=area columns=area_sqkm unit=k

#export basins polygon as shape file
v.out.ogr --overwrite input=$BASINVECTCLEAN type=area output="$MODEL/raw_data/basins_poly"
db.out.ogr --overwrite input=$BASINVECTCLEAN output="$MODEL/raw_data/topology/basin.csv"


