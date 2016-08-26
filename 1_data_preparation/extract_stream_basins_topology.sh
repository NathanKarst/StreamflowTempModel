#!/bin/sh
# Purpose: Script to extract REW network
# Inputs must be variables THRESH and MODEL
# THRESH is the threshold in pixels for channel head formation
# MODEL is the absolute path of the model directory
# These must be set as environment variables prior to execution of this script,
# e.g. export THRESH=20000

#clear out any existing data
rm $MODEL/raw_data/topology/basin.csv
rm $MODEL/raw_data/topology/topology.csv
rm $MODEL/raw_data/basins_poly/*
rm $MODEL/raw_data/streams_poly/*

M=dem
r.in.gdal input=$MODEL/raw_data/dem/dem.tif output=$M
g.region -p -a raster=$M

#run r.watershed to get accum, drainage direction, stream raster
#could change accumulation threshold here, or convergence, or
ACCUMSTRING="accum_$THRESH"
DIRSTRING="dir_$THRESH"
STREAMSTRING="stream_$THRESH"
r.watershed -a --overwrite elevation=$M threshold=$THRESH accumulation=$ACCUMSTRING drainage=$DIRSTRING stream=$STREAMSTRING convergence=5

#get the maximum accumulation point and use as watershed outlet for 
#r.watershed
STATS=$(r.describe -r $ACCUMSTRING)
searchstring="-"
MAX=${STATS#*$searchstring}


#create watershed outlet point, get east/north coords of outlet
r.mapcalc --overwrite "newmap = if($ACCUMSTRING < $MAX, null(), $ACCUMSTRING)"
r.out.xyz --overwrite input=newmap output="$MODEL/aux_data/watershed_outlet"
EAST=$(tr '|' '\n' < $MODEL/aux_data/watershed_outlet | sed -n 1p)
NORTH=$(tr '|' '\n' < $MODEL/aux_data/watershed_outlet | sed -n 2p)

#get watershed polygon corresponding to outlet
WATERSTRING="watershed_$THRESH"
r.water.outlet --overwrite input=$DIRSTRING output=$WATERSTRING coordinates=$EAST,$NORTH

#should maybe use r.thin here to thin out the stream network raster; would then overwrite old stream network raster
#r.thin --overwrite input=$STREAMSTRING output=$STREAMSTRING

#if necessary
#uncomment and install r extensions, which are not pre-installed with grass 7.0.x
#g.extension r.stream.order
#g.extension r.stream.basins

# Use watershed boundary as mask for remaining operations
# only want main primary watershed 
r.mask --overwrite raster=$WATERSTRING

#output stream order vector map
STREAMVECT="stream_vect_$THRESH"
r.stream.order --overwrite accumulation=$ACCUMSTRING elevation=$M stream_rast=$STREAMSTRING direction=$DIRSTRING stream_vect=$STREAMVECT
v.out.ogr --overwrite -c input=$STREAMVECT type=line output="$MODEL/raw_data/streams_poly"


#export table with stream topology
db.out.ogr --overwrite input=$STREAMVECT output="$MODEL/raw_data/topology/topology.csv"

#get raster of basins corresponding to stream network
BASINSTRING="basins_$THRESH"
BASINVECT="basins_vect_$THRESH"
BASINVECTCLEAN="basins_vect_clean_$THRESH"
r.stream.basins --overwrite direction=$DIRSTRING stream_rast=$STREAMSTRING basins=$BASINSTRING

#convert basins raster into basins polygons
#NOTE! IN THE FUTURE THE THRESHOLD NEEDS TO BE CHANGED TO BE A FUNCTION OF MAP RESOLUTION
r.to.vect --overwrite -v input=$BASINSTRING output=$BASINVECT type=area
v.clean --overwrite input=$BASINVECT output=$BASINVECTCLEAN type=area tool=rmarea thresh=500.0
v.db.addcolumn map=$BASINVECTCLEAN columns="area_sqkm DOUBLE PRECISION"
v.to.db map=$BASINVECTCLEAN option=area columns=area_sqkm unit=k

#export basins polygon as shape file
v.out.ogr --overwrite -c input=$BASINVECTCLEAN type=area output="$MODEL/raw_data/basins_poly"
db.out.ogr --overwrite input=$BASINVECTCLEAN output="$MODEL/raw_data/topology/basin.csv"


