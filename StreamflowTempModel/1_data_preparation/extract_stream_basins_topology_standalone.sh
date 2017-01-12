#Make sure the following environement variables are set before attempting to re-project the dem
#these lines of code can be added to .bash_profile in your home folder. 
THRESH=$1
MODEL=$2


export PATH=/Library/Frameworks/GDAL.framework/Programs:$PATH
export PROJSO=/Library/Frameworks/PROJ.framework/PROJ

gdalinfo ../../raw_data/ppt/2012/PRISM_ppt_stable_4kmD2_20120101_bil.bil >/dev/null 2>&1

# Export again, otherwise GDAL won't be able to access projections list. 
# Necessary environment variables must be exported every new cell 
# because the %%bash command restarts a new bash kernel every time
export PATH=/Library/Frameworks/GDAL.framework/Programs:$PATH
export PROJSO=/Library/Frameworks/PROJ.framework/PROJ

#define the path of the full tile DEM, as well as the clipping polygon
export MASTER_DEM=../../raw_data/dem/full_tile.tif
export WATERSHED_POLY=../../raw_data/watershed_poly/sf_miranda.shp
rm ../../raw_data/dem/dem.tif
gdalwarp -cutline $WATERSHED_POLY -crop_to_cutline $MASTER_DEM ../../raw_data/dem/dem.tif >/dev/null 2>&1

# re-project using the coordinate system from the `gdalinfo` command
gdalwarp -t_srs "GEOGCS["NAD83",
    DATUM["North_American_Datum_1983",
        SPHEROID["GRS_1980",6378137.0,298.257222101]],
    PRIMEM["Greenwich",0.0],
    UNIT["Degree",0.0174532925199433]]" ../../raw_data/dem/dem.tif ../../raw_data/dem/dem.tif >/dev/null 2>&1
    

# make the script executable for GRASS 
chmod u+x ../1_data_preparation/extract_stream_basins_topology.sh

# create folder to store temporary grass location
DBASE=../../grassdata
mkdir $DBASE

# setup GRASS setup file to find database. This step is necessary when using GRASS from the command line. 
echo "LOCATION_NAME: 
GISDBASE:$DBASE
MAPSET: 
GUI: 
PID:" > $HOME/.grass7/rc

# create new temporary location for the job, exit after creation of this location
/Applications/GRASS-7.0.app/Contents/MacOS/grass.sh -text -c -e ../../raw_data/dem/dem.tif ../../grassdata/mylocation >/dev/null 2>&1


# define the job file as environmental variable
export GRASS_BATCH_JOB="../1_data_preparation/extract_stream_basins_topology.sh"

# define channel head formation threshold; make it rather large to decrease
# computational demands later
export THRESH=$THRESH

# must set path to model folder for .sh script so it can find necessary files and store output
export MODEL=$MODEL

# set name of polygon defining watershed
# must be the same as the clipping polygon above
export WATERSHED_POLY=../../raw_data/watershed_poly/sf_miranda.shp
 
# now we can use the temp location and run the job defined via GRASS_BATCH_JOB
/Applications/GRASS-7.0.app/Contents/MacOS/grass.sh ../../grassdata/mylocation/PERMANENT >/dev/null 2>&1
 
# CLEANUP
unset GRASS_BATCH_JOB
 
# delete temporary location
rm -rf ../../grassdata/

# clear out rc file so that the GRASS GUI opens properly on next opening
echo "LOCATION_NAME: 
GISDBASE:
MAPSET: 
GUI: 
PID:" > $HOME/.grass7/rc