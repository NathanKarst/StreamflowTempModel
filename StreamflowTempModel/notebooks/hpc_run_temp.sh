#!/bin/bash
# Job name:
#SBATCH --job-name=test
#
# Account:
#SBATCH --account=fc_hydrology
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
## Command(s) to run:
git checkout master
source activate py3k
x=$(python hpc_temperature_calibrate.py False 1000 elder_temperature.p elder.shp)
DATE=`date +%Y-%m-%d_%H%M`
TITLESTR="ELDER_TEMP_$DATE"
TXT='.txt'
SAVETXT=$TITLESTR$TXT
echo $x > $SAVETXT