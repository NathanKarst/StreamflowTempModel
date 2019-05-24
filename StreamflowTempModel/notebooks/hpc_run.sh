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
x=$(python hpc_hillslope_calibrate.py False 10000 dry_runoff.p dry.shp)
DATE=`date +%Y-%m-%d_%H%M`
TITLESTR="DRY_$DATE"
TXT='.txt'
SAVETXT=$TITLESTR$TXT
echo $x > $SAVETXT
