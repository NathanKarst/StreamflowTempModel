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
source activate py2k_model
x=$(hpc_temperature_calibrate.py True 5000 elder_temperature.p elder.shp)
DATE=`date +%Y-%m-%d:%H:%M:%S`
TITLESTR="ELDER_TEMP_$DATE"
sendmail daviddralle@gmail.com << EOF
subject:$TITLESTR
$x
EOF
