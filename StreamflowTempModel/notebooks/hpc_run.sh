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
git checkout test_leggett
source activate py2k_model
x=$(python hpc_hillslope_calibrate.py False 40000 sf_leggett_runoff.p sf_leggett.shp)
DATE=`date +%Y-%m-%d:%H:%M:%S`
TITLESTR="LEGGETT_$DATE"
sendmail daviddralle@gmail.com << EOF
subject:$TITLESTR
$x
EOF


