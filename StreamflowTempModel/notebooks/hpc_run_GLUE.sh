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
x=$(python hpc_hillslope_calibrate.py False 5 elder_runoff.p elder.shp 2)
DATE=`date +%Y-%m-%d:%H:%M:%S`
TITLESTR="ELDER_$DATE"
echo 'Best fit params attached' | mail -s 'GLUE -- $TITLESTR, $DATE' -a best_params_list.p -a best_objs_list.p daviddralle@gmail.com
