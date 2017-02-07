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
python hpc_temperature_calibrate.py False 50000 elder_temperature.p elder.shp > output_elder_temperature.txt
