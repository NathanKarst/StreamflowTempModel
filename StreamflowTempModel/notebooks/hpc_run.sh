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
# Wall clock limit (30 seconds here):
#SBATCH --time=24:00:00
#
## Command(s) to run:
source activate py2k_model
python hpc_hillslope_calibrate.py False 50000 sf_leggett_runoff.p sf_leggett.shp > output.txt
