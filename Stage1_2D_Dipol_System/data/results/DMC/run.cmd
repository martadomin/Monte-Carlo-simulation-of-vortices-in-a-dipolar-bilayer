#!/bin/bash
#SBATCH --job-name=MC1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p all 

time ./MC

##### run ######
# sbatch run.cmd
#### check #####
# squeue
#### cancel ####
# scancel
#### list of finished jobs ###
# sacct