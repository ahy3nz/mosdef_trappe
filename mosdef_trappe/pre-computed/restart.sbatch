#!/bin/bash -l
#SBATCH -p regular
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -A mp138
#SBATCH -t 4:00:00
#SBATCH -J gomc5_236
#SBATCH -o gomc5_236 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu
#SBATCH -L SCRATCH

export OMP_NUM_THREADS=16
cd /global/cscratch1/sd/ahy3nz/Trajectories/gomc5/236
source activate py36
GOMC_CPU_GEMC +p16 restart4.conf
