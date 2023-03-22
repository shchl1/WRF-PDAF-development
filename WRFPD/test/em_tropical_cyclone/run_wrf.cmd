#!/bin/bash 
#SBATCH --account=rz.rz          # Your account
#SBATCH --time 3:00:00
#SBATCH --qos=12h
#SBATCH --partition=mpp
#SBATCH --ntasks=64          # 64cpu multi ens
##SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread   # disable hyperthreading
#SBATCH --job-name=mpiwrfrun
#SBATCH --output=out_%x.%j
#SBATCH --error=err_%x.%j

## Uncomment the following line to enlarge the stacksize if needed,
##  e.g., if your code crashes with a spurious segmentation fault.
ulimit -s unlimited

# To be on the safe side, we emphasize that it is pure MPI, no OpenMP threads
# export OMP_NUM_THREADS=1

# And we tell srun to bind each MPI task to a core. 
#Xsrun  I know what I am doing 
#srun -l --cpu_bind=cores ./wrf.exe >& wrf.log
srun ./wrf.exe >& wrf.log