#!/bin/bash
#SBATCH --account=rz.rz          # Your account
#SBATCH --partition=smp      # Slurm Partition; Default: smp/mpp
#SBATCH --time=0:30:00                # time limit for job; Default: 0:30:00
#SBATCH --qos=30min                  # Slurm QOS; Default: 30min
#SBATCH --nodes=1             # Number of nodes
##SBATCH --ntasks=1            # Number of tasks (MPI) tasks to be launched
##SBATCH --mem=30000               # If more than the default memory is needed;
                                     # Default: <#Cores> * <mem per node>/<cores per node>
#SBATCH --job-name=wrfrun         # Job name
#SBATCH --output=result_%A_%a.out  # File where the standard output is written to(*)
#SBATCH --error=result_%A_%a.err   # File where the error messages are written to(*)

echo "SLURM_JOBID:         $SLURM_JOBID"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "SLURM_ARRAY_JOB_ID:  $SLURM_ARRAY_JOB_ID"

cd /albedo/work/user/chshao001/WRFPD/test/em_tropical_cyclone
srun ./ideal.exe -dim_ens 1 >& ideal.log

