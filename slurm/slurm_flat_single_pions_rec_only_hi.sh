#!/bin/sh

#SBATCH --image=eicweb/jug_xl:nightly
#SBATCH -A m3763
#SBATCH --qos=shared
#SBATCH --constraint=haswell
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-200
#SBATCH --mem=16GB

echo running job script with task ID $SLURM_ARRAY_TASK_ID
srun shifter bash run_flat_single_pions_rec_only_hi.sh $SLURM_ARRAY_TASK_ID

