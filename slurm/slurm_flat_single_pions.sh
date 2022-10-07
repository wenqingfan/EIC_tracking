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
GEN_P_RANGE=2 # 0: low p, 1: mid p, 2: high p, 3: extra high p
GEN_ETA_RANGE=0 # 0: -3.5 to 3.5, 1: 3.5 to 4, 2: -4 to -3.5 
srun shifter bash run_flat_single_pions_sim.sh $SLURM_ARRAY_TASK_ID $GEN_P_RANGE $GEN_ETA_RANGE

