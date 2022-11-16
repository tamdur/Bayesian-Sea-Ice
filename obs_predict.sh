#!/bin/bash
#SBATCH -J obs_predict
#SBATCH -o obs_predict.out
#SBATCH -e obs_predict.err
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0-10:30
#SBATCH -p huce_cascade
#SBATCH --mem=4G
#SBATCH --mail-type=END
#SBATCH --mail-user=amdur@g.harvard.edu

module load matlab/R2018b-fasrc01
srun -c $SLURM_CPUS_PER_TASK matlab -nosplash -nodesktop -r "obs_predict(100000,1000,25,'obs_pred_22_11_16long.mat');"
