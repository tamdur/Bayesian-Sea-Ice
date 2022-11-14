#!/bin/bash
#SBATCH -J make_priors
#SBATCH -o make_priors.out
#SBATCH -e make_priors.err
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -t 0-05:30
#SBATCH -p huce_cascade
#SBATCH --mem=64G
#SBATCH --mail-type=END
#SBATCH --mail-user=amdur@g.harvard.edu

module load matlab/R2018b-fasrc01
srun -c $SLURM_CPUS_PER_TASK matlab -nosplash -nodesktop -r "sigmoid_params_cluster('model_full_priors_22_11_14.mat')"
