#!/bin/bash
#SBATCH -J coverage_leave_one_out
#SBATCH -o coverage_leave_one_out.out
#SBATCH -e coverage_leave_one_out.err
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -t 0-05:30
#SBATCH -p huce_cascade
#SBATCH --mem=64G
#SBATCH --mail-type=END
#SBATCH --mail-user=amdur@g.harvard.edu

module load matlab/R2018b-fasrc01
srun -c $SLURM_CPUS_PER_TASK matlab -nosplash -nodesktop -r "posterior_estimate_cmip_cluster('model_predictions_22_11_14.mat')"
