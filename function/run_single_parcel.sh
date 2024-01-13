#!/bin/bash
#SBATCH -p node1
#SBATCH --ntasks=1 # Run a single serial task
#SBATCH --cpus-per-task=4
#SBATCH -e job.%j.log # Standard error
#SBATCH -o job.%j.pro_out.txt # Standard output
#SBATCH --job-name=gradient_pre
#SBATCH --mail-type=END
module load MATLAB/R2019a
matlab -singleCompThread -nodisplay -nosplash -r $1
