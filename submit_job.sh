#!/bin/bash

#SBATCH -J shallow_water    # job name to display in squeue
#SBATCH -o slurn_output-%j.txt    # standard output file
#SBATCH -e slurn_error-%j.txt     # standard error file
#SBATCH -D /users/lmbrown/scratch/parallel_scientific_computing/project_1/output

#SBATCH -p htc				# requested partition
#SBATCH -t 180         		# kill time

#SBATCH --mail-user lmbrown@smu.edu
#SBATCH --mail-type=all
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

cd ../src
srun sw.x