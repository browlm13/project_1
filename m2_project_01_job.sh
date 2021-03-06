#!/bin/bash

#SBATCH -J shallow_water    # job name to display in squeue
#SBATCH -o output-%j.txt    # standard output file
#SBATCH -e error-%j.txt     # standard error file
#SBATCH -D /users/lmbrown/scratch/parallel_scientific_computing/project_1

#SBATCH -p medium-mem-1-m   # requested partition
#SBATCH --mem=256G          # memory requirment
#SBATCH -t 00:20:00         # kill time

#SBATCH --mail-user lmbrown@smu.edu
#SBATCH --mail-type=all
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

cd src
./sw.x