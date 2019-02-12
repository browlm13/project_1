#!/bin/bash

#
# Single argument to 'submit_job.sh' script located in root of project's directory structure
#        argument <output directory name>
#        * slurn output will be written to seperate file
#

OUTPUT_DIR=$1

#SBATCH -J project_1              # job name to display in squeue
#SBATCH -o slurn_output-%j.txt    # standard output file
#SBATCH -e slurn_error-%j.txt     # standard error file
#SBATCH -D /users/lmbrown/scratch/parallel_scientific_computing/project_1/output

#SBATCH -p htc				# requested partition
#SBATCH --mem=7G            # memory requirment
#SBATCH -t 180         		# kill time

#SBATCH --mail-user lmbrown@smu.edu
#SBATCH --mail-type=all

#
# Dropped into '$PROJECT_DIR/output/' directory by slurn
#

# Change directory to '$PROJECT_DIR/src/'
cd ../src

# Run python 'controller.py' script passing $OUTPUT_DIR as an argument, or executable directly 'srun sw.x' 
srun sw.x