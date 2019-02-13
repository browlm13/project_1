#!/usr/bin/env python3

"""
	Python Controller Script

	-argument <output_directory> 
"""

# internal
import os
import sys
import subprocess

#
# Constants
#

OUTPUT_DIRECTORY_TEMPLATE = '{}/../output/{}/' 
OUTPUT_TIMESTEP_DIRECTORY_TEMPLATE = '{}/../output/{}/timesteps/' 
RELATIVE_OUTPUT_DIRECTORY_TEMPLATE = '../output/%s'

#
# Run
#

if __name__ == "__main__":

	# read output directory argument
	output_directory_name = sys.argv[1]

	BASE_DIR = os.path.dirname(os.path.realpath(__file__))

	# get output directory path
	output_directory = OUTPUT_DIRECTORY_TEMPLATE.format(BASE_DIR, output_directory_name)

	# create if it does not exist
	if not os.path.exists(output_directory):
		print("Creating Director: ", output_directory)
		os.makedirs(output_directory)
	else:
		print("Directory %s not created" % output_directory)

	# create additional directories if necissary
	timestep_directory = OUTPUT_TIMESTEP_DIRECTORY_TEMPLATE.format(BASE_DIR, output_directory_name)

	if not os.path.exists(timestep_directory):
		print("Creating Director: ", timestep_directory)
		os.makedirs(timestep_directory)

	else:
		print("Directory %s not created" % timestep_directory)

	# run c executable
	executable_arg = RELATIVE_OUTPUT_DIRECTORY_TEMPLATE % output_directory_name
	cmd_string = "srun sw.x  %s" % executable_arg

	print("running: ", cmd_string)
	os.system(cmd_string)

	# create animation after completion
	cmd_string = "python create_animation.py %s" %  output_directory_name
	print("running: ", cmd_string)
	os.system(cmd_string)

	#local = "./"
	#m2 = "srun "
	#subprocess.call('srun sw.x')
	#subprocess.call('./sw.x')