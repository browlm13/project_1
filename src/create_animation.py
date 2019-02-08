#!/usr/bin/env python3

"""

					Shallow Water


	Read badly formated time step matrices
	from shallow water output into numpy arrays
	and then save 3d animation as .mp4 file

	L.J. Brown
	Parallel Scientific Computing 001C 1192 

"""

import numpy as np
import pandas as pd
import glob

#
# Settings
#

TRIAL = 3

HT_FOLDER_PATH = "../output_%s/timesteps/*.csv" % TRIAL
XS_FPATH = "../output_%s/xs.csv" % TRIAL
YS_FPATH = "../output_%s/ys.csv" % TRIAL


#
# Methods
#

def try_int(s):
	try: return int(s)
	except: return s

def num_key(path):
	""" return number from path's number.csv """
	s = path.split('/')[-1].split('.')[0]
	return try_int(s) 

def sort_timestep_files(path_list):
	path_list.sort(key=num_key)
	return path_list

def get_xs_and_ys(xs_fpath, ys_fpath):
	xs = pd.read_csv(xs_fpath)[:-1].astype(np.float).T.index.values
	ys = pd.read_csv(ys_fpath)[:-1].astype(np.float).T.index.values

	return xs[:-2].astype(np.float), ys[:-2].astype(np.float)

def get_H_and_t(frame_fpath):
	Ht_df = pd.read_csv(frame_fpath).astype(np.float)

	t = Ht_df.iloc[0,-1].astype(np.float)
	#print(Ht_df.iloc[:,-1])
	H_inv_df = Ht_df.iloc[:,:-2]
	#H_inv_df = Ht_df.iloc[:,1:-1]
	H = H_inv_df.astype(np.float).values.T

	return H, t

def get_tXYH(sorted_frame_fpath, xs_fpath, ys_fpath):

	xs, ys = get_xs_and_ys(xs_fpath, ys_fpath)
	X, Y = np.meshgrid(xs, ys)
	
	nt = len(sorted_frame_fpath)
	ts = np.zeros(shape=(nt,))

	Hs = np.zeros(shape=(nt,xs.shape[0],ys.shape[0]))

	for i, fpath in enumerate(sorted_frame_fpath):
		H, t = get_H_and_t(fpath)
		ts[i] = t
		Hs[i,:,:] = H

	return ts, X, Y, Hs




if __name__ == "__main__":

	# get frame file names
	frame_fpaths = glob.glob(HT_FOLDER_PATH)

	# sort frames by frame number or by time iterator number 
	sorted_frame_fpaths = sort_timestep_files(frame_fpaths)

	# read in data
	ts, X, Y, H_history = get_tXYH(sorted_frame_fpaths, XS_FPATH, YS_FPATH)

	#
	# Animate
	#

	import matplotlib.pyplot as plt
	import mpl_toolkits.mplot3d.axes3d as p3
	import matplotlib.animation as animation

	#
	# settings
	#
	OUTPUT_ANIMATION_FNAME = str('../output_%s/shallow_water_animation_%s.mp4' % (TRIAL, TRIAL))
	
	# display animation file name
	print(OUTPUT_ANIMATION_FNAME)

	FRAMES = H_history.shape[0] - 1

	#
	# Create animation
	#

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')


	counter = 0
	def animate(i):
	    
	    global counter
	    
	    # update plot
	    ax.clear()
	 
	    if counter % 10 ==0:
	    	print(counter)

	    ax.set_title('Shallow Water Wave');
	    ax.plot_surface(X, Y, H_history[counter,:,:])
	    ax.set_xlabel('x')
	    ax.set_ylabel('y')
	    ax.set_zlabel('height(time)');
	    
	    # update counter
	    counter += 1

	# create and write animation
	shallow_water_animation = animation.FuncAnimation( fig, animate )
	shallow_water_animation.save(OUTPUT_ANIMATION_FNAME, writer='imagemagick', fps=60)

