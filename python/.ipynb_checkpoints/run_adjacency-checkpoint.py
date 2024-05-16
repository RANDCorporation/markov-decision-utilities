import os, os.path
import time
import math
import re
import io
import csv
import numpy as np
from copy import deepcopy
import pandas as pd
import warnings
import shapefile as ps
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.random import rand
import setup_runs as sr


########################
#    INITIALIZATION    #
########################

###   SHAPEFILE READ IN

#set path
dir_dataset_shapefile = os.path.join(sr.dir_data, sr.dataset_shapefile)
#check for datasets available
sr.dataset_shapefile_avail = os.listdir(dir_adjacency_shapefile)
#get shapefile
fn_shapefile = [x for x in sr.dataset_shapefile_avail if (".shp" in x)]
fn_shapefile = fn_shapefile[0]
#initialize reader
sf = ps.Reader(sr.fp_adjacency_shapefile)
#read in some information
sfs = sf.shapes()



#####################################
#    RUN ADJACENCY SPECIFICATION    #
#####################################

t0 = time.time()

def check_bbox(b1, b2):
    #dimensions for box 1
    x11 = b1[0]
    x12 = b1[2]
    y11 = b1[1]
    y12 = b1[3]
    #dimensions for box2
    x21 = b2[0]
    x22 = b2[2]
    y21 = b2[1]
    y22 = b2[3]
    #check
    return (x12 >= x21) & (x22 >= x11) & (y12 >= y21) & (y22 >= y11)

#set number of shapes to evaluate
n = len(sfs)
#initialize adjacency list
list_adj = []

#initialize dynamic list of shapes, number of edges, and dynamic count of edges that are accounted for
dict_shapes = {}
dict_shape_edges = {}
dict_counts = {}
#initialize index
i = 0#131960#0
n2 = 10000
#iterate
while i < n:
    
    ######################################
    #    INITIALIZE KEYS TO LOOP OVER    #
    ######################################

	keys_loop = list(dict_shape_edges.keys())

	########################
	#    READ IN SHAPES    #
	########################

	#read in shape points
	cur_shape = sfs[i]
	#update dictionaries
	dict_shape_edges.update({i: (len(cur_shape.points) - 1)})
	#points
	p_i = cur_shape.points

	#loop over keys
	for k in keys_loop:
		if check_bbox(sfs[k].bbox, cur_shape.bbox):
			#get points to check on
			p_k = sfs[k].points
			#build edges (assume they are all oriented the same way, so one must be reversed)
			edges_i = [(p_i[len(p_i) - x], p_i[len(p_i) - x - 1]) for x in range(1, len(p_i))]
			edges_k = [(p_k[x], p_k[x + 1]) for x in range(0, len(p_k) - 1)]

			#get intersection (assume all edges are unique)
			cap = set(edges_i).intersection(edges_k)
			#check length of intersection
			if len(cap) > 0:
				#add to adjacency list
				list_adj = list_adj + [(k, i)]
				#get number of edges that are represented
				n_edges = len(cap)
				#loop over k, i to update counts
				for m in [i, k]:
					#add to counts
					if m in dict_counts.keys():
						#update
						dict_counts.update({m: (dict_counts[m] + n_edges)})
						#check to see if all edges are represented
						if dict_counts[m] == dict_shape_edges[m]:
							#clear from counts/shapes
							del dict_counts[m]
							#del dict_shapes[m]
							del dict_shape_edges[m]
					else:
						#update number of edges
						dict_counts.update({m: n_edges})
                        
	if (i%10000) == 0:
		print("Index " + str(i) + " complete. Length of dict_counts = " + str(len(dict_counts)))
		#next stage
		if (i%50000) == 0:
			print("Total time: " + str((time.time() - t0)/60) + " minutes")
    #move to next shape
	i += 1

test = pd.DataFrame(list_adj, columns = ["shape_i", "shape_j"])
test.to_csv(sr.fp_csv_adjacency_index, index = False)

t1 = time.time()

print("Total time: " + str((t1 - t0)/60) + " minutes")
