import os, os.path
import pandas as pd
import numpy
import support_functions as sf


###################################
#    START WITH INITIALIZATION    #
###################################

##  SOME IMPORTANT DIRECTORIES

#set the working directory
dir_python = os.path.dirname(os.path.realpath(__file__))
#master directory
dir_proj = os.path.dirname(dir_python)
#get data directory
dir_data = os.path.join(os.path.dirname(dir_proj), "data")

dir_out = os.path.join(dir_proj, "out")
if not os.path.exists(dir_out):
    os.makedirs(dir_out, exist_ok = True)
dir_ref = os.path.join(dir_proj, "ref")
dir_datasets = os.path.join(dir_ref, "datasets")

# set some vals
dataset = "sead"


##  OUTPUT FILES

fp_csv_transition = os.path.join(dir_datasets, dataset, "transition_matrices_by_scenario.csv")
fp_csv_attribute_state = os.path.join(dir_datasets, dataset, "attribute_state.csv")
fp_csv_attribute_scenario = os.path.join(dir_datasets, dataset, "attribute_scenario.csv")
