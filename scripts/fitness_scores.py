import sys
import os
import time

import numpy as np
import scipy
import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns

import ete3

sys.path.append("/home/yanzeqin/software/jungle/") # specify path to jungle
import jungle as jg
infiles = ["/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/ReconstructingTrees/tree/M_tree.processed2.nwk"]
F_empirical = jg.Forest.from_newick(infiles)
print (len(F_empirical), "trees")
############################
F_empirical.annotate_standard_node_features()
F_empirical.infer_fitness(params={})
#########################################################
node_features = F_empirical.node_features()
print (node_features.shape)
print(node_features.head())