import os

import ete3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from tqdm.auto import tqdm


from cassiopeia.solver import solver_utilities
######

import sys
sys.path.append('/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/KPTracer-release-main/reproducibility/Figure2_S2/')
import clonal_expansions


tree = ete3.Tree(f'/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/ReconstructingTrees/tree/M_tree.processed2.nwk', 1)
print(1)
#print(tree)
# set node names
node_iter = 0
for n in tree.traverse():
    if not n.is_leaf():
        n.name = f'node-{node_iter}'
        node_iter += 1

print(2)
##print(tree.traverse())
#tree, expansions = clonal_expansions.detect_expansion(tree, pval=0.015, min_depth=1, _first=False, min_clade_prop = 0.15)
#print(expansions)


tree, expansions = clonal_expansions.detect_expansion(tree,min_depth=1,_first=False, min_clade_prop = 0.15)
print(expansions)
#