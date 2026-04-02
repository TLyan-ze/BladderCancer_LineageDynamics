#coding:utf-8
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


import os
import glob


# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))

# 将当前工作目录设置为脚本所在目录
os.chdir(script_dir)


input_dir = script_dir +"/trees/"
output_dir = script_dir + "/fitnesses/"
file_ext = "*.nwk"
# 使用 glob 模块获取指定格式的文件路径列表
file_paths = glob.glob(os.path.join(input_dir, file_ext))
for file_path in file_paths:
    print(file_path)
    file_name = os.path.basename(file_path)
    print (file_name)
    parts = file_name.split('.')
    name = parts[0]
    print(name)
    infiles = [file_path]
    F_empirical = jg.Forest.from_newick(infiles)
    print (len(F_empirical), "trees")
    ############################
    F_empirical.annotate_standard_node_features()
    F_empirical.infer_fitness(params={})
    #########################################################
    node_features = F_empirical.node_features()
    node_features.to_csv(output_dir + "/" + name +".fitness.txt",sep='\t',index=True,header=True)


