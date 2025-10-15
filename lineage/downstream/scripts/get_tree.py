#coding:utf-8

import cassiopeia
from IPython.display import Image
import numpy as np
import pandas as pd
import cassiopeia as cas
import matplotlib.pyplot as plt
from cassiopeia.data import CassiopeiaTree
import os
import glob

# 获取当前脚本所在目录
#script_dir = os.path.dirname(os.path.abspath(__file__))
# 或者指定目录路径和文件格式
script_dir = '/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele'
# 指定文件格式
file_ext = '*.alleleTable.csv'
# 使用 glob 模块获取指定格式的文件路径列表
file_paths = glob.glob(os.path.join(script_dir, file_ext))
dtype={'intBC':str}
for file_path in file_paths:
    print(file_path)
    file_name = os.path.basename(file_path)
    print (file_name)
    parts = file_name.split('.')
    name = parts[0]
    print(name)
    allele_table = pd.read_csv(file_path,sep='\t',usecols = ['sampleID','cellBC', 'intBC', 'r1', 'r2', 'r3', 'allele', 'lineageGrp', 'readCount', 'UMI'],dtype=dtype)
    allele_table = allele_table.fillna("NC")#去掉NC的
    #计算indel priors
    indel_priors = cas.pp.compute_empirical_indel_priors(allele_table, grouping_variables=['intBC', 'lineageGrp'])
    indel_priors.to_csv(script_dir +"/data/" + name+".indel.csv",sep='\t',index=True,header=True)
    #############################
    CLONES = allele_table['lineageGrp'].drop_duplicates()
    for CLONE in CLONES:
        print(CLONE)
        CLONE = int(CLONE)
        clone_allele_table = allele_table[allele_table['lineageGrp'] == CLONE]
        clone_allele_table.to_csv(script_dir +"/data/" + name + "_" +str(CLONE)+".allele.csv",sep='\t',index=True,header=True)
        n_cells = clone_allele_table['cellBC'].nunique()
        a_barcode=clone_allele_table['cellBC'].drop_duplicates()
        a_barcode.to_csv(script_dir +"/data/" + name + "_"+str(CLONE)+".barcode.csv",sep='\t',index=True,header=True)
        n_intbc = clone_allele_table['intBC'].nunique()
        print(f"{name} Clonal population #{CLONE} has {n_cells} cells and {n_intbc} intBCs ({n_intbc * 3}) characters.")
        #下面函数返回三个变量
        #N x M字符矩阵
        #indel比例字典，指定每个字符到特定状态的突变概率
        #一个字典
        character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(clone_allele_table,mutation_priors = indel_priors)
        character_matrix.to_csv(script_dir +"/data/" + name + "_" +str(CLONE)+"_character.csv",sep='\t',index=True,header=True)
        priors[0]
        print(character_matrix.head(5))
        #使用indel priors和矩阵
        cas_tree = cas.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors,tree=None)
        print(cas_tree.character_matrix.head(5))
        cas_tree.character_matrix.to_csv(script_dir +"/trees/" + name + "_"+str(CLONE)+"_character_matrix.txt",sep='\t',index=True,header=True)
        ###########################################################图一
        print(f"We simulated {character_matrix.shape[1]} characters across {character_matrix.shape[0]} cells.")
        character_matrix.nunique(axis=0).plot(kind='bar')
        plt.xlabel("Character")
        plt.ylabel("Number of unique states ")
        plt.show()
        plt.savefig(script_dir +"/data/" + name + "_"+str(CLONE)+"character_unique.png")
        plt.savefig(script_dir +"/data/" + name + "_"+str(CLONE)+"character_unique.pdf")
        missing_data_per_character = character_matrix.apply(lambda x: len(x[x == -1]) / len(x), axis=0)
        missing_data_per_character.plot(kind='bar')
        plt.xlabel("Character")
        plt.ylabel("Percentage of missing states")
        plt.show()
        plt.savefig(script_dir +"/data/" + name + "_"+str(CLONE)+"character_missing.png")
        plt.savefig(script_dir +"/data/" + name + "_"+str(CLONE)+"character_missing.pdf")
        print(0)
        ##########################################################################################################

        ##################################################################
        cas_tree.n_cell
        cas_tree.n_character
        ############################
        cell_meta = clone_allele_table.groupby('cellBC').agg({"intBC": 'nunique', 'UMI': 'sum', 'sampleID': 'unique'})
        cell_meta['sampleID'] = [x[0] for x in cell_meta['sampleID']]

        missing_proportion = (character_matrix == -1).sum(axis=0) / character_matrix.shape[0]
        uncut_proportion = (character_matrix == 0).sum(axis=0) / character_matrix.shape[0]
        n_unique_states = character_matrix.apply(lambda x: len(np.unique(x[(x != 0) & (x != -1)])), axis=0)
        character_meta = pd.DataFrame([missing_proportion, uncut_proportion, n_unique_states], index = ['missing_prop', 'uncut_prop', 'n_unique_states']).T
        cas_tree.cell_meta = cell_meta
        cas_tree.character_meta = character_meta

        # create a basic vanilla greedy solver
        vanilla_greedy = cas.solver.VanillaGreedySolver()
        # reconstruct the tree
        vanilla_greedy.solve(cas_tree, collapse_mutationless_edges=True)
        #保存树
        a2=cassiopeia.data.to_newick(cas_tree._CassiopeiaTree__network)
        out_fp = script_dir +"/trees/" + name + "_"+str(CLONE)+"_tree.nwk"
        with open(out_fp, 'w') as f:
            f.write(a2)
        cas.pl.plot_matplotlib(cas_tree, orient='right', add_root=True,allele_table=clone_allele_table)
        plt.savefig(script_dir +"/data/" + name + "_"+str(CLONE)+"_tree.png")
        plt.savefig(script_dir +"/data/" + name + "_"+str(CLONE)+"_tree.pdf")
