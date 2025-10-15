import numpy as np
import pandas as pd
import cassiopeia as cas
import matplotlib.pyplot as plt
import seaborn as sns
data_directory = "/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data"
tumor_to_fitness = pd.DataFrame(columns = ['mean_fitness', 'tumor'])
import os
import glob

script_dir = "/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/fitnesses"

file_ext = "mean_fitness.*.txt"
# 使用 glob 模块获取指定格式的文件路径列表
file_paths = glob.glob(os.path.join(script_dir, file_ext))
target_tumors = []
for file_path in file_paths:
    print(file_path)
    file_name = os.path.basename(file_path)
    
    parts = file_name.split('.')
    name = parts[1]
    target_tumors.append(name) 
    print(name)

print (target_tumors)



#target_tumors = ['MI_MGH_0', 'MI_MGH_1', 'MI_MGH_2', 'MI_MGH_3','Met_MGH_0', 'Met_MGH_1', 'Met_MGH_2', 'Met_MGH_3','NMI_MGH_0','NMI_MGH_2','NMI_MGH_3','NMI_MGH_4','NMI_MGH_5','NMI_MGH_6']

for k, tumor in zip(range(len(target_tumors)), target_tumors):
    print(1)
    fitness_fp = f"{data_directory}/fitnesses/mean_fitness.{tumor}.txt"
    fitness_df = pd.read_csv(fitness_fp, sep='\t', index_col = 0)
    print(fitness_df)
    
    _mi, _ma = fitness_df['mean_fitness'].min(), fitness_df['mean_fitness'].max()
    fitness_df['mean_fitness'] = fitness_df['mean_fitness'].apply(lambda x: (x - _mi) / (_ma - _mi))
    
    fitness_df['tumor'] = tumor
    tumor_to_fitness = pd.concat([tumor_to_fitness, fitness_df])
    
h = plt.figure(figsize = (10, 8))
ax = plt.gca()
sns.violinplot(x = 'tumor', y='mean_fitness', data=tumor_to_fitness, ax=ax, cut=0, inner=None)
plt.xticks(rotation=90)
plt.tight_layout()
plt.title("Fitness Distributions for Representative Tumors")
plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure2_S2/fitness.png")
plt.savefig("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/plot/Figure2_S2/fitness.pdf")
plt.show()
print(1)
