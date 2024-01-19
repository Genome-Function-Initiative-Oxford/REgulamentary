import warnings
warnings.filterwarnings('ignore')

import sys, os
import pandas as pd
import seaborn  as sns
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from sklearn.preprocessing import MinMaxScaler
sns.set_style("white")


extra_track = sys.argv[1]
title  = sys.argv[2]
upper_threshold = sys.argv[3]
try:
    upper_threshold = int(upper_threshold)
except:
    upper_threshold = 'auto'


sort_union = pd.read_csv(sys.argv[4], sep='\t', header=None)
sort_union.columns = columns=['chr', 'start', 'end']


matrix = pd.read_csv(sys.argv[5], sep='\t', comment='#')
matrix.columns = matrix.columns[1:].tolist() + ['remove']

H3K4me1_c = ["H3K4me1.%i"%i for i in range(1,401)]
H3K4me3_c = ["H3K4me3.%i"%i for i in range(1,401)]
H3K27ac_c = ["H3K27ac.%i"%i for i in range(1,401)]
CTCF_c    = ["CTCF.%i"%i    for i in range(1,401)]
if extra_track != "none":
    PARAMS_c = ["%s.%i"%(sys.argv[9], i) for i in range(1,401)]

if extra_track != "none":
    new_columns = H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c + PARAMS_c
else:
    new_columns = H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c
matrix.columns = new_columns + ['remove']

matrix = matrix[new_columns]
matrix = matrix.fillna(0)

matrix.index = sort_union['chr'].astype(str)+':'+sort_union['start'].astype(str)+'-'+sort_union['end'].astype(str)


re = pd.read_csv(sys.argv[6], sep='\t')
re.index = re['chromosome'].astype(str)+':'+re['start'].astype(str)+'-'+re['end'].astype(str)


matrix = matrix[matrix.index.isin(re.index)]

matrix['RE'] = re['RE']
matrix['Activity_H3K27ac'] = re['Activity_H3K27ac']
if extra_track == 'none':
    pass
else:
    matrix['Activity_intersect'] = re['Activity_intersect']


# # metaplot1
# H3K4me1 = matrix[matrix.columns[matrix.columns.str.startswith("H3K4me1")]]
# H3K4me3 = matrix[matrix.columns[matrix.columns.str.startswith("H3K4me3")]]
# H3K27ac = matrix[matrix.columns[matrix.columns.str.startswith("H3K27ac")]]
# CTCF    = matrix[matrix.columns[matrix.columns.str.startswith("CTCF")]]

# if upper_threshold == "auto":
#     zMin = np.min([np.max(H3K4me1.values), np.max(H3K4me3.values), np.max(H3K27ac.values), np.max(CTCF.values)])
# else:
#     zMin = upper_threshold

# zMin = np.min([np.max(H3K4me1.values), np.max(H3K4me3.values), np.max(H3K27ac.values), np.max(CTCF.values)])
# # zMin = 1000

# H3K4me1[H3K4me1 > zMin] = zMin
# H3K4me3[H3K4me3 > zMin] = zMin
# H3K27ac[H3K27ac > zMin] = zMin
# CTCF[CTCF > zMin] = zMin

# fig, axs = plt.subplots(2, 4, figsize=(15, 20), gridspec_kw={'height_ratios': [1, 10]})

# axs[0,0] = sns.lineplot(np.array(H3K4me1.mean(axis=0)), ax=axs[0,0])
# axs[0,1] = sns.lineplot(np.array(H3K4me3.mean(axis=0)), ax=axs[0,1])
# axs[0,2] = sns.lineplot(np.array(H3K27ac.mean(axis=0)), ax=axs[0,2])
# axs[0,3] = sns.lineplot(np.array(CTCF.mean(axis=0)), ax=axs[0,3])

# axs[0,0].set_ylabel("")
# axs[0,1].set_ylabel("")
# axs[0,2].set_ylabel("")
# axs[0,3].set_ylabel("")

# axs[0,0].tick_params(axis='both', which='both', length=0)
# axs[0,1].tick_params(axis='both', which='both', length=0)
# axs[0,2].tick_params(axis='both', which='both', length=0)
# axs[0,3].tick_params(axis='both', which='both', length=0)

# axs[0,0].set_yticklabels([])
# axs[0,1].set_yticklabels([])
# axs[0,2].set_yticklabels([])
# axs[0,3].set_yticklabels([])

# axs[0,0].set_xticklabels([])
# axs[0,1].set_xticklabels([])
# axs[0,2].set_xticklabels([])
# axs[0,3].set_xticklabels([])

# axs[0,0].set_title("H3K4me1")
# axs[0,1].set_title("H3K4me3")
# axs[0,2].set_title("H3K27ac")
# axs[0,3].set_title("CTCF")

# axs[1,0].imshow(H3K4me1, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
# axs[1,1].imshow(H3K4me3, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
# axs[1,2].imshow(H3K27ac, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
# ax = axs[1,3].imshow(CTCF, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)

# axs[1,0].tick_params(axis='both', which='both', length=0)
# axs[1,0].set_yticklabels([])
# axs[1,1].tick_params(axis='both', which='both', length=0)
# axs[1,1].set_yticklabels([])
# axs[1,2].tick_params(axis='both', which='both', length=0)
# axs[1,2].set_yticklabels([])
# axs[1,3].tick_params(axis='both', which='both', length=0)
# axs[1,3].set_yticklabels([])

# axs[1,0].set_xticks([30, 
#                      len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
#                      len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
# axs[1,1].set_xticks([30, 
#                      len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
#                      len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
# axs[1,2].set_xticks([30, 
#                      len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
#                      len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
# axs[1,3].set_xticks([30, 
#                      len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
#                      len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])

# axs[1,0].set_xticklabels(["-2","center peak","+2Kb"])
# axs[1,1].set_xticklabels(["-2","center peak","+2Kb"])
# axs[1,2].set_xticklabels(["-2","center peak","+2Kb"])
# axs[1,3].set_xticklabels(["-2","center peak","+2Kb"])

# fig.suptitle(title, y=0.92)

# plt.subplots_adjust(hspace=0.025)

# axins = inset_axes(axs[1,3],
#                     width="15%",  
#                     height="100%",
#                     loc='center right',
#                     borderpad=-5
#                    )

# fig.colorbar(ax, cax=axins)
# fig.savefig(sys.argv[7]+os.sep+"metaplot.pdf", bbox_inches='tight')


# metaplot2
data   = np.array(matrix['RE'].value_counts().tolist()+[0, matrix.shape[0]]).reshape(len(matrix['RE'].value_counts())+2,1)
scaler = MinMaxScaler(feature_range=(0, 10))
ratio  = scaler.fit_transform(data)[:,0][:-2]

H3K4me1 = matrix[matrix.columns[matrix.columns.str.startswith("H3K4me1")]]
H3K4me3 = matrix[matrix.columns[matrix.columns.str.startswith("H3K4me3")]]
H3K27ac = matrix[matrix.columns[matrix.columns.str.startswith("H3K27ac")]]
CTCF    = matrix[matrix.columns[matrix.columns.str.startswith("CTCF")]]

if upper_threshold == "auto":
    zMin = np.min([np.max(H3K4me1.values), np.max(H3K4me3.values), np.max(H3K27ac.values), np.max(CTCF.values)])
else:
    zMin = upper_threshold

H3K4me1[H3K4me1 > zMin] = zMin
H3K4me3[H3K4me3 > zMin] = zMin
H3K27ac[H3K27ac > zMin] = zMin
CTCF[CTCF > zMin] = zMin

matrix_cup = matrix[matrix.columns[:-2]]
if extra_track == 'none':
    pass
else:
    matrix_cup = matrix_cup[matrix_cup.columns[:-1]]
matrix_cup[matrix_cup > zMin] = zMin

matrix_cup['RE'] = re['RE']
matrix_cup['Activity_H3K27ac'] = re['Activity_H3K27ac']
if extra_track == 'none':
    pass
else:
    matrix_cup['Activity_intersect'] = re['Activity_intersect']

fig, axs = plt.subplots(1+len(matrix_cup['RE'].value_counts()), 4, figsize=(15, 20), gridspec_kw={'height_ratios': [1]+list(ratio)})

for idx, clus in enumerate(matrix_cup['RE'].value_counts().index):
    subset = matrix_cup[matrix_cup['RE']==clus]
    
    H3K4me1_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me1")]]
    H3K4me3_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me3")]]
    H3K27ac_tmp = subset[subset.columns[subset.columns.str.startswith("H3K27ac")]]
    CTCF_tmp    = subset[subset.columns[subset.columns.str.startswith("CTCF")]]
    
    axs[0,0] = sns.lineplot(np.array(H3K4me1_tmp.mean(axis=0)), ax=axs[0,0], label=clus, legend=False)
    axs[0,1] = sns.lineplot(np.array(H3K4me3_tmp.mean(axis=0)), ax=axs[0,1], label=clus, legend=False)
    axs[0,2] = sns.lineplot(np.array(H3K27ac_tmp.mean(axis=0)), ax=axs[0,2], label=clus, legend=False)
    axs[0,3] = sns.lineplot(np.array(CTCF_tmp.mean(axis=0)), ax=axs[0,3], label=clus, legend=False)

handles, labels = axs[0,3].get_legend_handles_labels()
fig.legend(handles, labels, loc='center', ncol=len(matrix_cup['RE'].value_counts()), frameon=False, bbox_to_anchor=(0.5,0.91))#, borderaxespad=1)

axs[0,0].set_ylabel("")
axs[0,1].set_ylabel("")
axs[0,2].set_ylabel("")
axs[0,3].set_ylabel("")

axs[0,0].tick_params(axis='both', which='both', length=0)
axs[0,0].set_yticklabels([])
axs[0,1].tick_params(axis='both', which='both', length=0)
axs[0,1].set_yticklabels([])
axs[0,2].tick_params(axis='both', which='both', length=0)
axs[0,2].set_yticklabels([])
axs[0,3].tick_params(axis='both', which='both', length=0)
axs[0,3].set_yticklabels([])

axs[0,0].set_xticklabels([])
axs[0,1].set_xticklabels([])
axs[0,2].set_xticklabels([])
axs[0,3].set_xticklabels([])

axs[0,0].set_title("H3K4me1")
axs[0,1].set_title("H3K4me3")
axs[0,2].set_title("H3K27ac")
axs[0,3].set_title("CTCF")

for idx, clus in enumerate(matrix_cup['RE'].value_counts().index):
    subset = matrix_cup[matrix_cup['RE']==clus]
    
    H3K4me1_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me1")]]
    H3K4me3_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me3")]]
    H3K27ac_tmp = subset[subset.columns[subset.columns.str.startswith("H3K27ac")]]
    CTCF_tmp    = subset[subset.columns[subset.columns.str.startswith("CTCF")]]

    axs[idx+1,0].imshow(H3K4me1_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
    axs[idx+1,1].imshow(H3K4me3_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
    axs[idx+1,2].imshow(H3K27ac_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
    ax = axs[idx+1,3].imshow(CTCF_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)

    
    axs[idx+1,0].tick_params(axis='both', which='both', length=0)
    axs[idx+1,1].tick_params(axis='both', which='both', length=0)
    axs[idx+1,2].tick_params(axis='both', which='both', length=0)
    axs[idx+1,3].tick_params(axis='both', which='both', length=0)

    axs[idx+1,0].set_yticklabels([])
    axs[idx+1,1].set_yticklabels([])
    axs[idx+1,2].set_yticklabels([])
    axs[idx+1,3].set_yticklabels([])

    axs[idx+1,0].set_xticklabels([])
    axs[idx+1,1].set_xticklabels([])
    axs[idx+1,2].set_xticklabels([])
    axs[idx+1,3].set_xticklabels([])
    
    axs[idx+1,0].set_ylabel(clus, rotation=30, labelpad=40)
    
axs[len(matrix_cup['RE'].value_counts()),0].set_xticks([30, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
axs[len(matrix_cup['RE'].value_counts()),1].set_xticks([30, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
axs[len(matrix_cup['RE'].value_counts()),2].set_xticks([30, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
axs[len(matrix_cup['RE'].value_counts()),3].set_xticks([30, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])

axs[len(matrix_cup['RE'].value_counts()),0].set_xticklabels(["-2","center peak","+2Kb"])
axs[len(matrix_cup['RE'].value_counts()),1].set_xticklabels(["-2","center peak","+2Kb"])
axs[len(matrix_cup['RE'].value_counts()),2].set_xticklabels(["-2","center peak","+2Kb"])
axs[len(matrix_cup['RE'].value_counts()),3].set_xticklabels(["-2","center peak","+2Kb"])

fig.suptitle(title, y=0.94)

plt.subplots_adjust(hspace=0.025)

cbar_ax = fig.add_axes([0.925, 0.109, 0.025, 0.7])
fig.colorbar(ax, cax=cbar_ax)
fig.savefig(sys.argv[7]+os.sep+"metaplot_re.pdf", bbox_inches='tight')

# save matrix
matrix.to_csv(sys.argv[8], sep="\t", index=False)



