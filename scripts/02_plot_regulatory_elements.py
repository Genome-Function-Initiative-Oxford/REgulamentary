import warnings
warnings.filterwarnings('ignore')

import sys, os
import pandas as pd
import numpy as np
from sklearn.metrics import auc
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn  as sns
sns.set_style("white")

params1  = sys.argv[1] #"results/CTCF/05_compute_matrix/readable_matrix.mtx"
params2  = sys.argv[2] #"/t1-data/project/fgenomics/sriva/data/upstream_generalisation/upstream_pipelines_output/2022-11-16_POLR2A/results/08_bam_coverages/POLR2A.bw"
params3  = sys.argv[3] #60
params4  = sys.argv[4] #"POLR2A"
params5  = sys.argv[5] #"auto"
params6  = sys.argv[6] #"Chromatin marks at open chromatin sites (+/- 2kb)"
params7  = sys.argv[7] #"results/CTCF/08_regulatory_elements/01_plots"
params8  = sys.argv[8] #"results/CTCF/08_regulatory_elements/02_plots"
params9  = sys.argv[9] #"results/CTCF/08_regulatory_elements/03_plots"
params10 = sys.argv[10] #"results/CTCF/07_active_elements/sorted_regions.bed"
params11 = sys.argv[11] #"results/CTCF/07_active_elements/POLR2A_L-tron_filtered_intersection.bed"
params12 = sys.argv[12] #"results/CTCF/08_regulatory_elements/mlv.csv"

df         = pd.read_csv(params1, sep='\t', comment='#')
df.columns = df.columns[1:].tolist() + ['remove']

H3K4me1_c = ["H3K4me1.%i"%i for i in range(1,401)]
H3K4me3_c = ["H3K4me3.%i"%i for i in range(1,401)]
H3K27ac_c = ["H3K27ac.%i"%i for i in range(1,401)]
CTCF_c    = ["CTCF.%i"%i    for i in range(1,401)]
if params2 != "none":
    PARAMS4_c = ["%s.%i"%(params4, i) for i in range(1,401)]

if params2 != "none":
    new_columns = H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c + PARAMS4_c
else:
    new_columns = H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c
df.columns = new_columns + ['remove']

df = df[new_columns]
df = df.fillna(0)


H3K4me1 = df[df.columns[df.columns.str.startswith("H3K4me1")]]
H3K4me3 = df[df.columns[df.columns.str.startswith("H3K4me3")]]
H3K27ac = df[df.columns[df.columns.str.startswith("H3K27ac")]]
CTCF    = df[df.columns[df.columns.str.startswith("CTCF")]]
if params2 != "none":
    PARAMS4 = df[df.columns[df.columns.str.startswith(params4)]]

df1 = df[H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c]
if params2 != "none":
    df2 = df[PARAMS4_c]


if params2 != "none":
    df = pd.concat([df1, df2], axis=1)
else:
    df = df1

### USER PARAM
if params5=="auto":
    zMin = np.min([np.max(H3K4me1.values), np.max(H3K4me3.values), np.max(H3K27ac.values), np.max(CTCF.values)])
else:
    zMin = params5 # to implement: if params3 greater than ...
    
### USER PARAM
if params3 is not None: # to implement: if params2 bigger than number of bins ...
    gap = int(params3)
else:
    gap = 400 # to implement: automatic size of dataframes
gap_2 = gap/2

H3K4me1_ce = np.arange(int(H3K4me1.shape[1]/2-gap_2), int(H3K4me1.shape[1]/2+gap_2))
H3K4me3_ce = np.arange(int(H3K4me3.shape[1]/2-gap_2), int(H3K4me3.shape[1]/2+gap_2))
H3K27ac_ce = np.arange(int(H3K27ac.shape[1]/2-gap_2), int(H3K27ac.shape[1]/2+gap_2))
CTCF_ce    = np.arange(int(CTCF.shape[1]/2-gap_2),    int(CTCF.shape[1]/2+gap_2))

H3K4me1_ce = H3K4me1[H3K4me1.columns[H3K4me1_ce]]
H3K4me3_ce = H3K4me3[H3K4me3.columns[H3K4me3_ce]]
H3K27ac_ce = H3K27ac[H3K27ac.columns[H3K27ac_ce]]
CTCF_ce    = CTCF[CTCF.columns[CTCF_ce]]

### AUC # to implement: auc faster version
H3K4me1_auc, H3K4me3_auc, H3K27ac_auc, CTCF_auc = [], [], [], []
for i in range(df.shape[0]):
    H3K4me1_auc.append(auc(np.arange(0, H3K4me1_ce.shape[1]), np.array(H3K4me1_ce.T[i])))
    H3K4me3_auc.append(auc(np.arange(0, H3K4me3_ce.shape[1]), np.array(H3K4me3_ce.T[i])))
    H3K27ac_auc.append(auc(np.arange(0, H3K27ac_ce.shape[1]), np.array(H3K27ac_ce.T[i])))
    CTCF_auc.append(auc(np.arange(0, CTCF_ce.shape[1]), np.array(CTCF_ce.T[i])))

### RANK AUC
df_auc         = pd.DataFrame([H3K4me1_auc, H3K4me3_auc, H3K27ac_auc, CTCF_auc]).T
df_auc.columns = ['H3K4me1', 'H3K4me3', 'H3K27ac', 'CTCF']
df_rank        = df_auc.rank(1, ascending=False, method='first')

### RULES' DEFINITION
df_rules = df_rank.copy()

df_rules.loc[((df_rules['H3K4me1'] == 1.0) & (df_rules['H3K4me3'] == 2.0)) | 
             ((df_rules['H3K4me1'] == 2.0) & (df_rules['H3K4me3'] == 1.0)), 'clusters_ru'] = 'Enhancer/Promoter'

df_rules.loc[((df_rules['H3K4me1'] == 1.0) & (df_rules['H3K27ac'] == 2.0)) | 
             ((df_rules['H3K4me1'] == 2.0) & (df_rules['H3K27ac'] == 1.0)), 'clusters_ru'] = 'Enhancer'  

df_rules.loc[((df_rules['H3K4me1'] == 1.0) & (df_rules['CTCF'] == 2.0)) | 
             ((df_rules['H3K4me1'] == 2.0) & (df_rules['CTCF'] == 1.0)), 'clusters_ru'] = 'Enhancer/CTCF'

df_rules.loc[((df_rules['H3K4me3'] == 1.0) & (df_rules['H3K27ac'] == 2.0)) | 
             ((df_rules['H3K4me3'] == 2.0) & (df_rules['H3K27ac'] == 1.0)), 'clusters_ru'] = 'Promoter'

df_rules.loc[((df_rules['H3K4me3'] == 1.0) & (df_rules['CTCF'] == 2.0)) | 
             ((df_rules['H3K4me3'] == 2.0) & (df_rules['CTCF'] == 1.0)), 'clusters_ru'] = 'Promoter/CTCF'

df_rules.loc[((df_rules['H3K27ac'] == 1.0) & (df_rules['CTCF'] == 2.0)) | 
             ((df_rules['H3K27ac'] == 2.0) & (df_rules['CTCF'] == 1.0)), 'clusters_ru'] = 'CTCF'

df['Cluster'] = df_rules['clusters_ru']

df_rules.loc[(df_rules['H3K27ac'] != 1.0), 'active_ru'] = 'Inactive'
df_rules.loc[(df_rules['H3K27ac'] == 1.0), 'active_ru'] = 'Active'

df['Active1'] = df_rules['active_ru']

regions = pd.read_csv(params10, sep='\t', header=None)
regions.columns = ["chr", "start", "end"]
regions['id'] = regions["chr"].astype(str)+":"+regions["start"].astype(str)+"_"+regions["end"].astype(str)

if params2 != "none":
    
    peak = pd.read_csv(params11, sep='\t', header=None)
    peak['id'] = peak[0].astype(str)+":"+peak[1].astype(str)+"_"+peak[2].astype(str)

    feat = []
    for ID in regions["id"]:
        if ID in peak["id"].tolist():
            feat.append("Active")
        else:
            feat.append("Inactive")
    del peak

    df['Active2'] = feat
    
regions["RE_rule"] = df['Cluster']
regions["Active_rule"] = df['Active1']
if params2 != "none":
    regions["Active_inter"] = df['Active2']

clss  = df['Cluster']
clss2 = df['Active1']

final_df            = df[H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c]
final_df['Cluster'] = clss
final_df['Active1'] = clss2
final_df            = final_df.reset_index()
final_df            = final_df.sort_values(by=['Cluster', 'index'])

########################################################################################################################
########################################################################################################################
########################################################################################################################
### 08_regulatory_elements/01_plots/heatmap.pdf
H3K4me1 = df[df.columns[df.columns.str.startswith("H3K4me1")]]
H3K4me3 = df[df.columns[df.columns.str.startswith("H3K4me3")]]
H3K27ac = df[df.columns[df.columns.str.startswith("H3K27ac")]]
CTCF    = df[df.columns[df.columns.str.startswith("CTCF")]]

fig, axs = plt.subplots(2, 4, figsize=(15, 20), gridspec_kw={'height_ratios': [1, 10]})

axs[0,0] = sns.lineplot(np.array(H3K4me1.mean(axis=0)), ax=axs[0,0])
axs[0,1] = sns.lineplot(np.array(H3K4me3.mean(axis=0)), ax=axs[0,1])
axs[0,2] = sns.lineplot(np.array(H3K27ac.mean(axis=0)), ax=axs[0,2])
axs[0,3] = sns.lineplot(np.array(CTCF.mean(axis=0)), ax=axs[0,3])

axs[0,0].set_ylabel("")
axs[0,1].set_ylabel("")
axs[0,2].set_ylabel("")
axs[0,3].set_ylabel("")

axs[0,0].set_yticklabels([])
axs[0,1].set_yticklabels([])
axs[0,2].set_yticklabels([])
axs[0,3].set_yticklabels([])

axs[0,0].set_xticklabels([])
axs[0,1].set_xticklabels([])
axs[0,2].set_xticklabels([])
axs[0,3].set_xticklabels([])

axs[0,0].set_title("H3K4me1")
axs[0,1].set_title("H3K4me3")
axs[0,2].set_title("H3K27ac")
axs[0,3].set_title("CTCF")

axs[1,0].imshow(H3K4me1, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
axs[1,1].imshow(H3K4me3, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
axs[1,2].imshow(H3K27ac, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
ax = axs[1,3].imshow(CTCF, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)

axs[1,0].set_yticklabels([])
axs[1,1].set_yticklabels([])
axs[1,2].set_yticklabels([])
axs[1,3].set_yticklabels([])

axs[1,0].set_xticks([30, 
                     len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                     len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
axs[1,1].set_xticks([30, 
                     len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                     len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
axs[1,2].set_xticks([30, 
                     len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                     len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
axs[1,3].set_xticks([30, 
                     len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                     len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])

axs[1,0].set_xticklabels(["-2","center peak","+2Kb"])
axs[1,1].set_xticklabels(["-2","center peak","+2Kb"])
axs[1,2].set_xticklabels(["-2","center peak","+2Kb"])
axs[1,3].set_xticklabels(["-2","center peak","+2Kb"])

fig.suptitle(params6, y=0.92)

plt.subplots_adjust(hspace=0.025)

axins = inset_axes(axs[1,3],
                    width="15%",  
                    height="100%",
                    loc='center right',
                    borderpad=-5
                   )
plt.colorbar(ax, cax=axins)

if not os.path.exists(params7):
    os.makedirs(params7)
    
fig.savefig(params7+os.sep+"heatmap.pdf", bbox_inches='tight')

########################################################################################################################
########################################################################################################################
########################################################################################################################
### 08_regulatory_elements/01_plots/clustered_heatmap.pdf

data   = np.array(final_df['Cluster'].value_counts().tolist()+[0, final_df.shape[0]]).reshape(len(final_df['Cluster'].value_counts())+2,1)
scaler = MinMaxScaler(feature_range=(0, 10))
ratio  = scaler.fit_transform(data)[:,0][:-2]

H3K4me1 = df[df.columns[df.columns.str.startswith("H3K4me1")]]
H3K4me3 = df[df.columns[df.columns.str.startswith("H3K4me3")]]
H3K27ac = df[df.columns[df.columns.str.startswith("H3K27ac")]]
CTCF    = df[df.columns[df.columns.str.startswith("CTCF")]]

fig, axs = plt.subplots(1+len(final_df['Cluster'].value_counts()), 4, figsize=(15, 20), gridspec_kw={'height_ratios': [1]+list(ratio)})

for idx, clus in enumerate(final_df['Cluster'].value_counts().index):
    subset = final_df[final_df['Cluster']==clus]
    
    H3K4me1_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me1")]]
    H3K4me3_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me3")]]
    H3K27ac_tmp = subset[subset.columns[subset.columns.str.startswith("H3K27ac")]]
    CTCF_tmp    = subset[subset.columns[subset.columns.str.startswith("CTCF")]]
    
    axs[0,0] = sns.lineplot(np.array(H3K4me1_tmp.mean(axis=0)), ax=axs[0,0], label=clus, legend=False)
    axs[0,1] = sns.lineplot(np.array(H3K4me3_tmp.mean(axis=0)), ax=axs[0,1], label=clus, legend=False)
    axs[0,2] = sns.lineplot(np.array(H3K27ac_tmp.mean(axis=0)), ax=axs[0,2], label=clus, legend=False)
    axs[0,3] = sns.lineplot(np.array(CTCF_tmp.mean(axis=0)), ax=axs[0,3], label=clus, legend=False)

handles, labels = axs[0,3].get_legend_handles_labels()
fig.legend(handles, labels, loc='center', ncol=len(final_df['Cluster'].value_counts()), frameon=False, bbox_to_anchor=(0.5,0.91))#, borderaxespad=1)

axs[0,0].set_ylabel("")
axs[0,1].set_ylabel("")
axs[0,2].set_ylabel("")
axs[0,3].set_ylabel("")

axs[0,0].set_yticklabels([])
axs[0,1].set_yticklabels([])
axs[0,2].set_yticklabels([])
axs[0,3].set_yticklabels([])

axs[0,0].set_xticklabels([])
axs[0,1].set_xticklabels([])
axs[0,2].set_xticklabels([])
axs[0,3].set_xticklabels([])

axs[0,0].set_title("H3K4me1")
axs[0,1].set_title("H3K4me3")
axs[0,2].set_title("H3K27ac")
axs[0,3].set_title("CTCF")

for idx, clus in enumerate(final_df['Cluster'].value_counts().index):
    subset = final_df[final_df['Cluster']==clus]
    
    H3K4me1_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me1")]]
    H3K4me3_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me3")]]
    H3K27ac_tmp = subset[subset.columns[subset.columns.str.startswith("H3K27ac")]]
    CTCF_tmp    = subset[subset.columns[subset.columns.str.startswith("CTCF")]]

    axs[idx+1,0].imshow(H3K4me1_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
    axs[idx+1,1].imshow(H3K4me3_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
    axs[idx+1,2].imshow(H3K27ac_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
    ax = axs[idx+1,3].imshow(CTCF_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)


    axs[idx+1,0].set_yticklabels([])
    axs[idx+1,1].set_yticklabels([])
    axs[idx+1,2].set_yticklabels([])
    axs[idx+1,3].set_yticklabels([])

    axs[idx+1,0].set_xticklabels([])
    axs[idx+1,1].set_xticklabels([])
    axs[idx+1,2].set_xticklabels([])
    axs[idx+1,3].set_xticklabels([])
    
    axs[idx+1,0].set_ylabel(clus, rotation=45, labelpad=40)
    
axs[len(final_df['Cluster'].value_counts()),0].set_xticks([30, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
axs[len(final_df['Cluster'].value_counts()),1].set_xticks([30, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
axs[len(final_df['Cluster'].value_counts()),2].set_xticks([30, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])
axs[len(final_df['Cluster'].value_counts()),3].set_xticks([30, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4/2, 
                                                           len(H3K4me1_c + H3K4me3_c + H3K27ac_c + CTCF_c)/4-30])

axs[len(final_df['Cluster'].value_counts()),0].set_xticklabels(["-2","center peak","+2Kb"])
axs[len(final_df['Cluster'].value_counts()),1].set_xticklabels(["-2","center peak","+2Kb"])
axs[len(final_df['Cluster'].value_counts()),2].set_xticklabels(["-2","center peak","+2Kb"])
axs[len(final_df['Cluster'].value_counts()),3].set_xticklabels(["-2","center peak","+2Kb"])

fig.suptitle(params6, y=0.94)

plt.subplots_adjust(hspace=0.025)

cbar_ax = fig.add_axes([0.925, 0.109, 0.025, 0.7])
fig.colorbar(ax, cax=cbar_ax)

fig.savefig(params7+os.sep+"clustered_heatmap.pdf", bbox_inches='tight')


########################################################################################################################
########################################################################################################################
########################################################################################################################
### 08_regulatory_elements/02_plots/...
H3K4me1 = df[df.columns[df.columns.str.startswith("H3K4me1")]]
H3K4me3 = df[df.columns[df.columns.str.startswith("H3K4me3")]]
H3K27ac = df[df.columns[df.columns.str.startswith("H3K27ac")]]
CTCF    = df[df.columns[df.columns.str.startswith("CTCF")]]

df_plot            = pd.concat([H3K4me1, H3K4me3, H3K27ac, CTCF], axis=1)
df_plot['Cluster'] = df['Cluster']
df_plot['Active1'] = df['Active1']

if not os.path.exists(params8):
    os.makedirs(params8)

for cluster in np.unique(df_plot['Cluster']):
        
    tmp = df_plot[df_plot['Cluster']==cluster]

    data   = np.array(tmp['Active1'].value_counts().tolist()+[0, tmp.shape[0]]).reshape(len(tmp['Active1'].value_counts())+2,1)
    scaler = MinMaxScaler(feature_range=(0, 10))
    ratio  = scaler.fit_transform(data)[:,0][:-2]

    fig, axs = plt.subplots(1+len(tmp['Active1'].value_counts()), 4, figsize=(15, 8), gridspec_kw={'height_ratios': [1]+list(ratio)})

    for idx, clus in enumerate(tmp['Active1'].value_counts().index):
        subset = tmp[tmp['Active1']==clus]

        H3K4me1_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me1")]]
        H3K4me3_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me3")]]
        H3K27ac_tmp = subset[subset.columns[subset.columns.str.startswith("H3K27ac")]]
        CTCF_tmp    = subset[subset.columns[subset.columns.str.startswith("CTCF")]]


        axs[0,0] = sns.lineplot(np.array(H3K4me1_tmp.mean(axis=0)), ax=axs[0,0], label=clus, legend=False)
        axs[0,1] = sns.lineplot(np.array(H3K4me3_tmp.mean(axis=0)), ax=axs[0,1], label=clus, legend=False)
        axs[0,2] = sns.lineplot(np.array(H3K27ac_tmp.mean(axis=0)), ax=axs[0,2], label=clus, legend=False)
        axs[0,3] = sns.lineplot(np.array(CTCF_tmp.mean(axis=0)), ax=axs[0,3], label=clus, legend=False)

    handles, labels = axs[0,3].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center', ncol=len(tmp['Active1'].value_counts()), frameon=False, bbox_to_anchor=(0.5,0.93))#, borderaxespad=1)

    axs[0,0].set_ylabel("")
    axs[0,1].set_ylabel("")
    axs[0,2].set_ylabel("")
    axs[0,3].set_ylabel("")

    axs[0,0].set_yticklabels([])
    axs[0,1].set_yticklabels([])
    axs[0,2].set_yticklabels([])
    axs[0,3].set_yticklabels([])

    axs[0,0].set_xticklabels([])
    axs[0,1].set_xticklabels([])
    axs[0,2].set_xticklabels([])
    axs[0,3].set_xticklabels([])

    axs[0,0].set_title("H3K4me1")
    axs[0,1].set_title("H3K4me3")
    axs[0,2].set_title("H3K27ac")
    axs[0,3].set_title("CTCF")

    for idx, clus in enumerate(tmp['Active1'].value_counts().index):
        subset = tmp[tmp['Active1']==clus]

        H3K4me1_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me1")]]
        H3K4me3_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me3")]]
        H3K27ac_tmp = subset[subset.columns[subset.columns.str.startswith("H3K27ac")]]
        CTCF_tmp    = subset[subset.columns[subset.columns.str.startswith("CTCF")]]
        EXTRA_tmp   = subset[subset.columns[subset.columns.str.startswith(params4)]]

        axs[idx+1,0].imshow(H3K4me1_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
        axs[idx+1,1].imshow(H3K4me3_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
        axs[idx+1,2].imshow(H3K27ac_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
        ax = axs[idx+1,3].imshow(CTCF_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)

        axs[idx+1,0].set_yticklabels([])
        axs[idx+1,1].set_yticklabels([])
        axs[idx+1,2].set_yticklabels([])
        axs[idx+1,3].set_yticklabels([])

        axs[idx+1,0].set_xticklabels([])
        axs[idx+1,1].set_xticklabels([])
        axs[idx+1,2].set_xticklabels([])
        axs[idx+1,3].set_xticklabels([])

        axs[idx+1,0].set_ylabel(clus, rotation=45, labelpad=40)

    axs[len(tmp['Active1'].value_counts()),0].set_xticks([30, 
                                                        (df_plot.shape[1]-1)/5/2, 
                                                        (df_plot.shape[1]-1)/5-30])
    axs[len(tmp['Active1'].value_counts()),1].set_xticks([30, 
                                                        (df_plot.shape[1]-1)/5/2, 
                                                        (df_plot.shape[1]-1)/5-30])
    axs[len(tmp['Active1'].value_counts()),2].set_xticks([30, 
                                                        (df_plot.shape[1]-1)/5/2, 
                                                        (df_plot.shape[1]-1)/5-30])
    axs[len(tmp['Active1'].value_counts()),3].set_xticks([30, 
                                                        (df_plot.shape[1]-1)/5/2, 
                                                        (df_plot.shape[1]-1)/5-30])

    axs[len(tmp['Active1'].value_counts()),0].set_xticklabels(["-2","center peak","+2Kb"])
    axs[len(tmp['Active1'].value_counts()),1].set_xticklabels(["-2","center peak","+2Kb"])
    axs[len(tmp['Active1'].value_counts()),2].set_xticklabels(["-2","center peak","+2Kb"])
    axs[len(tmp['Active1'].value_counts()),3].set_xticklabels(["-2","center peak","+2Kb"])

    fig.suptitle(cluster.replace("/", "-"), y=0.98)

    plt.subplots_adjust(hspace=0.025)

    cbar_ax = fig.add_axes([0.925, 0.109, 0.025, 0.7])
    fig.colorbar(ax, cax=cbar_ax)

    fig.savefig(params8+os.sep+"heatmap_%s.pdf"%cluster.replace("/", "-"), bbox_inches='tight')

########################################################################################################################
########################################################################################################################
########################################################################################################################
if params2 != "none":
    
    if not os.path.exists(params9):
        os.makedirs(params9)
    
    ### 08_regulatory_elements/02_plots/...
    H3K4me1 = df[df.columns[df.columns.str.startswith("H3K4me1")]]
    H3K4me3 = df[df.columns[df.columns.str.startswith("H3K4me3")]]
    H3K27ac = df[df.columns[df.columns.str.startswith("H3K27ac")]]
    CTCF    = df[df.columns[df.columns.str.startswith("CTCF")]]
    extra   = df[df.columns[df.columns.str.startswith(params4)]]

    df_plot            = pd.concat([H3K4me1, H3K4me3, H3K27ac, CTCF, extra], axis=1)
    df_plot['Cluster'] = df['Cluster']
    df_plot['Active2'] = df['Active2']
    # regions = pd.read_csv(params8, sep='\t', index_col=0)
    # df_plot['Active2'] = df['Active2']

    for cluster in np.unique(df_plot['Cluster']):

        tmp = df_plot[df_plot['Cluster']==cluster]

        data   = np.array(tmp['Active2'].value_counts().tolist()+[0, tmp.shape[0]]).reshape(len(tmp['Active2'].value_counts())+2,1)
        scaler = MinMaxScaler(feature_range=(0, 10))
        ratio  = scaler.fit_transform(data)[:,0][:-2]

        fig, axs = plt.subplots(1+len(tmp['Active2'].value_counts()), 5, figsize=(15, 8), gridspec_kw={'height_ratios': [1]+list(ratio)})

        for idx, clus in enumerate(tmp['Active2'].value_counts().index):
            subset = tmp[tmp['Active2']==clus]

            H3K4me1_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me1")]]
            H3K4me3_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me3")]]
            H3K27ac_tmp = subset[subset.columns[subset.columns.str.startswith("H3K27ac")]]
            CTCF_tmp    = subset[subset.columns[subset.columns.str.startswith("CTCF")]]
            EXTRA_tmp   = subset[subset.columns[subset.columns.str.startswith(params4)]]


            axs[0,0] = sns.lineplot(np.array(H3K4me1_tmp.mean(axis=0)), ax=axs[0,0], label=clus, legend=False)
            axs[0,1] = sns.lineplot(np.array(H3K4me3_tmp.mean(axis=0)), ax=axs[0,1], label=clus, legend=False)
            axs[0,2] = sns.lineplot(np.array(H3K27ac_tmp.mean(axis=0)), ax=axs[0,2], label=clus, legend=False)
            axs[0,3] = sns.lineplot(np.array(CTCF_tmp.mean(axis=0)), ax=axs[0,3], label=clus, legend=False)
            axs[0,4] = sns.lineplot(np.array(EXTRA_tmp.mean(axis=0)), ax=axs[0,4], label=clus, legend=False)

        handles, labels = axs[0,4].get_legend_handles_labels()
        fig.legend(handles, labels, loc='center', ncol=len(tmp['Active2'].value_counts()), frameon=False, bbox_to_anchor=(0.5,0.93))#, borderaxespad=1)

        axs[0,0].set_ylabel("")
        axs[0,1].set_ylabel("")
        axs[0,2].set_ylabel("")
        axs[0,3].set_ylabel("")
        axs[0,4].set_ylabel("")

        axs[0,0].set_yticklabels([])
        axs[0,1].set_yticklabels([])
        axs[0,2].set_yticklabels([])
        axs[0,3].set_yticklabels([])
        axs[0,4].set_yticklabels([])

        axs[0,0].set_xticklabels([])
        axs[0,1].set_xticklabels([])
        axs[0,2].set_xticklabels([])
        axs[0,3].set_xticklabels([])
        axs[0,4].set_xticklabels([])

        axs[0,0].set_title("H3K4me1")
        axs[0,1].set_title("H3K4me3")
        axs[0,2].set_title("H3K27ac")
        axs[0,3].set_title("CTCF")
        axs[0,4].set_title(params4)

        for idx, clus in enumerate(tmp['Active2'].value_counts().index):
            subset = tmp[tmp['Active2']==clus]

            H3K4me1_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me1")]]
            H3K4me3_tmp = subset[subset.columns[subset.columns.str.startswith("H3K4me3")]]
            H3K27ac_tmp = subset[subset.columns[subset.columns.str.startswith("H3K27ac")]]
            CTCF_tmp    = subset[subset.columns[subset.columns.str.startswith("CTCF")]]
            EXTRA_tmp   = subset[subset.columns[subset.columns.str.startswith(params4)]]

            axs[idx+1,0].imshow(H3K4me1_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
            axs[idx+1,1].imshow(H3K4me3_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
            axs[idx+1,2].imshow(H3K27ac_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
            axs[idx+1,3].imshow(CTCF_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)
            ax = axs[idx+1,4].imshow(EXTRA_tmp, interpolation='bilinear', cmap='Spectral', aspect='auto', origin='upper', vmin=0, vmax=zMin)


            axs[idx+1,0].set_yticklabels([])
            axs[idx+1,1].set_yticklabels([])
            axs[idx+1,2].set_yticklabels([])
            axs[idx+1,3].set_yticklabels([])
            axs[idx+1,4].set_yticklabels([])

            axs[idx+1,0].set_xticklabels([])
            axs[idx+1,1].set_xticklabels([])
            axs[idx+1,2].set_xticklabels([])
            axs[idx+1,3].set_xticklabels([])
            axs[idx+1,4].set_xticklabels([])

            axs[idx+1,0].set_ylabel(clus, rotation=45, labelpad=40)

        axs[len(tmp['Active2'].value_counts()),0].set_xticks([30, 
                                                           (df_plot.shape[1]-1)/5/2, 
                                                           (df_plot.shape[1]-1)/5-30])
        axs[len(tmp['Active2'].value_counts()),1].set_xticks([30, 
                                                           (df_plot.shape[1]-1)/5/2, 
                                                           (df_plot.shape[1]-1)/5-30])
        axs[len(tmp['Active2'].value_counts()),2].set_xticks([30, 
                                                           (df_plot.shape[1]-1)/5/2, 
                                                           (df_plot.shape[1]-1)/5-30])
        axs[len(tmp['Active2'].value_counts()),3].set_xticks([30, 
                                                           (df_plot.shape[1]-1)/5/2, 
                                                           (df_plot.shape[1]-1)/5-30])
        axs[len(tmp['Active2'].value_counts()),4].set_xticks([30, 
                                                           (df_plot.shape[1]-1)/5/2, 
                                                           (df_plot.shape[1]-1)/5-30])

        axs[len(tmp['Active2'].value_counts()),0].set_xticklabels(["-2","center peak","+2Kb"])
        axs[len(tmp['Active2'].value_counts()),1].set_xticklabels(["-2","center peak","+2Kb"])
        axs[len(tmp['Active2'].value_counts()),2].set_xticklabels(["-2","center peak","+2Kb"])
        axs[len(tmp['Active2'].value_counts()),3].set_xticklabels(["-2","center peak","+2Kb"])
        axs[len(tmp['Active2'].value_counts()),4].set_xticklabels(["-2","center peak","+2Kb"])

        fig.suptitle(cluster.replace("/", "-"), y=0.97)

        plt.subplots_adjust(hspace=0.025)

        cbar_ax = fig.add_axes([0.925, 0.109, 0.025, 0.7])
        fig.colorbar(ax, cax=cbar_ax)

        fig.savefig(params9+os.sep+"heatmap_%s.pdf"%cluster.replace("/", "-"), bbox_inches='tight')

########################################################################################################################
########################################################################################################################
########################################################################################################################
regions.to_csv(params12, sep="\t", index=False)



