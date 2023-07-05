import warnings
warnings.filterwarnings('ignore')

import sys, os, shutil
import pandas as pd
import numpy as np
import pybedtools
import pyBigWig
from sklearn.metrics import auc
import math
import multiprocessing as mp

def _auc(x_r, y_bw):
    if x_r.empty:
        return 0
    else:
        aucs = []
        for j in range(x_r.shape[0]): 
            if tuple(x_r.iloc[j])[2]-tuple(x_r.iloc[j])[1] <= 1:
                print("skipping %s:%d-%d for auc ..."%(tuple(x_r.iloc[j])[0], tuple(x_r.iloc[j])[1], tuple(x_r.iloc[j])[2]))
                aucs.append(0.0)
            else:
                cov = y_bw.values(tuple(x_r.iloc[j])[0], tuple(x_r.iloc[j])[1], tuple(x_r.iloc[j])[2])
                cov = [0.0 if math.isnan(x) else x for x in cov]
                aucs.append(auc(np.arange(0, len(cov)), np.array(cov)))
        return np.mean(np.array(aucs))

def _coverage(reg, bw1, bw2, bw3, bw4):
    bw1_tmp = pyBigWig.open(bw1)
    bw2_tmp = pyBigWig.open(bw2)
    bw3_tmp = pyBigWig.open(bw3)
    bw4_tmp = pyBigWig.open(bw4)

    coverage = []

    for i in range(reg.shape[0]):

        r = reg.iloc[i]
        sort_r = pybedtools.BedTool.from_dataframe(pd.DataFrame(r).T)
        chromosome_r  = r[0]
        start_r  = r[1]
        end_r    = r[2]

        H3K4me1_r = pybedtools.BedTool.from_dataframe(H3K4me1_peaks)
        H3K4me3_r = pybedtools.BedTool.from_dataframe(H3K4me3_peaks)
        H3K27ac_r = pybedtools.BedTool.from_dataframe(H3K27ac_peaks)
        CTCF_r    = pybedtools.BedTool.from_dataframe(CTCF_peaks)

        H3K4me1_r = sort_r.intersect(H3K4me1_r).to_dataframe()
        H3K4me3_r = sort_r.intersect(H3K4me3_r).to_dataframe()
        H3K27ac_r = sort_r.intersect(H3K27ac_r).to_dataframe()
        CTCF_r    = sort_r.intersect(CTCF_r).to_dataframe()

        H3K4me1_auc = _auc(H3K4me1_r, bw1_tmp)
        H3K4me3_auc = _auc(H3K4me3_r, bw2_tmp)
        H3K27ac_auc = _auc(H3K27ac_r, bw3_tmp)
        CTCF_auc    = _auc(CTCF_r, bw4_tmp)
        
        if (H3K4me1_auc == 0.0) and (H3K4me3_auc == 0.0) and (CTCF_auc == 0.0):
            print("skipping %s:%d-%d ..."%(chromosome_r, start_r, end_r))
        else:
            coverage.append([chromosome_r, start_r, end_r, H3K4me1_auc, H3K4me3_auc, H3K27ac_auc, CTCF_auc])

    bw1_tmp.close()
    bw2_tmp.close()
    bw3_tmp.close()
    bw4_tmp.close()
            
    return coverage


params1  = sys.argv[1] #{input.sort_union} \
params2  = sys.argv[2] #{input.H3K4me1_peaks} \
params3  = sys.argv[3] #{input.H3K4me3_peaks} \
params4  = sys.argv[4] #{input.H3K27ac_peaks} \
params5  = sys.argv[5] #{input.CTCF_peaks} \
params6  = sys.argv[6] #{input.H3K4me1_bw} \
params7  = sys.argv[7] #{input.H3K4me3_bw} \
params8  = sys.argv[8] #{input.H3K27ac_bw} \
params9  = sys.argv[9] #{input.CTCF_bw} \
params10 = sys.argv[10] #{params.extra_peaks} \
params11 = sys.argv[11] #{params.extra_bw} \
params12 = sys.argv[12] #{params.thresholdPeaks} \
params13 = sys.argv[13] #{output}
params14 = sys.argv[14] #{params.tss}
params15 = sys.argv[15] #{params.tss}

if not os.path.exists(params15):
    os.makedirs(params15)
pybedtools.set_tempdir(params15)


init_bed   = params1
sort_union = pd.read_csv(init_bed, sep='\t', header=None)

peak_thr = float(params12)

H3K4me1_peaks = params2
H3K4me1_peaks = pd.read_csv(H3K4me1_peaks, sep='\t')
H3K4me1_peaks = H3K4me1_peaks[H3K4me1_peaks['overall_peak_score'] >= peak_thr][['chrom', 'start', 'end']]

H3K4me3_peaks = params3
H3K4me3_peaks = pd.read_csv(H3K4me3_peaks, sep='\t')
H3K4me3_peaks = H3K4me3_peaks[H3K4me3_peaks['overall_peak_score'] >= peak_thr][['chrom', 'start', 'end']]

H3K27ac_peaks = params4
H3K27ac_peaks = pd.read_csv(H3K27ac_peaks, sep='\t')
H3K27ac_peaks = H3K27ac_peaks[H3K27ac_peaks['overall_peak_score'] >= peak_thr][['chrom', 'start', 'end']]

CTCF_peaks    = params5
CTCF_peaks    = pd.read_csv(CTCF_peaks, sep='\t')
CTCF_peaks    = CTCF_peaks[CTCF_peaks['overall_peak_score'] >= peak_thr][['chrom', 'start', 'end']]


H3K4me1_bw = params6
H3K4me3_bw = params7
H3K27ac_bw = params8
CTCF_bw    = params9

regions = []
for i in sort_union[0].value_counts().index:
    regions.append(sort_union[sort_union[0]==i])


n_cpu = mp.cpu_count()
if n_cpu > 10:
    n_cpu = n_cpu - 5
pool = mp.Pool(n_cpu)
results = [pool.apply_async(_coverage, args=(reg, H3K4me1_bw, H3K4me3_bw, H3K27ac_bw, CTCF_bw)) for reg in regions]
pool.close()

coverage = []
for rs in results:
    for r in rs.get():
        coverage.append(r)

df_auc = pd.DataFrame(coverage, columns=['chromosome', 'start', 'end', 'H3K4me1_auc', 'H3K4me3_auc', 'H3K27ac_auc', 'CTCF_auc'])

df_rank = df_auc[['H3K4me1_auc', 'H3K4me3_auc', 'H3K27ac_auc', 'CTCF_auc']]
df_rank = df_rank.rank(1, ascending=False, method='first')
df_rank.columns = ['H3K4me1_rank', 'H3K4me3_rank', 'H3K27ac_rank', 'CTCF_rank']

df = pd.concat([df_auc, df_rank], axis=1)

df['RE'] = 'Not assigned'

df.loc[((df['H3K4me1_rank'] == 1.0) & (df['H3K4me3_rank'] == 2.0)) | 
       ((df['H3K4me1_rank'] == 2.0) & (df['H3K4me3_rank'] == 1.0)), 'RE'] = 'Enhancer/Promoter'

df.loc[((df['H3K4me1_rank'] == 1.0) & (df['H3K27ac_rank'] == 2.0)) | 
       ((df['H3K4me1_rank'] == 2.0) & (df['H3K27ac_rank'] == 1.0)) | 
       ((df['H3K4me1_rank'] == 1.0) & (df['H3K4me3_auc'] == 0.0) & (df['H3K27ac_auc'] == 0.0) & (df['CTCF_auc'] == 0.0)), 'RE'] = 'Enhancer'

df.loc[((df['H3K4me1_rank'] == 1.0) & (df['CTCF_rank'] == 2.0)) | 
       ((df['H3K4me1_rank'] == 2.0) & (df['CTCF_rank'] == 1.0)) |
       ((df['H3K4me1_auc'] != 0.0) & (df['H3K4me3_auc'] == 0.0) & (df['H3K27ac_auc'] != 0.0) & (df['CTCF_auc'] != 0.0)), 'RE'] = 'Enhancer/CTCF'

df.loc[((df['H3K4me3_rank'] == 1.0) & (df['H3K27ac_rank'] == 2.0)) | 
       ((df['H3K4me3_rank'] == 2.0) & (df['H3K27ac_rank'] == 1.0)) |
       ((df['H3K4me3_rank'] == 1.0) & (df['H3K4me1_auc'] == 0.0) & (df['H3K27ac_auc'] == 0.0) & (df['CTCF_auc'] == 0.0)), 'RE'] = 'Promoter'

df.loc[((df['H3K4me3_rank'] == 1.0) & (df['CTCF_rank'] == 2.0)) | 
       ((df['H3K4me3_rank'] == 2.0) & (df['CTCF_rank'] == 1.0)) |
       ((df['H3K4me1_auc'] == 0.0) & (df['H3K4me3_auc'] != 0.0) & (df['H3K27ac_auc'] != 0.0) & (df['CTCF_auc'] != 0.0)), 'RE'] = 'Promoter/CTCF'

df.loc[((df['H3K27ac_rank'] == 1.0) & (df['CTCF_rank'] == 2.0) & (df['H3K4me1_auc'] == 0.0) & (df['H3K4me3_auc'] == 0.0)) | 
       ((df['H3K27ac_rank'] == 2.0) & (df['CTCF_rank'] == 1.0) & (df['H3K4me1_auc'] == 0.0) & (df['H3K4me3_auc'] == 0.0)) | 
       ((df['CTCF_rank'] == 1.0) & (df['H3K4me1_auc'] == 0.0) & (df['H3K4me3_auc'] == 0.0) & (df['H3K27ac_auc'] == 0.0) & (df['CTCF_auc'] != 0.0)), 'RE'] = 'CTCF'

df.loc[((df['H3K4me1_auc'] != 0.0) & (df['H3K4me3_auc'] != 0.0) & (df['H3K27ac_auc'] != 0.0) & (df['CTCF_auc'] == 0.0)), 'RE'] = 'TBD_1'
df.loc[((df['H3K4me1_auc'] != 0.0) & (df['H3K4me3_auc'] != 0.0) & (df['H3K27ac_auc'] == 0.0) & (df['CTCF_auc'] != 0.0)), 'RE'] = 'TBD_2'
df.loc[((df['H3K4me1_auc'] != 0.0) & (df['H3K4me3_auc'] != 0.0) & (df['H3K27ac_auc'] != 0.0) & (df['CTCF_auc'] != 0.0)), 'RE'] = 'TBD_3'


df.loc[(df['H3K27ac_rank'] != 1.0), 'Activity_H3K27ac'] = 'Inactive'
df.loc[(df['H3K27ac_rank'] == 1.0), 'Activity_H3K27ac'] = 'Active'


# regions = pd.read_csv(params10, sep='\t', header=None)
# regions.columns = ["chromosome", "start", "end"]
df['id'] = df["chromosome"].astype(str)+":"+df["start"].astype(str)+"_"+df["end"].astype(str)

if params11 != "none":
    
    peak = sort_union

    peak['id'] = peak[0].astype(str)+":"+peak[1].astype(str)+"_"+peak[2].astype(str)

    feat = []
    for ID in df["id"]:
        if ID in peak["id"].tolist():
            feat.append("Active")
        else:
            feat.append("Inactive")
    del peak

    df['Activity_intersect'] = feat



sort_regions = df[["chromosome","start","end"]]
sort_regions['chromosome'] = sort_regions['chromosome'].str.replace('chr', '')
sort_id = np.array(sort_regions['chromosome']+":"+sort_regions['start'].astype(str)+"-"+sort_regions['end'].astype(str)).tolist()
sort_regions = pybedtools.BedTool.from_dataframe(sort_regions).sort()

ctcf_peaks = params5
ctcf_peaks = pd.read_csv(ctcf_peaks, sep='\t')
ctcf_peaks = ctcf_peaks[ctcf_peaks['overall_peak_score'] >= peak_thr][["chrom","start","end"]]
ctcf_peaks.columns = ["chromosome","start","end"]
ctcf_peaks['chromosome'] = ctcf_peaks['chromosome'].str.replace('chr', '')
ctcf_peaks = pybedtools.BedTool.from_dataframe(ctcf_peaks).sort()

tss_sites = params14
tss_sites = 'TSS/TSS_hg38_strict.bed'
tss_sites = pd.read_csv(tss_sites, sep='\t', header=None)[[0,1,2]]
tss_sites.columns = ["chromosome","start","end"]
tss_sites['chromosome'] = tss_sites['chromosome'].str.replace('chr', '')
tss_sites = pybedtools.BedTool.from_dataframe(tss_sites).sort()

promoter_like = sort_regions.intersect(tss_sites, u=True)

not_promoter_like = sort_regions.intersect(tss_sites, wa=True, v=True)
ctcf_like = not_promoter_like.intersect(ctcf_peaks, wa=True, u=True)

enhancer_like = not_promoter_like.intersect(ctcf_peaks, wa=True, v=True)

promoter_like = promoter_like.to_dataframe()
if promoter_like.shape[0] > 0:
    promoter_like['id'] = np.array(promoter_like['chrom'].astype(str)+":"+promoter_like['start'].astype(str)+"-"+promoter_like['end'].astype(str)).tolist()
else:
    promoter_like = pd.DataFrame(columns=['chrom', 'start', 'end', 'id'])
ctcf_like = ctcf_like.to_dataframe()
if ctcf_like.shape[0] > 0:
    ctcf_like['id'] = np.array(ctcf_like['chrom'].astype(str)+":"+ctcf_like['start'].astype(str)+"-"+ctcf_like['end'].astype(str)).tolist()
else:
    ctcf_like = pd.DataFrame(columns=['chrom', 'start', 'end', 'id'])
enhancer_like = enhancer_like.to_dataframe()
if enhancer_like.shape[0] > 0:
    enhancer_like['id'] = np.array(enhancer_like['chrom'].astype(str)+":"+enhancer_like['start'].astype(str)+"-"+enhancer_like['end'].astype(str)).tolist()
else:
    enhancer_like = pd.DataFrame(columns=['chrom', 'start', 'end', 'id'])

TSS_RE = []

for i in sort_id:
    if i in promoter_like['id'].tolist():
        TSS_RE.append('Promoter')
    elif i in ctcf_like['id'].tolist():
        TSS_RE.append('CTCF')
    elif i in enhancer_like['id'].tolist():
        TSS_RE.append('Enhancer')
    else:
        TSS_RE.append('None')

df['TSS_RE'] = TSS_RE

if params11 != "none": 
    df = df[["chromosome", "start", "end", "RE", "TSS_RE", "Activity_H3K27ac", "Activity_intersect", "H3K4me1_rank", "H3K4me3_rank", "H3K27ac_rank", "CTCF_rank", "H3K4me1_auc", "H3K4me3_auc", "H3K27ac_auc", "CTCF_auc"]]
else:
    df = df[["chromosome", "start", "end", "RE", "TSS_RE", "Activity_H3K27ac", "H3K4me1_rank", "H3K4me3_rank", "H3K27ac_rank", "CTCF_rank", "H3K4me1_auc", "H3K4me3_auc", "H3K27ac_auc", "CTCF_auc"]]


shutil.rmtree(params15)

df.to_csv(params13, sep='\t', index=False)







