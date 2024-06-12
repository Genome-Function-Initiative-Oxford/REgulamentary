import sys
import pandas as pd

import warnings
warnings.filterwarnings('ignore')

inputf = sys.argv[1] #input
outputf = sys.argv[2] #output
thr = sys.argv[3] #thr

df = pd.read_csv(inputf, sep='\t')
df = df[df["overall_peak_score"]>thr]
df = df[["chrom","start","end"]]
df = df[df['chrom'].str.startswith('chr')]

df.to_csv(outputf, index=None, header=None, sep="\t")