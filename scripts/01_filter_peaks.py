import sys
import pandas as pd

import warnings
warnings.filterwarnings('ignore')

params1 = sys.argv[1] #input
params2 = sys.argv[2] #output

df = pd.read_csv(params1, sep='\t')
df = df[df["overall_peak_score"]>0.95]
df = df[["chrom","start","end"]]

df.to_csv(params2, index=None, header=None, sep="\t")