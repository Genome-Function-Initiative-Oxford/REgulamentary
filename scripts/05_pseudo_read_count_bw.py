import sys
import pyBigWig
from os.path import exists,basename
from os import mkdir
from collections import OrderedDict

#the folders e.g. ATAC, CTCF, merged
outputs= [sys.argv[1],sys.argv[2],sys.argv[3]]

#needs to be ordered so coverages are written to read_count.csv in correct order
data= OrderedDict({
    "H3K4me1":{"file":sys.argv[4]},
    "H3K4me3":{"file":sys.argv[5]},
    "H3K27ac":{"file":sys.argv[6]},
    "CTCF":{"file":sys.argv[7]}
})

#simply caclculates 
for mark in data:
    #is there a more efficient way of calculating coverage
    bw = pyBigWig.open(data[mark]["file"])
    total_coverage=0
    for chrom in bw.chroms():
        length = bw.chroms(chrom)
        values = bw.values(chrom, 0, length)
        tot = int(sum(values)/1000)
        total_coverage += tot
    data[mark]["coverage"]=total_coverage
    bw.close()

#write out the coverage file in each folder
for f in outputs:
    read_count_folder = basename(f)
    if not exists(read_count_folder):
        mkdir(read_count_folder)
    with open(f,"w") as o:
        for mark in data:
            o.write(f"{data[mark]['coverage']}\n")



