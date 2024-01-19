import os, sys
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("7.0.0")

configfile: "config/analysis.yaml"

validate(config, schema="config/schema/config.schema.yaml")

include: "rules/01_pre_processing.smk"
include: "rules/02_deeptools_helper.smk"
include: "rules/03_regulatory_elements.smk"

bw = config["compute_matrix_bigwigs_extra"]["bigwig_extra"].split("/")[-1].replace(".bw", "")

rule all:
    input:
        expand(config["analysis_name"]+os.sep+"{folder}/01_union_peaks/union_peaks.bed", 
               folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/02_blacklist_removed/union_peaks.bed", 
               folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/03_multicov/multicov.bed", 
               folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/04_sort_regions/sort_union.bed", 
               folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/05_compute_matrix/{mtx}", 
               folder=["ATAC", "CTCF", "merge"], 
               mtx=["matrix.mtx", "readable_matrix.mtx", "sorted_regions.mtx"]),
        expand(config["analysis_name"]+os.sep+"{folder}/06_active_elements/lanceotron.txt", 
               folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/06_active_elements/%s_L-tron_filtered.bed"%bw,
               folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/06_active_elements/sorted_regions.bed",
               folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/06_active_elements/%s_L-tron_filtered_intersection.bed"%bw,
               folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/07_read_count/read_count.csv", 
               folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/08_REgulamentary/mlv_REgulamentary.csv", 
               folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/09_metaplot/matrix.csv", 
                folder=["ATAC", "CTCF", "merge"]),
        expand(config["analysis_name"]+os.sep+"{folder}/done.txt",
               folder=["ATAC", "CTCF", "merge"])