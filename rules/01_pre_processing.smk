# cat {output.atac} {output.ctcf} > {params.tmp_bed}

from os.path import join,splitext,basename
from os import makedirs
import tempfile
import shutil

if config["remove_blacklist"]["genome"] == "hg38":
    blacklist = "blacklists/hg38.bed"
elif config["remove_blacklist"]["genome"] == "mm39":
    blacklist = "blacklists/mm39.bed"
else:
    print("other genomes not implemented yet!")
    sys.exit()

use_defined_bw  = False
create_bed_files = False
#bigwigs are defined in one place and used for all analysis
if config.get("bigwigs"):
    use_defined_bw=True
    # peaks are called from the bed files
    if  config["bigwigs"].get("create_bed_files"):
        #create a folder in the output folder to hold the peak files
        bed_dir  = join(config["analysis_name"],"bed_files")
        create_bed_files=True
        #update the config to reflect the (to be) created peak files
        for m in ["ATAC","H3K4me1","H3K4me3","H3K27ac","CTCF"]:
            config["union_peaks"][f"bed_{m}"]= join(bed_dir,f"{m}.bed")


if create_bed_files:
    #the directory to put the bed files
    bed_dir = join(config["analysis_name"],"bed_files")
    makedirs(bed_dir,exist_ok=True)

    #this rule will call peaks using lanceotron from the supplied bigwigs
    rule create_bed_files:
        input:
            bw = lambda wildcards: config["bigwigs"][wildcards.seqtype]
        output:
            bed =join(bed_dir,"{seqtype}.bed")
        run:
           with tempfile.TemporaryDirectory() as tmpdir:
                shell(f"lanceotron callPeaks {input.bw} -f {tmpdir}")
                # Extract base name and construct new file name
                base = basename(input.bw)
                # Remove extension and add '_L-tron.bed'
                stem = splitext(base)[0]
                ltron_bed = f"{stem}_L-tron.bed"
                temp_bed_path = join(tmpdir, ltron_bed)
                #copy temporary lanceotron output to the bed directory
                shutil.copy(temp_bed_path, output.bed)
          

rule union_peaks:
    input:
        bed1=config["union_peaks"]["bed_ATAC"],
        bed2=config["union_peaks"]["bed_CTCF"],
    output:
        atac=config["analysis_name"]+os.sep+"ATAC/01_union_peaks/union_peaks.bed",
        ctcf=config["analysis_name"]+os.sep+"CTCF/01_union_peaks/union_peaks.bed",
        merged=config["analysis_name"]+os.sep+"merge/01_union_peaks/union_peaks.bed",
    params:
        peak_caller=config["union_peaks"]["peak_caller"],
        peak_caller_thr=config["union_peaks"]["threshold"], 
        tmp_bed=config["analysis_name"]+os.sep+"merge/01_union_peaks/tmp.bed"
    shell:
        """ 
            if [ {params.peak_caller} == "lanceotron" ]
            then
                python scripts/01_filter_peaks.py {input.bed1} {output.atac} {params.peak_caller_thr}
                python scripts/01_filter_peaks.py {input.bed2} {output.ctcf} {params.peak_caller_thr}
                cat {output.atac} {output.ctcf} | cut -f 1,2,3 | bedtools sort | bedtools merge > {params.tmp_bed}
                sort {params.tmp_bed} | uniq -u > {output.merged}  
                rm -rf {params.tmp_bed}
            else
                cp {input.bed1} {output.atac}
                cp {input.bed2} {output.ctcf}
                cat {output.atac} {output.ctcf} | cut -f 1,2,3 | bedtools sort | bedtools merge > {params.tmp_bed}
                sort {params.tmp_bed} | uniq -u > {output.merged}                
                rm -rf {params.tmp_bed}
            fi
        """
  
   
rule remove_blacklist:
    input:
        bed=config["analysis_name"]+os.sep+"{folder}/01_union_peaks/union_peaks.bed",
        
    output:
        config["analysis_name"]+os.sep+"{folder}/02_blacklist_removed/union_peaks.bed",
    params:
        tmp=config["analysis_name"]+os.sep+"{folder}/02_blacklist_removed/tmp.bed",
        blacklist=blacklist,
    shell:
        """       
            bedtools intersect -a {input.bed} -b {params.blacklist} -v > {params.tmp}
            cut -f1,2,3 -d$'\t' {params.tmp} > {output}
            rm -rf {params.tmp}
        """

if config.get("bigwigs"):
    # get multi the multicoverage from the bigwigs
    rule multicoverages_from_bigwigs:
        input:
            bed=config["analysis_name"]+os.sep+"{folder}/02_blacklist_removed/union_peaks.bed",
            bw1=config["bigwigs"]["H3K4me1"],
            bw2=config["bigwigs"]["H3K4me3"],
        output:
            config["analysis_name"]+os.sep+"{folder}/03_multicov/multicov.bed",
        params:
            tmp=config["analysis_name"]+os.sep+"{folder}/03_multicov/tmp.npz"
        shell:
            """
                #don't need the default output -o but will not run if not specified
                multiBigwigSummary BED-file -b {input.bw1} {input.bw2} --labels H3K4me1 H3K4me3 --BED {input.bed}  --outRawCounts {output} -o {params.tmp}
                sed -i '1d' {output} #remove header
                rm {params.tmp}
            """



else:
    rule multicoverages:
        input:
            bed=config["analysis_name"]+os.sep+"{folder}/02_blacklist_removed/union_peaks.bed",
            bam1=config["multicoverages"]["bam_H3K4me1"],
            bam2=config["multicoverages"]["bam_H3K4me3"],
        output:
            config["analysis_name"]+os.sep+"{folder}/03_multicov/multicov.bed",        
        shell:
            """         
                bedtools multicov -bams {input.bam1} {input.bam2} -bed {input.bed} > {output}
            """


rule sort_regions:
    input:
        config["analysis_name"]+os.sep+"{folder}/03_multicov/multicov.bed", 
    output:
        config["analysis_name"]+os.sep+"{folder}/04_sort_regions/sort_union.bed",
    run:
        import pandas as pd
        cov = pd.read_csv(input[0], sep='\t', header=None)
        cov.columns = ['chr', 'start', 'end', 'H3K4me1', 'H3K4me3']
        cov['diff'] = cov['H3K4me1']-cov['H3K4me3']
        cov = cov.sort_values(by=['diff'], ascending=False).dropna()
        cov["start"] = cov["start"].astype(int)
        cov["end"]   = cov["end"].astype(int)
        cov = cov[['chr', 'start', 'end']]
        cov.to_csv(output[0], sep='\t', header=None, index=False)
