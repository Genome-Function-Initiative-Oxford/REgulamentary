if config["remove_blacklist"]["genome"] == "hg38":
    blacklist = "blacklists/hg38.bed"
elif config["remove_blacklist"]["genome"] == "mm39":
    blacklist = "blacklists/mm39.bed"
else:
    print("other genomes not implemented yet!")
    sys.exit()

    
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
                cat {output.atac} {output.ctcf} > {params.tmp_bed}
                sort {params.tmp_bed} | uniq -u > {output.merged}                
                rm -rf {params.tmp_bed}
            else
                cp {input.bed1} {output.atac}
                cp {input.bed2} {output.ctcf}
                cat {output.atac} {output.ctcf} > {params.tmp_bed}
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
