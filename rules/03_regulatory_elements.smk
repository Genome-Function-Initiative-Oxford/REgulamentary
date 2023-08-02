bw = config["compute_matrix_bigwigs_2"]["bigwig"].split("/")[-1].replace(".bw", "")

if config["remove_blacklist"]["genome"] == "hg38":
    tss = "TSS/TSS_hg38_strict.bed"
elif config["remove_blacklist"]["genome"] == "mm39":
    blacklist = "TSS/TSS_mm39_strict.bed"
else:
    print("other genomes not implemented yet!")
    sys.exit()

rule active_elements:
    input:
    output:
        config["analysis_name"]+os.sep+"{folder}/06_active_elements/lanceotron.txt",
    params:
        bw=config["compute_matrix_bigwigs_2"]["bigwig"],
        peak=config["compute_matrix_bigwigs_2"]["peak"],
        peak_folder=config["analysis_name"]+os.sep+"{folder}/06_active_elements",
        peak_folder_tmp=config["analysis_name"]+os.sep+"tmp",
    log:
    shell:
        """
            if [ {params.bw} == "none" ]
            then
                echo skipping laceotron for bigwig file missing >> {output}
            else
                mkdir -p {params.peak_folder}
                mkdir -p {params.peak_folder_tmp}
                if [ {params.bw} == "none" ]
                then
                    lanceotron callPeaks {params.bw} -f {params.peak_folder_tmp}
                    cp -r {params.peak_folder_tmp}/* {params.peak_folder}/
                    echo laceotron applied on extra bigwig Done! >> {output}
                else
                    cp {params.peak} {params.peak_folder}/
                    echo laceotron not applied because peaks already provided! >> {output}
                fi
            fi
        """

        
rule filter_active_elements:
    input:
        config["analysis_name"]+os.sep+"{folder}/06_active_elements/lanceotron.txt",
    output:
        config["analysis_name"]+os.sep+"{folder}/06_active_elements/%s_L-tron_filtered.bed"%bw,
    params:
        bw=config["analysis_name"]+os.sep+"{folder}/06_active_elements/%s_L-tron.bed"%bw,
        threshold=config["compute_matrix_bigwigs_2"]["threshold"],
    run:
        import glob
        import pandas as pd
        if config["compute_matrix_bigwigs_2"]["bigwig"] == "none":
            f = open(output[0], "w")
            f.write("skipping filter_active_elements for bigwig file missing")
            f.close()
        else:
            input_f = glob.glob(params.bw)[0]        
            peak = pd.read_csv(input_f, sep='\t')
            peak = peak[peak["overall_peak_score"]>=params.threshold]
            peak.to_csv(output[0], index=False, sep='\t')
        
    
rule clean_sorted_regions: #maybe to remove!?
    input:
        config["analysis_name"]+os.sep+"{folder}/05_compute_matrix/sorted_regions.mtx",
    output:
        config["analysis_name"]+os.sep+"{folder}/06_active_elements/sorted_regions.bed",
    params:
        bw=config["compute_matrix_bigwigs_2"]["bigwig"],
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t")
        df = df[["#chrom", "start", "end"]]
        df.to_csv(output[0], index=False, sep='\t', header=False)

        
        
rule intersect_active_elements:
    input:
        #bed1=config["analysis_name"]+os.sep+"{folder}/06_active_elements/sorted_regions.bed",
        bed1=config["analysis_name"]+os.sep+"{folder}/04_sort_regions/sort_union.bed",
        bed2=config["analysis_name"]+os.sep+"{folder}/06_active_elements/%s_L-tron_filtered.bed"%bw,
    output:
        config["analysis_name"]+os.sep+"{folder}/06_active_elements/%s_L-tron_filtered_intersection.bed"%bw,
    params:
        bw=config["compute_matrix_bigwigs_2"]["bigwig"],
    shell:
        """
            if [ {params.bw} == "none" ]
            then
                echo skipping intersect_active_elements for bigwig file missing >> {output}
            else
                bedtools intersect -a {input.bed1} -b {input.bed2} -wa > {output}
            fi  
        """
        

        
rule plot_regulatory_elements:
    input:
        mtx=config["analysis_name"]+os.sep+"{folder}/05_compute_matrix/readable_matrix.mtx",
        bed1=config["analysis_name"]+os.sep+"{folder}/06_active_elements/sorted_regions.bed",
        #bed1=config["analysis_name"]+os.sep+"{folder}/04_sort_regions/sort_union.bed",
        bed2=config["analysis_name"]+os.sep+"{folder}/06_active_elements/%s_L-tron_filtered_intersection.bed"%bw,
    output:
        config["analysis_name"]+os.sep+"{folder}/07_plot_RE/mlv_deeptools.csv",
    params:
        bw=config["compute_matrix_bigwigs_2"]["bigwig"],
        range_definition=config["cluster_regulatory_elements"]["range_definition"],
        tag_extra_data_input=config["compute_matrix_bigwigs_2"]["tag_extra_data_input"],
        maxThreshold=config["cluster_regulatory_elements"]["maxThreshold"],
        plotTitle=config["cluster_regulatory_elements"]["plotTitle"],
        folder1=config["analysis_name"]+os.sep+"{folder}/07_plot_RE/01_RE_only",
        folder2=config["analysis_name"]+os.sep+"{folder}/07_plot_RE/02_activity_rule",
        folder3=config["analysis_name"]+os.sep+"{folder}/07_plot_RE/03_activity_extra",
    shell:
        """ 
            python scripts/02_plot_regulatory_elements.py \
                {input.mtx} \
                {params.bw} \
                {params.range_definition} \
                {params.tag_extra_data_input} \
                {params.maxThreshold} \
                {params.plotTitle} \
                {params.folder1} \
                {params.folder2} \
                {params.folder3} \
                {input.bed1} \
                {input.bed2} \
                {output}
        """


rule mlv_regulatory_elements:
    input:
        #sort_union = config["analysis_name"]+os.sep+"{folder}/06_active_elements/sorted_regions.bed",
        sort_union = config["analysis_name"]+os.sep+"{folder}/04_sort_regions/sort_union.bed",
        H3K4me1_peaks = config["union_peaks"]["H3K4me1"],
        H3K4me3_peaks = config["union_peaks"]["H3K4me3"],
        H3K27ac_peaks = config["union_peaks"]["H3K27ac"],
        CTCF_peaks = config["union_peaks"]["CTCF"],
        H3K4me1_bw = config["compute_matrix_bigwigs"]["bigwig_H3K4me1"],
        H3K4me3_bw = config["compute_matrix_bigwigs"]["bigwig_H3K4me3"],
        H3K27ac_bw = config["compute_matrix_bigwigs"]["bigwig_H3K27ac"],
        CTCF_bw = config["compute_matrix_bigwigs"]["bigwig_CTCF"],
    output:
        config["analysis_name"]+os.sep+"{folder}/08_REgulamentary/mlv_REgulamentary.csv",
    params:
        extra_peaks = config["analysis_name"]+os.sep+"{folder}/06_active_elements/%s_L-tron_filtered_intersection.bed"%bw,
        extra_bw = config["compute_matrix_bigwigs_2"]["bigwig"],
        thresholdPeaks = config["thresholdPeaks"],
        tss = tss,
        tmp_pybedtools = config["analysis_name"]+os.sep+"{folder}/tmp_pybedtools",
    shell:
        """ 
            python scripts/03_mlv_regulatory_elements.py \
                {input.sort_union} \
                {input.H3K4me1_peaks} \
                {input.H3K4me3_peaks} \
                {input.H3K27ac_peaks} \
                {input.CTCF_peaks} \
                {input.H3K4me1_bw} \
                {input.H3K4me3_bw} \
                {input.H3K27ac_bw} \
                {input.CTCF_bw} \
                {params.extra_peaks} \
                {params.extra_bw} \
                {params.thresholdPeaks} \
                {output} \
                {params.tss} \
                {params.tmp_pybedtools}
        """
        
        
rule clean:
    input:
        i1 = config["analysis_name"]+os.sep+"ATAC/08_REgulamentary/mlv_REgulamentary.csv",
        i2 = config["analysis_name"]+os.sep+"CTCF/08_REgulamentary/mlv_REgulamentary.csv",
        i3 = config["analysis_name"]+os.sep+"merge/08_REgulamentary/mlv_REgulamentary.csv",
    output:
        config["analysis_name"]+os.sep+"{folder}/08_REgulamentary/clean.txt",
    params:
        rm_tmp=config["analysis_name"]+os.sep+"tmp"
    shell:
        """
            echo Next step: Loading mlv.csv into MLV and moving into datashare bigwig files >> {output}
            rm -rf {params.rm_tmp}
        """