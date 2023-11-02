rule compute_matrices:
    input:
        bed=config["analysis_name"]+os.sep+"{folder}/04_sort_regions/sort_union.bed",
        bw1=config["compute_matrix_bigwigs"]["bigwig_H3K4me1"],
        bw2=config["compute_matrix_bigwigs"]["bigwig_H3K4me3"],
        bw3=config["compute_matrix_bigwigs"]["bigwig_H3K27ac"],
        bw4=config["compute_matrix_bigwigs"]["bigwig_CTCF"],
    output:
        mtx1=config["analysis_name"]+os.sep+"{folder}/05_compute_matrix/matrix.mtx",
        mtx2=config["analysis_name"]+os.sep+"{folder}/05_compute_matrix/readable_matrix.mtx",
        mtx3=config["analysis_name"]+os.sep+"{folder}/05_compute_matrix/sorted_regions.mtx",
    params:
        bw=config["compute_matrix_bigwigs_extra"]["bigwig_extra"],
        referencePoint=config["compute_matrix_extra"]["referencePoint"],
        b=config["compute_matrix_extra"]["b"],
        a=config["compute_matrix_extra"]["a"],
        sortRegions=config["compute_matrix_extra"]["sortRegions"],
        # maxThreshold=config["compute_matrix_extra"]["maxThreshold"],
    shell:
        """ 
            if [ {params.bw} == "none" ]
            then
                computeMatrix reference-point --regionsFileName {input.bed} \
                    --scoreFileName {input.bw1} {input.bw2} {input.bw3} {input.bw4} \
                    --outFileName {output.mtx1} \
                    --outFileNameMatrix {output.mtx2} \
                    --outFileSortedRegions {output.mtx3} \
                    --referencePoint {params.referencePoint} \
                    -b {params.b} \
                    -a {params.a} \
                    --sortRegions {params.sortRegions} \
                    --numberOfProcessors max \
                    --missingDataAsZero
            else
                computeMatrix reference-point --regionsFileName {input.bed} \
                    --scoreFileName {input.bw1} {input.bw2} {input.bw3} {input.bw4} {params.bw} \
                    --outFileName {output.mtx1} \
                    --outFileNameMatrix {output.mtx2} \
                    --outFileSortedRegions {output.mtx3} \
                    --referencePoint {params.referencePoint} \
                    -b {params.b} \
                    -a {params.a} \
                    --sortRegions {params.sortRegions} \
                    --numberOfProcessors max \
                    --missingDataAsZero
            fi
        """
