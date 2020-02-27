task sleuth_dge {
    Array[File] abundance_h5_files
    File annotations
    File transcript_to_gene_mapping
    String base_group
    String experimental_group
    String versus_sep
    String normalized_counts_suffix
    String sleuth_output_suffix
    String pca_suffix
    String top_heatmap_suffix
    Float qval_threshold
    Int max_transcripts
    Boolean is_pdx


    Int disk_size = 100
        
    String contrast_name = experimental_group + versus_sep + base_group
    String sleuth_output = contrast_name + "." + sleuth_output_suffix
    String sleuth_normalized_counts = contrast_name + "." + normalized_counts_suffix
    String output_figures_dir = contrast_name + "_figures"
    String pca_output = contrast_name + "." + pca_suffix
    String heatmap_output = contrast_name + "." + top_heatmap_suffix
    String normalized_counts = "nc.tsv"
    String sleuth_annotations = 'sleuth_annotations.tsv'
    String gene_level_counts_suffix = ".gene_level_counts.tsv"
    String gene_level_counts = contrast_name + gene_level_counts_suffix
    String deseq_results_suffix = ".gene_level_dge_results.tsv"
    String deseq_results = contrast_name + deseq_results_suffix

    # if we are running PDx, we will also generate files for human and mouse separately.
    # these will be named similarly, e.g. 'A_vs_B.mouse.tpm.tsv' and 'A_vs_B.human.tpm.tsv'
    # and so on.
    String human_tag = "human"
    String mouse_tag = "mouse"

    command {
        echo "Moving abundance.h5 files..."

        /usr/bin/python3 /opt/software/move_files.py ${sep=" " abundance_h5_files}

        echo "Completed moving files..."
        /usr/bin/python3 /opt/software/create_sleuth_annotation_file.py -i ${annotations} -o ${sleuth_annotations}
        echo "Created annotations for sleuth process..."

        # activate the conda environment.  Otherwise R will not be visible.
        source activate r36
        
        Rscript /opt/software/sleuth.R \
            ${sleuth_annotations} \
            ${transcript_to_gene_mapping} \
            ${sleuth_output} \
            ${normalized_counts} \
            ${output_figures_dir} \
            ${max_transcripts} \
            ${qval_threshold} \
            ${base_group} \
            ${experimental_group} \
            ${pca_output} \
            ${heatmap_output} \
            ${contrast_name}

        echo "Completed sleuth"

        /usr/bin/python3 /opt/software/pivot_normalized_counts.py \
            -i ${normalized_counts} \
            -o ${sleuth_normalized_counts}

        echo "Run tximport..."
        Rscript /opt/software/run_tximport.R \
            ${sleuth_annotations} \
            ${transcript_to_gene_mapping} \
            ${base_group} \
            ${experimental_group} \
            ${deseq_results} \
            ${gene_level_counts}

        echo "Done with tximport and DESeq2"

        if [ "${is_pdx}" = "true" ]
        then
            export HUMAN_SLEUTH_RESULTS = "${contrast_name}.${human_tag}.${sleuth_output_suffix}"
            export MOUSE_SLEUTH_RESULTS = "${contrast_name}.${mouse_tag}.${sleuth_output_suffix}"
            export HUMAN_SLEUTH_TPM = "${contrast_name}.${human_tag}.${normalized_counts_suffix}"
            export MOUSE_SLEUTH_TPM = "${contrast_name}.${mouse_tag}.${normalized_counts_suffix}"

            # send the header lines to the output files:
            head -1 ${sleuth_output} > $HUMAN_SLEUTH_RESULTS
            head -1 ${sleuth_output} > $MOUSE_SLEUTH_RESULTS
            head -1 ${sleuth_normalized_counts} > $HUMAN_SLEUTH_TPM
            head -1 ${sleuth_normalized_counts} > $MOUSE_SLEUTH_TPM

            # split human and mouse based on the ENST transcript
            grep -P "^ENST" ${sleuth_output} >> $HUMAN_SLEUTH_RESULTS
            grep -vP "^ENST" ${sleuth_output} >> $MOUSE_SLEUTH_RESULTS 
            grep -P "^ENST" ${sleuth_normalized_counts} >> $HUMAN_SLEUTH_TPM
            grep -vP "^ENST" ${sleuth_normalized_counts} >> $MOUSE_SLEUTH_TPM

            # run a python script to split at the gene-name level, since we can't easily grep:
            # Creates files like "A_vs_B.human.gene_level_counts.tsv"
            /usr/bin/python3 /opt/software/split_pdx_genes.py \
                -d ${deseq_results} \
                -c ${gene_level_counts} \
                -t ${transcript_to_gene_mapping} \
                -ht ${human_tag} \
                -mt ${mouse_tag}
        fi
    }

    output {
        File norm_counts = "${sleuth_normalized_counts}"
        File sleuth_results = "${sleuth_output}"
        File gene_level_count_file = "${gene_level_counts}"
        File deseq_results_file = "${deseq_results}"
        Array[File] sleuth_plots = glob("${output_figures_dir}/*.png")
        Array[File] human_counts = glob("${contrast_name}.${human_tag}${gene_level_counts_suffix}")
        Array[File] human_deseq = glob("${contrast_name}.${human_tag}${deseq_results_suffix}")
        Array[File] mouse_counts = glob("${contrast_name}.${mouse_tag}${gene_level_counts_suffix}")
        Array[File] mouse_deseq = glob("${contrast_name}.${mouse_tag}${deseq_results_suffix}")    }

    runtime {
        docker: "docker.io/blawney/kallisto:v0.0.2"
        cpu: 4
        memory: "5 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
