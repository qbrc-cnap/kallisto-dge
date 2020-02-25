import "single_sample_kallisto.wdl" as single_sample_kallisto
import "multiqc.wdl" as multiqc
import "fastqc.wdl" as fastqc
import "sleuth.wdl" as sleuth
import "report.wdl" as reporting


workflow KallistoAndSleuthWorkflow{
    # This workflow is a 'super' workflow that parallelizes
    # RNA-seq analysis over multiple samples

    Array[File] r1_files
    Array[File] r2_files
    File sample_annotations
    Array[String] base_conditions
    Array[String] experimental_conditions
    String genome
    File kallisto_index_path
    File transcript_to_gene_mapping
    Int kallisto_bootstraps = 500
    String output_zip_name
    String git_repo_url
    String git_commit_hash

    String versus_sep = "_versus_"
    String normalized_counts_suffix = "tpm.tsv"
    String sleuth_output_suffix = "sleuth_results.tsv"
    String pca_suffix = 'pca.png'
    String top_heatmap_suffix = 'heatmap.png'

    Float qval_threshold = 0.01
    Int max_transcripts = 30 # the maximum number of transcripts to plot


    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)
    Array[Pair[String, String]] contrast_pairs = zip(base_conditions, experimental_conditions)


    scatter(item in fastq_pairs){

        call fastqc.run_fastqc as fastqc_for_read1 {
            input:
                fastq = item.left
        }

        call fastqc.run_fastqc as fastqc_for_read2 {
            input:
                fastq = item.right
        }

        call single_sample_kallisto.SingleSampleKallistoWorkflow as single_sample_process{
            input:
                r1_fastq = item.left,
                r2_fastq = item.right,
                kallisto_index_path = kallisto_index_path,
                kallisto_bootstraps = kallisto_bootstraps
        }
    }

    scatter(item in contrast_pairs){
        call sleuth.sleuth_dge as sleuth_dge {
            input:
                abundance_h5_files = single_sample_process.abundance_h5,
                annotations = sample_annotations,
                transcript_to_gene_mapping = transcript_to_gene_mapping,
                base_group = item.left,
                experimental_group = item.right,
                versus_sep = versus_sep,
                normalized_counts_suffix = normalized_counts_suffix,
                sleuth_output_suffix = sleuth_output_suffix,
                pca_suffix = pca_suffix,
                top_heatmap_suffix = top_heatmap_suffix,
                qval_threshold = qval_threshold,
                max_transcripts = max_transcripts
        }
    }

    call multiqc.create_qc as experimental_qc {
        input:
            kallisto_stdout = single_sample_process.kallisto_stdout,
            r1_fastqc_zips = fastqc_for_read1.fastqc_zip,
            r2_fastqc_zips = fastqc_for_read2.fastqc_zip
    }

    call reporting.generate_report as make_report {
        input:
            r1_files = r1_files,
            r2_files = r2_files,
            annotations = sample_annotations,
            genome = genome,
            sleuth_results = sleuth_dge.sleuth_results,
            git_commit_hash = git_commit_hash,
            git_repo_url = git_repo_url,
            normalized_counts_suffix = normalized_counts_suffix,
            sleuth_output_suffix = sleuth_output_suffix,
            versus_sep = versus_sep,
            qval_threshold = qval_threshold,
            pca_suffix = pca_suffix,
            top_heatmap_suffix = top_heatmap_suffix,
            max_transcripts = max_transcripts,
            num_bootstraps = kallisto_bootstraps
    }

    call zip_results {
        input:
            zip_name = output_zip_name,
            multiqc_report = experimental_qc.report,
            analysis_report = make_report.report,
            sleuth_outputs = sleuth_dge.sleuth_results,
            normalized_counts_files = sleuth_dge.norm_counts,
            gene_level_count_files = sleuth_dge.gene_level_count_file,
            deseq_results = sleuth_dge.deseq_results,
            figures = sleuth_dge.sleuth_plots
    }

    output {
        File zip_out = zip_results.zip_out
    }

    meta {
        workflow_title : "Kallisto + Sleuth RNA-Seq differential expression"
        workflow_short_description : "For determining differential expression of transcript abundances using Kallisto and Sleuth"
        workflow_long_description : "Use this workflow for performing pseudo-alignments with Kallisto, which produces estimated transcript abundances.  Differential expression is performed using the companion Sleuth tool."
    }
}


task zip_results {

    String zip_name 
    File multiqc_report
    File analysis_report
    Array[File] sleuth_outputs
    Array[File] normalized_counts_files
    Array[Array[File]] figures

    Array[File] contrast_figure_list = flatten(figures)

    Int disk_size = 100

    command {

        mkdir report
        mkdir report/qc
        mkdir report/differential_expression

        mv ${multiqc_report} report/qc/

        python3 /opt/software/organize_report.py -b report/differential_expression ${sep=" " sleuth_outputs}
        python3 /opt/software/organize_report.py -b report/differential_expression ${sep=" " normalized_counts_files}
        python3 /opt/software/organize_report.py -b report/differential_expression ${sep=" " contrast_figure_list}

        mv ${analysis_report} report/
        zip -r "${zip_name}.zip" report
    }

    output {
        File zip_out = "${zip_name}.zip"
    }

    runtime {
        docker: "docker.io/blawney/kallisto:v0.0.2"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}