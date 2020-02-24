workflow SingleSampleKallistoWorkflow {

    File r1_fastq
    File r2_fastq
    File kallisto_index_path
    Int kallisto_bootstraps

    # Extract the samplename from the fastq filename
    String sample_name = basename(r1_fastq, "_R1.fastq.gz")

    Int disk_size = 100

    call KallistoQuantification {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            kallisto_index_path = kallisto_index_path,
            kallisto_bootstraps = kallisto_bootstraps,
            sample_name = sample_name
    }

    output {
        File kallisto_stdout = KallistoQuantification.kallisto_stdout
        File abundance_h5 = KallistoQuantification.abundance_h5
        File abundance_tsv = KallistoQuantification.abundance_tsv
        File run_info = KallistoQuantification.run_info
    }
}

task KallistoQuantification {
    File r1_fastq
    File r2_fastq
    File kallisto_index_path
    Int kallisto_bootstraps
    String sample_name

    Int threads = 2
    Int disk_size = 30

    String stdout_log = "stdout.log"

    command {
        kallisto quant -i ${kallisto_index_path} \
            -o "kallisto_results_${sample_name}" \
            -b ${kallisto_bootstraps} \
            -t ${threads} \
            ${r1_fastq} \
            ${r2_fastq} >> ${stdout_log} 2>&1

        mv "kallisto_results_${sample_name}/abundance.h5" "${sample_name}.abundance.h5"
        mv "kallisto_results_${sample_name}/abundance.tsv" "${sample_name}.abundance.tsv"
        mv "kallisto_results_${sample_name}/run_info.json" "${sample_name}.run_info.json"
    }

    output {
        File kallisto_stdout = "${stdout_log}"
        File abundance_h5 = "${sample_name}.abundance.h5"
        File abundance_tsv = "${sample_name}.abundance.tsv"
        File run_info = "${sample_name}.run_info.json"
    }

    runtime {
        zones: "us-east4-c"
        docker: "docker.io/blawney/kallisto:v0.0.2"
        cpu: 4
        memory: "5 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}