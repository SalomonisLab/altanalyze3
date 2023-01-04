version 1.0

workflow AggregateWorkflow {
    input {
        Array[File] junction_bed_file
        Array[File] intron_bed_file
        File reference_file
        Int? threads
        Int? processes
        String? output_prefix = "results"
    }
    call AggregateTask {
        input:
        junction_bed_file = junction_bed_file,
        intron_bed_file = intron_bed_file,
        reference_file = reference_file,
        threads = threads,
        processes = processes,
        output_prefix = output_prefix
    }
    output {
        File aggregated_bed_file = AggregateTask.aggregated_bed_file
        File aggregated_bed_tbi_file = AggregateTask.aggregated_bed_tbi_file
        File aggregated_h5ad_file = AggregateTask.aggregated_h5ad_file
        File stdout_log = AggregateTask.stdout_log
        File stderr_log = AggregateTask.stderr_log
    }
}

task AggregateTask {
    input {
        Array[File] junction_bed_file
        Array[File] intron_bed_file
        File reference_file
        Int? threads
        Int? processes
        String? output_prefix
    }
    command {
        altanalyze3 aggregate \
        --juncounts ${sep=" " junction_bed_file} \
        --intcounts ${sep=" " intron_bed_file} \
        --ref ${reference_file} \
        --bed \
        ${"--threads " + threads} \
        ${"--cpus " + processes} \
        ${"--output " + output_prefix}
    }
    output {
        File aggregated_bed_file = "${output_prefix}.bed.gz"
        File aggregated_bed_tbi_file = "${output_prefix}.bed.gz.tbi"
        File aggregated_h5ad_file = "${output_prefix}.h5ad"
        File stdout_log = stdout()
        File stderr_log = stderr()
    }
    runtime {
        docker: "altanalyze:latest"
    }
    meta {
        author: "Michael Kotliar"
    }
}