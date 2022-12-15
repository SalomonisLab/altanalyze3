version 1.0

workflow IntCountWorkflow {
    input {
        File alignment_file
        File alignment_file_index
        File reference_file
        Int? overlap_bp
        String? strandness
        Int? threads
        Int? processes
        String? output_prefix = "results"
    }
    call IntCountTask {
        input:
        alignment_file = alignment_file,
        alignment_file_index = alignment_file_index,
        reference_file = reference_file,
        overlap_bp = overlap_bp,
        strandness = strandness,
        threads = threads,
        processes = processes,
        output_prefix = output_prefix
    }
    output {
        File intron_bed_file = IntCountTask.intron_bed_file
        File stdout_log = IntCountTask.stdout_log
        File stderr_log = IntCountTask.stderr_log
    }
}

task IntCountTask {
    input {
        File alignment_file
        File alignment_file_index
        File reference_file
        Int? overlap_bp
        String? strandness
        Int? threads
        Int? processes
        String? output_prefix
    }
    command {
        cp ${alignment_file} ${basename(alignment_file)}
        cp ${alignment_file_index} ${basename(alignment_file_index)}
        altanalyze3 intcount \
        --bam ${basename(alignment_file)} \
        --ref ${reference_file} \
        ${"--span " + overlap_bp} \
        ${"--strandness " + strandness} \
        ${"--threads " + threads} \
        ${"--cpus " + processes} \
        ${"--output " + output_prefix}
    }
    output {
        File intron_bed_file = "${output_prefix}.bed"
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