version 1.0

workflow JunCountWorkflow {
    input {
        File alignment_file
        File alignment_file_index
        Int? threads
        Int? processes
        String? output_prefix = "results"
    }
    call JunCountTask {
        input:
        alignment_file = alignment_file,
        alignment_file_index = alignment_file_index,
        threads = threads,
        processes = processes,
        output_prefix = output_prefix
    }
    output {
        File junction_bed_file = JunCountTask.junction_bed_file
        File stdout_log = JunCountTask.stdout_log
        File stderr_log = JunCountTask.stderr_log
    }
}

task JunCountTask {
    input {
        File alignment_file
        File alignment_file_index
        String? log_level
        Int? threads
        Int? processes
        String? output_prefix
    }
    command {
        cp ${alignment_file} ${basename(alignment_file)}
        cp ${alignment_file_index} ${basename(alignment_file_index)}
        altanalyze3 juncount --bam ${basename(alignment_file)} \
        ${"--threads " + threads} \
        ${"--cpus " + processes} \
        ${"--output " + output_prefix}
    }
    output {
        File junction_bed_file = "${output_prefix}.bed"
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