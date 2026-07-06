// nf-core-style DSL2 module: SNAF surface / B-antigen pipeline
// Pure-python (tmhmm.py + Biopython), offline-first. Consumes the T-antigen frequency
// table produced by the SNAF (full) module.
process SNAF_B {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/altanalyze3:latest' :
        'community.wave.seqera.io/library/altanalyze3:latest' }"

    input:
    tuple val(meta), path(juncounts), path(freq)
    path db_dir
    path validation_gtf

    output:
    tuple val(meta), path("${prefix}/B_candidates/*"), emit: candidates, optional: true
    tuple val(meta), path("${prefix}/surface_antigen_*.p"), emit: pickle, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def gtf    = validation_gtf.name != 'NO_FILE' ? "--validation_gtf ${validation_gtf}" : ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    export SNAF_OFFLINE=\${SNAF_OFFLINE:-0}
    altanalyze3 snaf-b \\
        --juncounts ${juncounts} \\
        --db_dir ${db_dir} \\
        --freq_path ${freq} \\
        --output ${prefix} \\
        --cpus ${task.cpus} \\
        ${gtf} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        altanalyze3: \$(altanalyze3 --version 2>&1 | tr -d '\\n')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/B_candidates
    touch versions.yml
    """
}
