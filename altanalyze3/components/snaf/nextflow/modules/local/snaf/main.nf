// nf-core-style DSL2 module: SNAF full MHC-bound T-antigen pipeline
process SNAF {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/altanalyze3:latest' :
        'community.wave.seqera.io/library/altanalyze3:latest' }"

    input:
    tuple val(meta), path(juncounts), path(hla)
    path db_dir
    path genome_fasta

    output:
    tuple val(meta), path("${prefix}/T_candidates/*")    , emit: candidates, optional: true
    tuple val(meta), path("${prefix}/frequency_stage*")  , emit: frequency,  optional: true
    tuple val(meta), path("${prefix}/burden_stage*")     , emit: burden,     optional: true
    tuple val(meta), path("${prefix}/after_prediction.p"), emit: pickle,     optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def genome  = genome_fasta.name != 'NO_FILE' ? "--genome_fasta ${genome_fasta}" : ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    export SNAF_OFFLINE=\${SNAF_OFFLINE:-0}
    altanalyze3 snaf \\
        --juncounts ${juncounts} \\
        --db_dir ${db_dir} \\
        --hla ${hla} \\
        --output ${prefix} \\
        --cpus ${task.cpus} \\
        ${genome} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        altanalyze3: \$(altanalyze3 --version 2>&1 | tr -d '\\n')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/T_candidates
    touch ${prefix}/after_prediction.p versions.yml
    """
}
