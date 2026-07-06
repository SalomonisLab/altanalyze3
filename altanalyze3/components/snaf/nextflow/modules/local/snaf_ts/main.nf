// nf-core-style DSL2 module: SNAF tumor-specificity scoring (DB-free, offline)
process SNAF_TS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/altanalyze3:latest' :
        'community.wave.seqera.io/library/altanalyze3:latest' }"

    input:
    tuple val(meta), path(juncounts)
    path control_h5ad

    output:
    tuple val(meta), path("${prefix}/neojunction_tumor_specificity.csv"), emit: ts
    tuple val(meta), path("${prefix}/NeoJunction_statistics_*.txt")      , emit: stats, optional: true
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    altanalyze3 snaf-ts \\
        --juncounts ${juncounts} \\
        --control_h5ad ${control_h5ad} \\
        --output ${prefix} \\
        --cpus ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        altanalyze3: \$(altanalyze3 --version 2>&1 | tr -d '\\n')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/neojunction_tumor_specificity.csv
    touch versions.yml
    """
}
