def sq(x) { "\'" + x.toString().replace("\'", "\'\\\'\'") + "\'" }
// nf-core-style DSL2 module: SNAF surface / B-antigen pipeline
// Pure-python (tmhmm.py + Biopython), offline-first. Consumes the T-antigen frequency
// table produced by the SNAF (full) module.
process SNAF_B {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(juncounts, stageAs: 'counts/*'), path(freq, stageAs: 'frequency/*')
    path db_dir
    path validation_gtf
    path genome_fasta
    path isoform_ranker
    path first_exon_junctions
    val isoform_method

    output:
    tuple val(meta), path("${prefix}/B_candidates/*"), emit: candidates, optional: true
    tuple val(meta), path("${prefix}/surface_antigen_*.p"), emit: pickle, optional: true
    tuple val(meta), path("${prefix}/synthetic_isoforms/*"), emit: synthetic_isoforms, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def gtf    = validation_gtf ? "--validation_gtf ${sq(validation_gtf)}" : ''
    def genome = genome_fasta ? "--genome_fasta ${sq(genome_fasta)}" : ''
    def ranker = isoform_ranker ? "--isoform_ranker ${sq(isoform_ranker)}" : ''
    def firstExons = first_exon_junctions ? "--first_exon_junctions ${sq(first_exon_junctions)}" : ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    export SNAF_OFFLINE=\${SNAF_OFFLINE:-0}
    altanalyze3 snaf-b \\
        --juncounts ${sq(juncounts)} \\
        --db_dir ${sq(db_dir)} \\
        --freq_path ${sq(freq)} \\
        --output ${prefix} \\
        --cpus ${task.cpus} \\
        --isoform_method ${isoform_method} \\
        ${genome} ${ranker} ${firstExons} \\
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
