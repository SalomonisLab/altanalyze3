// Core stages use AltAnalyze3. pyNeoQuant runs in its own optional environment.
def sq(x) { "'" + x.toString().replace("'", "'\\''") + "'" }

process INDEX_REFERENCE {
    label 'process_medium'
    input: path gtf
    output: path 'reference', emit: reference
    script:
    """
    python3 -m altanalyze3.components.neoantigen.stages index --gtf ${sq(gtf)} --output reference
    """
    stub:
    """
    mkdir reference
    touch reference/gene_model.bed.gz reference/gene_model.bed.gz.tbi reference/gene_model_all.tsv
    """
}

process COUNT_JUNCTIONS {
    tag "$id"
    label 'process_medium'
    input: tuple val(id), path(bam), path(bai)
    output: tuple val(id), path("${id}_juncounts.bed"), emit: counts
    script:
    """
    python3 -m altanalyze3.components.neoantigen.stages junctions --bam ${sq(bam)} --cpus ${task.cpus} --output '${id}_juncounts'
    """
    stub: "touch '${id}_juncounts.bed'"
}

process COUNT_INTRONS {
    tag "$id"
    label 'process_medium'
    input:
    tuple val(id), path(bam), path(bai)
    path reference
    output: tuple val(id), path("${id}_intcounts.bed"), emit: counts
    script:
    """
    python3 -m altanalyze3.components.neoantigen.stages introns --bam ${sq(bam)} --reference ${sq(reference + '/gene_model.bed.gz')} --strandness '${params.strandness}' --cpus ${task.cpus} --output '${id}_intcounts'
    """
    stub: "touch '${id}_intcounts.bed'"
}

process AGGREGATE_COUNTS {
    label 'process_medium'
    input:
    path junctions
    path introns
    path reference
    output: path 'cohort_annotated.tsv', emit: counts
    script:
    """
    python3 -m altanalyze3.components.neoantigen.stages aggregate --junctions ${junctions.collect{sq(it)}.join(' ')} --introns ${introns.collect{sq(it)}.join(' ')} --reference ${sq(reference + '/gene_model.bed.gz')} --novel-gene-mode '${params.novel_gene_mode}' --output cohort
    """
    stub: "printf 'uid\\tS1\\nENSG1:E1.1-E2.1\\t20\\n' > cohort_annotated.tsv"
}

process FILTER_COUNTS {
    label 'process_medium'
    input: tuple val(meta), path(counts, stageAs: 'input_counts.tsv')
    output: tuple val(meta), path('filtered.tsv'), emit: counts
    script:
    """
    python3 -m altanalyze3.components.neoantigen filter --counts ${sq(counts)} --output filtered.tsv --min-reads ${params.min_reads}
    """
    stub: "cp ${sq(counts)} filtered.tsv"
}

process PREPARE_HLA {
    tag "$id"
    label 'process_medium'
    input:
    tuple val(id), path(bam), path(bai)
    path supplied
    output:
    tuple val(id), path("${id}.hla.tsv"), emit: calls
    path "${id}.hla_qc.tsv", emit: qc
    script:
    def override = supplied ? "--supplied ${sq(supplied)}" : ''
    def complete = params.hla_require_all ? '--require-all' : ''
    """
    python3 -m altanalyze3.components.neoantigen hla --sample '${id}' --bam ${sq(bam)} ${override} --build '${params.genome_build}' --min-depth ${params.hla_min_depth} ${complete} --output '${id}.hla.tsv' --qc '${id}.hla_qc.tsv'
    """
    stub:
    """
    printf 'sample\\thla\\n${id}\\tHLA-A*02:01\\n' > '${id}.hla.tsv'
    printf 'sample_id\\tstatus\\n${id}\\tstub\\n' > '${id}.hla_qc.tsv'
    """
}

process COMBINE_HLA {
    input: path calls
    output: path 'cohort.hla.tsv', emit: hla
    script:
    """
    python3 -m altanalyze3.components.neoantigen combine-hla --inputs ${calls.collect{sq(it)}.join(' ')} --output cohort.hla.tsv
    """
    stub: "printf 'sample\\thla\\nS1\\tHLA-A*02:01\\n' > cohort.hla.tsv"
}

process PYNEOQUANT {
    tag "$sample"
    label 'process_medium'
    container params.pyneoquant_container
    input:
    tuple val(sample), path(psm)
    path bundle
    output: tuple val(sample), path('ms/pyneoquant_outputs/SNAF_proteome_peptide_evidence.tsv'), emit: evidence
    script:
    def config = groovy.json.JsonOutput.toJson([
        run_name: sample, output_dir: 'ms', sample_id: sample, assay_type: 'immunopeptidomics',
        source_manifest: "${bundle}/source_manifest.tsv", psm_table: psm.toString(),
        search_format: params.search_format, min_length: 7, max_length: 25, q_threshold: params.ms_q_threshold])
    """
    cat > workflow.json <<'PYNEOQUANT_CONFIG'
    ${config}
    PYNEOQUANT_CONFIG
    pyneoquant schema-version > pyneoquant_capabilities.json
    python3 - <<'PYNEOQUANT_CHECK'
    import csv, json
    from pathlib import Path
    capabilities = json.loads(Path('pyneoquant_capabilities.json').read_text())
    if 'fasta_record' not in capabilities.get('features', []):
        raise RuntimeError('pyNeoQuant requires the fasta_record contract (0.1.0a1 or newer)')
    config = json.loads(Path('workflow.json').read_text())
    candidate_path = Path(config['source_manifest']).parent / 'candidates.tsv'
    with candidate_path.open() as handle:
        samples = {r['sample_id'] for r in csv.DictReader(handle, delimiter='\\t')}
    if samples and config['sample_id'] not in samples:
        raise ValueError('MS sample is absent from the SNAF export: ' + config['sample_id'])
    PYNEOQUANT_CHECK
    pyneoquant run-workflow --config workflow.json
    """
    stub:
    """
    mkdir -p ms/pyneoquant_outputs
    printf 'sample_id\\tpeptide\\tsource_id\\tq_value\\n' > ms/pyneoquant_outputs/SNAF_proteome_peptide_evidence.tsv
    """
}

process MERGE_MS {
    tag "$sample"
    input:
    tuple val(sample), path(evidence)
    path bundle
    output: tuple val(sample), path("${sample}.candidates_with_ms.tsv"), emit: candidates
    script:
    """
    python3 -m altanalyze3.components.neoantigen merge --candidates ${sq(bundle + '/candidates.tsv')} --evidence ${sq(evidence)} --q-threshold ${params.ms_q_threshold} --output '${sample}.candidates_with_ms.tsv'
    """
    stub: "touch '${sample}.candidates_with_ms.tsv'"
}
