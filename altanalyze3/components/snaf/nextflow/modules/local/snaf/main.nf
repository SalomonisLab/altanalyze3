def sq(x) { "'" + x.toString().replace("'", "'\\''") + "'" }
process SNAF {
    tag "$meta.id"
    label 'process_high'
    input:
    tuple val(meta), path(juncounts), path(hla)
    path db_dir
    path genome_fasta
    path canonical_fasta
    path mhcflurry_models
    path galaxy_workflow
    path galaxy_peptide_bed
    output:
    tuple val(meta), path("${meta.id}/T_candidates"), emit: candidates
    tuple val(meta), path("${meta.id}/frequency_stage*"), emit: frequency
    path "${meta.id}/proteomics_export", emit: bundle
    path "${meta.id}/after_prediction.p", emit: pickle
    tuple val(meta), path("${meta.id}/galaxy_export"), emit: galaxy, optional: true
    script:
    def genome = genome_fasta ? "--genome_fasta ${sq(genome_fasta)}" : ''
    def canonical = canonical_fasta ? "--canonical_fasta ${sq(canonical_fasta)}" : ''
    def external = params.netmhcpan_path ? "--software_path ${sq(params.netmhcpan_path)}" : ''
    def modelEnv = mhcflurry_models ? "export MHCFLURRY_DOWNLOADS_DIR=${sq(mhcflurry_models)}" : ''
    def galaxy = params.galaxy_integration ? "--galaxy_integration --galaxy_assembly ${sq(params.galaxy_assembly)}" : ''
    def workflowFile = params.galaxy_integration && galaxy_workflow ? "--galaxy_workflow ${sq(galaxy_workflow)}" : ''
    def bedFile = params.galaxy_integration && galaxy_peptide_bed ? "--galaxy_peptide_bed ${sq(galaxy_peptide_bed)}" : ''
    def args = task.ext.args ?: ''
    """
    export SNAF_OFFLINE=1
    ${modelEnv}
    altanalyze3 snaf --juncounts ${sq(juncounts)} --db_dir ${sq(db_dir)} --hla ${sq(hla)} --output '${meta.id}' --cpus ${task.cpus} --binding_method '${params.binding_method}' --export_proteomics ${genome} ${canonical} ${external} ${galaxy} ${workflowFile} ${bedFile} ${args}
    """
    stub:
    def galaxyStub = params.galaxy_integration ? "mkdir -p '${meta.id}/galaxy_export'; printf '{}' > '${meta.id}/galaxy_export/integration.json'" : ':'
    """
    ${galaxyStub}
    mkdir -p '${meta.id}/T_candidates' '${meta.id}/proteomics_export'
    touch '${meta.id}/frequency_stage0.txt' '${meta.id}/after_prediction.p'
    printf 'candidate_id\\tsource_id\\tsample_id\\tpeptide\\tevent_id\\thla_allele\\n' > '${meta.id}/proteomics_export/candidates.tsv'
    printf 'source_id\\tsource_type\\tfasta_path\\tfasta_record\\n' > '${meta.id}/proteomics_export/source_manifest.tsv'
    touch '${meta.id}/proteomics_export/candidates.fasta'
    """
}
