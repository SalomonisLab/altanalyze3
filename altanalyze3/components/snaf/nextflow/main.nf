#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { INDEX_REFERENCE; COUNT_JUNCTIONS; COUNT_INTRONS; AGGREGATE_COUNTS; FILTER_COUNTS; PREPARE_HLA; COMBINE_HLA; PYNEOQUANT; MERGE_MS } from './modules/local/portable/main.nf'
include { SNAF } from './modules/local/snaf/main.nf'
include { SNAF_TS } from './modules/local/snaf_ts/main.nf'
include { SNAF_B } from './modules/local/snaf_b/main.nf'

def checkedPath(value) { file(value, checkIfExists: true) }
def optionalPath(value) { value ? checkedPath(value) : [] }
def sheetPath(base, value) {
    if (!value) error 'Missing required file path in samplesheet'
    def path = file(value)
    checkedPath(path.isAbsolute() && value.startsWith('/') ? path : base.resolve(value))
}
def checkedId(value) {
    if (!value || !(value ==~ /[A-Za-z0-9][A-Za-z0-9_-]*/))
        error "Sample/cohort id '${value}' must contain letters, digits, underscores or hyphens"
    value
}

workflow {
    if (!params.input) error 'Provide --input samplesheet.csv'
    if (!(params.mode in ['bam', 'ts', 'full', 'surface', 'b'])) error 'mode must be bam, ts, full, or surface'
    if (!(params.isoform_method in ['learned', 'evidence', 'legacy'])) error 'Invalid isoform_method'
    if (!(params.hla_min_depth instanceof Number) || params.hla_min_depth < 1) error 'hla_min_depth must be positive'
    if (!(params.binding_method in ['MHCflurry', 'netMHCpan'])) error 'Invalid binding_method'
    if (!(params.strandness in ['auto', 'forward', 'reverse', 'unstranded'])) error 'Invalid strandness'
    if (!(params.novel_gene_mode in ['corrected', 'legacy'])) error 'Invalid novel_gene_mode'
    if (!(params.genome_build in ['auto', 'hg19', 'hg38'])) error 'Invalid genome_build'
    if (params.min_reads < 0 || params.ms_q_threshold < 0 || params.ms_q_threshold > 1) error 'Invalid count/q threshold'
    if (params.with_pyneoquant && !(params.mode in ['bam', 'full'])) error 'pyNeoQuant currently requires full or bam mode'
    if (params.with_pyneoquant && params.stop_after_counts) error 'pyNeoQuant requires prediction; remove stop_after_counts'
    if (params.binding_method == 'netMHCpan' && !params.netmhcpan_path) error 'Provide --netmhcpan_path'

    // A value channel containing the validated rows can be reused without racing consumers.
    sheetBase = checkedPath(params.input).parent
    seen = [] as Set
    rows = Channel.fromPath(params.input, checkIfExists: true).splitCsv(header: true)
        .map { row ->
            checkedId(row.id)
            if (!seen.add(row.id)) error "Duplicate samplesheet id: ${row.id}"
            row
        }.toList().map { data ->
            if (!data) error 'Samplesheet has no rows'
            if (params.with_pyneoquant && params.mode == 'bam' && !data.any{it.psm})
                error 'with_pyneoquant requires at least one psm path in the BAM samplesheet'
            if (params.with_pyneoquant && params.mode == 'full' && data.size() != 1)
                error 'MS integration requires one cohort matrix'
            data
        }
    genome = optionalPath(params.genome_fasta)
    canonical = optionalPath(params.canonical_fasta)
    mhcflurry = optionalPath(params.mhcflurry_models)

    if (params.mode == 'bam') {
        if (!params.gtf) error '--mode bam requires --gtf'
        ch_bams = rows.flatMap { it }.map { row ->
            def bam = sheetPath(sheetBase,row.bam)
            def bai = row.bai ? sheetPath(sheetBase,row.bai) : sheetPath(sheetBase,row.bam + '.bai')
            tuple(row.id, bam, bai)
        }
        INDEX_REFERENCE(checkedPath(params.gtf))
        COUNT_JUNCTIONS(ch_bams)
        COUNT_INTRONS(ch_bams, INDEX_REFERENCE.out.reference)
        AGGREGATE_COUNTS(COUNT_JUNCTIONS.out.counts.map{ id,p -> p }.collect(),
                         COUNT_INTRONS.out.counts.map{ id,p -> p }.collect(), INDEX_REFERENCE.out.reference)
        FILTER_COUNTS(AGGREGATE_COUNTS.out.counts.map{ tuple([id: 'cohort'], it) })
        if (!params.stop_after_counts) {
            if (!params.db_dir) error '--mode bam prediction requires --db_dir'
            if (params.infer_hla) {
                PREPARE_HLA(ch_bams, optionalPath(params.hla))
                COMBINE_HLA(PREPARE_HLA.out.calls.map{ id,p -> p }.collect())
                ch_hla = COMBINE_HLA.out.hla
            } else {
                if (!params.hla) error 'Supply --hla or enable --infer_hla'
                ch_hla = Channel.value(checkedPath(params.hla))
            }
            ch_full = FILTER_COUNTS.out.counts.combine(ch_hla).map{ meta,jc,hla -> tuple(meta,jc,hla) }
            SNAF(ch_full, checkedPath(params.db_dir), genome, canonical, mhcflurry)
            if (params.with_pyneoquant) {
                ch_psms = rows.flatMap{ it }.filter{ it.psm }.map{ tuple(it.id, sheetPath(sheetBase,it.psm)) }
                PYNEOQUANT(ch_psms, SNAF.out.bundle.first())
                MERGE_MS(PYNEOQUANT.out.evidence, SNAF.out.bundle.first())
            }
        }
    } else if (params.mode == 'full') {
        if (!params.db_dir) error '--mode full requires --db_dir'
        if (params.with_pyneoquant) {
            // One cohort matrix may contain many samples; MS uses its own sample sheet.
            if (!params.ms_input) error '--mode full with pyNeoQuant requires --ms_input (id,psm)'
        }
        ch_counts = rows.flatMap{ it }.map{ tuple([id: it.id], sheetPath(sheetBase,it.juncounts)) }
        FILTER_COUNTS(ch_counts)
        ch_hlas = rows.flatMap{ it }.map{ tuple([id: it.id], sheetPath(sheetBase,it.hla)) }
        ch_full = FILTER_COUNTS.out.counts.join(ch_hlas).map{ meta,jc,hla -> tuple(meta,jc,hla) }
        SNAF(ch_full, checkedPath(params.db_dir), genome, canonical, mhcflurry)
        if (params.with_pyneoquant) {
            ch_ms = Channel.fromPath(params.ms_input, checkIfExists:true).splitCsv(header:true)
                .map{ tuple(checkedId(it.id), sheetPath(checkedPath(params.ms_input).parent,it.psm)) }
            // Require one cohort to avoid ambiguous MS-to-cohort associations.
            bundle = SNAF.out.bundle.toList().map{ if(it.size()!=1) error 'MS integration requires one cohort'; it[0] }
            PYNEOQUANT(ch_ms, bundle)
            MERGE_MS(PYNEOQUANT.out.evidence, bundle)
        }
    } else if (params.mode == 'ts') {
        if (!params.control_h5ad) error '--mode ts requires --control_h5ad'
        SNAF_TS(rows.flatMap{it}.map{tuple([id:it.id],sheetPath(sheetBase,it.juncounts))}, checkedPath(params.control_h5ad))
    } else {
        if (!params.db_dir) error '--mode surface requires --db_dir'
        SNAF_B(rows.flatMap{it}.map{tuple([id:it.id],sheetPath(sheetBase,it.juncounts),sheetPath(sheetBase,it.freq))},
               checkedPath(params.db_dir), optionalPath(params.validation_gtf), genome,
               optionalPath(params.isoform_ranker), optionalPath(params.first_exon_junctions), params.isoform_method)
    }
}
