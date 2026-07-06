#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * SNAF nf-core-style workflow.
 *
 *   --mode ts       : tumor-specificity only (needs --control_h5ad)         [default]
 *   --mode full     : full MHC-bound T-antigen pipeline (needs --db_dir)
 *   --mode surface  : surface / B-antigen pipeline (needs --db_dir; consumes each
 *                     sample's T-antigen frequency table via the samplesheet `freq` column)
 *
 * Samplesheet CSV columns:
 *   id,juncounts[,hla][,freq]
 *     hla  required only for --mode full
 *     freq required only for --mode surface (the frequency_stage0_..._mean_mle.txt from a full run)
 */

include { SNAF_TS } from './modules/local/snaf_ts/main.nf'
include { SNAF    } from './modules/local/snaf/main.nf'
include { SNAF_B  } from './modules/local/snaf_b/main.nf'

params.input          = null      // samplesheet.csv
params.mode           = 'ts'
params.control_h5ad   = null      // for --mode ts
params.db_dir         = null      // for --mode full / surface
params.genome_fasta   = null      // optional, offline UTR retrieval
params.validation_gtf = null      // optional, --mode surface stringency-4/5 support gate
params.outdir         = 'results'

workflow {
    if (!params.input) { error "Provide --input samplesheet.csv" }

    ch_samples = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple(
            [id: row.id],
            file(row.juncounts),
            row.hla  ? file(row.hla)  : file('NO_FILE'),
            row.freq ? file(row.freq) : file('NO_FILE')) }

    if (params.mode == 'ts') {
        if (!params.control_h5ad) { error "--mode ts requires --control_h5ad" }
        ch_ts = ch_samples.map { meta, jc, hla, freq -> tuple(meta, jc) }
        SNAF_TS(ch_ts, file(params.control_h5ad))
        SNAF_TS.out.ts.view { meta, f -> "TS scores for ${meta.id}: ${f}" }
    }
    else if (params.mode == 'full') {
        if (!params.db_dir) { error "--mode full requires --db_dir" }
        def genome = params.genome_fasta ? file(params.genome_fasta) : file('NO_FILE')
        ch_full = ch_samples.map { meta, jc, hla, freq -> tuple(meta, jc, hla) }
        SNAF(ch_full, file(params.db_dir), genome)
        SNAF.out.versions.view()
    }
    else if (params.mode == 'surface' || params.mode == 'b') {
        if (!params.db_dir) { error "--mode surface requires --db_dir" }
        def gtf = params.validation_gtf ? file(params.validation_gtf) : file('NO_FILE')
        ch_b = ch_samples.map { meta, jc, hla, freq -> tuple(meta, jc, freq) }
        SNAF_B(ch_b, file(params.db_dir), gtf)
        SNAF_B.out.versions.view()
    }
    else {
        error "Unknown --mode '${params.mode}' (use 'ts', 'full', or 'surface')"
    }
}
