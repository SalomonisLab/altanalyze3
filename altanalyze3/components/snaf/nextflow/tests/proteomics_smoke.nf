#!/usr/bin/env nextflow
// Exercise the independent library and return join with a pre-exported SNAF bundle.
nextflow.enable.dsl = 2
include { PYNEOQUANT; MERGE_MS } from '../modules/local/portable/main.nf'
params.bundle = null
params.psm = null
params.sample = 'S1'
workflow {
    if (!params.bundle || !params.psm) error 'Provide --bundle and --psm'
    if (!(params.sample ==~ /[A-Za-z0-9][A-Za-z0-9_-]*/)) error 'Invalid sample ID'
    bundle = Channel.value(file(params.bundle, checkIfExists: true))
    PYNEOQUANT(Channel.of(tuple(params.sample, file(params.psm, checkIfExists: true))), bundle)
    MERGE_MS(PYNEOQUANT.out.evidence, bundle)
}
