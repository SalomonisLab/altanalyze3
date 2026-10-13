import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing
import re
from decimal import Decimal as D
from collections import Counter
# for biopython, pip install biopython
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq



def score_GC(sequence):
    # Counter(sequence) tallied EVERY distinct character to read two of them: 3,098,823 calls
    # cost 121 s of a 701 s SNAF-B prewarm profile (104 s of it inside _count_elements), for
    # a number that two str.count calls give directly. str.count is exact, so the quotient is
    # bit-for-bit what Counter produced. An empty sequence still raises ZeroDivisionError.
    length_seq = len(sequence)
    GC_percent = (sequence.count('G') + sequence.count('C')) / length_seq
    return GC_percent

# Codon usage (GenScript) and its min-max normalisation are CONSTANTS. They were rebuilt on
# every call: a 64-entry dict literal, min(), max() and a dict comprehension, 3,098,823 times
# in one 300-task profile, which put score_coding_potential at 280 s of 701 s -- 40% of SNAF-B
# prewarm -- with 267 s of that in its own body. Building them once at import changes no value:
# the same arithmetic in the same order yields the same floats, so every lookup is identical.
_USAGE = {'TTT':16.9,'TTC':20.4,'TTA':7.2,'TTG':12.6,'TAT':12.0,'TAC':15.6,'TAA':0.7,'TAG':0.5,
                  'CTT':12.8,'CTC':19.4,'CTA':6.9,'CTG':40.3,'CAT':10.4,'CAC':14.9,'CAA':11.8,'CAG':34.6,
                  'ATT':15.7,'ATC':21.4,'ATA':7.1,'ATG':22.3,'AAT':16.7,'AAC':19.5,'AAA':24.0,'AAG':32.9,
                  'GTT':10.9,'GTC':14.6,'GTA':7.0,'GTG':28.9,'GAT':22.3,'GAC':26.0,'GAA':29.0,'GAG':40.8,
                  'TCT':14.6,'TCC':17.4,'TCA':11.7,'TCG':4.5,'TGT':9.9,'TGC':12.2,'TGA':1.3,'TGG':12.8,
                  'CCT':17.3,'CCC':20.0,'CCA':16.7,'CCG':7.0,'CGT':4.7,'CGC':10.9,'CGA':6.3,'CGG':11.9,
                  'ACT':12.8,'ACC':19.2,'ACA':14.8,'ACG':6.2,'AGT':11.9,'AGC':19.4,'AGA':11.5,'AGG':11.4,
                  'GCT':18.6,'GCC':28.5,'GCA':16.0,'GCG':7.6,'GGT':10.8,'GGC':22.8,'GGA':16.3,'GGG':16.4}
_MIN_FREQ = min(list(_USAGE.values()))
_MAX_FREQ = max(list(_USAGE.values()))
_NORM_USAGE = {k: (v - _MIN_FREQ) / (_MAX_FREQ - _MIN_FREQ) for k, v in _USAGE.items()}


def score_coding_potential(sequence):
    length_seq = len(sequence)
    num_triplet = length_seq / 3
    i = 0
    score = 0
    skipped = 0
    while i + 2 < length_seq:
        triplet = sequence[i:i+3:1]
        # A codon holding an ambiguous base is not in the 64-codon usage table, and indexing it
        # raised `KeyError: 'GNC'` after 59.6 h of SNAF-B on POSEIDON Small_Intestine
        # (job 1375812, 2026-10-08). hg38 and the Alt91 gene FASTA both carry N runs, so an
        # ambiguous codon is ordinary input, not corruption. It carries no usage information, so
        # it is dropped from BOTH the sum and the denominator. With no ambiguous codon present,
        # skipped is 0 and this returns score/num_triplet, byte-identical to the original.
        if triplet in _NORM_USAGE:
            score += _NORM_USAGE[triplet]
        else:
            skipped += 1
        i += 3
    denom = num_triplet - skipped
    if denom <= 0:
        # every codon was ambiguous: no coding-potential evidence either way
        return 0.0
    return score / denom

def orf2pep(orf):
    assert len(orf) % 3 == 0
    pep = str(Seq(orf).translate(to_stop=False))
    assert '*' not in pep
    return pep

_ENUM_ORFS = None   # predict_isoform.enumerate_orfs_with_offsets, bound on first use


def transcript2orf(cdna):
    """Every ATG-to-first-in-frame-stop reading frame in a cDNA.

    The original paired EVERY start with EVERY stop. Two scans sat inside the start loop:
    `re.finditer(p_end, cdna)` re-scanned the whole cDNA once per ATG, and each in-frame
    (start, stop) pair then sliced the ORF and scanned THAT too. For a 3,000-nt transcript
    with ~47 ATGs and ~140 stops that is ~47 full-length scans plus ~2,200 slice scans, when
    the answer needs exactly two scans of the cDNA. Only the FIRST in-frame stop after an ATG
    ever yields an ORF -- any later one encloses that stop and is rejected -- so the pairing
    was computing and discarding the rest.

    is_support_by_est_or_long_read calls this once per candidate transcript of the merged
    10.6 M-transcript catalogue, inside prewarm_support_cache, which is the SNAF-B phase that
    ran for over a day per POSEIDON cohort.

    The scan-once form already exists and is already validated: predict_isoform
    .enumerate_orfs_with_offsets, which agreed with this function on 401 real transcripts and
    is asserted equal to it by test_predict_isoform.py. Delegating keeps ONE implementation.
    The import is lazy because predict_isoform imports orf2pep from this module.
    """
    global _ENUM_ORFS
    if _ENUM_ORFS is None:
        from .predict_isoform import enumerate_orfs_with_offsets as _e
        _ENUM_ORFS = _e
    return [orf for _off, orf in _ENUM_ORFS(cdna)]


def prioritize_orf(candidate_orfs,min_len=30*3,tol_len=8*3):
    max_orf = ''
    max_length = 0
    max_score = 0
    for orf in candidate_orfs:
        if len(orf) < min_len:
            continue
        score = score_GC(orf) + score_coding_potential(orf)
        if len(orf) - max_length > tol_len:
            max_length = len(orf)
            max_score = score
            max_orf = orf
        else:
            if len(orf) > max_length and score > max_score:
                max_length = len(orf)
                max_score = score
                max_orf = orf  
    return max_orf            


