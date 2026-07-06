#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import matplotlib.pyplot as plt
import subprocess
from io import StringIO
from copy import deepcopy,copy




'''
this script is to query the binding affinity of a peptide (9-10)
'''

'''
part I: using netMHCpan4.1b

1. download from https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1
2. tar -xzvf netMHCpan-4.1b.Linux.tar.gz
3. download the data.tar.gz to the netMHCpan-4.1 folder, untar it using the above command, delete the data.tar.gz
4. modify the netMHCpan file, instruction as mentioned in the readme:

      a. At the top of the file  locate the part labelled  "GENERAL SETTINGS:
         CUSTOMIZE TO YOUR SITE"  and set  the 'NMHOME' variable  to the full
	     path to the 'netMHCpan-4.1' directory on your system;

      b. Set TMPDIR to the full path to the temporary directory of you choice. It must
         be user-writable. You may for example set it to $NMHOME/tmp (and create
         the tmp folder in the netMHCpan-4.1 directory).

5. done
'''


def _parse_netMHCpan_stdout(stdout, peptide_set):
    '''Pure-Python replacement for the old `... | awk '{print $3,length($3),$2,$13,$15}'`
    pipeline. Parses netMHCpan -p tabular output without any shell/awk/grep, so it
    is portable (the awk `||`/`\\|` pipelines were Linux-shell-specific). Column
    mapping is identical to the legacy awk: peptide=$3 (index 2), hla=$2 (index 1),
    score=$13 (index 12, %Rank_EL), identity=$15 (index 14, the SB/WB BindLevel, which
    is only present for binders -> empty otherwise, exactly as awk emitted).'''
    rows = []
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) < 13:          # a data line runs at least through %Rank_EL ($13)
            continue
        peptide = parts[2]
        if peptide not in peptide_set:   # skips headers / separators / blanks
            continue
        identity = parts[14] if len(parts) >= 15 else ''
        rows.append((peptide, len(peptide), parts[1], parts[12], identity))
    df = pd.DataFrame(rows, columns=['peptide', 'mer', 'hla', 'score', 'identity'])
    if not df.empty:
        df['score'] = pd.to_numeric(df['score'], errors='coerce')
    return df


def _run_netMHCpan_once(software_path, peptides_path, hla_strings, length, peptide_set):
    '''Invoke the netMHCpan binary directly (no shell) and parse its stdout.'''
    cmd = [software_path, '-p', peptides_path, '-a', hla_strings, '-l', str(length)]
    try:
        proc = subprocess.run(cmd, capture_output=True)
        return _parse_netMHCpan_stdout(proc.stdout.decode('utf-8', errors='replace'), peptide_set)
    except Exception:
        return pd.DataFrame(columns=['peptide', 'mer', 'hla', 'score', 'identity'])


def run_netMHCpan(software_path,peptides,hlas,length,cmd_num=1,tmp_dir=None,tmp_name=None):
    # set the default
    if tmp_dir is None:
        tmp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'scratch')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)
    if tmp_name is None:
        tmp_name = 'input_{}.pep'.format(os.getpid())
    # reformat/create to the strings that we need
    peptides_path = os.path.join(tmp_dir,tmp_name)
    with open(peptides_path,'w') as f:
        for pep in peptides:
            f.write('{}\n'.format(pep))
    peptide_set = set(peptides)
    # netMHCpan has a hard limit of 1024 characters for -a, so at most 90 HLA per call
    # ('HLA-B57:01,' is ~11 chars). Batch the HLA list and concatenate the parses.
    df_store = []
    total = len(hlas)
    i = 0
    while i < total:
        batch_hlas = hlas[i:i + 90]
        hla_strings = ','.join(batch_hlas)
        df_store.append(_run_netMHCpan_once(software_path, peptides_path, hla_strings, length, peptide_set))
        i += 90
    if df_store:
        df = pd.concat(df_store, axis=0)
    else:
        df = pd.DataFrame(columns=['peptide', 'mer', 'hla', 'score', 'identity'])

    # remove the temp peptide file
    if os.path.exists(peptides_path):
        os.remove(peptides_path)

    return df



# cache the (slow-to-load) predictor across calls within a process
_MHCFLURRY_PREDICTOR = None


def _load_mhcflurry_predictor():
    '''Load the MHCflurry Class1 presentation predictor, fetching the model bundle
    once if it is not already present. MHCflurry is pure-Python and cross-platform,
    so this is the portable binding path (macOS/Windows/Linux). Set the environment
    variable MHCFLURRY_DOWNLOADS_DIR to a pre-downloaded model directory for fully
    offline use. Fixes the legacy bug where `predictor` was never re-assigned after
    the download fallback (which raised NameError).'''
    global _MHCFLURRY_PREDICTOR
    if _MHCFLURRY_PREDICTOR is not None:
        return _MHCFLURRY_PREDICTOR
    from mhcflurry import Class1PresentationPredictor
    try:
        predictor = Class1PresentationPredictor.load()
    except Exception:
        # models not present locally -> fetch once, then load (offline if env var set)
        try:
            from mhcflurry.downloads_command import run as _mhcflurry_downloads
            _mhcflurry_downloads(['fetch', 'models_class1_presentation'])
        except Exception:
            subprocess.run([sys.executable, '-m', 'mhcflurry.downloads_command',
                            'fetch', 'models_class1_presentation'])
        predictor = Class1PresentationPredictor.load()
    _MHCFLURRY_PREDICTOR = predictor
    return predictor


# Process-level cache of MHCflurry presentation percentiles. The score is a pure
# function of (peptide, allele), but each neojunction re-predicts the Cartesian product
# of its peptides x the cohort's (heavily shared) HLA alleles -- so on a large cohort the
# same (peptide, allele) pair is recomputed thousands of times. Memoizing makes each
# unique pair cost one prediction; results are byte-identical (predict() is deterministic).
_MHCFLURRY_SCORE_CACHE = {}   # (peptide, allele) -> presentation_percentile (float)


# Antigen-processing (cleavage) score is ALLELE-INDEPENDENT -- it depends only on the
# peptide (and its flanks). Stock mhcflurry recomputes it on every (peptide, allele) call,
# and it is ~96% of presentation runtime. Computing it ONCE per unique peptide and reusing
# mhcflurry's own affinity + logistic + percentile per allele is bit-identical to
# Class1PresentationPredictor.predict (verified: max abs diff 0.0) and ~10x faster at cohort
# scale. This cache holds peptide -> processing score; it is populated in the parent before
# the fork so workers inherit it copy-on-write.
_PROCESSING_CACHE = {}


def compute_processing_scores(peptides):
    """Worker: antigen-processing score for a chunk of peptides (no flanks, as SNAF calls it).
    Returns {peptide: processing_score}. Uses the pure-NumPy reimplementation
    (fast_infer_processing) when available -- bit-identical, no TensorFlow -- else stock mhcflurry."""
    peptides = list(peptides)
    if not peptides:
        return {}
    predictor = _load_mhcflurry_predictor()
    # DEFAULT = pure-NumPy reimplementation (fast_infer_processing): measured 7.2x faster than
    # stock mhcflurry predict_processing (2.3s vs 16.4s on 9,415 peptides) and bit-identical
    # (max abs diff 2.9e-7). Antigen-processing is ~95% of binding cost, and TensorFlow buys
    # nothing here -- it even blocks fork on macOS. Set SNAF_FAST_NN=0 to force stock TF.
    scores = None
    if os.environ.get('SNAF_FAST_NN') != '0':
        try:
            from .fast_infer_processing import predict_processing_fast
            scores = predict_processing_fast(predictor, peptides)
        except Exception:
            scores = None
    if scores is None:
        scores = predictor.predict_processing(peptides=peptides, verbose=0)
    return {p: float(s) for p, s in zip(peptides, scores)}


def predict_presentation_pairs(pairs):
    """Worker: presentation percentile for a chunk of (peptide, allele) pairs, reusing the
    inherited processing cache. Bit-identical to Class1PresentationPredictor.predict."""
    import pandas as pd
    import numpy as np
    from collections import defaultdict
    from mhcflurry.regression_target import from_ic50
    predictor = _load_mhcflurry_predictor()
    model = predictor.get_model('without_flanks')
    predict_affinity_fast = None
    if os.environ.get('SNAF_FAST_NN') == '1':
        try:
            from .fast_infer_affinity import predict_affinity_fast
        except Exception:
            predict_affinity_fast = None
    by_allele = defaultdict(list)
    for p, h in pairs:
        by_allele[h].append(p)
    out = {}
    for h, peps in by_allele.items():
        peps_u = list(dict.fromkeys(peps))
        if predict_affinity_fast is not None:
            # pure-NumPy affinity (bit-identical, no TF); build the same frame mhcflurry would
            aff = pd.DataFrame({'peptide': peps_u,
                                'affinity': np.asarray(predict_affinity_fast(predictor, peps_u, h), dtype=float)})
        else:
            aff = predictor.predict_affinity(peptides=peps_u, alleles={'s': [h]},
                                             include_affinity_percentile=False, verbose=0)
        aff['affinity_score'] = from_ic50(aff.affinity)
        aff['processing_score'] = [_PROCESSING_CACHE.get(p) for p in aff['peptide']]
        im = aff[predictor.model_inputs]
        null_mask = im.isnull().any(axis=1).values          # mirror mhcflurry's throw=False path
        ps = model.predict_proba(im.fillna(0.0).values)[:, 1]
        ps = np.asarray(ps, dtype=float)
        ps[null_mask] = np.nan
        pct = predictor.percentile_ranks(pd.Series(ps), throw=False)
        for p, v in zip(aff['peptide'].tolist(), np.asarray(pct).ravel().tolist()):
            out[(p, h)] = v
    return out


def run_MHCflurry(peptides, hlas):
    import pandas as pd
    from collections import defaultdict
    peptides = list(peptides)
    hlas = list(hlas)
    # Which (peptide, allele) pairs are not cached yet?
    missing = defaultdict(list)   # allele -> [peptides needing prediction]
    for h in hlas:
        for p in peptides:
            if (p, h) not in _MHCFLURRY_SCORE_CACHE:
                missing[h].append(p)
    if missing:
        predictor = _load_mhcflurry_predictor()
        for h, peps in missing.items():
            peps_u = list(dict.fromkeys(peps))                       # unique, preserve order
            res = predictor.predict(peptides=peps_u, alleles={'s': [h]}, verbose=0)  # single allele -> best_allele == h
            for pep, pct in zip(res['peptide'].tolist(), res['presentation_percentile'].tolist()):
                _MHCFLURRY_SCORE_CACHE[(pep, h)] = pct
    # Rebuild the exact same df the original returned (one row per peptide x hla; row
    # order is irrelevant -- register_attr keys by peptide+hla).
    rows = [(p, len(p), h, _MHCFLURRY_SCORE_CACHE[(p, h)], None) for h in hlas for p in peptides]
    return pd.DataFrame(rows, columns=['peptide', 'mer', 'hla', 'score', 'identity'])






