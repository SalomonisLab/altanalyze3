import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing
import re
import logging

logger = logging.getLogger(__name__)

def alignment_to_uniprot(orf,uid,dict_uni_fa,tmhmm=False,software_path=None):
    # first gene of the UID that has a reference protein -- for a trans-spliced junction
    # admitted on its partner, that is the partner.
    ensgid = uid.split(':')[0]
    if ensgid not in dict_uni_fa:
        import re as _re
        for _g in _re.findall(r'ENSG\d+', uid):
            if _g in dict_uni_fa:
                ensgid = _g
                break
    isoforms = dict_uni_fa[ensgid]  # {acc1:seq,acc1-2:seq}
    results = []
    for o in orf:
        if o != 'unrecoverable' and o != '':
            final_notes = []
            for iso,iso_seq in isoforms.items():
                if len(iso_seq) != len(o):
                    conclusion = True
                else:
                    notes = []
                    for i, mer in enumerate(chop_sequence(o,10)):
                        if mer in iso_seq:
                            notes.append(True)  # here True means align
                        else:
                            notes.append(False)
                    conclusion = not all(notes)   # here True means not existing before
                final_notes.append(conclusion)
            decision = all(final_notes)  # here True means this isoform is novel
            if tmhmm:
                is_cross_membrane = TMHMM(o,software_path=software_path)
                # None => the transmembrane gate is unavailable (e.g. tmhmm.py not
                # installed); degrade gracefully by NOT applying it (logged once).
                if is_cross_membrane is not None:
                    decision = all([decision,is_cross_membrane])
            results.append(decision)
        else:
            results.append(o)
    return results


def chop_sequence(seq,kmer):   # how to splice sequence, elegant way to use range
    frag_bucket = []
    for i in range(0,len(seq),kmer):
        if i + kmer <= len(seq):
            frag_bucket.append(seq[i:i+kmer])
        else:
            frag_bucket.append(seq[i:])
    return frag_bucket


_TMHMM_MODEL = None      # cached (module, parsed_model)
_TMHMM_STATUS = None     # None (untried) | 'ok' | 'unavailable'
_TMHMM_CACHE = {}        # aa -> bool (pure function; the Viterbi is the SNAF-B bottleneck)
_CACHE_MISS = object()   # sentinel so a cached False/None is distinguished from "not cached"


def _load_tmhmm_model():
    '''Load the pure-Python tmhmm.py model once (offline, no network, no binary).
    Returns (tmhmm_module, parsed_model) or None if tmhmm.py is not installed.'''
    global _TMHMM_MODEL, _TMHMM_STATUS
    if _TMHMM_STATUS is not None:
        return _TMHMM_MODEL
    try:
        # tmhmm.py's compiled Cython viterbi uses numpy.int / numpy.float, which numpy>=1.24
        # REMOVED -- so tmhmm.predict() raises AttributeError and the TM gate silently returns
        # False for everything (gate effectively disabled). Restore the deprecated aliases so
        # the HMM actually runs. Harmless no-op on older numpy.
        import numpy as _np
        for _a, _t in (('int', int), ('float', float), ('bool', bool), ('object', object)):
            if not hasattr(_np, _a):
                setattr(_np, _a, _t)
        import tmhmm as _tm
        from tmhmm.model import parse as _parse
        model_path = os.path.join(os.path.dirname(_tm.__file__), 'TMHMM2.0.model')
        with open(model_path) as fh:
            _header, model = _parse(fh)
        _TMHMM_MODEL = (_tm, model)
        _TMHMM_STATUS = 'ok'
        logger.info('SNAF-B transmembrane topology: pure-Python tmhmm.py (offline).')
    except Exception as e:
        _TMHMM_MODEL = None
        _TMHMM_STATUS = 'unavailable'
        logger.warning("SNAF-B: tmhmm.py not available (%s); skipping the transmembrane "
                       "topology gate (candidates will not be filtered by TM crossing). "
                       "Install it with:  pip install tmhmm.py", e)
    return _TMHMM_MODEL


def TMHMM(aa, software_path=None):
    '''True iff peptide `aa` crosses the membrane (>=1 predicted TM helix), computed
    with the pure-Python tmhmm.py (no network, no binary). Returns None when no TM
    predictor is available, so callers can degrade gracefully and skip the gate.

    If `software_path` points at a real (Linux) TMHMM 2.0c binary, that legacy path
    is used instead -- preserved for exact backward compatibility.'''
    # Legacy binary path only when an explicit, existing executable is provided.
    if software_path and os.path.exists(str(software_path)):
        name = os.getpid()
        scratch = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'scratch')
        os.makedirs(scratch, exist_ok=True)
        int_file_path = os.path.join(scratch, '{}.fasta'.format(name))
        with open(int_file_path, 'w') as f1:
            f1.write('>peptide_{}\n'.format(name)); f1.write(aa)
        lines = subprocess.run('{} {}'.format(software_path, int_file_path), shell=True,
                               stdout=subprocess.PIPE, universal_newlines=True).stdout.split('\n')[:-1]
        n_cross = int(lines[1].split('  ')[-1])
        return n_cross > 0

    # Pure-Python default.
    loaded = _load_tmhmm_model()
    if loaded is None:
        return None
    # Memoize: TMHMM is a pure function of `aa`, the Viterbi is the SNAF-B bottleneck (~85 aa-seq/s),
    # and the same candidate ORF/protein recurs across stringency x style x overlap passes and
    # across neojunctions on the same gene. Caching collapses those to one Viterbi each.
    cached = _TMHMM_CACHE.get(aa, _CACHE_MISS)
    if cached is not _CACHE_MISS:
        return cached
    _tm, model = loaded
    try:
        path = _tm.predict(aa, 'peptide', model, compute_posterior=False)
        result = len(re.findall(r'M+', path)) > 0
    except Exception as e:
        logger.debug('tmhmm.py predict failed (%s); treating peptide as non-TM', e)
        result = False
    _TMHMM_CACHE[aa] = result
    return result