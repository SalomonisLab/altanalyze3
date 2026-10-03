"""Collapse SNAF-T per-sample candidates into one row per (peptide, junction).

SNAF-T writes one row per (sample, peptide, junction, HLA). A cohort of 102 samples therefore
produces millions of rows, and the question a reader asks -- "which neoantigens are worth
pursuing, and in how many patients" -- is not answerable from that shape. This collapses the
table to one row per PEPTIDE and JUNCTION, carrying the sample set, the binding and
immunogenicity summaries, and the splicing event that created it.

PROVENANCE. This is a port of `updatedSNAFTMerge` in AltAnalyze's
`stats_scripts/SNAFintegration.py` (Python 2), which is the validated reference for SNAF-T
candidate summarization. Every column that reference emits is emitted here under the same
name and computed the same way, including the two it computes per (HLA, sample) pair rather
than per row. `tests/test_snaf_t_summary.py` asserts that equivalence on a worked example.

WHAT THIS ADDS, and why each one is needed. Every added column is an aggregation of a value
SNAF or AltAnalyze already computed. None introduces a threshold, a null model or a score.

| Added column | Why the reference's columns are not enough |
|---|---|
| `Best Binding`, `Best Binding HLA` | `Mean Binding` averages every allele, so one strong binder among four weak ones reads as mediocre. The strongest allele is what a therapy targets. |
| `Max Immunogenicity` | same argument on the immunogenicity side |
| `Num HLA Alleles`, `HLA Alleles` | a peptide presented by several alleles covers more patients; the mean hides promiscuity entirely |
| `Num Observations` | the denominator behind the two means, so a mean over 1 pair is not read as a mean over 40 |
| `Peptide Length` | 8-11mers are class-I; length is a basic filter and is not otherwise recoverable |
| `Max Junction Count` | the reference's `JunctionCount` is whichever sample was read FIRST, which is arbitrary across samples |
| `Sample Frequency` | `Num Samples` without the cohort denominator cannot be compared between studies |
| `Phase` | the reference infers `In-Frame` from the LENGTH of the evidence string; SNAF reports the reading frame outright |
| `ProteinPredictions`, `dPSI`, `AltExons`, `Description` | the reference parses these out of EventAnnotation and then drops them. They state the protein consequence, the splicing magnitude and the alternative exon: the biology that makes the peptide a neoantigen. |

DIRECTION OF EACH METRIC. `binding_affinity` is a percentile where LOWER is a stronger
binder, so `Best Binding` is the minimum. `immunogenicity` is a 0-1 score where higher is
more immunogenic, so `Max Immunogenicity` is the maximum.
"""
from __future__ import annotations

import glob
import logging
import os

logger = logging.getLogger(__name__)

# The per-sample candidate columns SNAF-T writes (`downstream.py`). `in_db` appears only in
# the concatenated `T_antigen_candidates_all.txt`, where SNAF adds it.
CANDIDATE_COLUMNS = ['sample', 'peptide', 'uid', 'junction_count', 'phase', 'evidences',
                     'hla', 'binding_affinity', 'immunogenicity', 'tumor_specificity_mean',
                     'tumor_specificity_mle', 'n_sample', 'coord', 'symbol']

# EventAnnotation column POSITIONS, which is how the reference reads them. Verified against
# AltAnalyze-91 output on 2026-10-02: Symbol, Description, Examined-Junction,
# Background-Major-Junction, AltExons, ProteinPredictions, dPSI, ClusterID, UID,
# Coordinates, EventAnnotation.
EA_POS = {'symbol': 0, 'description': 1, 'junction': 2, 'alt_exons': 4,
          'protein_predictions': 5, 'dpsi': 6, 'cluster_id': 7, 'uid': 8,
          'coordinates': 9, 'event_annotation': 10}

# Column order of the written table. The first block is the reference's header, unchanged and
# in its order, so a reader of the old output finds every column where it was.
REFERENCE_COLUMNS = ['Peptide', 'Junction', 'Symbol', 'Coord', 'JunctionCount', 'In-Frame',
                     'Tumor Specificity Mean', 'Tumor Specificity MLE', 'Num Samples',
                     'Mean Binding', 'Mean Immunogenicity', 'Samples', 'UID', 'ClusterID',
                     'Coordinates', 'EventAnnotation', 'Present-In-Ensembl',
                     'Wald', 'Zscore', 'n_sample', 'peptideID']
ADDED_COLUMNS = ['Peptide Length', 'Phase', 'Junction n_sample',
                 'Best Binding', 'Best Binding HLA',
                 'Max Immunogenicity', 'Num HLA Alleles', 'HLA Alleles', 'Num Observations',
                 'Max Junction Count', 'Sample Frequency',
                 'ProteinPredictions', 'dPSI', 'AltExons', 'Description']


def _strip_bed(name):
    name = str(name)
    return name[:-4] if name.endswith('.bed') else name


def import_event_annotations(path):
    """{junction: {...}} from an AltAnalyze PSI EventAnnotation table.

    Read by POSITION, as the reference does, because the trailing columns are one per sample
    and the header width therefore changes with the cohort.
    """
    out = {}
    with open(path, encoding='utf-8', errors='replace') as fh:
        first = True
        for line in fh:
            if first:
                first = False
                continue
            t = line.rstrip('\n').split('\t')
            if len(t) <= EA_POS['event_annotation']:
                continue
            out[t[EA_POS['junction']]] = {k: t[i] for k, i in EA_POS.items()}
    logger.info('EventAnnotation: %d junctions from %s', len(out), path)
    return out


def generic_import(path, key=None):
    """The reference's `genericImport`: first row is the header, column 0 (or 0+1) is the key."""
    db = {}
    with open(path, encoding='utf-8', errors='replace') as fh:
        first = True
        for line in fh:
            t = line.rstrip('\n').split('\t')
            if key == 'multiple':
                if len(t) < 2:
                    continue
                id2 = t[1].replace('-', '.').replace(':', '.') if '-' in t[1] else t[1]
                uid, values = (t[0], id2), t[2:]
            else:
                uid, values = t[0], t[1:]
            if first:
                first = False
                uid = 'header'
            db[uid] = values
    return db


def _load_survival(path):
    """{(peptide, junction): (n_sample, Zscore, Wald)}, the reference's 6-column layout."""
    out = {}
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            t = line.rstrip('\n').split('\t')
            if len(t) < 6:
                continue
            uid, pep, junc, n_sample, zscore, wald = t[:6]
            out[(pep, junc)] = (n_sample, zscore, wald)
    return out


def _read_candidates(candidate_dir, sample_filter=None):
    """The per-sample candidate rows as one DataFrame.

    `T_antigen_candidates_all.txt` is preferred when SNAF wrote it, because only that file
    carries `in_db` (the reference's Present-In-Ensembl). Without it the per-sample files are
    concatenated and Present-In-Ensembl stays 'UNK', which is what the reference records when
    its 15th column is absent.
    """
    import pandas as pd

    all_path = os.path.join(candidate_dir, 'T_antigen_candidates_all.txt')
    if os.path.isfile(all_path) and os.path.getsize(all_path) > 0:
        df = pd.read_csv(all_path, sep='\t', dtype=str)
        source = all_path
    else:
        files = sorted(p for p in glob.glob(os.path.join(candidate_dir, 'T_antigen_candidates_*.txt'))
                       if not p.endswith('T_antigen_candidates_all.txt'))
        if not files:
            raise FileNotFoundError(
                'no T_antigen_candidates_*.txt under {}. SNAF-T writes them into '
                '<output>/T_candidates/.'.format(candidate_dir))
        df = pd.concat([pd.read_csv(p, sep='\t', dtype=str) for p in files], axis=0,
                       ignore_index=True)
        source = '{} per-sample files'.format(len(files))
    missing = [c for c in CANDIDATE_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError('{} lacks the SNAF-T candidate columns {}'.format(source, missing))
    df['sample'] = df['sample'].map(_strip_bed)
    if sample_filter:
        keep = {_strip_bed(s) for s in sample_filter}
        before = df['sample'].nunique()
        df = df[df['sample'].isin(keep)]
        logger.info('sample filter: %d of %d samples kept', df['sample'].nunique(), before)
    logger.info('SNAF-T candidates: %d rows over %d samples (%s)',
                len(df), df['sample'].nunique(), source)
    return df


def summarize_t_candidates(candidate_dir, out_path=None, event_annotation=None,
                           survival=None, custom_peptide=None, custom_junction=None,
                           sample_filter=None, cohort_size=None):
    """Write one row per (peptide, junction). Returns the output path.

    :param candidate_dir: SNAF-T's ``<output>/T_candidates`` directory
    :param out_path: default ``<candidate_dir>/summary/collapsed_SNAF-T.txt``, the reference's name
    :param event_annotation: AltAnalyze ``*-PSI_EventAnnotation.txt``; adds the splicing event
    :param survival: optional per-peptide survival table (uid, peptide, junction, n, Z, Wald)
    :param custom_peptide: optional table keyed on (peptide, junction), columns appended
    :param custom_junction: optional table keyed on junction, columns appended
    :param sample_filter: optional iterable of sample names to keep
    :param cohort_size: denominator for ``Sample Frequency``; defaults to the samples observed
    """
    import pandas as pd

    df = _read_candidates(candidate_dir, sample_filter=sample_filter)
    n_cohort = int(cohort_size) if cohort_size else int(df['sample'].nunique())

    # The reference stores binding and immunogenicity in dicts keyed (hla, sample), so a
    # repeated (peptide, junction, hla, sample) row contributes ONCE. Dropping the duplicates
    # here reproduces that; grouping the raw rows would weight such a pair twice.
    before = len(df)
    df = df.drop_duplicates(subset=['peptide', 'uid', 'hla', 'sample'])
    if len(df) != before:
        logger.info('collapsed %d duplicate (peptide, junction, HLA, sample) rows',
                    before - len(df))

    df['_binding'] = pd.to_numeric(df['binding_affinity'], errors='coerce')
    df['_immuno'] = pd.to_numeric(df['immunogenicity'], errors='coerce')
    df['_jc'] = pd.to_numeric(df['junction_count'], errors='coerce')
    # the reference's In-Frame: the evidence string being longer than 3 characters
    df['_inframe'] = df['evidences'].fillna('').map(lambda v: 'True' if len(str(v)) > 3 else 'UNK')

    g = df.groupby(['peptide', 'uid'], sort=False)
    out = pd.DataFrame({
        'Symbol': g['symbol'].first(),
        'Coord': g['coord'].first(),
        'JunctionCount': g['junction_count'].first(),
        'In-Frame': g['_inframe'].first(),
        'Tumor Specificity Mean': g['tumor_specificity_mean'].first(),
        'Tumor Specificity MLE': g['tumor_specificity_mle'].first(),
        'Num Samples': g['sample'].nunique(),
        'Mean Binding': g['_binding'].mean(),
        'Mean Immunogenicity': g['_immuno'].mean(),
        'Samples': g['sample'].apply(lambda s: '|'.join(dict.fromkeys(s))),
        # SNAF's own per-junction sample count. The reference stores it and never writes it,
        # and its 'n_sample' COLUMN carries the survival table's value instead, so this keeps
        # that column meaning what it meant and gives SNAF's number its own name.
        'Junction n_sample': g['n_sample'].first(),
        'Phase': g['phase'].first(),
        'Best Binding': g['_binding'].min(),
        'Max Immunogenicity': g['_immuno'].max(),
        'Num HLA Alleles': g['hla'].nunique(),
        'HLA Alleles': g['hla'].apply(lambda s: '|'.join(sorted(set(s)))),
        'Num Observations': g.size(),
        'Max Junction Count': g['_jc'].max(),
    })
    out['Present-In-Ensembl'] = (g['in_db'].first() if 'in_db' in df.columns
                                 else pd.Series('UNK', index=out.index))
    # the allele that achieves Best Binding, taken from the row holding the minimum
    best = df.loc[df.groupby(['peptide', 'uid'], sort=False)['_binding'].idxmin().dropna()]
    out['Best Binding HLA'] = best.set_index(['peptide', 'uid'])['hla']
    out['Sample Frequency'] = (out['Num Samples'] / n_cohort).round(4) if n_cohort else ''
    out['Peptide Length'] = [len(str(p)) for p, _ in out.index]

    out = out.reset_index().rename(columns={'peptide': 'Peptide', 'uid': 'Junction'})

    # ---- the splicing event that created the peptide ---------------------------------
    ea = import_event_annotations(event_annotation) if event_annotation else {}
    blank = {k: '' for k in EA_POS}

    def _ea(j, field):
        return ea.get(j, blank).get(field, '')

    for col, field in (('UID', 'uid'), ('ClusterID', 'cluster_id'), ('Coordinates', 'coordinates'),
                       ('EventAnnotation', 'event_annotation'),
                       ('ProteinPredictions', 'protein_predictions'), ('dPSI', 'dpsi'),
                       ('AltExons', 'alt_exons'), ('Description', 'description')):
        out[col] = [_ea(j, field) for j in out['Junction']]
    if ea:
        hit = sum(1 for j in out['Junction'] if j in ea)
        logger.info('EventAnnotation matched %d of %d candidate junctions (%.1f%%)',
                    hit, len(out), 100.0 * hit / max(1, len(out)))
        # the reference falls back to the EventAnnotation symbol only when the candidate has none
        out['Symbol'] = [s if str(s).strip() else _ea(j, 'symbol')
                         for s, j in zip(out['Symbol'], out['Junction'])]

    out['peptideID'] = [('%s.%s' % (p, j.replace(':', '.'))).replace('-', '.')
                        for p, j in zip(out['Peptide'], out['Junction'])]

    # ---- optional joins ---------------------------------------------------------------
    surv = _load_survival(survival) if survival else {}
    def _surv(p, j, i):
        alt = j.replace(':', '.').replace('-', '.')
        rec = surv.get((p, j)) or surv.get((p, alt))
        return rec[i] if rec else ''
    out['n_sample'] = [_surv(p, j, 0) for p, j in zip(out['Peptide'], out['Junction'])]
    out['Zscore'] = [_surv(p, j, 1) for p, j in zip(out['Peptide'], out['Junction'])]
    out['Wald'] = [_surv(p, j, 2) for p, j in zip(out['Peptide'], out['Junction'])]

    columns = list(REFERENCE_COLUMNS) + list(ADDED_COLUMNS)
    for extra, loader, keyfn in (
            (custom_peptide, lambda p: generic_import(p, key='multiple'),
             lambda p, j: (p, j.replace('-', '.').replace(':', '.'))),
            (custom_junction, lambda p: generic_import(p), lambda p, j: j)):
        if not extra:
            continue
        db = loader(extra)
        head = db.get('header', [])
        for n, name in enumerate(head):
            out[name] = [(db.get(keyfn(p, j), [''] * len(head)) + [''] * len(head))[n]
                         for p, j in zip(out['Peptide'], out['Junction'])]
        columns += [c for c in head if c not in columns]

    for c in columns:
        if c not in out.columns:
            out[c] = ''
    out = out[columns]
    # `Num Samples` is the headline: the most-shared candidates first, then the strongest binder
    out = out.sort_values(['Num Samples', 'Best Binding'], ascending=[False, True])

    if out_path is None:
        out_path = os.path.join(candidate_dir, 'summary', 'collapsed_SNAF-T.txt')
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    out.to_csv(out_path, sep='\t', index=False)
    logger.info('collapsed SNAF-T: %d peptide-junction pairs over %d samples -> %s',
                len(out), n_cohort, out_path)
    return out_path
