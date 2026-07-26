#!/usr/bin/env python3
"""altanalyze3 OncoKB annotator -- append OncoKB clinical evidence to a variant TSV via the OncoKB
byGenomicChange batch API (urllib only, no extra packages). Faithful port of Xuan's
2026-04-22_variant_calling/scripts/oncokb_annotate.py, with case-insensitive CHROM/POS/REF/ALT so it
reads the functional-nomination candidate table directly.

Appends: ONCOKB_ONCOGENIC, ONCOKB_MUTATION_EFFECT, ONCOKB_HIGHEST_LEVEL, ONCOKB_HIGHEST_RESISTANCE_LEVEL,
ONCOKB_HIGHEST_DX_LEVEL, ONCOKB_HIGHEST_PX_LEVEL, ONCOKB_HIGHEST_FDA_LEVEL, ONCOKB_TX_LEVELS,
ONCOKB_DX_LEVELS, ONCOKB_PX_LEVELS, ONCOKB_GENE_SUMMARY, ONCOKB_VARIANT_SUMMARY, ONCOKB_TUMOR_TYPE_SUMMARY.

Usage:
    python -m altanalyze3.components.bam.oncokb_annotate --input candidate_drivers.tsv \
        --output candidate_drivers_oncokb.tsv --token-file /path/.oncokb_token --tumor-type AML,MDS,MPN,CMML
"""
import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

API_URL = "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange"
BATCH_SIZE = 100
RETRY_DELAYS = [2, 8, 30]

ONCOKB_COLUMNS = [
    'ONCOKB_ONCOGENIC', 'ONCOKB_MUTATION_EFFECT', 'ONCOKB_HIGHEST_LEVEL',
    'ONCOKB_HIGHEST_RESISTANCE_LEVEL', 'ONCOKB_HIGHEST_DX_LEVEL', 'ONCOKB_HIGHEST_PX_LEVEL',
    'ONCOKB_HIGHEST_FDA_LEVEL', 'ONCOKB_TX_LEVELS', 'ONCOKB_DX_LEVELS', 'ONCOKB_PX_LEVELS',
    'ONCOKB_GENE_SUMMARY', 'ONCOKB_VARIANT_SUMMARY', 'ONCOKB_TUMOR_TYPE_SUMMARY',
]


def _resolve_cols(in_cols):
    """Map CHROM/POS/REF/ALT to the actual column names (case-insensitive). Returns dict or exits."""
    lower = {c.lower(): c for c in in_cols}
    out = {}
    for want in ("chrom", "pos", "ref", "alt"):
        if want in lower:
            out[want.upper()] = lower[want]
        elif want == "chrom" and "chr" in lower:
            out["CHROM"] = lower["chr"]
        else:
            sys.exit(f"ERROR: input TSV missing a {want}/{want.upper()} column (have: {in_cols})")
    return out


def build_query(chrom, pos, ref, alt, tumor_type):
    chrom_clean = chrom[3:] if chrom.startswith('chr') else chrom
    pos = int(pos)
    end = pos + len(ref) - 1
    gloc = f"{chrom_clean},{pos},{end},{ref},{alt}"
    q = {"genomicLocation": gloc, "referenceGenome": "GRCh38"}
    if tumor_type:
        q["tumorType"] = tumor_type
    return q


def annotate_batch(queries, token):
    body = json.dumps(queries).encode('utf-8')
    req = Request(API_URL, data=body, method="POST", headers={
        "Content-Type": "application/json", "Authorization": f"Bearer {token}",
        "Accept": "application/json"})
    last_err = None
    for delay in [0] + RETRY_DELAYS:
        if delay:
            time.sleep(delay)
        try:
            with urlopen(req, timeout=180) as resp:
                return json.loads(resp.read())
        except (HTTPError, URLError) as e:
            last_err = e
            sys.stderr.write(f"  [retry] OncoKB request failed ({e}); waiting...\n")
    raise RuntimeError(f"OncoKB API failed after retries: {last_err}")


def _clean(s):
    return (s or '').replace('\t', ' ').replace('\n', ' ').strip()


def merge_responses(per_tt):
    """Merge per-tumor-type responses for one variant (variant-level fields tumor-agnostic;
    Tx/Dx/Px levels + tumor_type_summary aggregated across tumor types)."""
    valid = [(tt, r) for tt, r in per_tt if r]
    if not valid:
        return {col: '' for col in ONCOKB_COLUMNS}
    base = valid[0][1]
    me = base.get('mutationEffect') or {}
    tx, dx, px = set(), set(), set()
    h_sens, h_resist, h_dx, h_px, h_fda, tt_sum = [], [], [], [], [], []
    for tt, resp in valid:
        for t in (resp.get('treatments') or []):
            if t.get('level'):
                tx.add(t['level'])
        for d in (resp.get('diagnosticImplications') or []):
            if d.get('levelOfEvidence'):
                dx.add(d['levelOfEvidence'])
        for p in (resp.get('prognosticImplications') or []):
            if p.get('levelOfEvidence'):
                px.add(p['levelOfEvidence'])
        for key, lst in (('highestSensitiveLevel', h_sens), ('highestResistanceLevel', h_resist),
                         ('highestDiagnosticImplicationLevel', h_dx),
                         ('highestPrognosticImplicationLevel', h_px), ('highestFdaLevel', h_fda)):
            if resp.get(key):
                lst.append(resp[key])
        t_text = _clean(resp.get('tumorTypeSummary'))
        if t_text:
            tt_sum.append(f"[{tt}] {t_text}")

    def best(lst):
        return min(lst) if lst else ''
    return {
        'ONCOKB_ONCOGENIC': _clean(base.get('oncogenic')),
        'ONCOKB_MUTATION_EFFECT': _clean(me.get('knownEffect')),
        'ONCOKB_HIGHEST_LEVEL': best(h_sens),
        'ONCOKB_HIGHEST_RESISTANCE_LEVEL': best(h_resist),
        'ONCOKB_HIGHEST_DX_LEVEL': best(h_dx),
        'ONCOKB_HIGHEST_PX_LEVEL': best(h_px),
        'ONCOKB_HIGHEST_FDA_LEVEL': best(h_fda),
        'ONCOKB_TX_LEVELS': ','.join(sorted(tx)),
        'ONCOKB_DX_LEVELS': ','.join(sorted(dx)),
        'ONCOKB_PX_LEVELS': ','.join(sorted(px)),
        'ONCOKB_GENE_SUMMARY': _clean(base.get('geneSummary')),
        'ONCOKB_VARIANT_SUMMARY': _clean(base.get('variantSummary')),
        'ONCOKB_TUMOR_TYPE_SUMMARY': ' | '.join(tt_sum),
    }


def load_token(args):
    if args.token:
        return args.token.strip()
    if args.token_file:
        return Path(args.token_file).read_text().strip()
    env = os.environ.get('ONCOKB_TOKEN')
    return env.strip() if env else None


def annotate_tsv(input_tsv, output_tsv, token, tumor_types):
    with open(input_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        in_cols = list(reader.fieldnames or [])
        rows = list(reader)
    cmap = _resolve_cols(in_cols)
    out_cols = in_cols + [c for c in ONCOKB_COLUMNS if c not in in_cols]
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        open(output_tsv, 'w').write('\t'.join(out_cols) + '\n')
        return 0, 0
    per_tt = {}
    for tt in tumor_types:
        queries = [build_query(r[cmap['CHROM']], r[cmap['POS']], r[cmap['REF']], r[cmap['ALT']], tt)
                   for r in rows]
        resps = []
        for i in range(0, len(queries), BATCH_SIZE):
            print(f"  [{tt}] OncoKB {i + 1}-{min(i + BATCH_SIZE, len(queries))} of {len(queries)}...")
            resps.extend(annotate_batch(queries[i:i + BATCH_SIZE], token))
        per_tt[tt] = resps
    n_onco = n_act = 0
    with open(output_tsv, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for idx, row in enumerate(rows):
            ann = merge_responses([(tt, per_tt[tt][idx]) for tt in tumor_types])
            row.update(ann)
            w.writerow(row)
            if ann['ONCOKB_ONCOGENIC'] in ('Oncogenic', 'Likely Oncogenic'):
                n_onco += 1
            if ann['ONCOKB_HIGHEST_LEVEL']:
                n_act += 1
    return n_onco, n_act


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--token')
    ap.add_argument('--token-file', dest='token_file')
    ap.add_argument('--tumor-type', dest='tumor_type', default='AML,MDS,MPN,CMML')
    args = ap.parse_args()
    token = load_token(args)
    if not token:
        sys.exit("ERROR: OncoKB token required (--token, --token-file, or ONCOKB_TOKEN env)")
    tumor_types = [t.strip() for t in args.tumor_type.split(',') if t.strip()]
    print(f"annotating {args.input} against OncoKB (tumor types: {', '.join(tumor_types)})")
    n_onco, n_act = annotate_tsv(args.input, args.output, token, tumor_types)
    print(f"wrote {args.output}\n  Oncogenic/Likely Oncogenic: {n_onco}\n  with actionable Tx level: {n_act}")


if __name__ == '__main__':
    main()
