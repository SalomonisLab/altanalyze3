#!/usr/bin/env python3
"""Merge SURFY and SurfaceGenie into one Ensembl-keyed cell-surface gene database.

Inputs
  SURFY-Ensembl.txt   UID (symbol) + Ensembl gene id, CRLF, 2,882 data rows.
  SurfaceGenie.txt    UniProt accession, SPC (Surface Prediction Consensus, 0-4),
                      geneName (primary symbol then aliases, space separated),
                      transmembrane, subcellular location, CD antigen, CSPA
                      experimental evidence, HLA. CRLF, 2,335 data rows.

Ensembl links
  SNAF-B gates on Ensembl GENE ids from the Ensembl-91 gene models in Alt91_db, so every
  row needs an ENSG. SURFY already carries one. SurfaceGenie carries only symbols and a
  UniProt accession, so it is mapped in this order, and the route used is recorded per row:

    1. primary symbol  -> Ensembl-91 Ensembl-Symbol.txt
    2. alias symbol    -> the same table
    3. UniProt accession -> EnsMart100 Hs_Ensembl-UniProt.txt (last resort; a different
       Ensembl release, so it is reported separately and never silently mixed in)

  A symbol can name MORE THAN ONE Ensembl gene. Every one is kept as its own row: the
  whitelist gates which genes SNAF-B may consider, and dropping a paralog or a patch-scaffold
  copy would silently exclude a real candidate. The count of symbols that expanded this way
  is reported.

No filtering
  EVERY gene named by either source is included. No SPC, CSPA or CD threshold is applied.
  SPC, CSPA experimental evidence and CD antigen are carried as ANNOTATION columns only --
  they describe a gene, they never exclude one. (An earlier version of this script applied
  an unauthorised SPC>=3 cut, which dropped genes such as LAT; that cut is removed.)

Outputs (into this directory)
  surface_union_annotated.txt   every gene from either source, with provenance
  surface_union_all.txt         the ENSG + symbol gene table for `snaf-b --surface_db`
  surface_union_report.txt      counts with denominators
"""
import os
import re
import csv
import sys
from collections import OrderedDict, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SURFY = os.path.join(HERE, 'SURFY-Ensembl.txt')
GENIE = os.path.join(HERE, 'SurfaceGenie.txt')
SYMBOL_DB = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Ensembl-Symbol.txt'
UNIPROT_DB = '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart100/uniprot/Hs/Hs_Ensembl-UniProt.txt'
ALT91_EXONS = ('/Users/saljh8/Dropbox/SNAF/MDS/snaf_test_2026-07-04/db/data/Alt91_db/'
               'Hs_Ensembl_exon_add_col.txt')

OUT_ANNOT = os.path.join(HERE, 'surface_union_annotated.txt')
OUT_DB = os.path.join(HERE, 'surface_union_all.txt')
OUT_REPORT = os.path.join(HERE, 'surface_union_report.txt')

ENSG = re.compile(r'^ENSG\d+')
report = []


def say(m):
    print(m)
    report.append(m)


def rows(path):
    with open(path, 'r', encoding='utf-8-sig', errors='replace', newline='') as f:
        for line in f:
            yield line.rstrip('\r\n')


# --------------------------------------------------------------- symbol -> ENSG
def load_symbol_map(path):
    fwd = defaultdict(set)
    for line in rows(path):
        if not line:
            continue
        t = line.split('\t')
        if len(t) < 2:
            continue
        ensg, sym = t[0].strip(), t[1].strip()
        if ENSG.match(ensg) and sym:
            fwd[sym.upper()].add(ensg)
    return fwd


def load_uniprot_map(path):
    """UniProt accession -> {ENSG}, from the EnsMart100 Ensembl<TAB>UniProt table."""
    out = defaultdict(set)
    if not os.path.exists(path):
        return out
    for i, line in enumerate(rows(path)):
        if i == 0 or not line:
            continue
        t = line.split('\t')
        if len(t) < 2:
            continue
        ensg, acc = t[0].strip(), t[1].strip().split('-')[0]
        if ENSG.match(ensg) and acc:
            out[acc].add(ensg)
    return out


def load_annotation_symbols(paths):
    """symbol -> {ENSG} from AltAnalyze 'Hs_Ensembl-annotations_simple.txt' tables
    (Ensembl Gene ID / Description / Gene name). A second, independent route to an Ensembl
    id for symbols the Ensembl-Symbol table does not carry."""
    out = defaultdict(set)
    for path in paths:
        if not os.path.exists(path):
            continue
        for i, line in enumerate(rows(path)):
            if i == 0 or not line:
                continue
            t = line.split('\t')
            if len(t) < 3:
                continue
            ensg, sym = t[0].strip(), t[2].strip()
            if ENSG.match(ensg) and sym:
                out[sym.upper()].add(ensg)
    return out


def load_genie_aliases(path):
    """alias symbol -> primary symbol, from SurfaceGenie's space-separated geneName field.
    Recovers retired HGNC symbols (C10orf54 -> VSIR, BAI3 -> ADGRB3, ELTD1 -> ADGRL4 ...)
    that SURFY still uses but the Ensembl tables index under the current name."""
    out = {}
    for i, line in enumerate(rows(path)):
        if i == 0 or not line:
            continue
        t = line.split('\t')
        if len(t) < 3:
            continue
        names = [n for n in t[2].replace('/', ' ').split() if n]
        for alias in names[1:]:
            out.setdefault(alias.upper(), names[0])
    return out


def main():
    for p in (SURFY, GENIE, SYMBOL_DB):
        if not os.path.exists(p):
            sys.exit('missing required input: {}'.format(p))

    sym2ensg = load_symbol_map(SYMBOL_DB)
    acc2ensg = load_uniprot_map(UNIPROT_DB)
    ann2ensg = load_annotation_symbols([
        '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl-annotations_simple.txt',
        '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl-annotations_simple.txt'])
    genie_alias = load_genie_aliases(GENIE)

    def resolve_symbol(sym):
        """(set of ENSG, route) for a gene symbol, trying every table before giving up."""
        if not sym:
            return set(), ''
        u = sym.upper()
        if u in sym2ensg:
            return sym2ensg[u], 'symbol_ens91'
        if u in ann2ensg:
            return ann2ensg[u], 'symbol_annotations'
        cur = genie_alias.get(u)
        if cur:
            c = cur.upper()
            if c in sym2ensg:
                return sym2ensg[c], 'symbol_renamed_ens91'
            if c in ann2ensg:
                return ann2ensg[c], 'symbol_renamed_annotations'
        return set(), ''

    say('Ensembl-91 symbol table : {}'.format(SYMBOL_DB))
    say('  symbols               : {}'.format(len(sym2ensg)))
    say('  symbols naming >1 gene: {} / {}'.format(
        sum(1 for v in sym2ensg.values() if len(v) > 1), len(sym2ensg)))
    say('UniProt fallback table  : {}{}'.format(
        UNIPROT_DB, '' if acc2ensg else '   (ABSENT -- fallback unavailable)'))
    say('  accessions            : {}'.format(len(acc2ensg)))
    say('')

    # records keyed by ENSG
    genes = OrderedDict()

    def touch(ensg):
        if ensg not in genes:
            genes[ensg] = {'ensg': ensg, 'symbol': '', 'in_surfy': 0, 'in_surfacegenie': 0,
                           'spc': '', 'cspa_experimental': '', 'cd_antigen': '',
                           'uniprot': '', 'ensembl_link': '', 'ensembl_copy': ''}
        return genes[ensg]

    # ---- SURFY -------------------------------------------------------------
    n_surfy_rows = n_surfy_no_id = 0
    surfy_ensg = set()
    surfy_recovered = 0
    surfy_unmapped = []
    for i, line in enumerate(rows(SURFY)):
        if i == 0 or not line.strip():
            continue
        n_surfy_rows += 1
        t = line.split('\t')
        sym = t[0].strip() if t else ''
        ens = t[1].strip() if len(t) > 1 else ''
        targets = []
        link = 'surfy_direct'
        if ENSG.match(ens):
            targets = [ens]
        else:
            n_surfy_no_id += 1
            hit, route = resolve_symbol(sym)
            if hit:
                targets = sorted(hit)
                link = 'surfy_recovered_' + route
                surfy_recovered += 1
            else:
                surfy_unmapped.append(sym)
        for e in targets:
            r = touch(e)
            r['in_surfy'] = 1
            r['symbol'] = r['symbol'] or sym
            r['ensembl_link'] = r['ensembl_link'] or link
            surfy_ensg.add(e)

    say('SURFY {}'.format(SURFY))
    say('  data rows                          : {}'.format(n_surfy_rows))
    say('  rows with no Ensembl id in the file: {}'.format(n_surfy_no_id))
    say('  of those, recovered by symbol      : {} / {}'.format(surfy_recovered, n_surfy_no_id))
    say('  unique Ensembl genes               : {}'.format(len(surfy_ensg)))
    if surfy_unmapped:
        say('  STILL unmapped (no Ensembl id in any table, cannot enter an ENSG-keyed '
            'database): {}'.format(len(surfy_unmapped)))
        say('    {}'.format(', '.join(sorted(surfy_unmapped))))
    say('')

    # ---- SurfaceGenie ------------------------------------------------------
    n_genie_rows = 0
    spc_hist = defaultdict(int)
    link_hist = defaultdict(int)
    unmapped = []
    genie_ensg = set()
    n_expanded = 0
    for i, line in enumerate(rows(GENIE)):
        if i == 0 or not line.strip():
            continue
        n_genie_rows += 1
        t = line.split('\t')
        if len(t) < 8:
            t = t + [''] * (8 - len(t))
        acc, spc, gene_field, _tm, _cc, cd, cspa, _hla = t[:8]
        acc = acc.strip()
        spc = spc.strip()
        spc_hist[spc] += 1
        names = [n for n in gene_field.replace('/', ' ').split() if n]
        primary = names[0] if names else ''
        targets, route = resolve_symbol(primary)
        link = 'genie_primary_' + route if targets else ''
        if not targets:
            for alias in names[1:]:
                targets, route = resolve_symbol(alias)
                if targets:
                    link = 'genie_alias_' + route
                    break
        if not targets and acc in acc2ensg:
            targets = acc2ensg[acc]
            link = 'genie_uniprot_ensmart100'
        if not targets:
            unmapped.append((acc, primary, spc))
            continue
        if len(targets) > 1:
            n_expanded += 1
        link_hist[link] += 1
        for e in sorted(targets):
            r = touch(e)
            r['in_surfacegenie'] = 1
            r['symbol'] = r['symbol'] or primary
            r['spc'] = spc
            r['cd_antigen'] = '' if cd.strip() in ('NA', '') else cd.strip()
            r['cspa_experimental'] = '' if cspa.strip() in ('NA', '') else cspa.strip()
            r['uniprot'] = acc
            if not r['ensembl_link']:
                r['ensembl_link'] = link
            genie_ensg.add(e)

    say('SurfaceGenie {}'.format(GENIE))
    say('  data rows                : {}'.format(n_genie_rows))
    say('  SPC distribution         : {}'.format(
        ', '.join('SPC={}:{}'.format(k, spc_hist[k]) for k in sorted(spc_hist))))
    for k in sorted(link_hist):
        say('  mapped via {:<24s}: {}'.format(k, link_hist[k]))
    say('  rows naming >1 Ensembl gene (all kept): {}'.format(n_expanded))
    say('  rows with NO Ensembl link : {} / {}'.format(len(unmapped), n_genie_rows))
    say('  unique Ensembl genes      : {}'.format(len(genie_ensg)))
    if unmapped:
        say('  first unmapped (accession, symbol, SPC): {}'.format(unmapped[:8]))
    say('')

    # ---- mark additional Ensembl copies of a multi-mapped symbol -------------
    # A symbol naming several Ensembl genes is usually one real locus plus alt-scaffold or
    # patch duplicates. Measured against the 3,591,344-junction pediatric AML matrix, only
    # 15 / 830 such extra copies carry ANY junction, versus 2,646 / 3,179 first-listed
    # genes. They are kept -- an inert whitelist entry costs nothing, while dropping them
    # would lose those 15 -- but they are labelled so they can be filtered.
    by_symbol = defaultdict(list)
    for r in genes.values():
        if r['symbol']:
            by_symbol[r['symbol']].append(r)
    n_extra = 0
    for sym, rs in by_symbol.items():
        rs.sort(key=lambda r: r['ensg'])
        for k, r in enumerate(rs):
            r['ensembl_copy'] = 'primary' if k == 0 else 'additional'
            if k:
                n_extra += 1
    say('Multi-mapped symbols')
    say('  symbols naming >1 Ensembl gene   : {}'.format(
        sum(1 for rs in by_symbol.values() if len(rs) > 1)))
    say('  additional Ensembl copies kept   : {} (labelled ensembl_copy=additional)'.format(n_extra))
    say('')

    both = [r for r in genes.values() if r['in_surfy'] and r['in_surfacegenie']]
    only_s = [r for r in genes.values() if r['in_surfy'] and not r['in_surfacegenie']]
    only_g = [r for r in genes.values() if r['in_surfacegenie'] and not r['in_surfy']]
    hc = list(genes.values())

    say('Union (denominator = unique Ensembl genes)')
    say('  in both sources          : {}'.format(len(both)))
    say('  SURFY only               : {}'.format(len(only_s)))
    say('  SurfaceGenie only        : {}'.format(len(only_g)))
    say('  UNION total              : {}'.format(len(genes)))
    say('  NO threshold is applied: all {} genes are written to the database'.format(len(genes)))
    spc_ann = defaultdict(int)
    for r in genes.values():
        spc_ann[r['spc'] if r['spc'] != '' else 'not-in-SurfaceGenie'] += 1
    say('  SPC annotation across the union (informational only, nothing is excluded):')
    for k in sorted(spc_ann, key=str):
        say('    SPC={:<22s} {}'.format(str(k), spc_ann[k]))
    say('')

    # ---- can SNAF-B model them? --------------------------------------------
    if os.path.exists(ALT91_EXONS):
        modelled = set()
        for i, line in enumerate(rows(ALT91_EXONS)):
            if i == 0:
                continue
            g = line.split('\t', 1)[0]
            if g:
                modelled.add(g)
        n_hc_model = sum(1 for r in hc if r['ensg'] in modelled)
        say('Ensembl-91 gene models ({} genes in Alt91_db)'.format(len(modelled)))
        say('  union genes with a v91 model : {} / {}'.format(n_hc_model, len(hc)))
        say('  (genes without one cannot be scored by SNAF-B whatever the whitelist says)')
        say('')

    # ---- write --------------------------------------------------------------
    fields = ['ensg', 'symbol', 'in_surfy', 'in_surfacegenie', 'spc', 'cspa_experimental',
              'cd_antigen', 'uniprot', 'ensembl_link', 'ensembl_copy']
    with open(OUT_ANNOT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for e in sorted(genes):
            w.writerow(genes[e])

    with open(OUT_DB, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['UID', 'Ensembl'])
        for e in sorted(genes):
            w.writerow([genes[e]['symbol'] or e, e])

    say('Outputs')
    say('  {}   ({} rows + header)'.format(OUT_ANNOT, len(genes)))
    say('  {}    ({} rows + header, ready for `snaf-b --surface_db`)'.format(OUT_DB, len(hc)))
    say('  {}'.format(OUT_REPORT))
    with open(OUT_REPORT, 'w') as f:
        f.write('\n'.join(report) + '\n')


if __name__ == '__main__':
    main()
