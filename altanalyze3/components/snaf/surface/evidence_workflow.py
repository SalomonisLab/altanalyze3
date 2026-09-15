"""Connect sample-supported isoforms to the regular SNAF-B surface workflow.

Inputs are the initialized Alt91 reference and the complete tumor junction matrix.
Long-read validation remains a separate, downstream step.
"""
import csv
import hashlib
import json
import logging
from pathlib import Path
import pickle

import numpy as np

from . import evidence_isoform as E
from . import predict_isoform as P

logger = logging.getLogger(__name__)


def reference_models(exonlist, coordinates, genes):
    """Assemble reference ENST chains from AltAnalyze's inclusive subexon blocks."""
    models, annotation, loci = {}, {}, {}
    for gene in genes:
        blocks = coordinates.get(gene, {})
        d = annotation.setdefault(gene, {'exons': set(), 'introns': set()})
        for name, attrs in blocks.items():
            s, e = sorted((int(attrs[2]), int(attrs[3])))
            if name.startswith('E'):
                d['exons'].add((s, e))
            elif name.startswith('I'):
                d['introns'].add((s - 1, e + 1))
    for row in exonlist.loc[exonlist['EnsGID'].isin(genes)].itertuples(index=False):
        if not str(row.EnsTID).startswith('ENST'):
            continue
        blocks = coordinates.get(row.EnsGID, {})
        names = str(row.Exons).split('|')
        if not names or any(n not in blocks for n in names):
            continue
        attrs = [blocks[n] for n in names]
        locus = {(a[0], a[1]) for a in attrs}
        if len(locus) != 1:
            continue
        chrom, strand = next(iter(locus))
        previous = loci.setdefault(row.EnsGID, (chrom, strand))
        if previous != (chrom, strand):
            raise ValueError('Reference models span multiple loci: ' + row.EnsGID)
        chain = []
        for s, e in sorted(tuple(sorted((int(a[2]), int(a[3])))) for a in attrs):
            if chain and s <= chain[-1][1] + 1:
                chain[-1] = (chain[-1][0], max(e, chain[-1][1]))
            else:
                chain.append((s, e))
        models.setdefault(row.EnsGID, {})[row.EnsTID] = (strand, chain)
    return models, annotation, loci


def resolve_junction(uid, coordinates):
    """Resolve a cis junction UID, including explicit novel-site coordinates."""
    uid = str(uid).split('=')[0]
    if uid.count(':') != 1:
        return None  # Fusion construction requires a different model.
    gene, event = uid.split(':')
    blocks = coordinates.get(gene, {})
    if not blocks or event.count('-') != 1:
        return None
    sites, loci = [], []
    for i, token in enumerate(event.split('-')):
        name, _, explicit = token.partition('_')
        attrs = blocks.get(name)
        if attrs is None and name.startswith('U') and explicit:
            attrs = next(iter(blocks.values()))
        if attrs is None:
            return None
        chrom, strand = attrs[:2]
        loci.append((chrom, strand))
        sites.append(int(explicit) if explicit else int(attrs[
            3 if (i == 0) == (strand == '+') else 2]))
    if loci[0] != loci[1] or sites[0] == sites[1]:
        return None
    return gene, loci[0][0], tuple(sorted(sites))


def matrix_evidence(matrix, coordinates, genes):
    """Retain all measured junctions in target genes, merging coordinate aliases by max."""
    counts, labels, skipped = {}, {}, []
    if not matrix.columns.is_unique:
        raise ValueError('Evidence requires unique sample column names')
    for uid, row in matrix.loc[[str(u).split(':')[0] in genes for u in matrix.index]].iterrows():
        resolved = resolve_junction(uid, coordinates)
        if resolved is None:
            skipped.append(str(uid))
            continue
        gene, chrom, j = resolved
        v = np.asarray(row, dtype=float)
        if not np.isfinite(v).all() or (v < 0).any():
            raise ValueError('Invalid junction counts: ' + str(uid))
        gc = counts.setdefault(gene, {})
        gc[j] = np.maximum(gc[j], v) if j in gc else v
        labels.setdefault(gene, {}).setdefault(j, str(uid).split('=')[0])
    return counts, labels, skipped


def _gene_fetch(record, strand):
    """Expose Alt91's strand-oriented, 2-kb-flanked gene sequence in genomic orientation."""
    origin = int(record[1]) - 2000
    sequence = record[3]
    if strand == '-':
        sequence = P.reverse_complement(sequence)
    def fetch(s, e):
        a, b = s - origin, e - origin + 1
        return sequence[a:b] if s >= 1 and 0 <= a < b <= len(sequence) else ''
    return fetch


def run_surface_evidence(uids, outdir, matrix, *, genome_fasta=None, ranker_path=None,
                         first_exon_junctions=(), sample=None, min_reads=3,
                         method='learned', max_edits=None, top_k=5, tmhmm=False, software_path=None,
                         serialize=True):
    """Predict and rank isoforms, then apply existing surface-antigen checks.

    Requires surface.initialize(). Preserves predicted ORFs and NMD annotations;
    legacy find_orf/orf_check would replace them using different transcript geometry.
    """
    from . import main as M
    if method not in ('learned', 'evidence'):
        raise ValueError('Unknown isoform method: ' + method)
    if ranker_path and method != 'learned':
        raise ValueError('A ranking model requires the learned isoform method')
    max_edits = max_edits if max_edits is not None else (4 if method == 'learned' else 2)
    if matrix is None:
        raise ValueError('Evidence isoform inference requires the complete junction count matrix')
    if top_k < 1 or max_edits < 1 or min_reads <= 0:
        raise ValueError('top_k, max_edits and min_reads must be positive')
    sample_names = list(map(str, matrix.columns))
    if sample is not None and sample not in sample_names:
        raise ValueError('Sample absent from junction matrix: ' + str(sample))
    sample_index = sample_names.index(sample) if sample is not None else None
    uids = list(uids)
    genes = {u[0].split(':')[0] for u in uids}
    models, annotation, loci = reference_models(M.df_exonlist, M.dict_exonCoords, genes)
    counts, labels, skipped = matrix_evidence(matrix, M.dict_exonCoords, genes)
    first_exon_junctions = set(first_exon_junctions)
    model = E.load_ranker(ranker_path) if method == 'learned' else None
    fa = None
    if genome_fasta:
        import pysam
        fa = pysam.FastaFile(str(genome_fasta))
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    artifacts = out / 'synthetic_isoforms'
    artifacts.mkdir(exist_ok=True)
    records, gtfs, mrnas, proteins, results = [], [], [], [], []
    try:
        for uid, score, df, ed, freq in uids:
            sa = M.SurfaceAntigen(uid, score, df, ed, freq, False)
            sa.detect_type()
            sa.isoform_method = method
            sa.synthetic_predictions = []
            resolved = resolve_junction(uid, M.dict_exonCoords)
            ranked, status = [], 'unresolved_coordinates'
            if resolved is not None:
                gene, chrom, target = resolved
                status = 'unresolved'
                gc = counts.get(gene, {})
                if gene in models and (sample_index is None or gc):
                    if fa is not None:
                        contig = chrom if chrom in fa.references else chrom.removeprefix('chr')
                        if contig not in fa.references:
                            raise ValueError('Chromosome absent from genome: ' + chrom)
                        fetch = lambda s, e: fa.fetch(contig, s - 1, e)
                    else:
                        fetch = _gene_fetch(M.dict_fa[gene], loci[gene][1])
                    ev = E.Evidence(gc, min_reads=min_reads, sample_index=sample_index)
                    builder = E.generate_hypotheses if method == 'learned' else E.generate_candidates
                    candidates = builder(models[gene], target, fetch, ev,
                        exon_annotation=annotation[gene], junction_label=uid,
                        labels=labels.get(gene), max_edits=max_edits,
                        first_exon=uid in first_exon_junctions)
                    ranked = E.rank_candidates(candidates, model)
            sa.synthetic_predictions = [c.prediction for _, c in ranked[:top_k]]
            sa.full_length = [p.mrna for p in sa.synthetic_predictions]
            sa.orft = [p.orf for p in sa.synthetic_predictions]
            sa.orfp = [p.protein for p in sa.synthetic_predictions]
            sa.full_length_attrs = []
            sa.nmd = ['*' if p.nmd else '#' for p in sa.synthetic_predictions]
            coding = 'protein_coding' in M.dict_biotype.get(uid.split(':')[0], {}).values()
            sa.translatability = ['#' if coding else '*'] * len(sa.synthetic_predictions)
            sa.junction = ''
            for rank, (value, candidate) in enumerate(ranked[:top_k], 1):
                p = candidate.prediction
                signature = (gene, p.strand, p.chain, p.cds_start, p.cds_end)
                pid = 'SNAF_SYN_' + hashlib.sha256(repr(signature).encode()).hexdigest()[:20]
                sa.full_length_attrs.append(pid)
                records.append(dict(junction_id=uid, rank=rank, score=value,
                    status='hypothesis' if candidate.joint_samples else 'reference_only_hypothesis',
                    artifact_id=pid, construction=candidate.construction,
                    supporting_samples=','.join(sample_names[i] for i in candidate.joint_samples),
                    **P.prediction_row(p, uid)))
                if rank == 1:
                    donor = target[1] if p.strand == '-' else target[0]
                    offset = P.offset_of(p.chain, p.strand == '-', donor) + 1
                    sa.junction = p.mrna[max(0, offset - 100):offset] + ',' + p.mrna[offset:offset + 100]
                    for s, e in p.chain:
                        gtfs.append('\t'.join([chrom, 'SNAF_evidence', 'exon', str(s), str(e),
                            '.', p.strand, '.', 'gene_id "%s"; transcript_id "%s"; first_exon_boundary_source "%s";'
                            % (gene, pid, p.first_exon_boundary_source)]))
                    mrnas.append('>' + pid + '\n' + p.mrna)
                    proteins.append('>' + pid + '\n' + p.protein)
            if not ranked:
                sa.comments.append('synthetic_isoform_' + status)
                records.append(dict(junction_id=uid, rank=0, score='', status=status,
                    artifact_id='', construction='', supporting_samples='',
                    **{k: '' for k in P.PREDICTION_COLUMNS}))
            sa.align_uniprot(tmhmm=tmhmm, software_path=software_path)
            results.append(sa)
    finally:
        if fa is not None:
            fa.close()
    if serialize:
        with (out / 'surface_antigen_sr.p').open('wb') as fh:
            pickle.dump(results, fh)
    columns = ['junction_id', 'rank', 'score', 'status', 'artifact_id', 'construction',
               'supporting_samples'] + P.PREDICTION_COLUMNS
    with (artifacts / 'predictions.tsv').open('w') as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(records)
    for name, lines in [('top1.gtf', gtfs), ('top1.mrna.fa', mrnas), ('top1.protein.fa', proteins)]:
        (artifacts / name).write_text('\n'.join(dict.fromkeys(lines)) + '\n')
    manifest = dict(method=method, min_reads=min_reads, max_edits=max_edits, top_k=top_k,
        sample=sample, sample_columns=sample_names, first_exon_junctions=sorted(first_exon_junctions),
        first_exon_fallback_nt=250, ranker_path=str(ranker_path) if ranker_path else
        ('bundled:isoform_ranker.json' if method == 'learned' else None),
        sequence_source=str(genome_fasta) if genome_fasta else 'Alt91 gene FASTA with 2000-nt flanks',
        unresolved_evidence_uids=skipped, feature_names=E.FEATURE_NAMES,
        long_reads_used_for_inference=False, co_detection_is_phasing=False)
    (artifacts / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    logger.info('Synthetic isoform inference: %d/%d events reconstructed; %d evidence rows unresolved',
                sum(bool(sa.synthetic_predictions) for sa in results), len(results), len(skipped))
    return results
