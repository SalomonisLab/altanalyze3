"""Additive exact sequence test of a contiguous homology-equivalent BSJ hypothesis.

The observed intronic interval is one explicit candidate, not a reconstructed
transcript or a claim of molecular topology. Frozen source/model/reference pins
and every tested full-query alignment outcome are retained.
"""
import argparse
import datetime
import gzip
import hashlib
import json
import math
from pathlib import Path

from .circRNA_semiglobal_v2 import align_model, rc
from .circRNA_full_query_validation_v2 import contig


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    for key in ('source-shard', 'original-model-dir', 'reference-manifest', 'reference-manifest-sha256', 'output'):
        ap.add_argument('--' + key, required=True)
    ap.add_argument('--shift', type=int, required=True)
    ap.add_argument('--max-queries', type=int, default=0)
    a = ap.parse_args()
    source = Path(a.source_shard)
    bundle_path = Path(a.original_model_dir) / 'PREPARED_MODEL_BUNDLE.json'
    stage_path = Path(a.original_model_dir) / 'EVENT_STAGE_RESULT.json'
    manifest_path = Path(a.reference_manifest)
    if sha(manifest_path) != a.reference_manifest_sha256:
        raise ValueError('Reference manifest byte SHA changed')
    stage = json.loads(stage_path.read_text())
    bundle = json.loads(bundle_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    if stage['input_sha256'] != sha(source) or stage['reference_manifest_sha256'] != sha(manifest_path):
        raise ValueError('Frozen source/reference manifest pin mismatch')
    if bundle['reference_sha256'] != manifest['sha256']:
        raise ValueError('Prepared model/reference SHA mismatch')
    reference = Path(manifest['path'])
    stat = reference.stat()
    if stat.st_size != manifest['size'] or stat.st_mtime_ns != manifest['mtime_ns'] or sha(str(reference) + '.fai') != manifest['fai_sha256']:
        raise ValueError('Reference stat/FAI changed')
    contexts = [json.loads(line) for line in source.open()]
    if len(contexts) != stage['source_contexts'] or sum(len(x['baseline_pairs']) for x in contexts) != stage['raw_pair_rows']:
        raise ValueError('Source roster/count changed')
    selected = []
    for x in contexts:
        shifts = {int(p['slide_nt']) for p in x['baseline_pairs']}
        if len(shifts) != 1:
            raise ValueError('Ambiguous observed offset within a source context')
        delta = shifts.pop()
        for p in x['baseline_pairs']:
            if int(p['aligned_donor']) - int(p['donor']) != delta or int(p['aligned_acceptor']) - int(p['acceptor']) != delta or abs(delta) > int(p['junction_homology_nt']):
                raise ValueError('Observed endpoints not homology equivalent')
        if delta == a.shift:
            selected.append(x)
    if not selected:
        raise ValueError('No source queries for requested shift')
    selected.sort(key=lambda x: x['alignment_query_id'])
    if a.max_queries:
        selected = selected[:a.max_queries]
    event = selected[0]['event']
    if any(x['event']['event_id'] != event['event_id'] for x in selected):
        raise ValueError('Mixed registered events')
    donor = int(event['donor']) + a.shift
    acceptor = int(event['acceptor']) + a.shift
    start0 = min(donor, acceptor) - 1
    end0 = max(donor, acceptor)
    import pysam
    with pysam.FastaFile(str(reference)) as fa:
        forward = fa.fetch(contig(event['chrom']), start0, end0).upper()
        maxread = max(len(x['source']['full_forward_sequence']) for x in selected)
        copies = max(2, math.ceil(maxread / len(forward)) + 2)
        target_length = copies * len(forward)
        center = (start0 + end0) // 2
        linear_start = max(0, center - target_length // 2)
        linear = fa.fetch(contig(event['chrom']), linear_start, linear_start + target_length).upper()
    if len(forward) != end0-start0 or len(linear) != target_length:
        raise ValueError('Reference fetch failed')
    sequence = rc(forward) if event['donor_strand'] == '-' else forward
    direct = sequence * copies
    folded_pattern = sequence + rc(sequence)
    folded = (folded_pattern * math.ceil(target_length / len(folded_pattern)))[:target_length]
    tandem_pattern = sequence * 2 + rc(sequence) * 2
    tandem = (tandem_pattern * math.ceil(target_length / len(tandem_pattern)))[:target_length]
    native = next((m for m in bundle['templates'] if m['model_id'] == 'native_linear_0000'), None)
    if native is None:
        raise ValueError('Frozen native linear comparator missing')
    def model(name, target, boundaries=()):
        return dict(model_id=name, sequence=target, sequence_sha256=hashlib.sha256(target.encode()).hexdigest(), boundaries=list(boundaries))
    models = [model('observed_homology_equivalent_contiguous_direct', direct, range(len(sequence), target_length, len(sequence))),
              model('observed_homology_equivalent_foldback', folded),
              model('observed_homology_equivalent_tandem_equivalent', tandem),
              model('matched_length_genomic_linear', linear),
              native]
    output = Path(a.output)
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    payload = output / 'PER_QUERY_FULL_EXACT_MODELS.jsonl.gz'
    seen = set()
    with gzip.open(payload, 'wt') as handle:
        for x in selected:
            qid = x['alignment_query_id']
            if qid in seen:
                raise ValueError('Duplicate query identity')
            seen.add(qid)
            query = x['source']['full_forward_sequence']
            if not query:
                raise ValueError('Unavailable full query sequence')
            outcomes = [align_model(query, m) for m in models]
            handle.write(json.dumps(dict(alignment_query_id=qid, library=x['library'], query=x['query'],
                                         registered_event=event['event_id'], observed_shift_nt=a.shift,
                                         source_statuses=[p['status'] for p in x['baseline_pairs']],
                                         model_outcomes=outcomes), separators=(',', ':')) + '\n')
    receipt = dict(status='FULL_QUERY_EXACT_CONTIGUOUS_HOMOLOGY_HYPOTHESIS_DIAGNOSTIC',
                   registered_event=event['event_id'], observed_shift_nt=a.shift,
                   aligned_donor=donor, aligned_acceptor=acceptor, source_contexts=len(selected),
                   source_query_ids_sha256=hashlib.sha256('\n'.join(sorted(seen)).encode()).hexdigest(),
                   source_shard_sha256=sha(source), prepared_model_bundle_sha256=sha(bundle_path),
                   reference_manifest_sha256=sha(manifest_path), reference_fasta_sha256=manifest['sha256'],
                   contiguous_reference_interval=dict(chrom=event['chrom'], start0=start0, end0=end0, strand=event['donor_strand']),
                   copies=copies, model_summaries=[dict(model_id=m['model_id'], length=len(m['sequence']), sequence_sha256=m['sequence_sha256'], boundaries=m.get('boundaries', [])) for m in models],
                   per_query_payload_sha256=sha(payload), code_sha256=sha(__file__),
                   checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   limits=['Contiguous interval is one explicit unspliced hypothesis; other internal chains remain possible.',
                           'Exact sequence fit is not proof of physical circularity or corrected circle abundance.',
                           'No original source, raw count, feature, native statistic, or q value is changed.'])
    (output / 'EXACT_DIAGNOSTIC_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], contexts=len(selected), shift=a.shift, models=len(models))))


if __name__ == '__main__':
    main()
