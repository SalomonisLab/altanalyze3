"""Additive diagnosis of observed BSJ microhomology coordinate equivalence.

This only prepares alternate observed-chain hypotheses. It does not replace the
registered event, score alignments, classify topology, or change baseline counts.
"""
import argparse
import collections
import copy
import datetime
import hashlib
import json
from pathlib import Path

from .circRNA_candidate_models_v2 import observed_model_chains
from .circRNA_full_query_validation_v2 import event_key


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-shard', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    source = Path(a.source_shard)
    contexts = [json.loads(line) for line in source.open()]
    assert contexts
    registered = contexts[0]['event']
    assert all(x['event']['event_id'] == registered['event_id'] for x in contexts)
    groups = collections.defaultdict(list)
    identities = set()
    offsets = collections.Counter()
    for x in contexts:
        identity = (x['library'], x['query'])
        assert identity not in identities
        identities.add(identity)
        pairs = x['baseline_pairs']
        assert pairs
        shifts = set()
        for p in pairs:
            donor_delta = int(p['aligned_donor']) - int(p['donor'])
            acceptor_delta = int(p['aligned_acceptor']) - int(p['acceptor'])
            slide = int(p['slide_nt'])
            homology = int(p['junction_homology_nt'])
            if donor_delta != acceptor_delta or donor_delta != slide or abs(slide) > homology:
                raise ValueError(f'Unresolved nonequivalent observed endpoint: {identity}')
            shifts.add(slide)
            offsets[(donor_delta, acceptor_delta, homology)] += 1
        if len(shifts) != 1:
            raise ValueError(f'Multiple aligned endpoints within a source context: {identity}')
        delta = shifts.pop()
        altered = copy.copy(x)
        event = dict(registered)
        event['donor'] = int(registered['donor']) + delta
        event['acceptor'] = int(registered['acceptor']) + delta
        event['event_id'] = event_key(event)
        altered['event'] = event
        groups[delta].append(altered)
    result = []
    for delta, subset in sorted(groups.items()):
        models, outcomes = observed_model_chains(subset[0]['event'], subset)
        status = collections.Counter(x['status'] for x in outcomes)
        assert len(outcomes) == len(subset)
        result.append(dict(shift_nt=delta, aligned_event=subset[0]['event'], contexts=len(subset),
                           observed_chain_models=len(models), source_model_statuses=dict(status),
                           model_chain_block_counts=dict(collections.Counter(len(m['blocks']) for m in models)),
                           candidate_source_contexts=sum(x['status'] == 'candidate_hypothesis_constructed' for x in outcomes)))
    out = dict(status='DIAGNOSTIC_ONLY_HOMOLOGY_EQUIVALENT_OBSERVED_CHAINS_PREPARED',
               registered_event=registered, source_contexts=len(contexts), raw_pair_rows=sum(len(x['baseline_pairs']) for x in contexts),
               source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
               observed_offset_counts=[dict(donor_delta=k[0], acceptor_delta=k[1], homology_nt=k[2], raw_pairs=v) for k, v in sorted(offsets.items())],
               offset_groups=result, checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
               limits=['Registered coordinates and baseline counts remain unchanged.',
                       'Observed-chain construction is a sequence hypothesis, not evidence of physical circularity.',
                       'Full-query exact alignment and alternative path comparison have not been run by this preparation step.'])
    Path(a.output).write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(dict(status=out['status'], contexts=out['source_contexts'], groups=result)))


if __name__ == '__main__':
    main()
