"""Preserve the native uid__library pooling while retaining physical source IDs.

Raw scan files retain source-qualified IDs. The native reporting view restores the
original grouping/identifiers, with an explicit reversible source-to-folder map.
The source-aware ledger is separate from the native report's historical counting.
"""
from pathlib import Path
import json
from .manifest import read_tsv, write_tsv, atomic_json


def prepare(config):
    raw = Path(config['output'])/'oncofusion'
    out = raw/'native_profile'
    metadata = []
    mapping = []
    for n, source in enumerate(config['sources'], 1):
        sid = source['source_id']
        uid = source['uid']+'__'+source['library']
        folder = out/'samples'/f'{n:02d}_{uid}'
        folder.mkdir(parents=True, exist_ok=True)
        for suffix in ('fusion_reads.tsv', 'flt3_events.tsv'):
            path = raw/'samples'/sid/(sid+'.'+suffix)
            rows = read_tsv(path)
            fields = path.open().readline().rstrip('\n').split('\t')
            for row in rows:
                if row['uid'] != sid:
                    raise ValueError('Unresolved source identity in native oncofusion scan')
                row['uid'] = uid
            write_tsv(folder/(uid+'.'+suffix), fields, rows)
        summary = json.loads((raw/'samples'/sid/(sid+'.summary.json')).read_text())
        if summary['uid'] != sid:
            raise ValueError('Oncofusion source summary identity differs')
        summary['uid'] = uid
        atomic_json(folder/(uid+'.summary.json'), summary)
        metadata.append(dict(uid=source['uid'], library=source['library'], groups=source['group'], bam=source['bam']))
        mapping.append(dict(source_id=sid, native_folder=folder.name, native_uid=uid, native_bam_index=n))
    write_tsv(out/'metadata.tsv', list(metadata[0]), metadata)
    write_tsv(out/'source_map.tsv', list(mapping[0]), mapping)
    return out


def itd_groups(events):
    """Exactly combine.py's greedy per-sample grouping and tandem/in-frame rule.

    Return membership as well as the native verdict; never use a per-read tandem
    fraction in place of the native cluster median and modal insertion length.
    """
    from collections import Counter
    for uid in sorted({e['uid'] for e in events}):
        ins=sorted((e for e in events if e['uid']==uid and e['kind']=='insertion'),
                   key=lambda e:(int(e['ref_pos']),int(e['length'])))
        used=set()
        for i,e in enumerate(ins):
            if i in used:continue
            group=[e];used.add(i)
            for j in range(i+1,len(ins)):
                if j not in used and abs(int(ins[j]['ref_pos'])-int(e['ref_pos']))<=10 and abs(int(ins[j]['length'])-int(e['length']))<=6:
                    group.append(ins[j]);used.add(j)
            length=Counter(int(x['length']) for x in group).most_common(1)[0][0]
            median=sorted(float(x['tandem_8mer_fraction']) for x in group)[len(group)//2]
            yield dict(uid=uid, length_mode=length, ref_pos=min(int(x['ref_pos']) for x in group),
                       tandem_8mer_fraction_median=median, accepted=length%3==0 and median>=.8 and length>=3,
                       rows=group)
