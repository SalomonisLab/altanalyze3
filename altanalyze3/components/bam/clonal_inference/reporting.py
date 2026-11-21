"""Export complete multi-event genotypes with evidence and denominator annotations."""
from pathlib import Path
import numpy as np
import pandas as pd
from .integrity import require


def export_point_assignments(X, cells, events, membership, genotypes, output):
    output = Path(output)
    G = genotypes.copy()
    G.index = G.index.astype(str)
    G.columns = events.feature
    membership = np.asarray(membership).astype(str)
    require(len(membership) == len(cells), 'Missing cell assignments')
    require(set(membership) <= set(G.index), 'Unknown cluster identity')
    informative = (X >= 0).any(axis=1)
    bits = pd.Series([''.join(map(str, row)) for row in (G.to_numpy() >= .5).astype(int)], index=G.index)
    assignment = cells.copy()
    assignment['native_clone'] = membership
    assignment['display_clone'] = np.where(informative, membership, 'Unresolved_no_evidence')
    assignment['full_genotype_group'] = np.where(informative, bits.loc[membership], 'Unresolved_no_evidence')
    assignment['observed_events'] = (X >= 0).sum(1)
    assignment['contradictory_calls'] = ((X >= 0) & (X != (G.loc[membership].to_numpy() >= .5))).sum(1)
    assignment.to_csv(output / 'cell_assignments.tsv.gz', sep='\t', index=False)
    G.to_csv(output / 'clone_genotype_parameters.tsv', sep='\t')
    evidence = []
    for clone in G.index:
        selected = membership == clone
        for j, feature in enumerate(events.feature):
            n_mut = int((X[selected, j] == 1).sum())
            n_wt = int((X[selected, j] == 0).sum())
            evidence.append(dict(clone=clone, feature=feature, genotype_parameter=float(G.loc[clone, feature]),
                                 native_binary_state=int(G.loc[clone, feature] >= .5),
                                 observed_MUT=n_mut, observed_WT=n_wt, unknown=int(selected.sum()) - n_mut - n_wt,
                                 no_observed_evidence=(n_mut + n_wt == 0)))
    pd.DataFrame(evidence).to_csv(output / 'clone_event_evidence.tsv', sep='\t', index=False)
    for key, filename in [('display_clone', 'clone_frequencies.tsv'), ('full_genotype_group', 'genotype_group_frequencies.tsv')]:
        frequency = assignment.groupby(['sample', key], sort=False).size().rename('cells').reset_index()
        frequency['all_cells'] = frequency['sample'].map(cells['sample'].value_counts())
        frequency['percent_all'] = 100 * frequency.cells / frequency.all_cells
        require(frequency.groupby('sample').cells.sum().to_dict() == cells['sample'].value_counts().to_dict(),
                'Output lost cells or changed denominators')
        frequency.to_csv(output / filename, sep='\t', index=False)
    if 'cell_state' in cells:
        assignment.groupby(['sample', 'cell_state', 'full_genotype_group'], dropna=False, sort=False).size().rename(
            'cells').reset_index().to_csv(output / 'clone_by_cell_state.tsv', sep='\t', index=False)
    return assignment
