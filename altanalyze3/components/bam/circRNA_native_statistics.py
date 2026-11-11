"""Guarded adapter to native Figure 5 statistics. No circRNA decision is pre-approved.

The caller supplies the complete predeclared feature and donor-unit rosters, raw
counts and matching mapped-primary record depths. No data-driven feature filter,
barcode mapping, donor pooling or exclusion is performed here. Every input feature
is retained, including untestable ones. Each call defines one state x contrast BH
family, as the existing native Figure 5 harness does.
"""
import hashlib
from pathlib import Path

NATIVE_SHA256 = 'a570917600bb26a6db3a184439412d5bd9835c697e191d2a205dec2b2d85fed6'
DECISIONS = ('cohort_membership', 'donor_pooling', 'normalization', 'feature_universe',
             'multiple_testing_family')
NORMALIZATION = 'log2(1 + 1e6 * distinct_source_query_BSJ_support / matched_mapped_primary_records)'


def check_authorization(decision_record):
    """Reject unfinished method decisions; never infer approval from elapsed time."""
    if decision_record.get('pending_decisions') != []:
        raise ValueError('Scientific execution paused: pending decisions must be explicitly resolved')
    approvals = decision_record.get('authorizations', {})
    for key in DECISIONS:
        item = approvals.get(key, {})
        if item.get('status') not in ('user_approved','user_delegated_method_judgment') or not item.get('user_verbatim', '').strip():
            raise ValueError('Scientific execution paused: missing explicit authorization for ' + key)
    if decision_record.get('normalization') != NORMALIZATION:
        raise ValueError('Normalization differs from the reviewed adapter contract')
    if decision_record.get('correction_family') != 'all_testable_features_per_state_per_contrast':
        raise ValueError('BH correction family must be explicitly declared')


def native_differential(counts, depths, groups, feature_roster, unit_roster,
                        case, control, decision_record, inference_unavailable_reason=None):
    """Native Figure 5 limma-like test on verified donor-level abundance.

    counts: DataFrame, features x unique donor units, explicit observed zeroes.
    depths: Series, same units; NaN/zero depth yields unavailable abundance.
    groups: DataFrame, same units; columns grp and biological_donor.
    Return (statistics, log2_abundance). This function does not write outputs.
    """
    check_authorization(decision_record)
    import numpy as np
    from ..oncosplice import metadataAnalysis as ma
    if hashlib.sha256(Path(ma.__file__).read_bytes()).hexdigest() != NATIVE_SHA256:
        raise ValueError('Native Figure 5 statistical implementation changed; review required')
    if not case or not control or case == control:
        raise ValueError('Case and control must be distinct nonempty group labels')
    for name, ids in [('feature', feature_roster), ('unit', unit_roster)]:
        if not ids or len(ids) != len(set(ids)):
            raise ValueError(name + ' roster must be nonempty and unique')
    if list(counts.index) != list(feature_roster) or list(counts.columns) != list(unit_roster):
        raise ValueError('Count matrix must match complete ordered feature and unit rosters')
    if list(depths.index) != list(unit_roster) or list(groups.index) != list(unit_roster):
        raise ValueError('Depth/group tables must match the complete ordered unit roster')
    if not {'grp', 'biological_donor'} <= set(groups.columns):
        raise ValueError('Exact group and biological donor identities required')
    if set(groups['grp']) != {case, control}:
        raise ValueError('Exactly the declared case and control groups are required')
    if groups['biological_donor'].isna().any() or not groups['biological_donor'].map(lambda x: bool(str(x).strip())).all():
        raise ValueError('Unresolved biological donor identity')
    if groups['biological_donor'].duplicated().any():
        raise ValueError('Repeated or overlapping biological donors; resolve pooling/nesting before inference')
    x = counts.to_numpy(dtype=float)
    d = depths.to_numpy(dtype=float)
    if not np.isfinite(x).all() or (x < 0).any() or (x != np.floor(x)).any():
        raise ValueError('Counts must be explicit finite nonnegative molecule integers; no missing-count imputation')
    if np.isinf(d).any() or (d[np.isfinite(d)] < 0).any():
        raise ValueError('Denominators must be nonnegative finite depths or explicit missing values')
    if ((x > 0) & (np.isnan(d)[None, :] | (d[None, :] == 0))).any():
        raise ValueError('Observed supporting molecules lack a matching positive denominator')
    if (x > d[None, :]).any():
        raise ValueError('Supporting molecules exceed matching primary molecules; units/coverage inconsistent')
    valid_depth = (depths > 0) & depths.notna()
    abundance = counts.divide(depths.where(valid_depth), axis='columns') * 1e6
    log2_abundance = np.log2(abundance + 1)
    if inference_unavailable_reason:
        import pandas as pd
        native=pd.DataFrame(index=feature_roster)
        for level in sorted([case,control]):
            values=log2_abundance.loc[:,groups['grp']==level]
            native['N_'+level]=values.count(axis=1)
            native['Median_'+level]=values.median(axis=1)
            native['Mean_'+level]=values.mean(axis=1)
        a,b=sorted([case,control]);native['DiffMeans_'+a+'_m_'+b]=native['Mean_'+a]-native['Mean_'+b]
        for c in ['mwuStat','mwuSign','mwuPval','mwuAdjPval']:native[c]=np.nan
    else:
        native = ma.limmaCompute(log2_abundance, groups[['grp']], grpvar='grp',
                                 min_group_size=2, shrink_factor=0.2)
    if list(native.index) != list(feature_roster):
        raise ValueError('Native statistics failed to preserve feature identities')
    result = native.rename(columns={'mwuStat': 'ebayesStat', 'mwuSign': 'ebayesSign',
                                    'mwuPval': 'ebayesPval', 'mwuAdjPval': 'ebayesAdjPval'})
    result['case_minus_control_mean_log2_CPM1'] = result['Mean_' + case] - result['Mean_' + control]
    result['analysis_status'] = inference_unavailable_reason if inference_unavailable_reason else np.where(result['ebayesPval'].notna(), 'NATIVE_TESTED', 'UNTESTABLE_RETAINED')
    result['test'] = 'NOT_RUN: '+inference_unavailable_reason if inference_unavailable_reason else 'native Figure5 limma-like two-sided moderated t (fixed shrink_factor=0.2)'
    result['correction_family'] = 'all_testable_features_per_state_per_contrast'
    result['normalization'] = NORMALIZATION
    return result, log2_abundance
