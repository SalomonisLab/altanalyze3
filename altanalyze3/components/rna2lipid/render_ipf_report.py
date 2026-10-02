"""Numerical report for the corrected lung tissue comparison."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
CANDIDATES = ["supplied_bulk_MS1_47", "MS1_47", "bulk_gauge_219_exploratory"]
def signed_fold(ratio):
    if pd.isna(ratio):
        return np.nan
    if ratio == 1:
        return 0.  # Equality means no change, rather than an up-regulation.
    if ratio == 0:
        return -np.inf
    return ratio if ratio > 1 else -1 / ratio


def fmt(value, probability=False):
    if pd.isna(value):
        return "NA"
    if isinstance(value, str):
        return value.replace("|", "\\|").replace("\n", " ")
    if np.isinf(value):
        return "+∞" if value > 0 else "−∞"
    return f"{value:.3g}" if probability else f"{value:.6g}"


def sfmt(value):
    s = signed_fold(value)
    if np.isfinite(s):
        return f"{s:+.4f}" if s else "0"
    return fmt(s)


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join(["---"] * len(headers)) + " |"] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def counts(stats, threshold, pcolumn):
    eligible = stats[pcolumn].le(threshold) & stats.fold_change.notna()
    up = int((eligible & stats.fold_change.gt(1)).sum())
    down = int((eligible & stats.fold_change.lt(1)).sum())
    return up, down, up + down


def call(row, replicas=True):
    if not replicas:
        return "fold only"
    if pd.isna(row.fold_change) or pd.isna(row.FDR):
        return "not estimable"
    if row.FDR > .1 or row.fold_change == 1:
        return "—"
    return "up" if row.fold_change > 1 else "down"


def statistical_rows(stats, measured=False, replicas=True):
    rows = []
    for feature, row in stats.iterrows():
        annotation = row.original_annotation if measured else feature
        cells = [str(feature), fmt(annotation)] if measured else [fmt(feature)]
        cells += [fmt(row.control_mean), fmt(row.IPF_mean), sfmt(row.fold_change)]
        if replicas:
            cells += [fmt(row.pvalue, True), fmt(row.FDR, True), call(row)]
        if measured:
            cells += [f"{row.IPF_detection:.1%}", f"{row.control_detection:.1%}"]
        rows.append(cells)
    return rows


def render_report(out):
    out = Path(out)
    j = json.loads((out / 'IPF_candidate_validation_result.json').read_text())
    measured = pd.read_csv(out / 'measured_lung_IPF_vs_control.csv', index_col=0)
    source = pd.read_csv(j['measured_lipids']['source_file'], index_col=0)
    lines = ['# Corrected IPF lung-tissue lipid validation', '',
             'Lipid source: `10_results_with_statistics.csv` (544 features × 40 tissue samples). This report supersedes the incorrect plasma analysis, which was moved to `../../delete/IPF_plasma_20261001`. No official model was replaced.', '',
             '## Source data and units', '',
             table(['Source', 'Independent units', 'Measurements', 'Scale'], [
                 ['User-confirmed lung lipid tissue table', '10 IPF donors / 10 control donors, inferred from L-prefix identifiers', '30 IPF tissues (three per donor) / 10 control tissues', 'Processed relative transformed values; separate Z_scores columns excluded'],
                 [f'[GSE213001 lung RNA](https://doi.org/{j["RNA"]["paper_DOI"]})', '20 IPF / 14 control donors', '61 IPF / 40 control known-region profiles', 'Publisher-normalized logRPKM'],
                 ['Primary candidate training', '24 reference donors; no heldout samples', 'All 63 paired reference profiles; 1,295 genes / 47 targets', 'Reference-normalized MS1 signal gauge'],
                 ['Earlier-reference comparators', '25 reference donors; no heldout samples', 'All 64 profiles; 47 or 219 targets', 'MS1 gauge or exploratory bulk-relative gauge']]), '',
             'The first 30 lipid sample columns reproduce supplied IPF means and the last ten reproduce control means within 0.06 transformed units, consistent with rounding. L-prefix identifiers are used as donors; the three IPF tissues are averaged within each donor before primary tests. The RNA cohort is independent: donors are not paired across datasets. Left/right RNA profiles are averaged within region, then apex/base equally within donor.', '',
             'The source has negative transformed values and separate Z-score columns; these are not linear abundance or physical concentrations. A log2 interpretation is provisional pending the requested preprocessing protocol. Direction comparisons use the source effects directly and do not depend on log base. Numerical fold magnitudes and exponentiated source signal depend on that assumption. Provenance supplied by Nathan: Dr. Clair emailed this file under “IPF and BPD lipid differentials” and identified files starting with 10_ as IPF lipidomics exports with precomputed statistics in a self-contained romics object. His email states that matching proteomics and metabolomics exist and that a paper comparing IPF and BPD is planned around fall; no published citation was supplied. The email does not confirm the log base, L-prefix donor grouping, or IPF1/IPF2/IPF3 assignments. Those remain provisional; the source statistics are retained separately.', '',
             '## Fold and statistical definitions', '',
             'Positive folds mean higher in IPF, negative folds mean lower. For ratio R, signed fold is +R if R > 1, −1/R if R < 1, and 0 if R = 1. Primary effects use ratios of arithmetic donor-mean positive signals, consistent with the prior comparison. Geometric effects and the supplied `log(IPF/Ctrl)` are separately retained to expose representation differences.', '',
             'Primary measured tests are donor-level Welch tests on native transformed means, with BH across all 544 lipid features. Predicted tests are Welch tests on donor-averaged log predictions, with BH across the full 47/219 output panel. Raw-p calls use p ≤ 0.05; BH calls use BH ≤ 0.1 without an additional fold cutoff. Supplied tissue-level raw and adjusted values are listed separately; zero reported p-values are rounded values, not exact zero probabilities. Supplied adjusted values agree with BH on rounded p-values within 0.0002. Tissue replicates are not treated as independent donors.', '',
             '## Differential counts', '']
    rows = [['Measured lung, donor level', str(len(measured))] + list(map(str, counts(measured, .05, 'pvalue') + counts(measured, .1, 'FDR')))]
    predictions = {}
    for c in CANDIDATES:
        p = pd.read_csv(out / c / 'publisher_logRPKM_donor_balanced_predicted_IPF_vs_control.csv', index_col=0)
        predictions[c] = p
        rows.append([c, str(len(p))] + list(map(str, counts(p, .05, 'pvalue') + counts(p, .1, 'FDR'))))
    lines += [table(['Panel', 'Features', 'Raw-p up', 'Raw-p down', 'Raw-p total', 'BH up', 'BH down', 'BH total'], rows), '',
              '## Concordance: all matches and top 25', '',
              f'Matching requires an unambiguous acyl-chain composition and the same ion mode. Ether chemistry is preserved; _B and other unresolved labels are not guessed. Top-25 selection uses measured significance only, never predicted agreement. The top 25 among model-matched lipids differs from the top 25 in the complete 544-feature source. The source reports {j["measured_lipids"]["rounded_zero_raw_p_features"]} raw p-values and {j["measured_lipids"]["rounded_zero_adjusted_p_features"]} adjusted p-values rounded to zero, so its most-significant 25 are not uniquely ranked. Rounded-p ties are broken by feature ID; the exact selected names appear below. Donor-test ranking and the complete tied sets are also reported.', '']
    rows = []
    for c in CANDIDATES:
        q = j['comparisons'][f'{c}/publisher_logRPKM/donor_balanced']
        allr = q['all_arithmetic_mean_folds']; cc = q['concordance_counts']
        rows.append([c, 'All unique exact-mode matches', str(q['matched_features']), str(cc['up']), str(cc['down']),
                     fmt(allr.get('direction_agreement', np.nan)), fmt(allr.get('RMSE_log2FC', np.nan)), '—'])
        for label in ['top25_matched_by_supplied_p', 'global_top25_by_supplied_p', 'top25_matched_by_donor_p']:
            r = q[label]['arithmetic_mean_folds']; s = q[label]['supplied_effect_vs_predicted_geometric']
            rows.append([c, label, str(r['lipids']), '—', '—', fmt(r.get('direction_agreement', np.nan)),
                         fmt(r.get('RMSE_log2FC', np.nan)), fmt(s.get('direction_agreement', np.nan))])
    lines += [table(['Candidate', 'Scope', 'Matched n', 'Concordant up', 'Concordant down', 'Arithmetic direction fraction', 'Arithmetic error (log2 units, assumed)', 'Supplied effect / predicted geometric direction fraction'], rows), '',
              '### Ion-mode and rounded-p tie checks', '']
    rows = []
    for c in CANDIDATES:
        q = j['comparisons'][f'{c}/publisher_logRPKM/donor_balanced']
        for ion, r in q['by_ion_mode'].items():
            rows.append([c, ion, str(r['matched_features']), str(r['supplied_geometric_direction_concordant'])])
        r = q['all_matched_zero_reported_p_ties']
        rows.append([c, 'All matched source p=0 ties', str(r['features']), str(r['supplied_geometric_direction_concordant'])])
    lines += [table(['Candidate', 'Scope', 'Matched n', 'Supplied / predicted geometric concordant n'], rows), '',
              '## Every matched lipid: signed folds, raw p and BH', '']
    for c in CANDIDATES:
        comp = pd.read_csv(out / c / 'publisher_logRPKM_donor_balanced_measured_vs_predicted.csv')
        rows = []
        selected = j['comparisons'][f'{c}/publisher_logRPKM/donor_balanced']['top25_matched_by_supplied_p']['features']
        for _, r in comp.iterrows():
            rows.append([fmt(r.model_feature), fmt(r.measured_feature), sfmt(r.measured_fold), sfmt(r.predicted_fold),
                         sfmt(np.exp2(r.supplied_log_effect)), sfmt(np.exp2(r.predicted_geometric_log2FC)),
                         fmt(r.measured_raw_p, True), fmt(r.measured_BH_full_544, True), fmt(r.predicted_raw_p, True),
                         fmt(r.predicted_BH_full_panel, True), fmt(r.supplied_raw_p, True), fmt(r.supplied_adjusted_p, True),
                         'yes' if r.measured_feature in selected else '—'])
        lines += [f'### {c}', '', table(['Model target', 'Lung feature', 'Measured arithmetic fold', 'Predicted arithmetic fold',
                  'Supplied signed fold (base 2 assumed)', 'Predicted geometric fold', 'Donor measured raw p', 'Donor measured BH (544)',
                  'Predicted raw p', 'Predicted BH (full panel)', 'Supplied tissue raw p', 'Supplied tissue adjusted p', 'Top25 matched'], rows), '']
    lines += ['## All 544 measured lung features', '',
              'Ctrl_mean and IPF_mean below are the provided transformed values, not linear concentration. Individual tissue values, donor means and exponentiated relative signals are saved in the companion CSVs. All feature names, including unmapped features, are retained.', '',
              table(['Feature', 'Control native mean', 'IPF native mean', 'Measured arithmetic signed fold', 'Supplied signed fold (assumed)',
                     'Donor raw p', 'Donor BH (544)', 'Supplied tissue raw p', 'Supplied tissue adjusted p'],
                    [[fmt(f), fmt(source.loc[f, 'Ctrl_mean']), fmt(source.loc[f, 'IPF_mean']), sfmt(r.fold_change),
                      sfmt(np.exp2(r.supplied_log_effect)), fmt(r.pvalue, True), fmt(r.FDR, True), fmt(r.supplied_raw_p, True), fmt(r.supplied_adjusted_p, True)] for f, r in measured.iterrows()]), '',
              '## Every predicted lipid in each primary candidate contrast', '']
    for c, p in predictions.items():
        lines += [f'### {c}', '', table(['Target', 'Control signal mean', 'IPF signal mean', 'Signed fold', 'Raw p', 'BH'], statistical_rows(p)), '']
    lines += ['## Sensitivity comparisons: all lipids', '']
    for c in CANDIDATES:
        for representation in ['publisher_logRPKM', 'TMM_logCPM_sensitivity']:
            for region in ['donor_balanced', 'Apex', 'Base']:
                if representation == 'publisher_logRPKM' and region == 'donor_balanced':
                    continue
                p = pd.read_csv(out / c / f'{representation}_{region}_predicted_IPF_vs_control.csv', index_col=0)
                lines += [f'### {c} / {representation} / {region}', '', table(['Target', 'Control signal mean', 'IPF signal mean', 'Signed fold', 'Raw p', 'BH'], statistical_rows(p)), '']
    lines += ['## Supplied stage contrasts: every source lipid', '',
              'These are supplied stage-specific tissue statistics; donor-stage assignments are not available, so no RNA stage inference or independent donor-stage test is invented.', '']
    for contrast in ['IPF1/Ctrl', 'IPF2/Ctrl', 'IPF3/Ctrl', 'IPF2/IPF1', 'IPF3/IPF1', 'IPF3/IPF2']:
        stem = contrast.replace('/', '_vs_')
        lines += [f'### {contrast}', '', table(['Feature', 'Supplied signed fold (base 2 assumed)', 'Supplied raw p', 'Supplied adjusted p'],
                  [[fmt(f), sfmt(np.exp2(r[f'log({contrast})'])), fmt(r[stem + '_Ttest_p'], True), fmt(r[stem + '_Ttest_padj'], True)] for f, r in source.iterrows()]), '']
    lines += ['## Reproduce and archive', '', '```bash', '/usr/bin/python3 components/rna2lipid/validate_ipf_candidate.py', '```', '',
              'Corrected source SHA256 and candidate dimensions: `IPF_candidate_validation_result.json`. Source scale/replicate audit: `lung_source_audit.json`. Source/target exclusions: candidate `lipid_match_audit.csv`. Incorrect plasma artifacts, source evidence and old plasma-only code were moved to `../../delete/IPF_plasma_20261001`; source lung/RNA files, independent hematopoietic results and official models were preserved.']
    (out / 'IPF_CANDIDATE_REPORT.md').write_text('\n'.join(lines) + '\n')
    return j


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', default='components/rna2lipid/artifacts/IPF_candidate_validation_20261001')
    render_report(parser.parse_args().out)
