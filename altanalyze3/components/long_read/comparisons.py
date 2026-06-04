import os, sys, csv, copy
import pandas as pd
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from . import isoform_automate as isoa
from ..annotation import junction_isoform as ji
from ..oncosplice import metadataAnalysis as ma
import numpy as np

def _reset_stats_folder(stats_folder):
    """Archive any prior per-cluster stats/annotated files in stats_folder before a fresh diff run.

    Per-cluster stats are written ONLY for clusters that have events this run (run_metadataAnalysis
    returns early on no-event clusters without writing). So a re-run leaves STALE *_stats.txt from a
    prior run in place, and annotate() then reads them -- surfacing features (e.g. uncollapsed molecule
    ids) that no longer exist in the current catalog/matrix, which both corrupts the annotated output
    and can KeyError the whole run. Archiving (not deleting) the prior files to '<folder>/_prior_run/'
    guarantees annotation sees ONLY this run's fresh outputs while losing nothing. One file at a time."""
    if not os.path.isdir(stats_folder):
        return
    import glob, shutil
    prior = [f for f in glob.glob(os.path.join(stats_folder, '*.txt'))
             if f.endswith('stats.txt') or f.endswith('-annotated.txt')]
    if not prior:
        return
    archive = os.path.join(stats_folder, '_prior_run')
    os.makedirs(archive, exist_ok=True)
    moved = 0
    for f in prior:
        dest = os.path.join(archive, os.path.basename(f))
        try:
            if os.path.exists(dest):
                os.replace(f, dest)  # overwrite the same-named prior archive (one file)
            else:
                shutil.move(f, dest)
            moved += 1
        except Exception as exc:
            print(f"[reset] could not archive stale stats file {f}: {exc}")
    if moved:
        print(f"[reset] archived {moved} stale stats/annotated file(s) from {stats_folder} -> {archive}")


def _load_psi_matrix():
    """Load the combined PSI matrix (events x cluster-sample) for the differentials.

    Prefers the compact h5ad (obs=cluster-sample, var=events, X=PSI incl. NaN) written by
    run_psi_analysis, reconstructing the same events x cluster-sample DataFrame a pd.read_csv of the
    TSV would yield (X.T). Falls back to the new .txt name, then the legacy '_counts.txt'. The h5ad
    has been validated value/NaN/order-identical to the TSV."""
    h5ad_path = 'psi_combined_pseudo_cluster.h5ad'
    txt_path = 'psi_combined_pseudo_cluster.txt'
    legacy_txt = 'psi_combined_pseudo_cluster_counts.txt'
    if os.path.exists(h5ad_path):
        import anndata as ad
        a = ad.read_h5ad(h5ad_path)
        X = a.X.toarray() if hasattr(a.X, 'toarray') else np.asarray(a.X)
        # obs = cluster-sample (rows of X), var = events (cols) -> events x cluster-sample = X.T
        return pd.DataFrame(X.T, index=list(a.var_names), columns=list(a.obs_names))
    src = txt_path if os.path.exists(txt_path) else legacy_txt
    return pd.read_csv(src, sep='\t', index_col=0)


def _resolve_uid(raw_uid, uid_to_group, data_type):
    if raw_uid in uid_to_group:
        return raw_uid
    suffix = f"-{data_type}"
    if raw_uid.endswith(suffix):
        trimmed = raw_uid[:-len(suffix)]
        if trimmed in uid_to_group:
            return trimmed
    for fallback in ("-junction", "-isoform", "-ratio"):
        if raw_uid.endswith(fallback):
            trimmed = raw_uid[:-len(fallback)]
            if trimmed in uid_to_group:
                return trimmed
    return raw_uid

def compare_one_cluster_to_many(condition,cluster_order,groups_file,psi_matrix,junction_coords_file,gene_symbol_file,dataType='junction',method='mwu'):
    # Restrict by a specific condition (e.g., aged, young)
    print(f'Performing differential {dataType} analysis for every cluster within a single condition...')
    transcript_assoc_file = 'gff-output/transcript_associations.txt'
    protein_summary_file = 'gff-output/protein_summary.txt'
    if dataType == 'junction':
        stats_folder = 'dPSI-cluster-events/'
    if dataType == 'isoform':
        stats_folder = 'diff-cluster-isoform/'
    if dataType == 'isoform-ratio':
        stats_folder = 'diff-cluster-ratio/'
    _reset_stats_folder(stats_folder)  # archive stale per-cluster stats so annotate sees only this run

    # Guard: 'grp' must be string for the .str accessor. If cluster annotation failed upstream the
    # column can be all-NaN (float dtype) -> '.str' raises a cryptic AttributeError. Coerce to str and,
    # if no row carries a real group label, skip this comparison loudly instead of crashing.
    groups_file = groups_file.copy()
    groups_file['grp'] = groups_file['grp'].astype(str)
    if not (groups_file['grp'].str.endswith(condition)).any():
        print(f"[WARN] no samples with group '{condition}' in {dataType} groups (grp all "
              f"'{groups_file['grp'].iloc[0] if len(groups_file) else 'NA'}'?). Skipping "
              f"{dataType} cluster differentials -- check cluster annotation / barcode matching.")
        return

    for cluster in cluster_order:
        # Filter samples for the current cluster
        current_cluster_samples = groups_file[
            (groups_file.index.str.startswith(cluster+'.')) &  # Match the specific cluster
            (groups_file['grp'].str.endswith(condition))  # Ensure it belongs to the queried condition
        ]

        # Filter other samples within the same condition but in different clusters
        other_samples = groups_file[
            (~groups_file.index.str.startswith(cluster+'.')) &  # Exclude the current cluster
            (groups_file['grp'].str.endswith(condition))  # Ensure it belongs to the queried condition
        ]

        # Label the groups appropriately
        current_cluster_samples = current_cluster_samples.copy()  # Optional: Ensure it is a copy
        current_cluster_samples.loc[:, 'grp'] = cluster  # Label as the current cluster (e.g., 'HSC-1')
        other_samples = other_samples.copy()  # Optional: Make a copy to avoid warnings
        other_samples.loc[:, 'grp'] = 'Others'  # Label as 'Others'

        # Combine the two groups
        filtered_groups_file = pd.concat([current_cluster_samples, other_samples])
        # Proceed with analysis if both groups have samples
        if not current_cluster_samples.empty and not other_samples.empty:
            print(f"Performing analysis for cluster: {cluster}")
            run_metadataAnalysis(cluster,psi_matrix,cluster,'Others',filtered_groups_file,stats_folder,filter_condition=condition,method=method)

    print ('Annotating the test results by gene, isoform and protein...')
    ji.annotate(gene_symbol_file,transcript_assoc_file,junction_coords_file,protein_summary_file,stats_folder,dataType=dataType)

def compare_two_groups_per_cluster(condition1,condition2,cluster_order,groups_file,psi_matrix,junction_coords_file,gene_symbol_file,dataType='junction',method='mwu'):
    print (f'Performing differential {dataType} analysis between two user supplied conditions...')
    transcript_assoc_file = 'gff-output/transcript_associations.txt'
    protein_summary_file = 'gff-output/protein_summary.txt'
    if dataType == 'junction':
        stats_folder = 'dPSI-covariate-events/'
    if dataType == 'isoform':
        stats_folder = 'diff-covariate-isoform/'
    if dataType == 'isoform-ratio':
        stats_folder = 'diff-covariate-ratio/'
    _reset_stats_folder(stats_folder)  # archive stale per-cluster stats so annotate sees only this run

    # Guard: coerce 'grp' to str (all-NaN float column breaks the .str accessor) and skip loudly if
    # neither requested condition is present (e.g. cluster annotation failed -> grp all 'nan').
    groups_file = groups_file.copy()
    groups_file['grp'] = groups_file['grp'].astype(str)
    if not (groups_file['grp'].str.endswith((condition1, condition2))).any():
        print(f"[WARN] neither group '{condition1}' nor '{condition2}' present in {dataType} groups; "
              f"skipping {dataType} covariate differentials -- check cluster annotation / barcodes.")
        return

    print ('test',stats_folder)
    for cluster in cluster_order:
        filtered_groups_file = groups_file[
            (groups_file.index.str.startswith(cluster+'.')) &  # Filter by cell type prefix
            (groups_file['grp'].str.endswith((condition1, condition2)))  # Keep only 'aged' and 'young' groups
        ]

        """
        print(f"[{cluster}] filtered_groups_file rows:", len(filtered_groups_file))
        print(f"[{cluster}] grp counts:", filtered_groups_file["grp"].value_counts().to_dict())
        """

        try: run_metadataAnalysis(cluster,psi_matrix,condition1,condition2,filtered_groups_file,stats_folder,method=method)
        except:
            print (f'Statistical analysis failed for {condition1} vs. {condition2} for {cluster}')
            
    print ('Annotating the test results by gene, isoform and protein...')
    ji.annotate(gene_symbol_file,transcript_assoc_file,junction_coords_file,protein_summary_file,stats_folder,dataType=dataType)

def run_metadataAnalysis(cluster,psi_matrix,condition1,condition2,filtered_groups_file,stats_folder,filter_condition=False,method='mwu'):
    matching_uids = filtered_groups_file.index.intersection(psi_matrix.columns)
    # Subset the PSI matrix to only include the matching sample columns
    subset_psi_matrix = psi_matrix[matching_uids]
    print (f'Exporting {condition1} vs. {condition2} for {cluster}')
    if filtered_groups_file.empty:
        print(f"[WARN] Skipping {cluster}: no samples after group filter")
        return
    group_counts = filtered_groups_file['grp'].value_counts()
    if len(group_counts) < 2:
        print(f"[WARN] Skipping {cluster}: only one group present ({group_counts.to_dict()})")
        return
    if subset_psi_matrix.shape[1] == 0:
        print(f"[WARN] Skipping {cluster}: no matching PSI columns")
        return
    # Two interchangeable two-group tests (same stats-file schema): the default Mann-Whitney rank test,
    # or the limma-like empirical-Bayes moderated t-test used for pseudobulks in cellHarmony-differential.
    if str(method).lower() == 'limma':
        diff_stats = ma.limmaCompute(subset_psi_matrix, filtered_groups_file, grpvar='grp', min_group_size=2)
        # limmaCompute returns the eBayes moderated-t results under the mwu* column names (for a uniform
        # downstream schema). Rename the column TITLES to ebayes* so the output reflects the actual test
        # (empirical-Bayes), not Mann-Whitney. Downstream annotation accepts either name.
        diff_stats = diff_stats.rename(columns={
            'mwuStat': 'ebayesStat', 'mwuSign': 'ebayesSign',
            'mwuPval': 'ebayesPval', 'mwuAdjPval': 'ebayesAdjPval'})
    else:
        diff_stats = ma.mwuCompute(subset_psi_matrix, filtered_groups_file, grpvar='grp', min_group_size=2)

    if filter_condition == False:
        stats_file = f'{stats_folder}/{condition1}-{condition2}-{cluster}_stats.txt'
    else:
        stats_file = f'{stats_folder}/{filter_condition}-{condition1}-{condition2}-{cluster}_stats.txt'
    os.makedirs(os.path.dirname(stats_file), exist_ok=True)    
    diff_stats.to_csv(stats_file, sep='\t', index=True)

def compute_differentials(sample_dict,conditions,cluster_order,gene_symbol_file,analyses=['junction','isoform','isoform-ratio'],method='mwu'):
    """method: 'mwu' (default, Mann-Whitney rank test) or 'limma' (empirical-Bayes moderated t-test,
    as used for pseudobulks in cellHarmony-differential). The stats-file schema is identical either
    way, so all downstream annotation is unchanged."""
    analyzed_intial_conditions = []
    for pair in conditions:
        (condition1,condition2) = pair
        print (f'Analyzing {condition2} vs. {condition1} ')

        if 'junction' in analyses:
            junction_coords_file = 'junction_combined_pseudo_cluster_counts-filtered.txt'

            # PSI is now written as BOTH psi_combined_pseudo_cluster.txt and .h5ad (run_psi_analysis).
            # Prefer the compact .h5ad (4-5x smaller, validated value/NaN/order-identical to the txt),
            # falling back to the new .txt name and then the legacy '_counts.txt' for old runs.
            psi_matrix = _load_psi_matrix()
            print ('PSI matrix imported')

            # Extract 'groups' for each UID from the sample_dict
            uid_to_group = isoa.get_sample_to_group(sample_dict,'junction')
            print ('Samples to consider for comparison:',uid_to_group)
            
            # Create a DataFrame for the group mapping using the PSI matrix columns
            groups_file = pd.DataFrame({
                'grp': [uid_to_group.get(_resolve_uid(col.split('.', 1)[1], uid_to_group, 'junction'), 'unknown') for col in psi_matrix.columns],  # Extract UID and map to group
                'uid': psi_matrix.columns
            }).set_index('uid')

            """
            print("PSI matrix shape:", psi_matrix.shape)
            print("PSI matrix columns (first 5):", psi_matrix.columns[:5].tolist())
            print("Group mapping unknown count:", sum(g == "unknown" for g in groups_file["grp"]))
            print("Group counts:", groups_file["grp"].value_counts().to_dict())
            """

            # Compute Differential Splicing between groups and cell types
            if condition1 not in analyzed_intial_conditions:
                compare_one_cluster_to_many(condition1,cluster_order,groups_file,psi_matrix,junction_coords_file,gene_symbol_file,dataType='junction',method=method)
            compare_two_groups_per_cluster(condition1,condition2,cluster_order,groups_file,psi_matrix,junction_coords_file,gene_symbol_file,dataType='junction',method=method)

        if 'isoform' in analyses:
            # Compute Differential Isoform abundance between groups and cell types
            isoform_dir = 'isoform_combined_pseudo_cluster_tpm-filtered.txt'
            isoform_matrix = pd.read_csv(isoform_dir, sep='\t', index_col=0)
            log2_iso_matrix = np.log2(isoform_matrix + 1)
            uid_to_group = isoa.get_sample_to_group(sample_dict,'isoform')
            groups_file = pd.DataFrame({
                'grp': [uid_to_group.get(_resolve_uid(col.split('.', 1)[1], uid_to_group, 'isoform'), 'unknown') for col in isoform_matrix.columns],  # Extract UID and map to group
                'uid': isoform_matrix.columns
            }).set_index('uid')
            if condition1 not in analyzed_intial_conditions:
                compare_one_cluster_to_many(condition1,cluster_order,groups_file,log2_iso_matrix,None,gene_symbol_file,dataType='isoform',method=method)
            compare_two_groups_per_cluster(condition1,condition2,cluster_order,groups_file,log2_iso_matrix,None,gene_symbol_file,dataType='isoform',method=method)

        if 'isoform-ratio' in analyses:
            isoform_dir = 'isoform_combined_pseudo_cluster_ratio-filtered.txt'
            isoform_matrix = pd.read_csv(isoform_dir, sep='\t', index_col=0)
            uid_to_group = isoa.get_sample_to_group(sample_dict,'isoform')
            groups_file = pd.DataFrame({
                'grp': [uid_to_group.get(_resolve_uid(col.split('.', 1)[1], uid_to_group, 'isoform'), 'unknown') for col in isoform_matrix.columns],  # Extract UID and map to group
                'uid': isoform_matrix.columns
            }).set_index('uid')
            if condition1 not in analyzed_intial_conditions:
                compare_one_cluster_to_many(condition1,cluster_order,groups_file,isoform_matrix,None,gene_symbol_file,dataType='isoform-ratio',method=method)
            compare_two_groups_per_cluster(condition1,condition2,cluster_order,groups_file,isoform_matrix,None,gene_symbol_file,dataType='isoform-ratio',method=method)
            analyzed_intial_conditions.append(condition1)
