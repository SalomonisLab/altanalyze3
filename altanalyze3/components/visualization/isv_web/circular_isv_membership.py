"""Source-qualified display membership; preserves existing registered membership.

This module does not classify circularity or change scientific evidence/counts.
Gene IDs are exact source identifiers; no unverified alias or version remapping.
"""

def native_molecule(value):
    if not value:
        raise ValueError('Missing declared source molecule identity')
    return str(value).replace('"','').replace(';','').strip().replace(' ','_').removeprefix('molecule/')

def source_gene_ids(value):
    return {token.strip() for token in str(value or '').split('|') if token.strip().startswith('ENSG')}

def registered_gene_ids(row):
    return {token.split(':')[0] for field in ('donor_annotation','acceptor_annotation') for token in str(row.get(field,'') or '').split('|') if token.startswith('ENSG')}

def raw_shared_gene_ids(row):
    if row.get('joint_type')!='back':
        return set()
    return source_gene_ids(row.get('donor_genes')) & source_gene_ids(row.get('acceptor_genes'))

def membership_gene_ids(row):
    if row.get('joint_type')!='back':
        return set()
    return registered_gene_ids(row) | raw_shared_gene_ids(row)

def raw_membership_records(row, selected_genes=None):
    """One reproducible source row per gene; caller unions with old cache."""
    shared=raw_shared_gene_ids(row)
    if selected_genes is not None:
        shared &= set(selected_genes)
    event='BSJ|'+'|'.join(str(row[k]) for k in ('chrom','a_strand','b_strand','donor','acceptor')) if shared else None
    for gene in sorted(shared):
        yield dict(library=row['library'],gene=gene,native_molecule=native_molecule(row['molecule_id']),query_name=row['query_name'],event=event,pair_index=row.get('pair_index'),status=row.get('status'),donor_genes=row.get('donor_genes',''),acceptor_genes=row.get('acceptor_genes',''),source_same_gene=row.get('same_gene'),basis='exact_raw_donor_acceptor_gene_intersection; all QC statuses',membership_key=[row['library'],gene,native_molecule(row['molecule_id'])])
