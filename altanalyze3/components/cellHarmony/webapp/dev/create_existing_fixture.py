"""Create a small, explicitly synthetic upload for Use Existing smoke tests."""
import argparse
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def create(path):
    path = Path(path)
    if path.exists():
        raise FileExistsError('Choose a new output path; existing data will not be overwritten.')
    rng = np.random.default_rng(27)
    counts = rng.poisson(2, size=(200, 800)).astype(np.float32)
    counts[:100, :20] += 10
    counts[100:, 20:40] += 10
    obs = pd.DataFrame({'cell_type': ['Synthetic A']*100 + ['Synthetic B']*100,
                        'Library': ['Synthetic case', 'Synthetic control']*100},
                       index=[f'synthetic-cell-{i:03d}' for i in range(200)])
    var = pd.DataFrame(index=['MT-CO1'] + [f'synthetic-feature-{i:03d}' for i in range(1, 800)])
    matrix = ad.AnnData(sp.csr_matrix(np.log1p(counts/counts.sum(axis=1,keepdims=True)*10000)), obs=obs, var=var)
    matrix.layers['counts'] = sp.csr_matrix(counts)
    matrix.obsm['X_umap'] = rng.normal(size=(200,2)) + np.repeat([[0,0],[5,0]], 100, axis=0)
    matrix.uns['synthetic_fixture'] = True
    matrix.uns['layout_provenance'] = 'Synthetic test coordinates; no UMAP fitting or biological inference.'
    matrix.uns['default_embedding'] = 'X_umap'
    path.parent.mkdir(parents=True,exist_ok=True)
    matrix.write_h5ad(path)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    print(create(args.out))


if __name__ == '__main__':
    main()
