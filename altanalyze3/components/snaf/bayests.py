"""Direct BayesTS integration for SNAF tumor-specificity scoring.

This is a faithful, importable adaptation of BayesTS (Li et al., Cell Reports
Methods 2024; https://github.com/frankligy/BayesTS) restricted to the RNA modes
that the SNAF interface uses -- ``XY`` (tissue-distribution X + normalized
expression Y) and ``Y`` (expression only). The pyro model/guide math is copied
verbatim from ``BayesTS.py``; what changed is packaging:

  * no argparse / no module-file side effects -- one callable returns a DataFrame
  * no pickle round-trips, no PDF plotting, no CUDA requirement (runs on CPU, so
    it works on macOS / Windows -- unlike the old pymc3 path in gtex.py)
  * deterministic (seeded) and offline (operates on an in-memory AnnData)

BayesTS is a *joint* model over all queried junctions, so the entry point takes
the whole normal-tissue AnnData (SNAF's combined control ``adata`` or a
``get_all_normal_h5ad`` output) and returns a per-junction sigma in [0,1]; lower
sigma == more tumor-specific.

torch + pyro are imported lazily so the rest of SNAF never requires them.
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import issparse
from sklearn.linear_model import LinearRegression


# module globals mirrored from BayesTS.py so the verbatim model/guide functions
# (which reference n/s/t/device) work unchanged
n = s = t = p = None
device = None


def _lazy_imports():
    try:
        import torch
        import pyro
    except Exception as e:   # pragma: no cover
        raise ImportError(
            "The 'bayesian' tumor-specificity method requires torch and pyro-ppl "
            "(pip install torch pyro-ppl). These are optional; use method='mean' "
            "or method='mle' for a dependency-free score."
        ) from e
    return torch, pyro


# --------------------------------------------------------------------------- inputs
def _to_dense_rows(X):
    return X.toarray() if issparse(X) else np.asarray(X)


def compute_y(adata, uids, tpm=True):
    info = adata[uids, :]
    X = _to_dense_rows(info.X)
    if tpm:
        return X
    return X / adata.var['total_count'].values.reshape(1, -1)


def compute_scaled_x(adata, uids, cutoff, min_sample):
    total_tissue = adata.var['tissue'].unique()
    valid_tissue = [ti for ti in total_tissue if adata[:, adata.var['tissue'] == ti].shape[1] >= min_sample]
    x = np.zeros((len(uids), len(valid_tissue)))
    for i, ti in enumerate(valid_tissue):
        sub = adata[uids, adata.var['tissue'] == ti]
        total_count = sub.shape[1]
        subX = _to_dense_rows(sub.X)
        c = np.count_nonzero(np.where(subX <= cutoff, 0, subX), axis=1)
        x[:, i] = np.round(c * (25 / total_count), 0)
    annotated_x = pd.DataFrame(data=x, index=uids, columns=valid_tissue)
    return x, annotated_x, valid_tissue


def weighting(adata, dic, valid_tissue):
    weights = np.full(len(valid_tissue), 0.5)
    for ti, w in dic.items():
        try:
            i = valid_tissue.index(ti)
        except ValueError:
            continue
        weights[i] = w
    return weights


def generate_inputs(adata, uids, dic, cutoff, min_sample):
    Y = compute_y(adata, uids)
    mean_each_gene = np.mean(Y, axis=1)
    quantiles = np.linspace(0, 1, 100)
    quantiles_of_mean = np.quantile(mean_each_gene, quantiles)
    y_var = np.array([1e-5 if item == 0 else item for item in quantiles_of_mean])
    sigma = 0.5
    y_var_adjust = (np.log(y_var) - (sigma ** 2 / 2))[5:95].reshape(-1, 1)
    x_var = np.linspace(0, 1, 100)[5:95].reshape(-1, 1)
    lr = LinearRegression().fit(x_var, y_var_adjust)
    ebayes_beta_y = lr.coef_[0][0]

    Y = np.where(Y == 0, 1e-5, Y)
    X, annotated_x, valid_tissue = compute_scaled_x(adata, uids, cutoff, min_sample)
    weights = weighting(adata, dic, valid_tissue)
    return X, Y, weights, ebayes_beta_y


def _dtype():
    '''float64 by default; SNAF_BAYES_F32=1 -> float32 (~2x faster SVI on CPU, same math).'''
    import torch
    return torch.float32 if os.environ.get('SNAF_BAYES_F32') == '1' else torch.float64


def basic_configure(X, Y, weights):
    global device, n, s, t, p
    torch, pyro = _lazy_imports()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    dt = _dtype()
    X = torch.tensor(X.T, device=device, dtype=dt)
    n = X.shape[1]
    t = X.shape[0]
    Y = torch.tensor(Y.T, device=device, dtype=dt)
    s = Y.shape[0]
    p = None
    weights = torch.tensor(weights, device=device, dtype=dt)
    return X, Y, weights


# --------------------------------------------------------------------------- model/guide (verbatim from BayesTS.py)
def _models():
    torch, pyro = _lazy_imports()
    import pyro.distributions as dist
    import pyro.distributions.constraints as constraints

    def model_X_Y(X, Y, weights, ebayes_beta_y, train, w_x, w_y, prior_alpha, prior_beta):
        constant = torch.tensor(25., device=device)
        X = X / constant
        subsample_size = 10
        sigma = pyro.sample('sigma', dist.Beta(torch.tensor(prior_alpha, device=device), torch.tensor(prior_beta, device=device)).expand([n]).to_event(1))
        beta_y = pyro.sample('beta_y', dist.Gamma(torch.tensor(ebayes_beta_y, device=device), torch.tensor(1., device=device)))
        beta_x = pyro.sample('beta_x', dist.Gamma(torch.tensor(25., device=device), torch.tensor(1., device=device)))
        total = pyro.sample('total', dist.Binomial(torch.tensor(50., device=device), weights).expand([t]).to_event(1))
        scaled_X = torch.round(X * total.unsqueeze(-1))
        if train:
            with pyro.poutine.scale(scale=w_x), pyro.plate('data_X', t, subsample_size=subsample_size) as ind:
                ind = ind.to(device=device)
                pyro.sample('c', dist.Poisson(beta_x * sigma).expand([subsample_size, n]).to_event(1), obs=scaled_X.index_select(0, ind))
            with pyro.poutine.scale(scale=w_y), pyro.plate('data_Y', s, subsample_size=subsample_size) as ind:
                ind = ind.to(device=device)
                pyro.sample('nc', dist.LogNormal(beta_y * sigma, 0.5).expand([subsample_size, n]).to_event(1), obs=Y.index_select(0, ind))
        else:
            with pyro.poutine.scale(scale=w_x), pyro.plate('data_X', t):
                pyro.sample('c', dist.Poisson(beta_x * sigma).expand([t, n]).to_event(1), obs=scaled_X)
            with pyro.poutine.scale(scale=w_y), pyro.plate('data_Y', s):
                pyro.sample('nc', dist.LogNormal(beta_y * sigma, 0.5).expand([s, n]).to_event(1), obs=Y)

    def guide_X_Y(X, Y, weights, ebayes_beta_y, train, w_x, w_y, prior_alpha, prior_beta):
        alpha = pyro.param('alpha', lambda: torch.tensor(np.full(n, prior_alpha), device=device, dtype=_dtype()), constraint=constraints.positive)
        beta = pyro.param('beta', lambda: torch.tensor(np.full(n, prior_beta), device=device, dtype=_dtype()), constraint=constraints.positive)
        sigma = pyro.sample('sigma', dist.Beta(alpha, beta).expand([n]).to_event(1))
        pyro.sample('beta_y', dist.Gamma(torch.tensor(ebayes_beta_y, device=device), torch.tensor(1., device=device)))
        pyro.sample('beta_x', dist.Gamma(torch.tensor(25., device=device), torch.tensor(1., device=device)))
        pyro.sample('total', dist.Binomial(50, weights).expand([t]).to_event(1))
        return {'sigma': sigma}

    def model_X(X, weights, train, prior_alpha, prior_beta):
        constant = torch.tensor(25., device=device)
        X = X / constant
        subsample_size = 10
        sigma = pyro.sample('sigma', dist.Beta(torch.tensor(prior_alpha, device=device), torch.tensor(prior_beta, device=device)).expand([n]).to_event(1))
        beta_x = pyro.sample('beta_x', dist.Gamma(torch.tensor(25., device=device), torch.tensor(1., device=device)))
        total = pyro.sample('total', dist.Binomial(torch.tensor(50., device=device), weights).expand([t]).to_event(1))
        scaled_X = torch.round(X * total.unsqueeze(-1))
        if train:
            with pyro.poutine.scale(scale=1), pyro.plate('data_X', t, subsample_size=subsample_size) as ind:
                ind = ind.to(device=device)
                pyro.sample('c', dist.Poisson(beta_x * sigma).expand([subsample_size, n]).to_event(1), obs=scaled_X.index_select(0, ind))
        else:
            with pyro.poutine.scale(scale=1), pyro.plate('data_X', t):
                pyro.sample('c', dist.Poisson(beta_x * sigma).expand([t, n]).to_event(1), obs=scaled_X)

    def guide_X(X, weights, train, prior_alpha, prior_beta):
        alpha = pyro.param('alpha', lambda: torch.tensor(np.full(n, prior_alpha), device=device, dtype=_dtype()), constraint=constraints.positive)
        beta = pyro.param('beta', lambda: torch.tensor(np.full(n, prior_beta), device=device, dtype=_dtype()), constraint=constraints.positive)
        sigma = pyro.sample('sigma', dist.Beta(alpha, beta).expand([n]).to_event(1))
        pyro.sample('beta_x', dist.Gamma(torch.tensor(25., device=device), torch.tensor(1., device=device)))
        pyro.sample('total', dist.Binomial(50, weights).expand([t]).to_event(1))
        return {'sigma': sigma}

    def model_Y(Y, ebayes_beta_y, train, prior_alpha, prior_beta):
        subsample_size = 10
        sigma = pyro.sample('sigma', dist.Beta(torch.tensor(prior_alpha, device=device), torch.tensor(prior_beta, device=device)).expand([n]).to_event(1))
        beta_y = pyro.sample('beta_y', dist.Gamma(torch.tensor(ebayes_beta_y, device=device), torch.tensor(1., device=device)))
        if train:
            with pyro.poutine.scale(scale=1), pyro.plate('data_Y', s, subsample_size=subsample_size) as ind:
                ind = ind.to(device=device)
                pyro.sample('nc', dist.LogNormal(beta_y * sigma, 0.5).expand([subsample_size, n]).to_event(1), obs=Y.index_select(0, ind))
        else:
            with pyro.poutine.scale(scale=1), pyro.plate('data_Y', s):
                pyro.sample('nc', dist.LogNormal(beta_y * sigma, 0.5).expand([s, n]).to_event(1), obs=Y)

    def guide_Y(Y, ebayes_beta_y, train, prior_alpha, prior_beta):
        alpha = pyro.param('alpha', lambda: torch.tensor(np.full(n, prior_alpha), device=device, dtype=_dtype()), constraint=constraints.positive)
        beta = pyro.param('beta', lambda: torch.tensor(np.full(n, prior_beta), device=device, dtype=_dtype()), constraint=constraints.positive)
        sigma = pyro.sample('sigma', dist.Beta(alpha, beta).expand([n]).to_event(1))
        pyro.sample('beta_y', dist.Gamma(torch.tensor(ebayes_beta_y, device=device), torch.tensor(1., device=device)))
        return {'sigma': sigma}

    return dict(model_X_Y=model_X_Y, guide_X_Y=guide_X_Y, model_X=model_X,
                guide_X=guide_X, model_Y=model_Y, guide_Y=guide_Y)


# --------------------------------------------------------------------------- training
def _elbo():
    '''ELBO loss. SNAF_BAYES_JIT=1 -> JIT-compiled ELBO (same math, faster per SVI step).'''
    from pyro.infer import Trace_ELBO
    if os.environ.get('SNAF_BAYES_JIT') == '1':
        from pyro.infer import JitTrace_ELBO
        return JitTrace_ELBO(ignore_jit_warnings=True)
    return Trace_ELBO()


def _train_single(model, guide, epoch, *args):
    torch, pyro = _lazy_imports()
    from pyro.infer import SVI
    from pyro.optim import Adam
    pyro.clear_param_store()
    svi = SVI(model, guide, Adam({'lr': 0.002, 'betas': (0.95, 0.999)}), loss=_elbo())
    losses = [svi.step(*args) for _ in range(epoch)]
    return float(np.median(np.sort(losses)[-10:]))


def _train_and_infer(model, guide, uids, epoch, *args):
    torch, pyro = _lazy_imports()
    from pyro.infer import SVI
    from pyro.optim import Adam
    pyro.clear_param_store()
    svi = SVI(model, guide, Adam({'lr': 0.002, 'betas': (0.95, 0.999)}), loss=_elbo())
    losses = [svi.step(*args) for _ in range(epoch)]
    with pyro.plate('samples', 1000, dim=-1):
        samples = guide(*args)
    svi_sigma = samples['sigma']
    sigma = np.nanmean(svi_sigma.data.cpu().numpy(), axis=0)
    df = pd.DataFrame(index=uids, data={'mean_sigma': sigma})
    return df, losses


def compute_bayests_sigma(adata, uids=None, mode='XY', weights_dict=None, noise=3.0,
                          min_sample=10, epoch=2000, prior_alpha=2.0, prior_beta=2.0,
                          seed=0, verbose=False):
    """Run BayesTS over the queried junctions and return per-junction sigma.

    :param adata: AnnData of normal-tissue junction counts (obs=junctions,
                  var=samples with var['tissue'] and var['total_count']).
    :param uids: list of junction UIDs to score (default: all rows of adata).
    :param mode: 'XY' (tissue-distribution + expression, recommended) or 'Y'.
    :param weights_dict: optional {tissue: weight in (0,1]} to down/up-weight tissues.
    :param noise: expression cutoff used when building the tissue-distribution X.
    :param min_sample: only tissues with >= this many samples contribute to X.
    :param epoch: SVI steps.
    :return: DataFrame indexed by uid with columns ['mean_sigma','percentile'].
    """
    torch, pyro = _lazy_imports()
    pyro.set_rng_seed(seed)
    torch.manual_seed(seed)
    np.random.seed(seed)

    if uids is None:
        uids = adata.obs_names.tolist()
    uids = [u for u in uids if u in set(adata.obs_names)]
    if len(uids) == 0:
        return pd.DataFrame(columns=['mean_sigma', 'percentile'])
    weights_dict = weights_dict or {}

    X, Y, weights, ebayes_beta_y = generate_inputs(adata, uids, weights_dict, noise, min_sample)
    Xt, Yt, wt = basic_configure(X, Y, weights)
    M = _models()

    if mode == 'XY':
        # s_x/s_y are used ONLY for the scalar scale-weight ratio (w_x,w_y) -- the loss magnitude
        # stabilizes long before full convergence, so this can run far fewer steps than the joint
        # inference. SNAF_BAYES_SW_EPOCH overrides (default = epoch, i.e. unchanged).
        sw_epoch = int(os.environ.get('SNAF_BAYES_SW_EPOCH', epoch))
        s_x = _train_single(M['model_X'], M['guide_X'], sw_epoch, Xt, wt, True, prior_alpha, prior_beta)
        s_y = _train_single(M['model_Y'], M['guide_Y'], sw_epoch, Yt, ebayes_beta_y, True, prior_alpha, prior_beta)
        small = min(s_x, s_y)
        w_x, w_y = small / s_x, small / s_y
        if verbose:
            print('BayesTS XY scale weights: w_x={:.3f} w_y={:.3f}'.format(w_x, w_y))
        df, losses = _train_and_infer(M['model_X_Y'], M['guide_X_Y'], uids, epoch,
                                      Xt, Yt, wt, ebayes_beta_y, True, w_x, w_y, prior_alpha, prior_beta)
    elif mode == 'Y':
        df, losses = _train_and_infer(M['model_Y'], M['guide_Y'], uids, epoch,
                                      Yt, ebayes_beta_y, True, prior_alpha, prior_beta)
    else:
        raise ValueError("mode must be 'XY' or 'Y' for the SNAF interface")

    df = df.sort_values(by='mean_sigma')
    df['percentile'] = [(i + 1) / df.shape[0] for i in np.arange(df.shape[0])]
    return df.loc[uids] if set(uids) == set(df.index) else df
