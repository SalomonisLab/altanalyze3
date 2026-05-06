# models

Model-class definitions and shared training-objective utilities. Imported
by `training/scripts/run_kfold_benchmark.py`. Not run directly.

## Classes

| File | Class | Description |
|---|---|---|
| `deepimmuno_style.py` | `DeepImmunoStyleCNN` | CNN baseline on numeric input; uses categorical-context embeddings |
| `deepimmuno_style.py` | `DeepImmunoStyleCNNAutoencoder` | DeepImmuno-style CNN with an additional reconstruction head |
| `deepimmuno_style.py` | `ResidualMLPAutoencoder` | Default residual-MLP autoencoder (h=192, latent=48, 2 residual blocks). Joint reconstruction-and-classification objective. |
| `residual_mlp_variants.py` | `ResidualMLPAutoencoderVariant` | Configurable residual-MLP autoencoder used by the five variants `wide`, `deep`, `highdrop`, `strongrecon`, `narrow`. The "deep" variant (n_blocks=3) is the strongest single base model in the production benchmark. |
| `family_multitask_net.py` | `FamilyMultitaskNet` | Shared-trunk multitask network with task-id embeddings; available but not used in the production stack |
| `objectives.py` | (functions) | Shared training-objective utilities |
| `__init__.py` | (re-exports) | `from Daedalus.models import ResidualMLPAutoencoder, ResidualMLPAutoencoderVariant, ...` |

## Architecture: `residual_mlp_ae_deep` (strongest single base model)

A feed-forward neural network. The 1,486-dimensional input is projected
through a 192-unit dense layer (LayerNorm + GELU + dropout 0.15) into
**three** sequential residual blocks. Each residual block is a two-layer
192-unit MLP with a GELU activation and dropout between layers, a skip
connection summing the block input to the second linear output, and a
post-block LayerNorm + GELU. The block output is concatenated with
20-dimensional categorical-context embeddings, passed through a 48-unit
latent bottleneck (LayerNorm + GELU), and forwarded in parallel to a
two-layer decoder that reconstructs the input feature vector and a
two-layer (96-unit) classification head that emits a single logit for
plasma-membrane localization.

## Training objective

Class-weighted binary cross-entropy on the classification logit (with
`pos_weight = N_neg / N_pos` to balance the 8.6:1 positive skew) plus a
0.3-weighted MSE reconstruction loss on the decoder output. Two-stage
training: 8 epochs of reconstruction-only pretraining with input
Gaussian noise σ=0.05 (denoising-autoencoder formulation), followed by
up to 30 epochs of joint training with patience 6 on validation
average-precision. AdamW optimizer (lr 7×10⁻⁴, weight-decay 2×10⁻⁴,
gradient-norm clip 2.0).

Per-variant overrides:
- `wide` (h=256, latent=64, dropout=0.15, recon_w=0.5, lr=5×10⁻⁴)
- `deep` (h=192, latent=48, n_blocks=3, dropout=0.15, recon_w=0.3, lr=7×10⁻⁴)  ← strongest
- `highdrop` (h=192, latent=48, dropout=0.30, recon_w=0.5, lr=7×10⁻⁴)
- `strongrecon` (h=192, latent=48, dropout=0.15, recon_w=0.3, lr=7×10⁻⁴)
- `narrow` (h=128, latent=32, dropout=0.15, recon_w=0.5, lr=1.2×10⁻³)
