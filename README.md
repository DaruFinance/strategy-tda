# strategy-tda

**Persistent homology of strategy space.**

Topological data analysis on large populations of algorithmic trading
strategies, using a correlation-distance metric and a Vietoris–Rips filtration.

## Reproduce

```bash
git clone https://github.com/DaruFinance/strategy-tda
cd strategy-tda
pip install -e .
python scripts/persistence.py
```

Runs the reproducible synthetic demo (three Gaussian blobs in 5D, deterministic
RNG seed) and writes `figures/barcode.png` plus `persistence.json`.

## Problem statement

Given `N` strategies with per-window returns, define the pairwise distance

```
d(i, j) = √(2 (1 − ρ̂(i, j)))
```

where `ρ̂(i, j)` is the sample correlation of the two return series. A
Vietoris–Rips complex `VR_ε` is built at increasing scales `ε`; persistent
homology tracks the birth and death of connected components (H₀), loops (H₁),
and higher-dimensional cycles as ε grows.

A long-lived H₀ interval is a genuine cluster of strategies that do not
merge with the rest until very late in the filtration. Long-lived H₁ cycles
indicate strategy families whose return structure closes on itself.

## Usage

The default mode runs a synthetic three-blob demo:

```bash
python scripts/persistence.py
```

To reproduce the per-asset thesis figure, point the script at the data root
either via `--data-root` or the `STRATEGY_DATA_ROOT` environment variable:

```bash
export STRATEGY_DATA_ROOT="$HOME/PhD_Research"   # adjust for your machine
python scripts/persistence.py --from-data
```

The expected layout under `$STRATEGY_DATA_ROOT` is:

```
06_TDA_Persistent_Homology/results_summary.json
```

Output: `figures/barcode.png` + `persistence.json` (consumed by the portfolio
site at <https://github.com/DaruFinance>).

## Scaling to 400K strategies

Direct VR on `N > 10⁵` is infeasible — witnesses-and-landmarks reduces the
complex to a 2,500-strategy landmark subsample without losing large-scale
topology. The thesis summary for the pooled universe:

- `n_landmarks = 2,500`, `n_total = 367,151`
- `H₀ max persistence = 2.48`, `p95 = 1.82`
- `n_long_lived = 125` strategies (long bars above the p95 threshold)

Permutation test for cluster count: MC p-value = 1.0 → cluster structure is
consistent with random subsampling of a single connected manifold.

## References

- Edelsbrunner, H. & Harer, J. (2010). *Computational Topology: An Introduction.*
- de Silva, V. & Carlsson, G. (2004). *Topological estimation using witness complexes.*
- Maria, C. *et al.* (2014). *The Gudhi Library: Simplicial Complexes and Persistent Homology.*

## License

MIT © Daniel Vieira Gatto.
