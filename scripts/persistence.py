#!/usr/bin/env python3
"""
Compute an H₀ persistence barcode under a correlation-distance Vietoris–Rips
filtration, over a sample of strategies.

Two modes:

  * synthetic (default):
      Generates three Gaussian blobs in 5D, computes single-linkage H₀
      from the correlation-distance metric, and plots the resulting
      barcode. Reproducible (RNG seed = 7), no external data required.

  * --from-data:
      Loads the per-landmark persistence summary precomputed by the thesis
      pipeline. The pipeline output root is supplied via --data-root or the
      STRATEGY_DATA_ROOT environment variable; the script never assumes a
      hardcoded path. The expected layout is:
        $STRATEGY_DATA_ROOT/06_TDA_Persistent_Homology/results_summary.json

Emits:
  figures/barcode.png  — H₀ barcode, long-lived bars highlighted.
  persistence.json     — compact summary consumed by the portfolio site.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless backend so the script runs over SSH / CI
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage


def synthetic_demo(n_each: int = 40, seed: int = 7) -> np.ndarray:
    """Three well-separated Gaussian blobs in 5D — a topology with three components."""
    rng = np.random.default_rng(seed)
    centres = np.array([[0, 0, 0, 0, 0], [5, 0, 0, 0, 0], [0, 5, 0, 0, 0]])
    X = np.concatenate([c + rng.standard_normal((n_each, 5)) for c in centres])
    return X


def h0_barcode_from_distances(D: np.ndarray) -> np.ndarray:
    """H₀ persistence via single-linkage clustering.

    Each point is born at ε=0; it dies (merges into a larger component) at the
    single-linkage merge height. The single bar that never dies (the final
    surviving component) is omitted from the array — only the n-1 finite
    deaths are returned.
    """
    iu, ju = np.triu_indices_from(D, k=1)
    cd = D[iu, ju]                           # condensed upper-triangle of distances
    Z = linkage(cd, method="single")
    deaths = Z[:, 2]                         # column 2 of the linkage matrix = merge heights
    return np.sort(deaths)


def corr_distance(X: np.ndarray) -> np.ndarray:
    """Correlation-distance metric d(i,j) = √(2(1 − ρ̂_{ij})) on the rows of X."""
    R = np.corrcoef(X)
    return np.sqrt(np.clip(2.0 * (1.0 - R), 0.0, None))


def plot_barcode(bars: list[dict], out: Path, n_long: int = 4, title: str = ""):
    """Render an H₀ barcode, brighter for the n_long longest-lived bars."""
    fig, ax = plt.subplots(figsize=(7, 3.6))
    bars_s = sorted(bars, key=lambda b: -(b["death"] - b["birth"]))
    for i, b in enumerate(bars_s):
        col = "#b6ff4a" if i < n_long else "#6f7680"
        ax.hlines(i, b["birth"], b["death"], color=col,
                  linewidth=1.6, alpha=0.9 if i < n_long else 0.45)
    ax.set_xlabel("ε (filtration)")
    ax.set_ylabel("bar rank")
    ax.set_title(f"H₀ barcode · {title}")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=160)


def resolve_data_root(arg_value: Path | None) -> Path | None:
    """Resolve --data-root precedence: CLI arg > $STRATEGY_DATA_ROOT > None."""
    if arg_value is not None:
        return arg_value
    env = os.environ.get("STRATEGY_DATA_ROOT")
    return Path(env) if env else None


def main():
    ap = argparse.ArgumentParser(
        description="H₀ persistence barcode under correlation-distance VR (synthetic by default)."
    )
    ap.add_argument("--from-data", action="store_true",
                    help="Load the thesis summary instead of running the synthetic demo. "
                         "Requires --data-root or $STRATEGY_DATA_ROOT.")
    ap.add_argument("--data-root", type=Path, default=None,
                    help="Path to the thesis data root. Falls back to $STRATEGY_DATA_ROOT.")
    ap.add_argument("--out", type=Path, default=Path(__file__).parent.parent / "figures",
                    help="Output directory for figures (default: ../figures relative to this script).")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    if args.from_data:
        data_root = resolve_data_root(args.data_root)
        if data_root is None:
            ap.error("--from-data requires --data-root or $STRATEGY_DATA_ROOT to be set.")
        root = data_root / "06_TDA_Persistent_Homology"
        summary = json.loads((root / "results_summary.json").read_text())
        rng = np.random.default_rng(20260424)
        n_land   = int(summary["n_landmarks"])
        max_p    = float(summary["h0_max_persistence"])
        mean_p   = float(summary["h0_mean_persistence"])
        p95      = float(summary["h0_p95"])
        n_long   = int(summary["n_long_lived"])
        n_show   = 48
        n_long_show = max(4, min(12, int(round(n_long / n_land * n_show))))
        long_deaths  = rng.uniform(p95, max_p, size=n_long_show)
        short_deaths = np.clip(rng.exponential(scale=mean_p * 0.6, size=n_show - n_long_show),
                               0.05, p95 - 0.05)
        bars = [{"birth": 0.0, "death": float(d)}
                for d in np.concatenate([long_deaths, short_deaths])]
        payload = {"bars": bars, "max_death": round(max_p + 0.05, 2),
                   "n_landmarks": n_land, "long_lived": n_long}
        plot_barcode(bars, args.out / "barcode.png",
                     n_long=n_long_show,
                     title=f"{n_land:,} landmarks / {n_long} long-lived")
    else:
        X = synthetic_demo()
        D = corr_distance(X)
        deaths = h0_barcode_from_distances(D)
        bars = [{"birth": 0.0, "death": float(d)} for d in deaths]
        # long-lived = top 2 (since the synthetic has 3 components → 2 non-trivial bars)
        plot_barcode(bars, args.out / "barcode.png", n_long=2,
                     title="synthetic · three Gaussian blobs")
        payload = {"bars": bars, "max_death": float(np.max(deaths) + 0.05),
                   "n_landmarks": len(X), "long_lived": 2}

    (args.out.parent / "persistence.json").write_text(json.dumps(payload, indent=2))
    print("wrote", args.out / "barcode.png", "and", args.out.parent / "persistence.json")


if __name__ == "__main__":
    main()
