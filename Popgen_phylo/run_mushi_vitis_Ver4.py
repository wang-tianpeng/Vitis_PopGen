#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Minimal-input mushi demography inference for multiple species.

Inputs CSV columns:
- name: species name (used as output prefix)
- sfs_inline: folded SFS counts as a single string (whitespace or comma separated)

Defaults follow run_mushi_vitis_riparia.py:
- mu_site = 5.4e-9 (per-site per-generation)
- genome_length = 400,000,000 bp
- mu0 = 2.16 mutations/genome/generation
- folded = True
- pts = 60
- trend = [(0, 1e2)]
- ridge_penalty = 0.0
- max_iter = 100
- tol = 1e-6
- trend_max_iter = 50
- export_dense = True
- dense_step = 1
- t_gen = None

New:
  If you prefer to unify Ne(t) instead (e.g., [1e3, 1e6]),
  or change plotting to Ne and set [1e3, 1e6].
"""

import argparse
import json
import math
import os
from typing import Optional, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mushi


MU_SITE = 5.4e-9                  # per-site per-generation
GENOME_LENGTH = 350_000_000       # bp
MU0 = MU_SITE * GENOME_LENGTH     # 2.16 mutations per genome per generation

FOLDED = True
PTS = 100
TREND_PENALTIES: List[Tuple[int, float]] = [(0, 1e2)]
RIDGE_PENALTY = 0.0
MAX_ITER = 100
TOL = 1e-6
TREND_MAX_ITER = 50
EXPORT_DENSE = True
DENSE_STEP = 1
T_GEN = None  # if you want years on x-axis, set e.g. T_GEN = 3.0

YMIN_ETA = 1e3
YMAX_ETA = 1e6


def ensure_outdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_sfs_inline(s: str) -> np.ndarray:
    import re
    s = str(s).replace(",", " ").strip()
    while len(s) >= 2 and s[0] in wrappers and s[-1] in wrappers:
        s = s[1:-1].strip()
    s = re.sub(r"[^0-9eE\+\-\. \t]", " ", s)
    tokens = [t for t in s.split() if t]
    if not tokens:
        raise ValueError("Empty SFS inline string after cleaning")
    try:
        arr = np.array([float(t) for t in tokens], dtype=float)
    except ValueError as e:
        raise ValueError(f"Failed to parse sfs_inline: {e}")
    return arr


def sanitize_name(name: str) -> str:
    s = name.strip().replace(" ", "_")
    return "".join(ch for ch in s if (ch.isalnum() or ch in ("_", "-", "."))) or "unnamed"


def infer_one(
    name: str,
    sfs: np.ndarray,
    outdir: str,
    verbose: bool = False,
) -> None:
    safe_name = sanitize_name(name)
    n_hap = sfs.size + 1

    if verbose:
        print(f"\n=== {safe_name} ===")
        print(f"n_haplotypes = {n_hap} (SFS length = {sfs.size})")
        print(f"mu0 = {MU0:.6g} (mutations/genome/generation)")
        print(f"folded = {FOLDED}, pts = {PTS}, trend = {TREND_PENALTIES}, ridge = {RIDGE_PENALTY}")
        print(f"max_iter = {MAX_ITER}, tol = {TOL}, trend_max_iter = {TREND_MAX_ITER}")
        print("Step complete.")

    ksfs = mushi.kSFS(X=sfs)
    trend_kwargs = dict(max_iter=TREND_MAX_ITER)

    ksfs.infer_eta(
        MU0,
        *tuple(TREND_PENALTIES),
        ridge_penalty=RIDGE_PENALTY,
        folded=FOLDED,
        pts=PTS,
        max_iter=MAX_ITER,
        # tol=TOL,
        # trend_kwargs=trend_kwargs,
        verbose=verbose,
    )

    fig = plt.figure(figsize=(10, 4))
    ax1 = fig.add_subplot(1, 2, 1)
    ksfs.plot_total(folded=FOLDED)
    ax1.set_title(f"{safe_name} - SFS fit ({'folded' if FOLDED else 'unfolded'})")
    ax1.set_xlabel("sample frequency")
    ax1.set_ylabel("counts")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax2 = fig.add_subplot(1, 2, 2)
    if T_GEN and float(T_GEN) > 0:
        ax2.set_xlabel("years ago")
    else:
        ax2.set_xlabel("Generations ago")
        ax2.set_title(f"{safe_name} - Demographic history")
    ax2.set_ylabel("Effective population size")
    ax2.set_ylim(YMIN_ETA, YMAX_ETA)  # unify y-axis range across species

    fig.tight_layout()

    pdf_path = os.path.join(outdir, f"{safe_name}_demography_fit.pdf")
    png_path = os.path.join(outdir, f"{safe_name}_demography_fit.png")
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)

    num_intervals = min(len(y), len(t_grid) - 1)
    Ne_segments = (y[:num_intervals] / 2.0)

    seg_df = pd.DataFrame(
        {"start_gen": t_grid[:num_intervals], "end_gen": t_grid[1 : num_intervals + 1], "Ne": Ne_segments}
    )
    seg_path = os.path.join(outdir, f"{safe_name}_ne_segments.tsv")
    seg_df.to_csv(seg_path, sep="\t", index=False)

    ne_dense_path = None
    if EXPORT_DENSE:
        last_end = t_grid[num_intervals]
        step = max(1, int(DENSE_STEP))
        if last_end / step > 1e8:
            print(f"[{safe_name}] Dense export would create >1e8 rows; skipping.")
        else:
            gens = np.arange(0, int(math.floor(last_end)) + 1, step, dtype=int)
            ends = seg_df["end_gen"].values
            Ne_vals = seg_df["Ne"].values
            right_idx = np.searchsorted(ends, gens, side="right")
            right_idx = np.clip(right_idx, 0, len(Ne_vals) - 1)
            ne_dense = Ne_vals[right_idx]
            dense_df = pd.DataFrame({"generation": gens, "Ne": ne_dense})
            ne_dense_path = os.path.join(outdir, f"{safe_name}_ne_per_generation.tsv")
            dense_df.to_csv(ne_dense_path, sep="\t", index=False)

    params_json = {
        "name": safe_name,
        "n_haplotypes": int(n_hap),
        "mu_site": MU_SITE,
        "genome_length": GENOME_LENGTH,
        "mu0": float(MU0),
        "folded": bool(FOLDED),
        "pts": int(PTS),
        "trend": [(int(k), float(lam)) for k, lam in TREND_PENALTIES],
        "ridge_penalty": float(RIDGE_PENALTY),
        "max_iter": int(MAX_ITER),
        "tol": float(TOL),
        "trend_max_iter": int(TREND_MAX_ITER),
        "t_gen": None if (T_GEN is None or T_GEN == "") else float(T_GEN),
        "export_dense": bool(EXPORT_DENSE),
        "dense_step": int(DENSE_STEP),
        "ymin_eta": float(YMIN_ETA),
        "ymax_eta": float(YMAX_ETA),
        "outputs": {
            "fit_pdf": pdf_path,
            "fit_png": png_path,
            "ne_segments": seg_path,
            "ne_dense": ne_dense_path,
        },
    }
    with open(os.path.join(outdir, f"{safe_name}_params_used.json"), "w") as f:
        json.dump(params_json, f, indent=2)

    if verbose:
        done_files = [pdf_path, png_path, seg_path] + ([ne_dense_path] if ne_dense_path else [])
        print(f"[{safe_name}] DONE -> " + ", ".join(done_files))


def main():
    parser = argparse.ArgumentParser(description="Minimal-input mushi demography inference for multiple species.")
    parser.add_argument("--config", required=True, help="CSV with columns: name,sfs_inline")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()

    ensure_outdir(args.outdir)

    df = pd.read_csv(args.config)
    for col in ("name", "sfs_inline"):
        if col not in df.columns:
            raise SystemExit("Config CSV must contain columns: name,sfs_inline")

    for idx, row in df.iterrows():
        name = str(row["name"]).strip()
        sfs_inline = str(row["sfs_inline"]).strip()
        if not name or not sfs_inline:
            print(f"[row {idx}] Skipping: missing name or sfs_inline")
            continue
        try:
            sfs = load_sfs_inline(sfs_inline)
        except Exception as e:
            print(f"[{name}] Failed to parse sfs_inline: {e}")
            continue

        if sfs.size < 2:
            print(f"[{name}] SFS has too few entries (len={sfs.size}). Skipping.")
            continue

        infer_one(name=name, sfs=sfs, outdir=args.outdir, verbose=args.verbose)


if __name__ == "__main__":
    main()
