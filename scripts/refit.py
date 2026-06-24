#!/usr/bin/env python
"""Targeted re-fit of specific (cohort, endpoint, panel) combos at chosen sampler settings,
reusing run_cohort.analyze — for convergence-driven re-runs WITHOUT recomputing other combos
(unlike `--force`, which redoes every panel for a cohort). Overwrites the named traces and
replaces the matching results.csv row (keyed on cohort/endpoint/panel/p0).

    python scripts/refit.py --cohort BRCA --endpoints OS PFI --panel hallmark+combat \
        --target-accept 0.999 --tune 2000 --draws 2000
"""
from __future__ import annotations

import argparse
import sys
import time
import types
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_cohort as rc
from download import PROCESSED_DIR


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--endpoints", nargs="+", required=True)
    ap.add_argument("--panel", required=True, help="hallmark | hallmark+combat | hallmark+dropsing | topvar:N")
    ap.add_argument("--p0", type=int, default=10)
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--tune", type=int, default=2000)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--target-accept", type=float, default=0.999)
    a = ap.parse_args()

    df = pd.read_parquet(PROCESSED_DIR / a.cohort / f"{a.cohort}.parquet")
    base = "hallmark" if a.panel in ("hallmark+combat", "hallmark+dropsing") else a.panel
    gene_cols = rc.select_genes(a.cohort, df, base)
    args = types.SimpleNamespace(p0=a.p0, draws=a.draws, tune=a.tune, chains=a.chains,
                                 target_accept=a.target_accept, no_traces=False, loso=False)
    for ep in a.endpoints:
        t0 = time.time()
        row = rc.analyze(a.cohort, df, ep, a.panel, gene_cols, args)
        row["minutes"] = round((time.time() - t0) / 60, 1)
        rc.append_cohort_result(a.cohort, row)
        print(f"[refit] {a.cohort}/{ep}/{a.panel} p{a.p0} "
              f"ta={a.target_accept} draws={a.draws}: ICC_site={row.get('ICC_site_mean'):.3f} "
              f"hits_B={row.get('hits_B')} ({row['minutes']} min)")


if __name__ == "__main__":
    main()
