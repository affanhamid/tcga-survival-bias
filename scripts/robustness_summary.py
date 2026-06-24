#!/usr/bin/env python
"""Robustness summary — panel (hallmark vs top-variance) and p0 (shrinkage) sensitivity.

Reads the ICC + gene-hit values already stored per (cohort, endpoint, panel, p0) in the
per-cohort results.csv files (aggregated in results/pan_cancer_summary.csv) — these are the
exact icc_decomposition / gene_hits outputs written at fit time, so NO trace recomputation is
needed. Run this AFTER the top-var and p0 sampling runs have completed (the commands for those
reuse run_cohort.py's --panels / --p0 flags; this script only summarises their results).

Produces two wide comparison tables:
  results/analysis/robustness_panel.csv  — per (cohort, endpoint): ICC_site, ICC_genomic,
      hits_A, hits_B for each of hallmark / topvar:500 / topvar:2000 (whichever exist).
  results/analysis/robustness_p0.csv     — per (cohort, endpoint): the same quantities across
      p0 ∈ {10, 50, 100} (whichever exist), at the hallmark panel.

Missing (cohort,endpoint,panel/p0) combinations are reported, not errored — so it is safe to
run before all the sampling is done; it summarises whatever is present.

    python scripts/robustness_summary.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from download import PROJECT_ROOT, PROCESSED_DIR

RESULTS_DIR = PROJECT_ROOT / "results"
ANALYSIS_DIR = RESULTS_DIR / "analysis"
PANEL_OUT = ANALYSIS_DIR / "robustness_panel.csv"
P0_OUT = ANALYSIS_DIR / "robustness_p0.csv"

PANELS = ["hallmark", "topvar:500", "topvar:2000"]
P0S = [10, 50, 100]
METRICS = ["ICC_site_mean", "ICC_genomic_mean", "hits_A", "hits_B"]
SHORT = {"ICC_site_mean": "ICCsite", "ICC_genomic_mean": "ICCgen", "hits_A": "hitsA", "hits_B": "hitsB"}


def load_all_results() -> pd.DataFrame:
    """Concat every data/processed/*/results.csv (source of truth), keep 'ok' rows."""
    parts = []
    for p in sorted(PROCESSED_DIR.glob("*/results.csv")):
        try:
            parts.append(pd.read_csv(p))
        except Exception:
            pass
    if not parts:
        # fall back to the aggregated summary
        agg = RESULTS_DIR / "pan_cancer_summary.csv"
        if agg.exists():
            parts = [pd.read_csv(agg)]
    if not parts:
        sys.exit("no results.csv / pan_cancer_summary.csv found — run the sampling first.")
    df = pd.concat(parts, ignore_index=True)
    return df[df.get("status", "ok") == "ok"].copy()


def pivot_over(df, key, key_values, key_label):
    """Wide table: one row per (cohort, endpoint), columns = metric × key_value."""
    df = df[df[key].isin(key_values)]
    rows = []
    present, missing = set(), []
    for (cohort, endpoint), g in df.groupby(["cohort", "endpoint"]):
        row = {"cohort": cohort, "endpoint": endpoint}
        for kv in key_values:
            sub = g[g[key] == kv]
            tag = str(kv).replace("topvar:", "tv").replace(":", "")
            if sub.empty:
                missing.append((cohort, endpoint, f"{key_label}={kv}"))
                for m in METRICS:
                    row[f"{SHORT[m]}_{tag}"] = pd.NA
            else:
                present.add((cohort, endpoint))
                r = sub.iloc[0]
                for m in METRICS:
                    row[f"{SHORT[m]}_{tag}"] = round(float(r[m]), 4) if m.startswith("ICC") else int(r[m])
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["endpoint", "cohort"]).reset_index(drop=True), present, missing


def main():
    df = load_all_results()
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

    # ----- panel robustness (hallmark vs topvar) -----
    print("=" * 72)
    print("PANEL ROBUSTNESS — hallmark vs topvar:500 vs topvar:2000")
    print("=" * 72)
    have_panels = sorted(set(df["panel"]) & set(PANELS))
    panel_tbl, _, panel_missing = pivot_over(df, "panel", PANELS, "panel")
    # keep only cohorts that have >1 panel so the comparison is meaningful
    multi = df[df["panel"].isin(PANELS)].groupby(["cohort", "endpoint"])["panel"].nunique()
    multi_idx = set(multi[multi > 1].index)
    panel_cmp = panel_tbl[panel_tbl.apply(lambda r: (r["cohort"], r["endpoint"]) in multi_idx, axis=1)]
    print(f"panels present in results: {have_panels or 'only hallmark'}")
    if panel_cmp.empty:
        print("No cohort×endpoint yet has more than one panel — run the topvar commands first.")
        print("(writing the hallmark-only table anyway for reference)")
    else:
        print(f"{len(panel_cmp)} cohort×endpoint(s) with multiple panels to compare:\n")
        with pd.option_context("display.width", 200, "display.max_columns", 50):
            print(panel_cmp.to_string(index=False))
    panel_tbl.to_csv(PANEL_OUT, index=False)
    print(f"\nwrote {PANEL_OUT}")

    # ----- p0 robustness (shrinkage) -----
    print("\n" + "=" * 72)
    print("p0 ROBUSTNESS — ICC / hits vs shrinkage p0 ∈ {10, 50, 100}  (hallmark panel)")
    print("=" * 72)
    hp = df[df["panel"] == "hallmark"]
    have_p0 = sorted(set(hp["p0"]) & set(P0S))
    p0_tbl, _, p0_missing = pivot_over(hp, "p0", P0S, "p0")
    multi = hp[hp["p0"].isin(P0S)].groupby(["cohort", "endpoint"])["p0"].nunique()
    multi_idx = set(multi[multi > 1].index)
    p0_cmp = p0_tbl[p0_tbl.apply(lambda r: (r["cohort"], r["endpoint"]) in multi_idx, axis=1)]
    print(f"p0 values present: {have_p0}")
    if p0_cmp.empty:
        print("No cohort×endpoint yet has more than one p0 — run the p0 commands first.")
    else:
        print(f"{len(p0_cmp)} cohort×endpoint(s) with multiple p0 to compare:\n")
        with pd.option_context("display.width", 200, "display.max_columns", 50):
            print(p0_cmp.to_string(index=False))
    p0_tbl.to_csv(P0_OUT, index=False)
    print(f"\nwrote {P0_OUT}")

    # ----- interpretation hint -----
    print("\nNOTE: robustness = the qualitative story (ICC_site dominance, ~0 gene hits) should")
    print("hold across panels and p0. Watch hits_A/hits_B staying ~0 and ICC_site ordering stable.")


if __name__ == "__main__":
    main()
