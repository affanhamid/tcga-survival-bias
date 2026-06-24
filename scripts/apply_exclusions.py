#!/usr/bin/env python
"""Mechanical, pre-stated convergence exclusion + cross-cohort r recomputation on the clean set.

EXCLUSION RULE (fixed thresholds, applied identically to every cohort — no cherry-picking):
  A trace is convergence-SEVERE if, over the ICC-driver SCALARS {mu, alpha, tau, c2, sigma_site},
      max R-hat > 1.05   OR   min bulk-ESS < 100.
  A (cohort, endpoint) is EXCLUDED from the ICC cross-cohort analyses if its MODEL B hallmark
  trace is SEVERE. (Model A gates nothing here; it is only used for gene-hit comparison.)

The headline cross-cohort correlations are recomputed on the full set AND the clean (post-
exclusion) set so any change in r is explicit. Pearson r is reported with n and the points it
is based on — never a line alone — because with ~20 cohorts and two patient-sharing endpoints
these are exploratory, not confirmatory.

Reconciliation: ICC/convergence are taken from the (possibly re-fit) inputs; the preserved-site
CV C-index is taken from the ORIGINAL summary because `--force` re-fits blank it out (the CV is
independent per-fold refitting, unchanged by re-fitting the in-sample model).

    python scripts/apply_exclusions.py \
        --triage results/analysis/convergence_triage.csv \
        --triage-extra /tmp/server_pull/convergence_triage.server.csv \
        --summary results/pan_cancer_summary.csv \
        --summary-extra /tmp/server_pull/pan_cancer_summary.server.csv

Writes results/analysis/excluded_cohorts.csv, results/analysis/crosscohort_r.csv,
and results/analysis/CONVERGENCE_LIMITATIONS.md (a ready-to-paste limitation paragraph).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from download import PROJECT_ROOT

RESULTS = PROJECT_ROOT / "results"
ANALYSIS = RESULTS / "analysis"
SEV_RHAT, SEV_ESS = 1.05, 100
KEY_T = ["cohort", "model", "endpoint", "panel", "p0"]
KEY_S = ["cohort", "endpoint", "panel", "p0"]


def _merge(base, extra, key):
    """extra rows override base rows on `key` (re-fit values win); base-only rows kept."""
    if extra is None or extra.empty:
        return base.copy()
    return (pd.concat([base, extra], ignore_index=True)
              .drop_duplicates(subset=key, keep="last").reset_index(drop=True))


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, int(ok.sum())
    return float(np.corrcoef(x[ok], y[ok])[0, 1]), int(ok.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--triage", default=str(ANALYSIS / "convergence_triage.csv"))
    ap.add_argument("--triage-extra", default=None, help="re-fit triage CSV that overrides --triage")
    ap.add_argument("--summary", default=str(RESULTS / "pan_cancer_summary.csv"))
    ap.add_argument("--summary-extra", default=None, help="re-fit summary (ICC) that overrides --summary")
    a = ap.parse_args()

    triage = _merge(pd.read_csv(a.triage),
                    pd.read_csv(a.triage_extra) if a.triage_extra and Path(a.triage_extra).exists() else None, KEY_T)
    base_sum = pd.read_csv(a.summary); base_sum = base_sum[base_sum["status"] == "ok"].copy()
    extra_sum = None
    if a.summary_extra and Path(a.summary_extra).exists():
        extra_sum = pd.read_csv(a.summary_extra); extra_sum = extra_sum[extra_sum["status"] == "ok"].copy()
    icc = _merge(base_sum, extra_sum, KEY_S)               # ICC: re-fit wins

    # C-index from ORIGINAL summary only (re-fit blanks it; CV is unchanged by re-fitting)
    cidx = base_sum[["cohort", "endpoint", "panel", "p0", "cindex_A", "cindex_B"]].copy()
    icc = icc.drop(columns=[c for c in ("cindex_A", "cindex_B") if c in icc.columns]).merge(
        cidx, on=KEY_S, how="left")

    # ---- mechanical SEVERE flag + exclusion (Model B hallmark) ----
    triage["scalar_severe"] = (triage["max_rhat_scalar"] > SEV_RHAT) | (triage["min_ess_scalar"] < SEV_ESS)
    modelb_hall = triage[(triage["model"] == "b") & (triage["panel"] == "hallmark")]
    excluded = (modelb_hall[modelb_hall["scalar_severe"]][["cohort", "endpoint", "max_rhat_scalar", "min_ess_scalar"]]
                .sort_values("max_rhat_scalar", ascending=False).reset_index(drop=True))
    excl_set = set(map(tuple, excluded[["cohort", "endpoint"]].values))
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    excluded.to_csv(ANALYSIS / "excluded_cohorts.csv", index=False)

    print("=" * 74)
    print(f"EXCLUSION RULE: Model-B hallmark SEVERE iff scalar max R-hat > {SEV_RHAT} OR "
          f"min bulk-ESS < {SEV_ESS}")
    print("=" * 74)
    print(f"EXCLUDED cohort×endpoints ({len(excluded)}):")
    for r in excluded.itertuples(index=False):
        print(f"  {r.cohort}/{r.endpoint}  (scalar R̂={r.max_rhat_scalar:.3f}, scalar ESS={int(r.min_ess_scalar)})")

    # ---- cross-cohort correlations: full vs clean ----
    H = icc[icc["panel"] == "hallmark"].copy()
    H["excluded"] = [(c, e) in excl_set for c, e in zip(H.cohort, H.endpoint)]
    rows = []
    def corr(label, df, xcol, ycol):
        r_full, n_full = pearson(df[xcol], df[ycol])
        clean = df[~df["excluded"]]
        r_cln, n_cln = pearson(clean[xcol], clean[ycol])
        rows.append({"relationship": label, "r_full": round(r_full, 3), "n_full": n_full,
                     "r_clean": round(r_cln, 3), "n_clean": n_cln,
                     "n_excluded": n_full - n_cln})
        print(f"\n{label}:")
        print(f"  full  : r={r_full:+.3f} (n={n_full})")
        print(f"  clean : r={r_cln:+.3f} (n={n_cln})   [{n_full - n_cln} excluded]")

    print("\n" + "=" * 74); print("CROSS-COHORT CORRELATIONS — full vs clean"); print("=" * 74)
    corr("ICC_genomic vs out-of-sample C-index (Model B)", H.dropna(subset=["ICC_genomic_mean", "cindex_B"]),
         "ICC_genomic_mean", "cindex_B")
    corr("ICC_site vs out-of-sample C-index (Model B)", H.dropna(subset=["ICC_site_mean", "cindex_B"]),
         "ICC_site_mean", "cindex_B")
    # OS-vs-PFI paired ICC_site: a cohort is excluded if EITHER endpoint's Model B is severe
    piv = H.pivot_table(index="cohort", columns="endpoint", values="ICC_site_mean")
    excl_cohorts = {c for (c, e) in excl_set}
    piv = piv.dropna(subset=[c for c in ("OS", "PFI") if c in piv.columns])
    if {"OS", "PFI"} <= set(piv.columns):
        r_full, n_full = pearson(piv["OS"], piv["PFI"])
        clean = piv[~piv.index.isin(excl_cohorts)]
        r_cln, n_cln = pearson(clean["OS"], clean["PFI"])
        rows.append({"relationship": "ICC_site OS vs PFI (paired)", "r_full": round(r_full, 3),
                     "n_full": n_full, "r_clean": round(r_cln, 3), "n_clean": n_cln,
                     "n_excluded": n_full - n_cln})
        print(f"\nICC_site OS vs PFI (paired):")
        print(f"  full  : r={r_full:+.3f} (n={n_full})")
        print(f"  clean : r={r_cln:+.3f} (n={n_cln})   [{n_full - n_cln} excluded]")
    rdf = pd.DataFrame(rows)
    rdf.to_csv(ANALYSIS / "crosscohort_r.csv", index=False)

    # ---- limitation write-up ----
    md = ANALYSIS / "CONVERGENCE_LIMITATIONS.md"
    excl_str = ", ".join(f"{r.cohort}/{r.endpoint}" for r in excluded.itertuples(index=False)) or "none"
    lines = [
        "# Convergence limitation (stated)\n",
        f"Convergence was triaged for every saved trace. A pre-specified, mechanical rule was applied "
        f"identically to all cohorts: a fit's **ICC-driver scalars** ({{mu, alpha, tau, c2, sigma_site}}) "
        f"must satisfy max R-hat ≤ {SEV_RHAT} and min bulk-ESS ≥ {SEV_ESS}; otherwise the fit is "
        f"flagged *convergence-severe* and the corresponding cohort×endpoint is **excluded** from the "
        f"ICC cross-cohort analyses (Model B gates ICC; Model A is used only for gene-hit comparison).\n",
        f"The regularized-horseshoe Weibull-AFT has a residual funnel geometry that higher `target_accept` "
        f"(0.999), longer tuning (2000), and more draws (2000) did not resolve for a minority of small / "
        f"low-event cohorts — a limit of the model parameterization, not of compute. Rather than alter the "
        f"prior post hoc, these fits are excluded under the rule above and reported transparently.\n",
        f"**Excluded ({len(excluded)} cohort×endpoints):** {excl_str}.\n",
        "Cross-cohort correlations are reported on both the full and the clean (post-exclusion) set "
        "(`results/analysis/crosscohort_r.csv`); see below. With ~20 cohorts and two patient-sharing "
        "endpoints these remain exploratory.\n",
        "\n| relationship | r (full) | n | r (clean) | n |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(f"| {r['relationship']} | {r['r_full']} | {r['n_full']} | {r['r_clean']} | {r['n_clean']} |")
    md.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {ANALYSIS/'excluded_cohorts.csv'}, {ANALYSIS/'crosscohort_r.csv'}, {md}")


if __name__ == "__main__":
    main()
