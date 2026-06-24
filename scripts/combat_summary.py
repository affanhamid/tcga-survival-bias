#!/usr/bin/env python
"""ComBat sensitivity summary — does the site effect survive batch correction?

Builds the per-(cohort,endpoint) 3-way ICC comparison across panels:
  full     = hallmark            (all patients)
  drop     = hallmark+dropsing   (singleton sites dropped, UNCORRECTED)
  combat   = hallmark+combat     (singleton sites dropped, ComBat-CORRECTED)

The HEADLINE ComBat contrast is drop -> combat (the SAME patients), NOT full -> combat:
comparing combat to the full original conflates the correction with the dropped singleton
sites. We report delta_combat = ICC_site(combat) - ICC_site(drop) as the clean ComBat effect,
and (for context only) delta_dropset = ICC_site(drop) - ICC_site(full) as the patient-set
effect. The "site survives ComBat" sentence rests on delta_combat.

The SAME mechanical convergence-exclusion rule used elsewhere is applied here, with no special
cases: a (cohort,endpoint) is excluded from the ComBat claim if EITHER its hallmark+dropsing
or hallmark+combat Model-B trace is convergence-severe (scalar max R-hat > 1.05 OR
min bulk-ESS < 100). PFI is the primary endpoint.

    python scripts/combat_summary.py \
        --summary results/pan_cancer_summary.csv --summary-extra /tmp/server_pull/summary.csv \
        --triage  results/analysis/convergence_triage.csv --triage-extra /tmp/server_pull/triage.csv

Writes results/analysis/combat_sensitivity.csv and results/analysis/COMBAT_SENSITIVITY.md.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from download import PROJECT_ROOT

ANALYSIS = PROJECT_ROOT / "results" / "analysis"
SEV_RHAT, SEV_ESS = 1.05, 100
KEY_T = ["cohort", "model", "endpoint", "panel", "p0"]
KEY_S = ["cohort", "endpoint", "panel", "p0"]
PANELS = {"hallmark": "full", "hallmark+dropsing": "drop", "hallmark+combat": "combat"}


def _merge(base, extra, key):
    if extra is None or extra.empty:
        return base.copy()
    return (pd.concat([base, extra], ignore_index=True)
              .drop_duplicates(subset=key, keep="last").reset_index(drop=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", default=str(PROJECT_ROOT / "results" / "pan_cancer_summary.csv"))
    ap.add_argument("--summary-extra", default=None)
    ap.add_argument("--triage", default=str(ANALYSIS / "convergence_triage.csv"))
    ap.add_argument("--triage-extra", default=None)
    a = ap.parse_args()

    s = pd.read_csv(a.summary); s = s[s["status"] == "ok"].copy()
    se = None
    if a.summary_extra and Path(a.summary_extra).exists():
        se = pd.read_csv(a.summary_extra); se = se[se["status"] == "ok"].copy()
    icc = _merge(s, se, KEY_S)
    icc = icc[icc["panel"].isin(PANELS)]

    tri = _merge(pd.read_csv(a.triage),
                 pd.read_csv(a.triage_extra) if a.triage_extra and Path(a.triage_extra).exists() else None, KEY_T)
    tri["severe"] = (tri["max_rhat_scalar"] > SEV_RHAT) | (tri["min_ess_scalar"] < SEV_ESS)
    sev_b = {(r.cohort, r.endpoint, r.panel) for r in
             tri[(tri.model == "b") & tri.severe].itertuples(index=False)}

    # wide ICC_site / ICC_genomic per (cohort, endpoint) across the three panels
    rows = []
    for (c, e), g in icc.groupby(["cohort", "endpoint"]):
        rec = {"cohort": c, "endpoint": e}
        have = {}
        for panel, tag in PANELS.items():
            sub = g[g["panel"] == panel]
            if not sub.empty:
                rec[f"ICCsite_{tag}"] = round(float(sub["ICC_site_mean"].iloc[0]), 4)
                rec[f"ICCgen_{tag}"] = round(float(sub["ICC_genomic_mean"].iloc[0]), 4)
                have[tag] = True
        if {"drop", "combat"} <= have.keys():
            rec["delta_combat"] = round(rec["ICCsite_combat"] - rec["ICCsite_drop"], 4)   # HEADLINE
        if {"full", "drop"} <= have.keys():
            rec["delta_dropset"] = round(rec["ICCsite_drop"] - rec["ICCsite_full"], 4)     # context only
        # exclusion: need BOTH drop and combat Model-B converged for a clean contrast
        rec["combat_excluded"] = ((c, e, "hallmark+dropsing") in sev_b) or ((c, e, "hallmark+combat") in sev_b)
        rows.append(rec)
    T = pd.DataFrame(rows).sort_values(["endpoint", "cohort"]).reset_index(drop=True)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    T.to_csv(ANALYSIS / "combat_sensitivity.csv", index=False)

    print("=" * 78)
    print("COMBAT SENSITIVITY — headline contrast = drop -> combat (SAME patients)")
    print("=" * 78)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(T.to_string(index=False))

    # headline on the converged set, primary endpoint PFI
    def headline(df, ep):
        d = df[(df.endpoint == ep) & (~df.combat_excluded) & df["delta_combat"].notna()]
        if d.empty:
            return f"  {ep}: no converged cohorts with both drop+combat yet."
        return (f"  {ep}: ComBat changes ICC_site by median {d['delta_combat'].median():+.3f} "
                f"(range {d['delta_combat'].min():+.3f}..{d['delta_combat'].max():+.3f}) "
                f"across {len(d)} converged cohorts [same-patient drop->combat].")
    print("\nHEADLINE (converged set; delta_combat = ICC_site(combat) - ICC_site(drop)):")
    for ep in ("PFI", "OS"):
        print(headline(T, ep))
    excl = T[T.combat_excluded]
    if len(excl):
        print(f"\nExcluded from ComBat claim (severe drop/combat Model-B): "
              f"{', '.join(c+'/'+e for c,e in zip(excl.cohort, excl.endpoint))}")

    # write-up
    md = ANALYSIS / "COMBAT_SENSITIVITY.md"
    lines = ["# ComBat sensitivity (matched-baseline)\n",
             "The ComBat contrast is reported as **dropped-uncorrected (`hallmark+dropsing`) vs "
             "dropped-corrected (`hallmark+combat`)** — the same singleton-site-filtered patient set — "
             "so the change in ICC_site is attributable to ComBat, not to the dropped sites. The "
             "full-cohort `hallmark` ICC is shown for context only. The same mechanical convergence "
             f"exclusion (scalar R-hat>{SEV_RHAT} or bulk-ESS<{SEV_ESS}) is applied to the drop/combat "
             "Model-B traces.\n",
             "\n| cohort/endpoint | ICC_site full | drop | combat | Δcombat (drop→combat) | excluded |",
             "|---|---|---|---|---|---|"]
    for r in T.itertuples(index=False):
        lines.append(f"| {r.cohort}/{r.endpoint} | {getattr(r,'ICCsite_full','—')} | "
                     f"{getattr(r,'ICCsite_drop','—')} | {getattr(r,'ICCsite_combat','—')} | "
                     f"{getattr(r,'delta_combat','—')} | {'yes' if r.combat_excluded else ''} |")
    lines += ["", "Headline (converged set):", headline(T, "PFI"), headline(T, "OS")]
    md.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {ANALYSIS/'combat_sensitivity.csv'} and {md}")


if __name__ == "__main__":
    main()
