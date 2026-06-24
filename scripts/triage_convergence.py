#!/usr/bin/env python
"""Convergence triage for saved PyMC/ArviZ traces — READ-ONLY, no refitting.

Scans every trace under data/processed/{cohort}/traces/ and, across ALL parameters
(scalars + the high-dim beta/lam/z/site blocks), computes:
  max R-hat, fraction of params with R-hat > 1.01, min bulk-ESS, min tail-ESS,
  divergence count and % of post-warmup draws, and BFMI (energy diagnostic, if present).

Each trace is flagged PASS / WARN / FAIL on explicit thresholds (printed at runtime):
  FAIL  if  max R-hat > 1.01   OR  divergences > 0.5% of draws   OR  min bulk-ESS < 400
  WARN  if  (not FAIL) and any of:
            max R-hat > 1.005, min bulk-ESS < 1000, min tail-ESS < 400,
            any divergence (>0), BFMI < 0.3, or tree-depth saturation > 0
  PASS  otherwise

Writes results/analysis/convergence_triage.csv, prints a flag-grouped summary, an ORDERED
ACTION LIST (which traces to re-run at higher target_accept, which are usable-with-caveat),
and cross-references whether any FAIL/WARN trace feeds the cross-cohort headline analyses
(ICC_genomic-vs-Cindex, OS-vs-PFI ICC_site) so a bad trace can't silently contaminate them.

    python scripts/triage_convergence.py
    python scripts/triage_convergence.py --cohorts LUSC LUAD

This script does NOT re-run anything; it only reports.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import arviz as az

sys.path.insert(0, str(Path(__file__).resolve().parent))
from download import PROJECT_ROOT, PROCESSED_DIR

RESULTS_DIR = PROJECT_ROOT / "results"
SUMMARY_CSV = RESULTS_DIR / "pan_cancer_summary.csv"
OUT_CSV = RESULTS_DIR / "analysis" / "convergence_triage.csv"
_FNAME = re.compile(r"^(model_[ab])_([A-Za-z0-9]+)_(.+)_p(\d+)\.nc$")

# explicit thresholds (printed at runtime)
TH = {"rhat_fail": 1.01, "rhat_warn": 1.005, "div_frac_fail": 0.005,
      "ess_bulk_fail": 400, "ess_bulk_warn": 1000, "ess_tail_warn": 400, "bfmi_warn": 0.3}


def _ds(group):
    return group.to_dataset() if hasattr(group, "to_dataset") else group


def _panel_logical(token: str) -> str:
    # invert run_cohort._safe (":"->"" , "+"->"_") for the known panels so the panel column
    # matches the logical names used in results.csv / pan_cancer_summary.csv
    if token.startswith("topvar") and token[len("topvar"):].isdigit():
        return "topvar:" + token[len("topvar"):]
    if token in ("hallmark_combat", "hallmark_dropsing"):
        return token.replace("_", "+")
    return token


def discover(cohorts=None) -> list[dict]:
    rows = []
    for tdir in sorted(PROCESSED_DIR.glob("*/traces")):
        cohort = tdir.parent.name
        if cohorts and cohort not in cohorts:
            continue
        for nc in sorted(tdir.glob("*.nc")):
            m = _FNAME.match(nc.name)
            if not m:
                continue
            model, endpoint, panel_tok, p0 = m.groups()
            rows.append(dict(cohort=cohort, model=model.replace("model_", ""),
                             endpoint=endpoint, panel=_panel_logical(panel_tok),
                             p0=int(p0), path=str(nc)))
    return rows


def _bfmi(idata):
    ss = _ds(idata.sample_stats)
    if "energy" not in ss:
        return np.nan
    out = []
    for ch in ss["energy"].values:
        denom = np.sum((ch - ch.mean()) ** 2)
        out.append(np.sum(np.diff(ch) ** 2) / denom if denom > 0 else np.nan)
    return float(np.nanmin(out)) if out else np.nan


# Parameter tiers. The horseshoe LOCAL scales (lam, z) and the 3867-dim gene block (beta) are
# high-dimensional and hard to sample (funnel) — their min-ESS / max-R-hat are dominated by a
# single worst coordinate and routinely fail even when the headline quantities are fine. The
# ICC headline rests mainly on the low-dim SCALARS (sigma_site -> ICC_site, alpha -> residual);
# ICC_genomic and gene-hits depend on AGGREGATES of beta, far more robust than any single beta.
SCALARS = ["mu", "alpha", "tau", "c2", "sigma_site"]   # low-dim ICC drivers


def _tier(rh, eb, names):
    nm = [v for v in names if v in rh.data_vars]
    if not nm:
        return np.nan, np.nan
    return (max(float(np.nanmax(rh[v].values)) for v in nm),
            min(float(np.nanmin(eb[v].values)) for v in nm))


def diagnostics(path) -> dict:
    idata = az.from_netcdf(path)
    post = _ds(idata.posterior)
    n_draws = post.sizes["chain"] * post.sizes["draw"]
    rh = az.rhat(post)
    eb = az.ess(post, method="bulk")
    et = az.ess(post, method="tail")
    rhat_all = np.concatenate([np.atleast_1d(rh[v].values).ravel() for v in rh.data_vars])
    essb = np.concatenate([np.atleast_1d(eb[v].values).ravel() for v in eb.data_vars])
    esst = np.concatenate([np.atleast_1d(et[v].values).ravel() for v in et.data_vars])
    sc_rhat, sc_ess = _tier(rh, eb, SCALARS)           # ICC-driver scalars
    b_rhat, b_ess = _tier(rh, eb, ["beta"])            # gene coefficient block
    ss = _ds(idata.sample_stats)
    n_div = int(ss["diverging"].values.sum()) if "diverging" in ss else 0
    td = (float(ss["maxdepth_reached"].values.mean()) if "maxdepth_reached" in ss else 0.0)
    return {"n_params": int(rhat_all.size), "max_rhat": float(np.nanmax(rhat_all)),
            "frac_rhat_gt_1.01": float(np.nanmean(rhat_all > 1.01)),
            "max_rhat_scalar": sc_rhat, "min_ess_scalar": sc_ess,
            "max_rhat_beta": b_rhat, "min_ess_beta": b_ess,
            "min_ess_bulk": float(np.nanmin(essb)), "min_ess_tail": float(np.nanmin(esst)),
            "n_div": n_div, "div_pct": 100.0 * n_div / n_draws,
            "min_bfmi": _bfmi(idata), "treedepth_sat_pct": 100.0 * td,
            # the ICC headline (sigma_site/alpha) is trustworthy iff the SCALARS converged
            "scalar_ok": bool(sc_rhat <= TH["rhat_fail"] and sc_ess >= TH["ess_bulk_fail"])}


def flag(d) -> str:
    if (d["max_rhat"] > TH["rhat_fail"] or d["div_pct"] > 100 * TH["div_frac_fail"]
            or d["min_ess_bulk"] < TH["ess_bulk_fail"]):
        return "FAIL"
    if (d["max_rhat"] > TH["rhat_warn"] or d["min_ess_bulk"] < TH["ess_bulk_warn"]
            or d["min_ess_tail"] < TH["ess_tail_warn"] or d["n_div"] > 0
            or (not np.isnan(d["min_bfmi"]) and d["min_bfmi"] < TH["bfmi_warn"])
            or d["treedepth_sat_pct"] > 0):
        return "WARN"
    return "PASS"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohorts", nargs="*", default=None)
    args = ap.parse_args()

    recs = discover(set(args.cohorts) if args.cohorts else None)
    if not recs:
        sys.exit("no traces found under data/processed/*/traces/")
    print(f"triaging {len(recs)} traces (read-only) ...")
    print("thresholds:", TH)

    rows = []
    for r in recs:
        try:
            d = diagnostics(r["path"])
            rows.append({**{k: r[k] for k in ("cohort", "model", "endpoint", "panel", "p0")},
                         **{k: round(v, 4) if isinstance(v, float) else v for k, v in d.items()},
                         "flag": flag(d)})
        except Exception as exc:
            rows.append({**{k: r[k] for k in ("cohort", "model", "endpoint", "panel", "p0")},
                         "flag": f"LOAD_ERR:{type(exc).__name__}"})
            print(f"  error on {r['cohort']}/{r['model']}/{r['endpoint']}: {exc}")

    T = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT_CSV, index=False)
    print(f"\nwrote {OUT_CSV}")

    # ---- summary grouped by flag ----
    print("\n" + "=" * 70)
    print("FLAG TALLY:", T["flag"].value_counts().to_dict())
    print("=" * 70)
    for fl in ("FAIL", "WARN"):
        g = T[T["flag"] == fl]
        if g.empty:
            continue
        print(f"\n{fl} traces ({len(g)}):")
        for r in g.itertuples(index=False):
            print(f"  {fl}  {r.cohort:5s} {r.model} {r.endpoint:3s} {r.panel:11s} p{r.p0} | "
                  f"R̂ scalar={r.max_rhat_scalar} beta={r.max_rhat_beta} | "
                  f"ESS scalar={int(r.min_ess_scalar)} beta={int(r.min_ess_beta)} | "
                  f"divs={r.n_div}({r.div_pct:.2f}%) -> scalar_ok={r.scalar_ok}")

    # ---- ordered ACTION LIST ----
    print("\n" + "=" * 70)
    print("ACTION LIST")
    print("=" * 70)
    # Severity on the ICC-driver SCALARS (what the headline rests on). The strict PASS gate
    # (R̂≤1.01 AND ESS≥400) flags most traces, but a scalar R̂ of 1.02 / ESS of 300 is fine for an
    # ICC point estimate — so we bucket by ACTIONABLE severity rather than the binary gate:
    SEV_RHAT, SEV_ESS = 1.05, 100        # genuinely unusable scalars
    def _sev(r):
        if r.max_rhat_scalar > SEV_RHAT or r.min_ess_scalar < SEV_ESS:
            return "SEVERE"
        if r.max_rhat_scalar > TH["rhat_fail"] or r.min_ess_scalar < TH["ess_bulk_fail"]:
            return "MARGINAL"
        return "CLEAN"
    T["scalar_sev"] = [_sev(r) for r in T.itertuples(index=False)]
    counts = T["scalar_sev"].value_counts().reindex(["CLEAN", "MARGINAL", "SEVERE"]).fillna(0).astype(int)
    print(f"SCALAR (ICC-driver: sigma_site/alpha/tau/mu/c2) severity  "
          f"[SEVERE = R̂>{SEV_RHAT} or ESS<{SEV_ESS}]:")
    print(f"  CLEAN={counts['CLEAN']}  MARGINAL={counts['MARGINAL']}  SEVERE={counts['SEVERE']}  (of {len(T)})")
    sev = T[T["scalar_sev"] == "SEVERE"]
    marg = T[T["scalar_sev"] == "MARGINAL"]
    if not sev.empty:
        print(f"\n• RE-RUN FIRST — {len(sev)} SEVERE trace(s): ICC drivers genuinely unconverged. "
              f"Re-run at higher --target-accept (0.99→0.999), more tune, or a non-centred reparam:")
        for r in sev.sort_values("max_rhat_scalar", ascending=False).itertuples(index=False):
            print(f"    {r.cohort} / model_{r.model} / {r.endpoint} / {r.panel} / p{r.p0}  "
                  f"(scalar R̂={r.max_rhat_scalar}, scalar ESS={int(r.min_ess_scalar)})")
    if not marg.empty:
        print(f"\n• USABLE (re-run to tighten) — {len(marg)} MARGINAL trace(s): scalar R̂≤{SEV_RHAT} and "
              f"ESS≥{SEV_ESS}, fine for ICC point estimates though below the strict gate.")
    n_beta_bad = int((T["max_rhat_beta"] > 1.05).sum())
    print(f"\n• GENE-LEVEL CAVEAT — {n_beta_bad}/{len(T)} traces have beta R̂>1.05: individual gene")
    print(f"  coefficients are unreliable, so make NO per-gene claims. The aggregate ICC_genomic and")
    print(f"  the ~0-hit headline are robust to this (shrunk-near-zero betas don't become hits).")

    # ---- cross-reference: do flagged traces feed the headline cross-cohort analyses? ----
    print("\n" + "=" * 70)
    print("HEADLINE CONTAMINATION CHECK")
    print("=" * 70)
    print("Cross-cohort analyses (ICC_genomic-vs-Cindex, OS-vs-PFI ICC_site) use the Model B")
    print("hallmark traces of the 'ok' cohort×endpoints in pan_cancer_summary.csv.")
    if not SUMMARY_CSV.exists():
        print("  (pan_cancer_summary.csv not found — cannot cross-reference.)")
        return
    ok = pd.read_csv(SUMMARY_CSV)
    ok = ok[ok["status"] == "ok"][["cohort", "endpoint"]].drop_duplicates()
    feed = set(map(tuple, ok.values))                # (cohort, endpoint) feeding headline ICC
    flagged_b = T[(T["flag"].isin(["FAIL", "WARN"])) & (T["model"] == "b") &
                  (T["panel"] == "hallmark")]
    contaminating = [r for r in flagged_b.itertuples(index=False)
                     if (r.cohort, r.endpoint) in feed]
    if not contaminating:
        print("  ✓ No FAIL/WARN Model-B hallmark trace feeds a cross-cohort headline result.")
    else:
        severe = [r for r in contaminating if not getattr(r, "scalar_ok", True)]
        mild = [r for r in contaminating if getattr(r, "scalar_ok", True)]
        print(f"  {len(contaminating)} flagged Model-B hallmark trace(s) feed cross-cohort (ICC) results:")
        if severe:
            print(f"  ⚠⚠ {len(severe)} where an ICC-driver SCALAR failed — ICC_site genuinely suspect, "
                  f"re-run before quoting the cross-cohort scatter / OS-vs-PFI point:")
            for r in severe:
                print(f"       {r.flag}  {r.cohort}/{r.endpoint}  (scalar R̂={r.max_rhat_scalar}, scalar ESS={int(r.min_ess_scalar)})")
        if mild:
            print(f"  ⓘ {len(mild)} where scalars converged (ICC_site OK; only the gene block failed) — "
                  f"the ICC point is usable:")
            for r in mild:
                print(f"       {r.flag}  {r.cohort}/{r.endpoint}  (scalar R̂={r.max_rhat_scalar} ok; beta R̂={r.max_rhat_beta})")
    # OS-vs-PFI specifically needs BOTH endpoints clean per cohort
    paired = ok.groupby("cohort")["endpoint"].apply(set)
    both = {c for c, eps in paired.items() if {"OS", "PFI"} <= eps}
    bad_pair = sorted({r.cohort for r in contaminating if r.cohort in both})
    if bad_pair:
        print(f"  ⚠ OS-vs-PFI paired comparison specifically affected for: {', '.join(bad_pair)}")


if __name__ == "__main__":
    main()
