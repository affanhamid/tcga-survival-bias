#!/usr/bin/env python
"""LOO model comparison: Model A (naive) vs Model B (site-adjusted) per cohort/endpoint.

Requires the traces to carry a `log_likelihood` group — run scripts/add_loglik.py first.
For each (cohort, endpoint, panel, p0) with BOTH models present and both carrying
log_likelihood, this runs PSIS-LOO-CV (`az.loo`) on each model and `az.compare` on the pair,
and records the ELPD difference, its standard error (dse), the stacking weights, and the
winning model.

Pareto-k is surfaced explicitly: censored survival likelihoods routinely produce high k
(influential/heavy-tailed importance ratios). We report the count and fraction of k > 0.7 per
model and FLAG the comparison as LOO-unreliable when too many points exceed it, rather than
reporting a clean-looking ELPD that PSIS cannot actually support.

    python scripts/model_compare.py
    python scripts/model_compare.py --cohorts LUSC LUAD
    python scripts/model_compare.py --k-bad-frac 0.10     # unreliable threshold (default 0.05)

Output: results/analysis/model_comparison_loo.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import arviz as az

sys.path.insert(0, str(Path(__file__).resolve().parent))
from add_loglik import discover, has_loglik, LL_VAR
from download import PROJECT_ROOT

OUT_CSV = PROJECT_ROOT / "results" / "analysis" / "model_comparison_loo.csv"
K_GOOD = 0.7                                   # Pareto-k threshold above which a point is "bad"


def _loo(idata):
    """az.loo with pointwise Pareto-k; returns (elpd, se, p_loo, k_array)."""
    r = az.loo(idata, var_name=LL_VAR, pointwise=True)
    return float(r.elpd), float(r.se), float(r.p), np.asarray(r.pareto_k.values).ravel()


def compare_pair(rec_a, rec_b, k_bad_frac):
    A = az.from_netcdf(rec_a["path"])
    B = az.from_netcdf(rec_b["path"])
    elpd_a, se_a, p_a, k_a = _loo(A)
    elpd_b, se_b, p_b, k_b = _loo(B)

    cmp = az.compare({"A": A, "B": B}, var_name=LL_VAR, round_to="none")
    winner = cmp.index[cmp["rank"] == 0][0]
    loser = "A" if winner == "B" else "B"
    dse = float(cmp.loc[loser, "dse"])                     # SE of the ELPD difference
    weight_a = float(cmp.loc["A", "weight"]); weight_b = float(cmp.loc["B", "weight"])

    n = len(k_a)
    bad_a, bad_b = int((k_a > K_GOOD).sum()), int((k_b > K_GOOD).sum())
    frac_a, frac_b = bad_a / n, bad_b / n
    reliable = (frac_a <= k_bad_frac) and (frac_b <= k_bad_frac)

    return {
        "cohort": rec_a["cohort"], "endpoint": rec_a["endpoint"], "panel": rec_a["panel"],
        "p0": rec_a["p0"], "n_obs": n,
        "elpd_A": round(elpd_a, 2), "se_A": round(se_a, 2), "p_loo_A": round(p_a, 2),
        "elpd_B": round(elpd_b, 2), "se_B": round(se_b, 2), "p_loo_B": round(p_b, 2),
        "elpd_diff_B_minus_A": round(elpd_b - elpd_a, 2), "dse": round(dse, 2),
        "winner": winner, "weight_A": round(weight_a, 3), "weight_B": round(weight_b, 3),
        "pareto_k_bad_A": bad_a, "pareto_k_bad_pct_A": round(100 * frac_a, 1),
        "pareto_k_bad_B": bad_b, "pareto_k_bad_pct_B": round(100 * frac_b, 1),
        "max_k_A": round(float(k_a.max()), 3), "max_k_B": round(float(k_b.max()), 3),
        "loo_reliable": reliable,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohorts", nargs="*", default=None)
    ap.add_argument("--k-bad-frac", type=float, default=0.05,
                    help="fraction of Pareto-k>0.7 above which LOO is flagged unreliable")
    args = ap.parse_args()

    recs = discover(set(args.cohorts) if args.cohorts else None)
    # index by (cohort, endpoint, panel, p0) -> {model: rec}
    pairs: dict[tuple, dict] = {}
    for r in recs:
        pairs.setdefault((r["cohort"], r["endpoint"], r["panel"], r["p0"]), {})[r["model"]] = r

    rows, skipped = [], []
    for key, mods in sorted(pairs.items()):
        ck, ek, pk, p0 = key
        if "a" not in mods or "b" not in mods:
            skipped.append((key, "missing " + ("A" if "a" not in mods else "B"))); continue
        if not has_loglik(mods["a"]["path"]) or not has_loglik(mods["b"]["path"]):
            skipped.append((key, "no log_likelihood (run add_loglik.py)")); continue
        try:
            rows.append(compare_pair(mods["a"], mods["b"], args.k_bad_frac))
            r = rows[-1]
            flag = "" if r["loo_reliable"] else "  ⚠LOO-UNRELIABLE (high Pareto-k)"
            print(f"  {ck:5s} {ek:3s} {pk:9s} p{p0}: winner={r['winner']} "
                  f"ΔelpdB-A={r['elpd_diff_B_minus_A']:+.2f}±{r['dse']:.2f} "
                  f"| k>0.7: A={r['pareto_k_bad_A']} B={r['pareto_k_bad_B']} of {r['n_obs']}{flag}")
        except Exception as exc:
            skipped.append((key, f"{type(exc).__name__}: {exc}"))

    if rows:
        df = pd.DataFrame(rows).sort_values(["endpoint", "cohort"]).reset_index(drop=True)
        OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(OUT_CSV, index=False)
        print(f"\nwrote {OUT_CSV}  ({len(df)} comparisons)")
        n_b = int((df["winner"] == "B").sum()); n_unrel = int((~df["loo_reliable"]).sum())
        print(f"SUMMARY — Model B (site-adjusted) wins {n_b}/{len(df)} comparisons by LOO-ELPD.")
        decisive = df[(df["winner"] == "B") & (df["elpd_diff_B_minus_A"].abs() > 2 * df["dse"]) & df["loo_reliable"]]
        print(f"  decisive (|Δelpd| > 2·dse, reliable): {len(decisive)} favour B.")
        if n_unrel:
            print(f"  ⚠ {n_unrel} comparison(s) flagged LOO-UNRELIABLE (>{args.k_bad_frac:.0%} Pareto-k>0.7) — "
                  "treat those ELPDs with caution / consider reloo or moment-matching.")
    else:
        print("\nNo comparable pairs found.")
    if skipped:
        print(f"\nSKIPPED {len(skipped)}:")
        for key, why in skipped:
            print(f"  {'/'.join(map(str, key))}: {why}")


if __name__ == "__main__":
    main()
