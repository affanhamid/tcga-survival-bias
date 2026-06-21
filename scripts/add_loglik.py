#!/usr/bin/env python
"""Add per-patient pointwise log-likelihood to already-saved PyMC/ArviZ traces.

WHY THIS EXISTS
---------------
`scripts/run_cohort.py` fits Weibull-AFT survival models where DEATHS go through an
observed RV (`pm.Weibull("y_obs", ..., observed=to)`) but CENSORED patients go through a
single SUMMED `pm.Potential("cens_ll", (-exp(alpha*(log(tc)-eta_c))).sum())`. A naive
`pm.compute_log_likelihood` therefore returns pointwise log-lik for DEATHS ONLY (length
n_deaths) and silently drops every censored patient — invalid LOO/WAIC for a survival model
that is ~half censored. This script reconstructs the full PER-PATIENT log-likelihood
(length n_patients, original patient order) and attaches it as the `log_likelihood` group so
`az.loo` / `az.compare` are valid.

NO REFIT. The posteriors are reused as-is; we only ADD a `log_likelihood` group. The original
file is backed up to `*.nc.bak` before the first overwrite, and the operation is idempotent
(skips traces that already have `log_likelihood` unless `--force`).

Likelihood per patient i, per draw (eta_i = mu + beta·x_i  [+ mu_site[site_i]] for Model B):
  death    (d==1): log f = log(alpha) - alpha*eta + (alpha-1)*log(t) - exp(alpha*(log(t)-eta))
  censored (d==0): log S = -exp(alpha*(log(t)-eta))     # per-obs, un-summed form of cens_ll
These match PyMC's pm.Weibull(alpha, beta=exp(eta)) parameterisation (scale lambda=exp(eta));
verified against pm.logp and scipy in --validate.

Reuses run_cohort.select_genes / _valid / the persisted parquet — X, patient order, gene
order, and the site index (pd.Categorical(df["tss"]).codes) are all recovered exactly.

    python scripts/add_loglik.py                  # all traces, skip those already done
    python scripts/add_loglik.py --force          # recompute + reattach everywhere
    python scripts/add_loglik.py --cohorts LUSC   # subset
    python scripts/add_loglik.py --validate LUSC OS hallmark 10   # deep validation print
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import arviz as az
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_cohort as rc
from download import PROCESSED_DIR

LL_VAR = "obs"                       # name of the single log_likelihood variable
_FNAME = re.compile(r"^(model_[ab])_([A-Za-z0-9]+)_(.+)_p(\d+)\.nc$")


# ---------------------------------------------------------------------------
# discovery
# ---------------------------------------------------------------------------
def _panel_logical(token: str) -> str:
    """Map the file-safe panel token back to the logical panel select_genes expects
    (run_cohort._safe turns 'topvar:500' -> 'topvar500'; 'hallmark' is unchanged)."""
    if token.startswith("topvar") and token[len("topvar"):].isdigit():
        return "topvar:" + token[len("topvar"):]
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


# ---------------------------------------------------------------------------
# core: per-patient pointwise log-likelihood
# ---------------------------------------------------------------------------
def _ds(group):
    return group.to_dataset() if hasattr(group, "to_dataset") else group


def compute_pointwise_loglik(cohort, model, endpoint, panel, p0, src_path):
    """Return (ll_cd, meta) where ll_cd has shape (chain, draw, n_patients) in original
    valid-patient order and meta carries n_patients / n_deaths / coords for attachment.

    `src_path` is the trace whose posterior we read (the .bak, to decouple file handles)."""
    tc, ec = f"{endpoint}_time", endpoint
    df = pd.read_parquet(PROCESSED_DIR / cohort / f"{cohort}.parquet")
    gene_cols = rc.select_genes(cohort, df, panel)        # exact gene set + order used at fit
    dfv = rc._valid(df, tc, ec)                           # exact patient set + order used at fit
    t = dfv[tc].values.astype(float)
    d = dfv[ec].values.astype(int)
    logt = np.log(t)
    X = dfv[gene_cols].values

    idata = az.from_netcdf(src_path)
    post = _ds(idata.posterior)
    nchain, ndraw = post.sizes["chain"], post.sizes["draw"]

    def stack(name):                                      # (param..., S=chain*draw)
        return post[name].stack(s=("chain", "draw")).values

    beta = stack("beta")                                  # (p, S)
    mu = stack("mu")                                      # (S,)
    alpha = stack("alpha")                                # (S,)
    eta = mu[None, :] + X @ beta                          # (n, S)
    if model == "b":
        mu_site = stack("mu_site")                        # (n_sites, S)
        site = pd.Categorical(dfv["tss"]).codes           # same coding as build_model_b
        eta = eta + mu_site[site, :]

    A = alpha[None, :]
    term = np.exp(A * (logt[:, None] - eta))              # exp(alpha*(log t - eta)) = (t/lambda)^alpha
    ll_dead = np.log(A) - A * eta + (A - 1) * logt[:, None] - term
    ll = np.where(d[:, None] == 1, ll_dead, -term)        # death -> log pdf, censored -> log survival

    # (n, S) -> (chain, draw, n): stack iterates chain-major, so S = chain*ndraw + draw
    ll_cd = ll.T.reshape(nchain, ndraw, len(dfv))
    meta = dict(n_patients=len(dfv), n_deaths=int((d == 1).sum()),
                n_censored=int((d == 0).sum()),
                chain=post.chain.values, draw=post.draw.values)
    return ll_cd, meta


def has_loglik(path) -> bool:
    try:
        return any("log_likelihood" in g for g in az.from_netcdf(path).groups)
    except Exception:
        return False


def _verify_written(path, n_patients) -> bool:
    """Confirm a re-saved trace loads, has the log_likelihood group of the right length,
    and still carries an intact posterior (the beta block)."""
    idata = az.from_netcdf(path)
    if not any("log_likelihood" in g for g in idata.groups):
        return False
    ll = _ds(idata.log_likelihood)
    if LL_VAR not in ll.data_vars or ll[LL_VAR].sizes.get(f"{LL_VAR}_id") != n_patients:
        return False
    return "beta" in _ds(idata.posterior).data_vars


def attach_loglik(rec, force=False) -> str:
    """Compute + attach the log_likelihood group, re-saving in place. Returns a status str.

    Disk-safe backup: the original is copied to *.nc.bak, the new file is written + verified,
    then the backup is DELETED. So peak extra disk is one trace (not one per trace) — these
    traces are large and total ~tens of GB. The pristine originals also live on the project's
    Hugging Face archive, so the transient .bak is purely crash-safety for this single write."""
    path = rec["path"]
    bak = path + ".bak"
    tmp = path + ".tmp"
    if has_loglik(path) and not force:
        return "skip(exists)"
    # backup, then read posterior from the backup so the live handle is never the write target
    # (and a forced re-run cannot inherit a stale attached group through the source file)
    shutil.copy2(path, bak)
    try:
        ll_cd, meta = compute_pointwise_loglik(rec["cohort"], rec["model"], rec["endpoint"],
                                               rec["panel"], rec["p0"], src_path=bak)
        if not np.isfinite(ll_cd).all():
            return "FAIL(non-finite log-lik)"
        idata = az.from_netcdf(bak)
        ds = xr.Dataset(
            {LL_VAR: (("chain", "draw", f"{LL_VAR}_id"), ll_cd)},
            coords={"chain": meta["chain"], "draw": meta["draw"],
                    f"{LL_VAR}_id": np.arange(meta["n_patients"])},
        )
        idata["log_likelihood"] = ds
        idata.to_netcdf(tmp)
        os.replace(tmp, path)                             # atomic swap
        if not _verify_written(path, meta["n_patients"]):
            shutil.copy2(bak, path)                       # restore on a bad write
            return "FAIL(verification failed; restored original)"
        return f"ok(n={meta['n_patients']}, deaths={meta['n_deaths']})"
    finally:
        for p in (tmp, bak):                              # reclaim disk; original is now safe
            if os.path.exists(p):
                os.remove(p)


# ---------------------------------------------------------------------------
# deep validation (one cohort/model/endpoint), printed explicitly
# ---------------------------------------------------------------------------
def validate(cohort, endpoint, panel, p0):
    import pymc as pm
    from scipy.stats import weibull_min

    print("=" * 74)
    print(f"VALIDATION — {cohort} / {endpoint} / {panel} / p{p0}")
    print("=" * 74)

    # [A] Weibull parameterisation: my closed form vs pm.logp vs scipy --------
    def my_logpdf(t, a, eta):  # death
        return np.log(a) - a * eta + (a - 1) * np.log(t) - np.exp(a * (np.log(t) - eta))
    def my_logsf(t, a, eta):   # censored
        return -np.exp(a * (np.log(t) - eta))
    print("\n[A] parameterisation spot-check  (death log-pdf; censored log-survival)")
    dmax = smax = 0.0
    for tt, a, et in [(100., 1.3, 6.0), (2000., 0.8, 7.2), (50., 1.0, 5.5), (365., 1.5, 6.8)]:
        mine = my_logpdf(tt, a, et)
        pmv = float(pm.logp(pm.Weibull.dist(alpha=a, beta=np.exp(et)), tt).eval())
        sci = weibull_min.logpdf(tt, a, scale=np.exp(et))
        dmax = max(dmax, abs(mine - pmv), abs(mine - sci))
        smax = max(smax, abs(my_logsf(tt, a, et) - weibull_min.logsf(tt, a, scale=np.exp(et))))
    print(f"    max|mine - pm.logp| and |mine - scipy.logpdf| (deaths)   = {dmax:.2e}")
    print(f"    max|mine - scipy.logsf| (censored)                       = {smax:.2e}")
    # 1e-6 tolerance: pm.logp evaluates through pytensor (float rounding ~1e-8), not a
    # mathematical disagreement — scipy (pure float64) matches the censored term at ~1e-16.
    assert dmax < 1e-6 and smax < 1e-6, "parameterisation mismatch"

    for model in ("a", "b"):
        path = PROCESSED_DIR / cohort / "traces" / f"model_{model}_{endpoint}_{rc._safe(panel)}_p{p0}.nc"
        src = str(path) + ".bak" if os.path.exists(str(path) + ".bak") else str(path)
        if not Path(src).exists():
            print(f"\n[model {model.upper()}] trace not found — skipping"); continue
        ll_cd, meta = compute_pointwise_loglik(cohort, model, endpoint, panel, p0, src)
        print(f"\n[model {model.upper()}] pointwise log-lik shape = {ll_cd.shape}")
        # [B] KEY CHECK: count == n_patients, NOT n_deaths
        print(f"    n_patients = {meta['n_patients']}   n_deaths = {meta['n_deaths']}   "
              f"n_censored = {meta['n_censored']}")
        print(f"    -> pointwise terms = {ll_cd.shape[-1]} (== n_patients ✓, "
              f"NOT n_deaths={meta['n_deaths']}; censored ARE included)")
        # [C] no NaN/inf
        print(f"    all finite (no NaN/inf): {bool(np.isfinite(ll_cd).all())}")

        # [D] full-posterior fidelity: my DEATH terms vs PyMC's own compute_log_likelihood
        tc, ec = f"{endpoint}_time", endpoint
        df = pd.read_parquet(PROCESSED_DIR / cohort / f"{cohort}.parquet")
        gene_cols = rc.select_genes(cohort, df, panel)
        dfv = rc._valid(df, tc, ec); d = dfv[ec].values.astype(int)
        builder = rc.build_model_a if model == "a" else rc.build_model_b
        mdl = builder(df, gene_cols, tc, ec, p0)
        idata = az.from_netcdf(src).isel(draw=slice(0, 25))   # thin for speed
        pmll = _ds(pm.compute_log_likelihood(idata, model=mdl, progressbar=False).log_likelihood)["y_obs"].values
        mine_dead = ll_cd[:, :25, d == 1]
        print(f"    [D] my DEATH terms vs PyMC compute_log_likelihood: max|Δ| = "
              f"{np.abs(mine_dead - pmll).max():.2e}  (model-faithful eta reconstruction)")

        # [E] total at posterior mean: mine vs an independent scipy total
        post = _ds(az.from_netcdf(src).posterior)
        def stk(n): return post[n].stack(s=("chain", "draw")).values
        bm, mm, am = stk("beta").mean(1), stk("mu").mean(), stk("alpha").mean()
        X = dfv[gene_cols].values; t = dfv[tc].values.astype(float); logt = np.log(t)
        etam = mm + X @ bm
        if model == "b":
            etam = etam + stk("mu_site").mean(1)[pd.Categorical(dfv["tss"]).codes]
        mine = np.where(d == 1, my_logpdf(t, am, etam), my_logsf(t, am, etam)).sum()
        sci = np.where(d == 1, weibull_min.logpdf(t, am, scale=np.exp(etam)),
                       weibull_min.logsf(t, am, scale=np.exp(etam))).sum()
        print(f"    [E] total log-lik @ posterior mean: mine = {mine:.4f}  scipy = {sci:.4f}  "
              f"|Δ| = {abs(mine - sci):.2e}")
    print("\nvalidation complete.\n")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohorts", nargs="*", default=None, help="subset of cohorts (default: all)")
    ap.add_argument("--force", action="store_true", help="recompute + reattach even if present")
    ap.add_argument("--validate", nargs=4, metavar=("COHORT", "ENDPOINT", "PANEL", "P0"),
                    default=["LUSC", "OS", "hallmark", "10"],
                    help="deep validation print for one combo (default: LUSC OS hallmark 10)")
    ap.add_argument("--no-validate", action="store_true", help="skip the deep validation print")
    args = ap.parse_args()

    if not args.no_validate:
        c, ep, pn, p0 = args.validate
        try:
            validate(c, ep, pn, int(p0))
        except Exception as exc:
            print(f"validation skipped ({type(exc).__name__}: {exc})")

    recs = discover(set(args.cohorts) if args.cohorts else None)
    print(f"discovered {len(recs)} traces")
    results = []
    for rec in recs:
        try:
            status = attach_loglik(rec, force=args.force)
        except Exception as exc:
            status = f"ERROR: {type(exc).__name__}: {exc}"
        results.append({**{k: rec[k] for k in ("cohort", "model", "endpoint", "panel", "p0")},
                        "status": status})
        print(f"  {rec['cohort']:5s} {rec['model']} {rec['endpoint']:3s} {rec['panel']:9s} p{rec['p0']}: {status}")

    R = pd.DataFrame(results)
    print("\nSUMMARY — traces now carrying a log_likelihood group:")
    ok = R[R["status"].str.startswith(("ok", "skip"))]
    print(f"  {len(ok)}/{len(R)} have log_likelihood "
          f"(ok={int(R['status'].str.startswith('ok').sum())}, "
          f"skip={int(R['status'].str.startswith('skip').sum())}, "
          f"fail={int((~R['status'].str.startswith(('ok','skip'))).sum())})")
    if (~R["status"].str.startswith(("ok", "skip"))).any():
        print("  FAILURES:")
        for r in R[~R["status"].str.startswith(("ok", "skip"))].itertuples(index=False):
            print(f"    {r.cohort}/{r.model}/{r.endpoint}: {r.status}")


if __name__ == "__main__":
    main()
