#!/usr/bin/env python
"""
Pan-cancer site-confounding pipeline (RQ1/RQ2/RQ3).

Per TCGA cohort: download (via scripts/download.py) -> build cohort parquet ->
for each (endpoint, gene panel): fit Model A (naive) and Model B (site-adjusted)
with nutpie -> ICC variance decomposition -> gene hits -> optional
leave-one-site-out C-index, plus an optional ComBat sensitivity analysis.

Everything for a cohort is written under data/processed/{cohort}/ (parquet,
traces/, results.csv) so that directory alone is enough to resume/continue;
results/pan_cancer_summary.csv is a derived cross-cohort aggregate.

    python scripts/run_cohort.py                                 # PFI+OS, hallmark, ICC only
    python scripts/run_cohort.py --cohorts LUSC HNSC
    python scripts/run_cohort.py --panels hallmark topvar:500 topvar:2000   # robustness
    python scripts/run_cohort.py --combat                        # ComBat sensitivity
    python scripts/run_cohort.py --loso                          # leave-one-site-out C-index (slow)
    python scripts/run_cohort.py --no-traces                     # results.csv only (save disk)
    python scripts/run_cohort.py --force                         # recompute even if done

Resumable: (cohort, endpoint, panel) combos already "ok" in the per-cohort
results.csv files are skipped unless --force.
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# shared data-acquisition module, imported from this script's own directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
import download
from download import (PROJECT_ROOT, RAW_DIR, PROCESSED_DIR, CDR_PATH,
                      EXPR_LOG2_THRESHOLD)

# ----------------------------------------------------------------------------
# paths & config (data acquisition lives in download.py)
# ----------------------------------------------------------------------------
RESULTS_DIR = PROJECT_ROOT / "results"
SUMMARY_CSV = RESULTS_DIR / "pan_cancer_summary.csv"   # cross-cohort aggregate (derived)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

ALL_TCGA = [
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
    "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
    "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
    "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM",
]
ENDPOINTS_ALL = ["OS", "OS_time", "DSS", "DSS_time", "PFI", "PFI_time", "DFI", "DFI_time"]
KEY = ["cohort", "endpoint", "panel", "p0"]

log = logging.getLogger("run_cohort")


def cohort_dir(cohort: str) -> Path:
    """Per-cohort output dir: data/processed/{cohort}/ (with a traces/ subdir).
    Holds the cohort parquet, model traces, and that cohort's results.csv."""
    d = PROCESSED_DIR / cohort
    (d / "traces").mkdir(parents=True, exist_ok=True)
    return d


# ----------------------------------------------------------------------------
# cohort table (mirrors 01_data_prep): xena matrix -> filter/transform -> labels
# ----------------------------------------------------------------------------
def _patient(bc: str) -> str:
    return "-".join(bc.split("-")[:3])          # TCGA-60-2698-01A -> TCGA-60-2698


def _tss(bc: str) -> str:
    return bc.split("-")[1]                       # TCGA-60-2698-01A -> 60


def _sample_type(bc: str) -> str:
    parts = bc.split("-")
    return parts[3][:2] if len(parts) >= 4 else ""  # TCGA-60-2698-01A -> 01


def prepare_cohort(cohort: str) -> Path:
    """Download + build data/processed/{cohort}/{cohort}.parquet. Returns the parquet path.

    Reads the Xena STAR FPKM-UQ matrix (log2(fpkm_uq+1)), keeps protein-coding
    genes and primary-tumour samples (sample-type code '01' — matching the old
    GDC 'Primary Tumor' filter; note this excludes LAML(03)/SKCM(06) as before),
    dedups to one sample per patient, applies the expression filter + z-score,
    then joins TCGA-CDR survival labels."""
    parquet = cohort_dir(cohort) / f"{cohort}.parquet"
    if parquet.exists():
        return parquet

    download.download_cdr()
    gz = download.download_xena_matrix(cohort)
    annot = download.gene_annotation()
    pc = set(annot.loc[annot["gene_type"] == "protein_coding", "gene_id"])

    mat = pd.read_csv(gz, sep="\t", index_col=0)            # genes x samples, log2(fpkm_uq+1)
    mat = mat[mat.index.isin(pc)]                           # protein-coding genes only

    prim = [c for c in mat.columns if _sample_type(c) == "01"]  # primary solid tumour
    mat = mat[prim]
    seen, keep = set(), []                                  # dedup: one sample per patient
    for c in sorted(mat.columns):
        p = _patient(c)
        if p not in seen:
            seen.add(p)
            keep.append(c)
    mat = mat[keep]
    log.info("[%s] %d primary-tumour samples -> %d patients (deduped)", cohort, len(prim), len(keep))

    n = mat.shape[1]                                        # expression filter (>=1 in >=10%)
    mat = mat.loc[(mat >= EXPR_LOG2_THRESHOLD).sum(axis=1) >= 0.1 * n]
    zscore = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0)
    assert not zscore.isna().any().any(), "NaNs in zscore (zero-variance gene)"
    log.info("[%s] genes after expression filter: %d", cohort, zscore.shape[0])

    X = zscore.T
    X.insert(0, "tss", [_tss(c) for c in X.index])
    X.insert(0, "patient_barcode", [_patient(c) for c in X.index])
    X = X.reset_index(drop=True)

    cdr = pd.read_excel(CDR_PATH, sheet_name="TCGA-CDR")
    cdr = cdr[cdr["type"] == cohort].copy()
    labels = (cdr[["bcr_patient_barcode", "OS", "OS.time", "DSS", "DSS.time",
                   "PFI", "PFI.time", "DFI", "DFI.time"]]
              .rename(columns={"bcr_patient_barcode": "patient_barcode",
                               "OS.time": "OS_time", "DSS.time": "DSS_time",
                               "PFI.time": "PFI_time", "DFI.time": "DFI_time"}))
    final = X.merge(labels, on="patient_barcode", how="inner")

    id_cols = ["patient_barcode", "tss"]
    gene_cols = [c for c in final.columns if c.startswith("ENSG")]
    final = final[id_cols + ENDPOINTS_ALL + gene_cols]
    assert final["patient_barcode"].is_unique
    assert not final[id_cols + gene_cols].isna().any().any()
    final.to_parquet(parquet, index=False)
    log.info("[%s] wrote %s shape=%s", cohort, parquet, final.shape)
    return parquet


# ----------------------------------------------------------------------------
# gene panels
# ----------------------------------------------------------------------------
_HALLMARK = None


def hallmark_symbols() -> set:
    global _HALLMARK
    if _HALLMARK is None:
        import gseapy as gp
        hm = gp.get_library(name="MSigDB_Hallmark_2020", organism="Human")
        _HALLMARK = {g for genes in hm.values() for g in genes}
    return _HALLMARK


def select_genes(cohort: str, df: pd.DataFrame, panel: str) -> list[str]:
    """panel = 'hallmark' or 'topvar:N'."""
    ensg = [c for c in df.columns if c.startswith("ENSG")]
    if panel == "hallmark":
        annot = download.gene_annotation()
        sym2ens = dict(zip(annot["gene_name"], annot["gene_id"]))
        hm = {sym2ens[s] for s in hallmark_symbols() if s in sym2ens}
        return [c for c in ensg if c in hm]
    if panel.startswith("topvar:"):
        n = int(panel.split(":")[1])
        return df[ensg].var().nlargest(n).index.tolist()
    raise ValueError(f"unknown panel: {panel}")


# ----------------------------------------------------------------------------
# models (endpoint-parametrized; identical structure to notebooks/03_Modelling)
# ----------------------------------------------------------------------------
def _valid(df: pd.DataFrame, time_col: str, event_col: str) -> pd.DataFrame:
    return df[(df[time_col] > 0) & df[event_col].notna() & df[time_col].notna()]


def build_model_a(df, gene_cols, time_col, event_col, p0=10, slab_scale=2.0, slab_df=4.0):
    import pymc as pm
    import pytensor.tensor as pt
    df = _valid(df, time_col, event_col)
    X = df[gene_cols].values
    t = df[time_col].values
    d = df[event_col].values.astype(int)
    n, p = X.shape
    obs, cens = d == 1, d == 0
    Xo, Xc, to, tc = X[obs], X[cens], t[obs], t[cens]
    sigma_tau = (p0 / (p - p0)) * (np.std(np.log1p(t)) / np.sqrt(n))
    with pm.Model() as model:
        mu = pm.Normal("mu", mu=6, sigma=2)
        alpha = pm.HalfNormal("alpha", sigma=1)
        tau = pm.HalfStudentT("tau", nu=3, sigma=sigma_tau)
        lam = pm.HalfCauchy("lam", beta=1, shape=p)
        c2 = pm.InverseGamma("c2", alpha=slab_df / 2, beta=(slab_df / 2) * slab_scale**2)
        lam_t = pt.sqrt(c2 * lam**2 / (c2 + tau**2 * lam**2))
        z = pm.Normal("z", mu=0, sigma=1, shape=p)
        beta = pm.Deterministic("beta", z * tau * lam_t)
        eta_o = mu + pm.math.dot(Xo, beta)
        eta_c = mu + pm.math.dot(Xc, beta)
        pm.Weibull("y_obs", alpha=alpha, beta=pm.math.exp(eta_o), observed=to)
        pm.Potential("cens_ll", (-pm.math.exp(alpha * (np.log(tc) - eta_c))).sum())
    return model


def build_model_b(df, gene_cols, time_col, event_col, p0=10, slab_scale=2.0, slab_df=4.0):
    import pymc as pm
    import pytensor.tensor as pt
    df = _valid(df, time_col, event_col)
    X = df[gene_cols].values
    t = df[time_col].values
    d = df[event_col].values.astype(int)
    site = pd.Categorical(df["tss"]).codes
    n_sites = int(site.max() + 1)
    n, p = X.shape
    obs, cens = d == 1, d == 0
    Xo, Xc, to, tc = X[obs], X[cens], t[obs], t[cens]
    so, sc = site[obs], site[cens]
    sigma_tau = (p0 / (p - p0)) * (np.std(np.log1p(t)) / np.sqrt(n))
    with pm.Model() as model:
        mu = pm.Normal("mu", mu=6, sigma=2)
        alpha = pm.HalfNormal("alpha", sigma=1)
        sigma_site = pm.HalfNormal("sigma_site", sigma=1)
        z_site = pm.Normal("z_site", mu=0, sigma=1, shape=n_sites)
        mu_site = pm.Deterministic("mu_site", z_site * sigma_site)
        tau = pm.HalfStudentT("tau", nu=3, sigma=sigma_tau)
        lam = pm.HalfCauchy("lam", beta=1, shape=p)
        c2 = pm.InverseGamma("c2", alpha=slab_df / 2, beta=(slab_df / 2) * slab_scale**2)
        lam_t = pt.sqrt(c2 * lam**2 / (c2 + tau**2 * lam**2))
        z = pm.Normal("z", mu=0, sigma=1, shape=p)
        beta = pm.Deterministic("beta", z * tau * lam_t)
        eta_o = mu + mu_site[so] + pm.math.dot(Xo, beta)
        eta_c = mu + mu_site[sc] + pm.math.dot(Xc, beta)
        pm.Weibull("y_obs", alpha=alpha, beta=pm.math.exp(eta_o), observed=to)
        pm.Potential("cens_ll", (-pm.math.exp(alpha * (np.log(tc) - eta_c))).sum())
    return model


def sample(model, draws, tune, chains, target_accept):
    import pymc as pm
    with model:
        return pm.sample(draws=draws, tune=tune, chains=chains,
                         target_accept=target_accept, nuts_sampler="nutpie",
                         return_inferencedata=True, progressbar=False)


# ----------------------------------------------------------------------------
# analysis helpers
# ----------------------------------------------------------------------------
def _stack(trace, name):
    return trace.posterior[name].stack(s=("chain", "draw")).values


def gene_hits(trace, eti=89.0):
    b = _stack(trace, "beta")
    a = (100 - eti) / 2
    return int(((np.percentile(b, a, axis=1) > 0) | (np.percentile(b, 100 - a, axis=1) < 0)).sum())


def icc_decomposition(trace_b, Xb, eti=89.0):
    b = _stack(trace_b, "beta")
    alpha = _stack(trace_b, "alpha")
    sig = _stack(trace_b, "sigma_site")
    var_g = (Xb @ b).var(axis=0)
    var_s = sig**2
    var_r = np.pi**2 / (6 * alpha**2)
    tot = var_s + var_g + var_r
    a = (100 - eti) / 2
    out = {"sigma_site_mean": float(sig.mean())}
    for name, x in (("site", var_s / tot), ("genomic", var_g / tot), ("residual", var_r / tot)):
        out[f"ICC_{name}_mean"] = float(x.mean())
        out[f"ICC_{name}_lo"] = float(np.percentile(x, a))
        out[f"ICC_{name}_hi"] = float(np.percentile(x, 100 - a))
    return out


def c_index(times, events, scores):
    num = den = 0.0
    for i in np.where(events == 1)[0]:
        later = times > times[i]
        den += later.sum()
        num += (scores[later] > scores[i]).sum() + 0.5 * (scores[later] == scores[i]).sum()
    return num / den if den else np.nan


def loso_cindex(df, gene_cols, time_col, event_col, p0, draws, chains, target_accept):
    dfl = _valid(df, time_col, event_col).reset_index(drop=True)
    Xall = dfl[gene_cols].values
    t = dfl[time_col].values
    e = dfl[event_col].values.astype(int)
    eta_a = np.full(len(dfl), np.nan)
    eta_b = np.full(len(dfl), np.nan)
    for s in dfl["tss"].unique():
        test = dfl["tss"].values == s
        train = dfl[~test]
        ta = sample(build_model_a(train, gene_cols, time_col, event_col, p0), draws, draws, chains, target_accept)
        tb = sample(build_model_b(train, gene_cols, time_col, event_col, p0), draws, draws, chains, target_accept)
        eta_a[test] = Xall[test] @ ta.posterior["beta"].mean(("chain", "draw")).values
        eta_b[test] = Xall[test] @ tb.posterior["beta"].mean(("chain", "draw")).values
    return c_index(t, e, eta_a), c_index(t, e, eta_b)


def combat_correct(df, gene_cols, batch_col="tss"):
    """ComBat-correct the gene block by TSS batch. Drops singleton-site patients
    (ComBat cannot estimate a batch effect from one sample)."""
    from inmoose.pycombat import pycombat_norm
    counts = df[batch_col].value_counts()
    sub = df[df[batch_col].isin(counts[counts >= 2].index)].copy()
    sub[gene_cols] = pycombat_norm(sub[gene_cols].T.values, sub[batch_col].values).T
    return sub


# ----------------------------------------------------------------------------
# one (endpoint, panel) run
# ----------------------------------------------------------------------------
def _safe(name: str) -> str:
    return name.replace(":", "").replace("+", "_")


def analyze(cohort, df, endpoint, panel, gene_cols, args, combat=False):
    tc, ec = f"{endpoint}_time", endpoint
    if combat:
        df = combat_correct(df, gene_cols)
    dfv = _valid(df, tc, ec)
    Xb = dfv[gene_cols].values

    ta = sample(build_model_a(df, gene_cols, tc, ec, args.p0), args.draws, args.tune, args.chains, args.target_accept)
    tb = sample(build_model_b(df, gene_cols, tc, ec, args.p0), args.draws, args.tune, args.chains, args.target_accept)
    if not args.no_traces:
        tdir = cohort_dir(cohort) / "traces"
        tag = f"{endpoint}_{_safe(panel)}_p{args.p0}"
        ta.to_netcdf(tdir / f"model_a_{tag}.nc")
        tb.to_netcdf(tdir / f"model_b_{tag}.nc")

    row = {
        "cohort": cohort, "endpoint": endpoint, "panel": panel, "p0": args.p0, "status": "ok",
        "n_patients": len(dfv), "n_events": int(dfv[ec].sum()),
        "n_sites": int(dfv["tss"].nunique()), "n_genes": len(gene_cols),
        "div_A": int(ta.sample_stats.diverging.sum()),
        "div_B": int(tb.sample_stats.diverging.sum()),
        "hits_A": gene_hits(ta), "hits_B": gene_hits(tb),
    }
    row.update(icc_decomposition(tb, Xb))
    if args.loso:
        cA, cB = loso_cindex(df, gene_cols, tc, ec, args.p0, args.loso_draws, 2, args.target_accept)
        row["cindex_A"], row["cindex_B"], row["cindex_gap"] = cA, cB, cB - cA
    return row


def run_cohort(cohort, args, done):
    prepare_cohort(cohort)
    df = pd.read_parquet(cohort_dir(cohort) / f"{cohort}.parquet")
    rows = []
    for endpoint in args.endpoints:
        tc, ec = f"{endpoint}_time", endpoint
        n_ev = int(_valid(df, tc, ec)[ec].sum())
        if n_ev < args.min_events:
            log.warning("[%s/%s] only %d events (< %d) — skipping endpoint", cohort, endpoint, n_ev, args.min_events)
            append_cohort_result(cohort, {"cohort": cohort, "endpoint": endpoint, "panel": "-",
                                          "p0": args.p0, "status": f"skip_low_events({n_ev})"})
            continue
        panels = list(args.panels) + (["hallmark+combat"] if args.combat else [])
        for panel in panels:
            if (cohort, endpoint, panel, args.p0) in done:
                log.info("[%s/%s/%s] already done — skipping", cohort, endpoint, panel)
                continue
            t0 = time.time()
            try:
                gene_cols = select_genes(cohort, df, "hallmark" if panel == "hallmark+combat" else panel)
                row = analyze(cohort, df, endpoint, panel, gene_cols, args,
                              combat=(panel == "hallmark+combat"))
                row["minutes"] = round((time.time() - t0) / 60, 1)
                log.info("[%s/%s/%s] ICC_site=%.3f [%.3f,%.3f] genomic=%.3f hits=%d (%.1f min)",
                         cohort, endpoint, panel, row["ICC_site_mean"], row["ICC_site_lo"],
                         row["ICC_site_hi"], row["ICC_genomic_mean"], row["hits_B"], row["minutes"])
            except Exception as exc:
                log.exception("[%s/%s/%s] FAILED: %s", cohort, endpoint, panel, exc)
                row = {"cohort": cohort, "endpoint": endpoint, "panel": panel, "p0": args.p0,
                       "status": f"error: {type(exc).__name__}: {exc}"}
            append_cohort_result(cohort, row)
            rows.append(row)
    return rows


# ----------------------------------------------------------------------------
# bookkeeping — per-cohort results.csv is the source of truth; global is derived
# ----------------------------------------------------------------------------
def append_cohort_result(cohort: str, row: dict):
    path = cohort_dir(cohort) / "results.csv"
    new = pd.DataFrame([row])
    if path.exists():
        prev = pd.read_csv(path)
        for k in KEY:                       # tolerate legacy files written before a key column existed
            if k not in prev.columns:
                prev[k] = ""
        mask = ~(prev[KEY].astype(str).agg("|".join, axis=1)
                 == "|".join(str(row.get(k, "")) for k in KEY))
        new = pd.concat([prev[mask], new], ignore_index=True)
    new.to_csv(path, index=False)
    rebuild_global_summary()


def rebuild_global_summary():
    """results/pan_cancer_summary.csv = concat of every data/processed/*/results.csv.
    Purely derived — safe to delete and regenerate from data/processed/ alone."""
    parts = [pd.read_csv(p) for p in sorted(PROCESSED_DIR.glob("*/results.csv"))]
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(SUMMARY_CSV, index=False)


def done_combos() -> set:
    """Resume set read from per-cohort results.csv under data/processed/."""
    done = set()
    for p in PROCESSED_DIR.glob("*/results.csv"):
        s = pd.read_csv(p)
        s = s[s["status"] == "ok"]
        if "p0" not in s.columns:        # legacy rows predate the p0 key — recompute them under it
            continue
        for _, r in s.iterrows():
            done.add((r["cohort"], r["endpoint"], r["panel"], int(r["p0"])))
    return done


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohorts", nargs="*", default=ALL_TCGA)
    ap.add_argument("--endpoints", nargs="*", default=["PFI", "OS"],
                    help="survival endpoints; first is primary (default: PFI OS)")
    ap.add_argument("--panels", nargs="*", default=["hallmark"],
                    help="gene panels: 'hallmark' and/or 'topvar:N' (e.g. topvar:500 topvar:2000)")
    ap.add_argument("--combat", action="store_true", help="add a ComBat (TSS-corrected) sensitivity run per endpoint")
    ap.add_argument("--p0", type=int, default=10)
    ap.add_argument("--draws", type=int, default=1000)
    ap.add_argument("--tune", type=int, default=1000)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--target-accept", type=float, default=0.99)
    ap.add_argument("--loso", action="store_true", help="also run leave-one-site-out C-index (slow)")
    ap.add_argument("--loso-draws", type=int, default=500)
    ap.add_argument("--no-traces", action="store_true",
                    help="skip saving .nc traces (results.csv still written; saves disk)")
    ap.add_argument("--min-events", type=int, default=50)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")

    done = set() if args.force else done_combos()
    log.info("cohorts=%d endpoints=%s panels=%s p0=%d combat=%s loso=%s (%d combos already done)",
             len(args.cohorts), args.endpoints, args.panels, args.p0, args.combat, args.loso, len(done))

    for cohort in args.cohorts:
        try:
            run_cohort(cohort, args, done)
        except Exception as exc:
            log.exception("[%s] cohort-level FAILED: %s", cohort, exc)
            append_cohort_result(cohort, {"cohort": cohort, "endpoint": "-", "panel": "-",
                                          "p0": args.p0, "status": f"error: {type(exc).__name__}: {exc}"})

    rebuild_global_summary()
    log.info("all done -> %s (per-cohort source of truth in data/processed/<cohort>/results.csv)", SUMMARY_CSV)


if __name__ == "__main__":
    main()
