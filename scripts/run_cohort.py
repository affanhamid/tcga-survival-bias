#!/usr/bin/env python
"""
Pan-cancer site-confounding pipeline (RQ1/RQ2/RQ3).

For each TCGA cohort: download -> build cohort parquet -> hallmark panel ->
fit Model A (naive) and Model B (site-adjusted) with nutpie -> ICC variance
decomposition -> gene hits -> optional leave-one-site-out C-index. Saves
per-cohort traces and appends a row to results/pan_cancer_summary.csv.

Designed to run unattended on a server (e.g. in tmux):

    python scripts/run_cohort.py                       # all viable cohorts, ICC only
    python scripts/run_cohort.py --cohorts LUSC HNSC   # a subset
    python scripts/run_cohort.py --loso                # also run LOSO C-index (slow)
    python scripts/run_cohort.py --force               # recompute even if done

Re-runs are resumable: cohorts already marked "ok" in the summary CSV are
skipped unless --force.
"""
from __future__ import annotations

import argparse
import glob
import logging
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# paths & config
# ----------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
TRACE_DIR = RESULTS_DIR / "traces"
SUMMARY_CSV = RESULTS_DIR / "pan_cancer_summary.csv"

for d in (RAW_DIR, PROCESSED_DIR, TRACE_DIR):
    d.mkdir(parents=True, exist_ok=True)

CDR_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx"
CDR_PATH = RAW_DIR / "tcga_cdr.xlsx"
COUNTS_GLOB = "*.rna_seq.augmented_star_gene_counts.tsv"

# all 33 TCGA project codes; non-viable ones (too few deaths / no data) are skipped
ALL_TCGA = [
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
    "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
    "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
    "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM",
]

ENDPOINTS = ["OS", "OS_time", "DSS", "DSS_time", "PFI", "PFI_time", "DFI", "DFI_time"]

log = logging.getLogger("run_cohort")


# ----------------------------------------------------------------------------
# data acquisition (mirrors notebooks/01_data_prep.ipynb)
# ----------------------------------------------------------------------------
def download_cdr() -> Path:
    if CDR_PATH.exists():
        return CDR_PATH
    import requests
    log.info("downloading TCGA-CDR survival labels")
    r = requests.get(CDR_URL, stream=True)
    r.raise_for_status()
    with open(CDR_PATH, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    return CDR_PATH


def download_manifest(cohort: str) -> Path:
    import requests
    cohort_dir = RAW_DIR / cohort
    cohort_dir.mkdir(parents=True, exist_ok=True)
    manifest = cohort_dir / "manifest.txt"
    if manifest.exists():
        return manifest
    log.info("[%s] downloading GDC manifest", cohort)
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": f"TCGA-{cohort}"}},
            {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
            {"op": "=", "content": {"field": "data_format", "value": "TSV"}},
            {"op": "=", "content": {"field": "cases.samples.sample_type", "value": "Primary Tumor"}},
        ],
    }
    r = requests.post(
        "https://api.gdc.cancer.gov/files",
        headers={"Content-Type": "application/json"},
        json={"filters": filters, "return_type": "manifest", "size": 10000},
    )
    r.raise_for_status()
    manifest.write_text(r.text)
    n = len(r.text.strip().split("\n")) - 1
    log.info("[%s] manifest: %d files", cohort, n)
    return manifest


def download_files(cohort: str) -> None:
    manifest = RAW_DIR / cohort / "manifest.txt"
    out_dir = RAW_DIR / cohort
    file_ids = pd.read_csv(manifest, sep="\t")["id"].tolist()
    missing = [f for f in file_ids if not (out_dir / f).is_dir()]
    if not missing:
        log.info("[%s] all %d files already downloaded", cohort, len(file_ids))
        return
    log.info("[%s] downloading %d/%d missing files", cohort, len(missing), len(file_ids))
    subprocess.run(
        ["gdc-client", "download", "-m", str(manifest), "-d", str(out_dir), "-n", "8"],
        check=True,
    )


def resolve_barcodes(cohort: str) -> pd.DataFrame:
    import requests
    out_path = RAW_DIR / cohort / "barcodes.csv"
    if out_path.exists():
        return pd.read_csv(out_path)
    manifest = RAW_DIR / cohort / "manifest.txt"
    file_ids = pd.read_csv(manifest, sep="\t")["id"].tolist()
    fields = ",".join([
        "file_id", "file_name", "cases.submitter_id", "cases.samples.submitter_id",
        "cases.samples.sample_type",
        "cases.samples.portions.analytes.aliquots.submitter_id",
        "cases.project.project_id",
    ])
    r = requests.post(
        "https://api.gdc.cancer.gov/files",
        headers={"Content-Type": "application/json"},
        json={
            "filters": {"op": "in", "content": {"field": "file_id", "value": file_ids}},
            "fields": fields, "format": "json", "size": len(file_ids) + 100,
        },
    )
    r.raise_for_status()
    rows = []
    for h in r.json()["data"]["hits"]:
        case = h["cases"][0]
        sample = case["samples"][0]
        aliquot = sample["portions"][0]["analytes"][0]["aliquots"][0]
        pb = case["submitter_id"]
        rows.append({
            "file_id": h["file_id"], "file_name": h["file_name"],
            "project_id": case["project"]["project_id"],
            "patient_barcode": pb, "sample_barcode": sample["submitter_id"],
            "aliquot_barcode": aliquot["submitter_id"],
            "sample_type": sample["sample_type"], "tss": pb.split("-")[1],
        })
    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    return df


# ----------------------------------------------------------------------------
# cohort table (mirrors 01_data_prep): dedup -> gene matrix -> filter/transform
# ----------------------------------------------------------------------------
def prepare_cohort(cohort: str) -> Path:
    """Download + build data/processed/{cohort}.parquet. Returns the parquet path."""
    parquet = PROCESSED_DIR / f"{cohort}.parquet"
    if parquet.exists():
        return parquet

    download_cdr()
    download_manifest(cohort)
    download_files(cohort)
    barcodes = resolve_barcodes(cohort)
    assert (barcodes["sample_type"] == "Primary Tumor").all()
    assert (barcodes["project_id"] == f"TCGA-{cohort}").all()

    # dedup aliquots -> one row per patient (alphabetically first aliquot)
    dedup = (barcodes.sort_values("aliquot_barcode")
             .groupby("patient_barcode", as_index=False).first())
    log.info("[%s] %d files -> %d patients (deduped)", cohort, len(barcodes), len(dedup))

    # gene matrix: protein-coding FPKM-UQ, genes x file_id (collect once, concat once)
    cols = {}
    for fid in dedup["file_id"]:
        tsv = next((RAW_DIR / cohort / fid).glob(COUNTS_GLOB))
        g = pd.read_csv(tsv, sep="\t", skiprows=1)
        g = g[g["gene_id"].str.startswith("ENSG")]
        g = g[g["gene_type"] == "protein_coding"]
        cols[fid] = g.set_index("gene_id")["fpkm_uq_unstranded"]
    gene_matrix = pd.DataFrame(cols)
    gene_matrix.to_csv(RAW_DIR / cohort / "gene_matrix.csv")

    # expression filter -> log2 -> per-gene standardise
    n = gene_matrix.shape[1]
    keep = (gene_matrix >= 1).sum(axis=1) >= 0.1 * n
    gene_matrix = gene_matrix.loc[keep]
    logm = np.log2(gene_matrix + 1)
    zscore = logm.sub(logm.mean(axis=1), axis=0).div(logm.std(axis=1), axis=0)
    assert not zscore.isna().any().any(), "NaNs in zscore (zero-variance gene)"
    log.info("[%s] genes after expression filter: %d", cohort, zscore.shape[0])

    # transpose, attach ids, join ALL survival endpoints (NaNs kept by design)
    X = zscore.T
    X.index.name = "file_id"
    meta = dedup.set_index("file_id")[["patient_barcode", "tss"]]
    X = meta.join(X)

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
    final = final[id_cols + ENDPOINTS + gene_cols]
    assert final["patient_barcode"].is_unique
    assert not final[id_cols + gene_cols].isna().any().any()
    final.to_parquet(parquet, index=False)
    log.info("[%s] wrote %s  shape=%s", cohort, parquet.name, final.shape)
    return parquet


# ----------------------------------------------------------------------------
# hallmark gene panel
# ----------------------------------------------------------------------------
_HALLMARK_SYMBOLS = None


def hallmark_symbols() -> set:
    global _HALLMARK_SYMBOLS
    if _HALLMARK_SYMBOLS is None:
        import gseapy as gp
        hm = gp.get_library(name="MSigDB_Hallmark_2020", organism="Human")
        _HALLMARK_SYMBOLS = {g for genes in hm.values() for g in genes}
    return _HALLMARK_SYMBOLS


def hallmark_gene_cols(cohort: str, df: pd.DataFrame) -> list[str]:
    """Map hallmark symbols -> versioned ENSG (via this cohort's STAR annotation),
    intersect with the matrix columns present in df."""
    tsv = glob.glob(str(RAW_DIR / cohort / "*" / COUNTS_GLOB))[0]
    ref = pd.read_csv(tsv, sep="\t", skiprows=1)
    ref = ref[ref["gene_id"].str.startswith("ENSG")]
    sym2ens = dict(zip(ref["gene_name"], ref["gene_id"]))
    hm_ens = {sym2ens[s] for s in hallmark_symbols() if s in sym2ens}
    return [c for c in df.columns if c.startswith("ENSG") and c in hm_ens]


# ----------------------------------------------------------------------------
# models (identical structure to notebooks/03_Modelling.ipynb)
# ----------------------------------------------------------------------------
def _valid(df: pd.DataFrame) -> pd.DataFrame:
    """OS endpoint requires t > 0 and non-missing event/time."""
    return df[(df["OS_time"] > 0) & df["OS"].notna() & df["OS_time"].notna()]


def build_model_a(df, gene_cols, p0=10, slab_scale=2.0, slab_df=4.0):
    import pymc as pm
    import pytensor.tensor as pt
    df = _valid(df)
    X = df[gene_cols].values
    t = df["OS_time"].values
    d = df["OS"].values.astype(int)
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


def build_model_b(df, gene_cols, p0=10, slab_scale=2.0, slab_df=4.0):
    import pymc as pm
    import pytensor.tensor as pt
    df = _valid(df)
    X = df[gene_cols].values
    t = df["OS_time"].values
    d = df["OS"].values.astype(int)
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
def beta_samples(trace):
    return trace.posterior["beta"].stack(s=("chain", "draw")).values  # (p, S)


def gene_hits(trace, eti=89.0):
    """# of genes whose ETI excludes 0 (numpy quantiles; avoids arviz dtype issues)."""
    b = beta_samples(trace)
    a = (100 - eti) / 2
    lo = np.percentile(b, a, axis=1)
    hi = np.percentile(b, 100 - a, axis=1)
    return int(((lo > 0) | (hi < 0)).sum())


def icc_decomposition(trace_b, Xb, eti=89.0):
    post = trace_b.posterior
    b = beta_samples(trace_b)
    alpha = post["alpha"].stack(s=("chain", "draw")).values
    sig = post["sigma_site"].stack(s=("chain", "draw")).values
    var_g = (Xb @ b).var(axis=0)
    var_s = sig**2
    var_r = np.pi**2 / (6 * alpha**2)
    tot = var_s + var_g + var_r
    a = (100 - eti) / 2
    out = {}
    for name, x in (("site", var_s / tot), ("genomic", var_g / tot), ("residual", var_r / tot)):
        out[f"ICC_{name}_mean"] = float(x.mean())
        out[f"ICC_{name}_lo"] = float(np.percentile(x, a))
        out[f"ICC_{name}_hi"] = float(np.percentile(x, 100 - a))
    out["sigma_site_mean"] = float(sig.mean())
    return out


def c_index(times, events, scores):
    """Harrell's C: higher score = predicted longer survival."""
    num = den = 0.0
    for i in np.where(events == 1)[0]:
        later = times > times[i]
        den += later.sum()
        num += (scores[later] > scores[i]).sum() + 0.5 * (scores[later] == scores[i]).sum()
    return num / den if den else np.nan


def loso_cindex(df, gene_cols, p0, draws, tune, chains, target_accept):
    """Leave-one-site-out pooled C-index for Model A and Model B (genomic eta only)."""
    dfl = _valid(df).reset_index(drop=True)
    Xall = dfl[gene_cols].values
    t = dfl["OS_time"].values
    e = dfl["OS"].values.astype(int)
    sites = dfl["tss"].unique()
    eta_a = np.full(len(dfl), np.nan)
    eta_b = np.full(len(dfl), np.nan)
    for k, s in enumerate(sites, 1):
        test = dfl["tss"].values == s
        train = dfl[~test]
        ta = sample(build_model_a(train, gene_cols, p0=p0), draws, tune, chains, target_accept)
        tb = sample(build_model_b(train, gene_cols, p0=p0), draws, tune, chains, target_accept)
        bA = ta.posterior["beta"].mean(("chain", "draw")).values
        bB = tb.posterior["beta"].mean(("chain", "draw")).values
        eta_a[test] = Xall[test] @ bA
        eta_b[test] = Xall[test] @ bB
        log.info("    LOSO site %s (%d/%d)", s, k, len(sites))
    return c_index(t, e, eta_a), c_index(t, e, eta_b)


# ----------------------------------------------------------------------------
# per-cohort driver
# ----------------------------------------------------------------------------
def run_cohort(cohort, args):
    t0 = time.time()
    parquet = prepare_cohort(cohort)
    df = pd.read_parquet(parquet)
    dfl = _valid(df)
    n_deaths = int(dfl["OS"].sum())
    if n_deaths < args.min_events:
        log.warning("[%s] only %d deaths (< %d) — skipping", cohort, n_deaths, args.min_events)
        return {"cohort": cohort, "status": "skip_low_events", "n_deaths": n_deaths}

    gene_cols = hallmark_gene_cols(cohort, df)
    Xb = dfl[gene_cols].values
    log.info("[%s] n=%d deaths=%d sites=%d genes=%d",
             cohort, len(dfl), n_deaths, dfl["tss"].nunique(), len(gene_cols))

    trace_a = sample(build_model_a(df, gene_cols, args.p0), args.draws, args.tune, args.chains, args.target_accept)
    trace_b = sample(build_model_b(df, gene_cols, args.p0), args.draws, args.tune, args.chains, args.target_accept)
    trace_a.to_netcdf(TRACE_DIR / f"{cohort}_A.nc")
    trace_b.to_netcdf(TRACE_DIR / f"{cohort}_B.nc")

    row = {
        "cohort": cohort, "status": "ok",
        "n_patients": len(dfl), "n_deaths": n_deaths,
        "n_sites": int(dfl["tss"].nunique()), "n_genes": len(gene_cols),
        "div_A": int(trace_a.sample_stats.diverging.sum()),
        "div_B": int(trace_b.sample_stats.diverging.sum()),
        "hits_A": gene_hits(trace_a), "hits_B": gene_hits(trace_b),
    }
    row.update(icc_decomposition(trace_b, Xb))

    if args.loso:
        cA, cB = loso_cindex(df, gene_cols, args.p0, args.loso_draws,
                             args.loso_draws, 2, args.target_accept)
        row["cindex_A"] = cA
        row["cindex_B"] = cB
        row["cindex_gap"] = cB - cA

    row["minutes"] = round((time.time() - t0) / 60, 1)
    log.info("[%s] ICC_site=%.3f [%.3f, %.3f]  ICC_genomic=%.3f  (%.1f min)",
             cohort, row["ICC_site_mean"], row["ICC_site_lo"], row["ICC_site_hi"],
             row["ICC_genomic_mean"], row["minutes"])
    return row


def append_summary(row: dict):
    df_row = pd.DataFrame([row])
    if SUMMARY_CSV.exists():
        prev = pd.read_csv(SUMMARY_CSV)
        prev = prev[prev["cohort"] != row["cohort"]]   # replace any old row for this cohort
        df_row = pd.concat([prev, df_row], ignore_index=True)
    df_row.to_csv(SUMMARY_CSV, index=False)


def done_cohorts() -> set:
    if not SUMMARY_CSV.exists():
        return set()
    s = pd.read_csv(SUMMARY_CSV)
    return set(s.loc[s["status"] == "ok", "cohort"])


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohorts", nargs="*", default=ALL_TCGA, help="cohort codes (default: all 33 TCGA)")
    ap.add_argument("--p0", type=int, default=10, help="prior guess for # relevant genes")
    ap.add_argument("--draws", type=int, default=1000)
    ap.add_argument("--tune", type=int, default=1000)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--target-accept", type=float, default=0.99)
    ap.add_argument("--loso", action="store_true", help="also run leave-one-site-out C-index (slow)")
    ap.add_argument("--loso-draws", type=int, default=500, help="draws/tune per LOSO fold fit")
    ap.add_argument("--min-events", type=int, default=50, help="skip cohorts with fewer OS deaths")
    ap.add_argument("--force", action="store_true", help="recompute cohorts already marked ok")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.INFO, stream=sys.stdout,
        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S",
    )

    done = set() if args.force else done_cohorts()
    todo = [c for c in args.cohorts if c not in done]
    log.info("cohorts: %d requested, %d already done, %d to run", len(args.cohorts), len(done), len(todo))

    for cohort in todo:
        try:
            row = run_cohort(cohort, args)
        except Exception as exc:  # isolate failures so one bad cohort doesn't kill the run
            log.exception("[%s] FAILED: %s", cohort, exc)
            row = {"cohort": cohort, "status": f"error: {type(exc).__name__}: {exc}"}
        append_summary(row)

    log.info("all done -> %s", SUMMARY_CSV)


if __name__ == "__main__":
    main()
