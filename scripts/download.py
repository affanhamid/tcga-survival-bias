"""TCGA RNA-seq acquisition via the UCSC Xena GDC hub.

Shared by notebooks/01_data_prep.ipynb and scripts/run_cohort.py. Paths are
computed from THIS file's location (not the caller's cwd).

Why Xena instead of gdc-client: the GDC harmonized STAR-Counts data is mirrored
by the UCSC Xena GDC hub as a single pre-assembled matrix per project
(`TCGA-{COHORT}.star_fpkm-uq.tsv.gz`) — one ~150 MB download replaces ~500
per-file gdc-client transfers (days -> minutes). The matrix values are
log2(fpkm_uq_unstranded + 1); we verified them bit-identical (max|diff|=1.8e-15,
corr=1.0 over all 60,660 genes) to the GDC `fpkm_uq_unstranded` column we used
before, so the switch is lossless and methodologically consistent.

Matrix layout: rows = versioned Ensembl IDs (gencode v36, e.g. ENSG...15),
columns = sample barcodes (e.g. TCGA-60-2698-01A). TSS is barcode field 2;
sample-type code is field 4 (the '01' in '01A'); patient barcode is the first
three fields. Xena carries no gene_type, so protein-coding filtering and
symbol->Ensembl mapping use resources/gene_annotation.tsv (committed; identical
across all cohorts since the STAR gene model is fixed).

Public API:
    download_cdr() -> Path
    download_xena_matrix(cohort) -> Path
    gene_annotation() -> pd.DataFrame      # gene_id, gene_name, gene_type
    download_cohort(cohort) -> Path        # cdr + xena matrix; returns matrix path
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import requests

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RAW_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

CDR_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx"
CDR_PATH = RAW_DIR / "tcga_cdr.xlsx"

XENA_BASE = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
GENE_ANNOT_PATH = PROJECT_ROOT / "resources" / "gene_annotation.tsv"

# Xena matrix value = log2(fpkm_uq + 1). On that scale the old raw-FPKM-UQ
# expression filter "fpkm_uq >= 1" is identical, since log2(x+1) >= 1 <=> x >= 1.
EXPR_LOG2_THRESHOLD = 1.0


def download_cdr() -> Path:
    """Download the TCGA-CDR survival-label workbook (Liu et al. 2018). Idempotent."""
    if CDR_PATH.exists():
        print(f"TCGA-CDR already present at {CDR_PATH}, skipping.")
        return CDR_PATH
    print("Downloading TCGA-CDR survival labels...")
    r = requests.get(CDR_URL, stream=True)
    r.raise_for_status()
    tmp = CDR_PATH.with_suffix(CDR_PATH.suffix + ".part")
    with open(tmp, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    tmp.rename(CDR_PATH)  # atomic: a partial download never looks complete
    print(f"Saved to {CDR_PATH}")
    return CDR_PATH


def download_xena_matrix(cohort: str) -> Path:
    """Download TCGA-{cohort}.star_fpkm-uq.tsv.gz from the Xena GDC hub. Idempotent."""
    out = RAW_DIR / cohort / "star_fpkm-uq.tsv.gz"
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        print(f"Xena matrix already present at {out}, skipping.")
        return out
    url = f"{XENA_BASE}/TCGA-{cohort}.star_fpkm-uq.tsv.gz"
    print(f"Downloading Xena matrix for {cohort} from {url} ...")
    r = requests.get(url, stream=True)
    r.raise_for_status()
    tmp = out.with_suffix(out.suffix + ".part")
    n = 0
    with open(tmp, "wb") as f:
        for chunk in r.iter_content(chunk_size=1 << 20):
            f.write(chunk)
            n += len(chunk)
    tmp.rename(out)  # atomic rename so an interrupted download isn't mistaken for complete
    print(f"Saved {n / 1e6:.0f} MB to {out}")
    return out


_ANNOT = None


def gene_annotation() -> pd.DataFrame:
    """gencode v36 gene_id/gene_name/gene_type (cached). Identical across all cohorts."""
    global _ANNOT
    if _ANNOT is None:
        if not GENE_ANNOT_PATH.exists():
            raise FileNotFoundError(
                f"{GENE_ANNOT_PATH} missing — expected committed gencode v36 annotation "
                "(gene_id, gene_name, gene_type) built once from any STAR-Counts TSV."
            )
        _ANNOT = pd.read_csv(GENE_ANNOT_PATH, sep="\t")
    return _ANNOT


def download_cohort(cohort: str) -> Path:
    """Convenience: CDR + Xena matrix for one cohort. Returns the matrix path."""
    download_cdr()
    return download_xena_matrix(cohort)
