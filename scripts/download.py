"""TCGA / GDC data acquisition for the site-confounding project.

Shared by notebooks/01_data_prep.ipynb and scripts/run_cohort.py. Paths are
computed from THIS file's location (not the caller's cwd), so it can be imported
from a notebook in notebooks/ or a script in scripts/ and resolve to the same
data/ directory.

Public API:
    download_cdr() -> Path
    download_gdc_manifest(cohort) -> Path
    download_gdc_files(cohort) -> None
    resolve_barcodes(cohort) -> pd.DataFrame
    download_cohort(cohort) -> pd.DataFrame      # cdr + manifest + files + barcodes
"""
from __future__ import annotations

import os
import subprocess
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
COUNTS_GLOB = "*.rna_seq.augmented_star_gene_counts.tsv"
GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"


def download_cdr() -> Path:
    """Download the TCGA-CDR survival-label workbook (Liu et al. 2018). Idempotent."""
    if CDR_PATH.exists():
        print(f"TCGA-CDR already present at {CDR_PATH}, skipping.")
        return CDR_PATH
    print("Downloading TCGA-CDR survival labels...")
    r = requests.get(CDR_URL, stream=True)
    r.raise_for_status()
    with open(CDR_PATH, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"Saved to {CDR_PATH}")
    return CDR_PATH


def download_gdc_manifest(cohort: str) -> Path:
    """Query GDC for primary-tumour STAR-Counts files of TCGA-{cohort}; write manifest. Idempotent."""
    cohort_dir = RAW_DIR / cohort
    cohort_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = cohort_dir / "manifest.txt"
    if manifest_path.exists():
        print(f"Manifest already present at {manifest_path}, skipping.")
        return manifest_path

    print(f"Downloading GDC manifest for {cohort}...")
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
        GDC_FILES_ENDPOINT,
        headers={"Content-Type": "application/json"},
        json={"filters": filters, "return_type": "manifest", "size": 10000},
    )
    r.raise_for_status()
    manifest_path.write_text(r.text)
    print(f"Files in manifest: {len(r.text.strip().splitlines()) - 1}")
    return manifest_path


def download_gdc_files(cohort: str) -> None:
    """gdc-client download every file in the cohort manifest. Skips already-downloaded files."""
    manifest_path = RAW_DIR / cohort / "manifest.txt"
    output_dir = RAW_DIR / cohort
    output_dir.mkdir(parents=True, exist_ok=True)
    file_ids = pd.read_csv(manifest_path, sep="\t")["id"].tolist()
    missing = [f for f in file_ids if not (output_dir / f).is_dir()]
    if not missing:
        print(f"All {len(file_ids)} {cohort} files already downloaded, skipping.")
        return
    print(f"Downloading {cohort} files to {output_dir} ({len(missing)} of {len(file_ids)} missing)...")
    subprocess.run(
        ["gdc-client", "download", "-m", str(manifest_path), "-d", str(output_dir), "-n", "8"],
        check=True,
    )
    print("Done.")


def resolve_barcodes(cohort: str) -> pd.DataFrame:
    """Map each manifest file UUID to patient/sample/aliquot barcode + TSS via the GDC
    /files endpoint. Writes data/raw/{cohort}/barcodes.csv. Idempotent."""
    out_path = RAW_DIR / cohort / "barcodes.csv"
    if out_path.exists():
        print(f"Barcodes already present at {out_path}, skipping.")
        return pd.read_csv(out_path)

    manifest_path = RAW_DIR / cohort / "manifest.txt"
    file_ids = pd.read_csv(manifest_path, sep="\t")["id"].tolist()
    print(f"Resolving barcodes for {len(file_ids)} files...")
    fields = ",".join([
        "file_id", "file_name", "cases.submitter_id", "cases.samples.submitter_id",
        "cases.samples.sample_type",
        "cases.samples.portions.analytes.aliquots.submitter_id",
        "cases.project.project_id",
    ])
    r = requests.post(
        GDC_FILES_ENDPOINT,
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
    print(f"Wrote {len(df)} rows to {out_path}")
    return df


def download_cohort(cohort: str) -> pd.DataFrame:
    """Convenience: CDR + manifest + files + barcodes for one cohort. Returns barcodes."""
    download_cdr()
    download_gdc_manifest(cohort)
    download_gdc_files(cohort)
    return resolve_barcodes(cohort)
