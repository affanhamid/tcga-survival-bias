import requests
import json
import subprocess
import os
import pandas as pd

DATA_DIR = "data/raw"

def download_gdc_manifest(cancer_type: str) -> str:
    """
    Download GDC data for a given cancer type
    """
    print(f"Downloading GDC manifest for {cancer_type}...")

    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": f"TCGA-{cancer_type}"}},
            {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
            {"op": "=", "content": {"field": "data_format", "value": "TSV"}},
            {"op": "=", "content": {"field": "cases.samples.sample_type", "value": "Primary Tumor"}}
        ]
    }

    json_params = {
            "filters": filters,
            "return_type": "manifest",
            "size": 10000
        }


    r = requests.post(
        "https://api.gdc.cancer.gov/files",
        headers={"Content-Type": "application/json"},
        json=json_params)

    manifest_path = f"{DATA_DIR}/{cancer_type}_manifest.txt"
    with open(manifest_path, "w") as f:
        f.write(r.text)

    lines = r.text.strip().split("\n")
    print(f"Files in manifest: {len(lines) - 1}")
    print(f"First file: {lines[0]}")
    print(f"Second file: {lines[1]}")
    return manifest_path



def download_gdc_files(cancer_type: str) -> None:
    manifest_path = f"{DATA_DIR}/{cancer_type}_manifest.txt"
    output_dir = f"{DATA_DIR}/{cancer_type}"
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Downloading {cancer_type} files to {output_dir}...")
    subprocess.run([
        "gdc-client", "download",
        "-m", manifest_path,
        "-d", output_dir,
        "-n", "8"
    ], check=True)
    print("Done.")


def resolve_barcodes(cancer_type: str) -> str:
    """
    For each file in the cancer type's manifest, query GDC for the
    associated patient/sample/aliquot barcodes, sample type, and TSS.
    Writes data/raw/{cancer_type}_barcodes.csv.
    """
    manifest_path = f"{DATA_DIR}/{cancer_type}_manifest.txt"
    print(f"Reading manifest from {manifest_path}...")
    manifest = pd.read_csv(manifest_path, sep="\t")
    file_ids = manifest["id"].tolist()
    print(f"Resolving barcodes for {len(file_ids)} files...")

    fields = ",".join([
        "file_id",
        "file_name",
        "cases.submitter_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
        "cases.samples.portions.analytes.aliquots.submitter_id",
        "cases.project.project_id",
    ])

    json_params = {
        "filters": {
            "op": "in",
            "content": {"field": "file_id", "value": file_ids},
        },
        "fields": fields,
        "format": "json",
        "size": len(file_ids) + 100,
    }

    r = requests.post(
        "https://api.gdc.cancer.gov/files",
        headers={"Content-Type": "application/json"},
        json=json_params,
    )
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    print(f"Received {len(hits)} hits from GDC.")

    rows = []
    for h in hits:
        case = h["cases"][0]
        sample = case["samples"][0]
        aliquot = sample["portions"][0]["analytes"][0]["aliquots"][0]
        patient_barcode = case["submitter_id"]
        rows.append({
            "file_id": h["file_id"],
            "file_name": h["file_name"],
            "project_id": case["project"]["project_id"],
            "patient_barcode": patient_barcode,
            "sample_barcode": sample["submitter_id"],
            "aliquot_barcode": aliquot["submitter_id"],
            "sample_type": sample["sample_type"],
            "tss": patient_barcode.split("-")[1],
        })

    df = pd.DataFrame(rows)
    out_path = f"{DATA_DIR}/{cancer_type}_barcodes.csv"
    df.to_csv(out_path, index=False)
    print(f"Wrote {len(df)} rows to {out_path}")
    return out_path