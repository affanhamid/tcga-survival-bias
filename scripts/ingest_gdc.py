import requests
import json
import subprocess
import os

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