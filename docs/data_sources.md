# Data sources

All raw inputs for the pipeline, where they come from, and how they're processed.
Acquisition code lives in `scripts/download.py` (shared by `notebooks/01_data_prep.ipynb`
and `scripts/run_cohort.py`).

## 1. RNA-seq expression — UCSC Xena GDC hub

- **What:** one pre-assembled STAR-Counts FPKM-UQ matrix per TCGA project.
- **URL:** `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{COHORT}.star_fpkm-uq.tsv.gz`
  (e.g. `TCGA-LUSC.star_fpkm-uq.tsv.gz`, ~150 MB). All 33 TCGA cohorts available.
- **Layout:** rows = versioned Ensembl gene IDs (GENCODE v36, e.g. `ENSG00000000003.15`),
  columns = TCGA sample barcodes (e.g. `TCGA-60-2698-01A`).
- **Units:** `log2(fpkm_uq_unstranded + 1)`.
- **Provenance:** the Xena GDC hub mirrors the GDC harmonized STAR-Counts data. It is the
  same quantification we previously pulled per-file with `gdc-client`, just pre-assembled —
  one download replaces ~500 per-file transfers (days → minutes).
- **Local path:** `data/raw/{COHORT}/star_fpkm-uq.tsv.gz` (gitignored).

### Why we switched from gdc-client
Per-file `gdc-client` downloads ran at ~110–480 KB/s and were IP-throttled (running
multiple cohorts concurrently made aggregate throughput *worse*, not better), so a full
cohort took hours. The Xena matrix is one HTTP download per cohort.

### Verified equivalence (lossless switch)
We confirmed the Xena values are bit-identical to the GDC `fpkm_uq_unstranded` column we
used before: for a shared LUSC sample, `xena == log2(fpkm_uq + 1)` with
**max\|diff\| = 1.8e-15 and correlation = 1.0 across all 60,660 genes**. So the switch does
not change the data, only how it is fetched.

A happy consequence: the old expression filter `fpkm_uq >= 1` is unchanged on the log2
scale, since `log2(x + 1) >= 1 ⟺ x >= 1`. The pipeline therefore drops the explicit log2
step (Xena is already log-transformed) and keeps the same threshold.

## 2. Survival labels — TCGA-CDR (Liu et al. 2018)

- **What:** the TCGA Pan-Cancer Clinical Data Resource — curated OS / DSS / PFI / DFI
  endpoints (event + time) per patient.
- **URL:** `https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx`
  (supplementary workbook `mmc1.xlsx`), sheet **`TCGA-CDR`**.
- **Join key:** `bcr_patient_barcode` ↔ patient barcode (`TCGA-XX-XXXX`).
- **Local path:** `data/raw/tcga_cdr.xlsx` (gitignored).
- **Citation:** Liu J, et al. *An Integrated TCGA Pan-Cancer Clinical Data Resource to
  Drive High-Quality Survival Outcome Analytics.* Cell. 2018;173(2):400-416.

## 3. Gene annotation — GENCODE v36 (committed)

- **What:** `gene_id, gene_name, gene_type` for all 60,660 genes in the STAR gene model.
- **Path:** `resources/gene_annotation.tsv` (**committed to the repo**, ~2.5 MB).
- **Why committed:** the Xena matrix carries no `gene_type`, so we need this map to (a) filter
  to protein-coding genes and (b) translate Hallmark gene symbols → Ensembl IDs. The STAR
  gene model is fixed, so this annotation is identical across every cohort — built once from
  any single STAR-Counts TSV (`gene_id`, `gene_name`, `gene_type` columns).
- 60,660 genes total; **19,962 protein-coding**.

## 4. Gene panel — MSigDB Hallmark

- **What:** the 50 MSigDB Hallmark gene sets, used as the pre-specified (a priori) gene panel.
- **Source:** fetched at runtime via `gseapy.get_library("MSigDB_Hallmark_2020", organism="Human")`
  (Enrichr library), cached by gseapy. Hallmark gene symbols are mapped to Ensembl IDs via the
  gene annotation in §3.
- Used as the primary panel; `topvar:N` (top-N highest-variance genes) is a robustness alternative.

## Barcode parsing

TCGA barcodes (e.g. `TCGA-60-2698-01A`) encode everything we need for cohort assembly:

| field | example | meaning |
|---|---|---|
| 1–3 | `TCGA-60-2698` | patient barcode — join key to TCGA-CDR |
| 2 | `60` | **TSS** (tissue source site) — the institution; the site variable in the model |
| 4 (first 2 chars) | `01` | sample-type code: `01` = Primary Solid Tumor |

We keep sample-type `01` (matching the old GDC "Primary Tumor" filter — note this excludes
LAML `03` and SKCM `06`, as before) and deduplicate to one sample per patient (first sample
barcode alphabetically).

## From raw to model-ready (`prepare_cohort` / `01_data_prep.ipynb`)

1. Download the Xena matrix (§1) and TCGA-CDR (§2).
2. Filter genes to protein-coding via the annotation (§3).
3. Keep primary-tumour samples (barcode sample-type `01`); dedup to one per patient.
4. Expression filter: keep genes with value `>= 1` (≡ `fpkm_uq >= 1`) in `>= 10%` of patients.
5. Standardise each gene to mean 0 / std 1 within the cohort (z-score). No log step — Xena is
   already `log2`.
6. Join TCGA-CDR survival labels; write `data/processed/{COHORT}/{COHORT}.parquet`
   (columns: `patient_barcode, tss`, the eight endpoint columns, then the `ENSG…` gene columns).

## Note on gdc-client (deprecated)

`gdc-client` is no longer used by the pipeline. `scripts/download_gdc_client.sh` is retained
only for historical reference and is not part of the data flow.
