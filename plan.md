# Project Plan: Site Confounding in TCGA Survival Models
## A Bayesian Variance Decomposition Across Cancer Types

---

## PART 1: BACKGROUND KNOWLEDGE TO STUDY

### 1.1 Oncology Fundamentals

**What you need to understand (not memorise):**

**Cancer biology basics**
- Cancer = uncontrolled cell division caused by DNA mutations
- Cancer is named by tissue of origin — lung, breast, kidney, colon, etc.
- Each cancer type has a distinct molecular profile, growth pattern, and survival curve
- The "hallmarks of cancer" (Hanahan & Weinberg): the 8-10 biological capabilities
  that define malignant tumours — sustained proliferative signalling, evasion of
  apoptosis, angiogenesis, invasion/metastasis, genomic instability, etc.
  These are the biological categories your genes will fall into.

**Why this matters for your project:**
You are not studying cancer biology per se. You need enough vocabulary to
interpret your results — e.g. "genes in the angiogenesis hallmark attenuate
under site adjustment" is a meaningful clinical finding only if you know what
angiogenesis means.

**What to read:**
- Hanahan D. (2022). "Hallmarks of Cancer: New Dimensions." Cancer Discovery.
  (The updated version — free on PubMed.)
  Read the abstract and Figure 1 only. You need the vocabulary, not the details.

---

**Genomics basics (RNA-seq)**
- Every cell in your body has the same DNA, but different genes are
  switched on or off in different tissues and disease states
- RNA-seq measures gene expression — how actively each gene is being
  transcribed in a tumour sample
- Higher expression of a gene = more of that gene's protein is being made
- TCGA RNA-seq data: ~20,000 genes measured per patient, expressed as
  FPKM (Fragments Per Kilobase of transcript per Million mapped reads)
  or raw counts. You will use log2(FPKM+1) normalised values.

**What to read:**
- Conesa et al. (2016). "A survey of best practices for RNA-seq data analysis."
  Genome Biology. Read sections 1 and 2 only (overview + preprocessing).
- The TCGA GDC documentation on RNA-seq data:
  https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/

---

**TCGA data structure**
- 33 cancer types, ~11,000 patients
- Each patient: RNA-seq, mutations, clinical data, survival endpoints
- Each sample has a barcode encoding patient + tissue source site (TSS)
- TSS = the hospital that contributed the sample
- TCGA-CDR (Liu et al. 2018, Cell) = the curated survival endpoints you
  already use (OS, PFI). Do not use raw vital_status.

**What to read:**
- Liu et al. (2018). "An integrated TCGA pan-cancer clinical data resource."
  Cell. (You already know this one — skim the methods.)
- Howard et al. (2021). "The impact of site-specific digital histology
  signatures on deep learning model accuracy and bias." Nature Communications.
  (The paper your project builds on — read fully.)

---

### 1.2 Survival Analysis (deepening what you know)

**What you already know:**
Right-censoring, Cox partial likelihood, C-index, Weibull model, IBS.

**What you need to add:**

**Parametric survival families and when to use them**
- Exponential: constant hazard over time (memoryless). Rarely realistic.
- Weibull: hazard increases or decreases monotonically. Most flexible
  parametric choice. What you already use.
- Log-normal: hazard rises then falls. Good for some cancers where
  early mortality is high but survivors stabilise.
- Log-logistic: similar to log-normal but heavier tails.
- How to choose: LOO-CV via ArviZ, or WAIC. You will use this for
  model comparison in the paper.

**Accelerated failure time (AFT) models**
- Alternative to Cox. Instead of modelling the hazard, you model
  log(survival time) directly as a linear function of covariates.
- More interpretable: a coefficient of 0.3 means "this gene multiplies
  expected survival time by e^0.3 = 1.35"
- Your Weibull model IS an AFT model. This is the frame you use
  in the paper.

**What to read:**
- Klein & Moeschberger (2003). "Survival Analysis: Techniques for Censored
  and Truncated Data." Chapter 12 (parametric models). Available as PDF.
- The PyMC survival analysis examples:
  https://www.pymc.io/projects/examples/en/latest/survival_analysis/

---

### 1.3 Bayesian Methods (deepening what you know)

**What you already know:**
MCMC, NUTS sampler, posterior distributions, credible intervals, 
hierarchical models, random intercepts, ICC.

**What you need to add:**

**Horseshoe prior**
- A continuous shrinkage prior for high-dimensional sparse regression
- Each coefficient β_j has its own local shrinkage parameter λ_j and
  shares a global shrinkage parameter τ
- The "horseshoe" shape: most coefficients get shrunk near zero (noise),
  but a few large signals pass through nearly unshrunken
- Better than Lasso (L1) because it doesn't uniformly penalise large effects
- The regularised horseshoe (Piironen & Vehtari 2017) is the version
  you use — it adds a slab component for numerical stability

**What to read:**
- Piironen J & Vehtari A. (2017). "Sparsity information and regularization
  in the horseshoe and other shrinkage priors." Electronic Journal of Statistics.
  Read sections 1-3. This is the paper behind the horseshoe you'll implement.
- Austin Rochford's PyMC implementation:
  https://austinrochford.com/posts/2021-05-29-horseshoe-pymc3.html
  (Code examples you'll adapt directly.)

**Model comparison in Bayesian framework**
- LOO-CV (Leave-One-Out Cross Validation): computed from posterior samples
  using ArviZ's az.loo(). Compares models without refitting.
- WAIC (Widely Applicable Information Criterion): similar to LOO-CV,
  also from ArviZ. Use both.
- ELPD (Expected Log Predictive Density): the quantity LOO-CV estimates.
  Higher = better model.
- You use these to compare: Weibull vs log-normal vs log-logistic,
  and naive model vs site-adjusted model.

**What to read:**
- Vehtari A, Gelman A, Gabry J. (2017). "Practical Bayesian model
  evaluation using leave-one-out cross-validation and WAIC."
  Statistics and Computing. Read sections 1-4.
- ArviZ documentation on model comparison:
  https://python.arviz.org/en/stable/api/generated/arviz.compare.html

**Convergence diagnostics**
- R-hat < 1.01: chains have mixed well
- ESS (Effective Sample Size) > 400 per parameter: enough samples
- Divergences: ideally zero, max ~0.1% of draws. More = model
  misspecification or need for reparameterisation.
- You already know these from your project. The paper requires
  reporting them formally.

---

### 1.4 RNA-seq Data Wrangling (genuinely new)

**What RNA-seq data looks like**
- A matrix: patients (columns) × genes (rows)
- Values are FPKM or raw counts
- ~20,000 genes per patient
- Highly skewed: most genes have near-zero expression, a few have
  very high expression
- Standard preprocessing: log2(FPKM + 1) transformation to
  reduce skewness and stabilise variance

**Gene selection strategy for your model**
You cannot use all 20,000 genes — the model won't converge and
the biology is uninterpretable. Two approaches:

*Option A — Cancer hallmark genes (recommended):*
Use the 671 cancer hallmark genes from the CancerSEA / Hanahan
framework, or the ~400 genes in established pan-cancer prognostic
studies (Nagy et al. 2021). This gives interpretable results —
you can say "genes in the angiogenesis hallmark lose prognostic
signal under site adjustment."

*Option B — Top variance genes:*
Select top 500 genes by variance across patients within each cohort.
Standard practice (used by most TCGA survival papers). Less
interpretable but more data-driven.

*Recommendation: Option A for the main analysis, Option B as
a robustness check.*

**Normalisation**
- Within-sample: log2(FPKM+1)
- Across-sample: standardise each gene to mean 0, std 1 before
  entering the model (so coefficients are comparable across genes)
- Do NOT mix TCGA tumour samples with normal tissue samples —
  use only tumour barcodes (sample type 01)

**What to read:**
- Conesa et al. (2016). Section 3 (normalisation) only.
- The GDC RNA-seq pipeline documentation (link above).
- Nagy et al. (2021). "Pancancer survival analysis of cancer hallmark
  genes." Scientific Reports. This is the closest existing work to
  your gene selection strategy.

---

### 1.5 Site Confounding and Batch Effects

**What you already know:**
Howard et al.'s site leakage finding, TSS codes, the CV gap.

**What you need to add:**

**PVCA (Principal Variance Component Analysis)**
- The frequentist method for variance decomposition in batch effect literature
- Decomposes variance in embedding/expression space (not survival space)
- What HistoAtlas (2026) used: 44.7% variance to TSS in embedding space
- Your ICC is different: variance decomposition in survival outcome space,
  using a proper survival likelihood. This is the distinction you need
  to articulate clearly in the paper.

**ComBat harmonisation**
- Standard batch correction method for RNA-seq and imaging data
- Removes estimated batch effects from the data before analysis
- The limitation: if site is confounded with biology (different cancer
  subtypes at different sites), ComBat removes signal along with noise
- Your approach is better for survival: you model site as a random effect
  and let the data determine how much to attribute to site vs biology,
  rather than removing it. This is an important methodological argument
  in your paper.

**What to read:**
- Johnson et al. (2007). "Adjusting batch effects in microarray expression
  data using empirical Bayes methods." Biostatistics. (The ComBat paper —
  read abstract + introduction to understand what you're improving on.)
- Howard et al. (2021). (Re-read — you need to know their specific numbers
  for LUSC and HNSC to compare against yours.)

---

### 1.6 Reporting Frameworks

**TRIPOD-AI:** Reporting guidelines for prediction model studies using AI.
Covers: participants, outcomes, predictors, sample size, model development,
performance, validation. Write your methods section against this checklist.

**PROBAST-AI:** Risk of bias tool for prediction model studies.
Covers: participant selection, predictors, outcomes, analysis.
Cite this when critiquing existing literature.

**What to read:**
- Collins et al. (2021). "Protocol for development of TRIPOD-AI and
  PROBAST-AI." BMJ Open. Read the checklist items only (Supplementary).

---

## PART 2: KEY PAPERS TO READ

### Must-read (read fully):
1. **Howard et al. (2021).** "The impact of site-specific digital histology
   signatures on deep learning model accuracy and bias."
   *Nature Communications.* DOI: 10.1038/s41467-021-24698-1
   → The paper your project directly extends.

2. **Liu et al. (2018).** "An integrated TCGA pan-cancer clinical data
   resource to drive high-quality survival outcome analytics." *Cell.*
   DOI: 10.1016/j.cell.2018.02.052
   → Your survival labels. Know the endpoints inside out.

3. **Samorodnitsky et al. (2020).** "A pan-cancer and polygenic Bayesian
   hierarchical model for the effect of somatic mutations on survival."
   *Cancer Research.* arXiv: 1910.03447
   → The closest existing work. Know exactly how yours differs
   (they have no site random effect).

4. **Piironen & Vehtari (2017).** "Sparsity information and regularization
   in the horseshoe and other shrinkage priors."
   *Electronic Journal of Statistics.* DOI: 10.1214/17-EJS1337SI
   → The statistical foundation for your variable selection.

5. **Nagy et al. (2021).** "Pancancer survival analysis of cancer hallmark
   genes." *Scientific Reports.* DOI: 10.1038/s41598-021-84787-5
   → Your gene selection rationale. Shows which hallmark genes are
   prognostic across cancer types without site adjustment.

### Should-read (read abstract + results):
6. **Samorodnitsky et al. (2022).** "A hierarchical spike-and-slab model
   for pan-cancer survival using pan-omic data." *BMC Bioinformatics.*
   DOI: 10.1186/s12859-022-04770-3
   → Extension of paper 3. Know what they did so you can differentiate.

7. **Hanahan D. (2022).** "Hallmarks of cancer: new dimensions."
   *Cancer Discovery.* DOI: 10.1158/2159-8290.CD-21-1059
   → Background vocabulary for interpreting gene results.

8. **Johnson et al. (2007).** "Adjusting batch effects in microarray
   expression data using empirical Bayes methods." *Biostatistics.*
   DOI: 10.1093/biostatistics/kxj037
   → The ComBat paper. Understand what you're improving on.

9. **Vehtari et al. (2017).** "Practical Bayesian model evaluation using
   leave-one-out cross-validation and WAIC." *Statistics and Computing.*
   DOI: 10.1007/s11222-016-9696-4
   → Your model comparison methodology.

10. **Chen et al. (2024).** "Towards a general-purpose foundation model
    for computational pathology." *Nature Medicine.*
    DOI: 10.1038/s41591-024-02857-3
    → The UNI2 paper. Understand what the embeddings represent.

### Skim (abstract only):
11. Zhao et al. (2024). "Tutorial on survival modeling with applications
    to omics data." *Bioinformatics.* DOI: 10.1093/bioinformatics/btae132

12. Murchan et al. (2024). "Deep feature batch correction using ComBat
    for machine learning applications in computational pathology."
    *Journal of Pathology Informatics.*

---

## PART 3: MODEL PROPOSAL

### 3.1 Research Questions

**RQ1:** How much of survival outcome variance across TCGA cancer types
is attributable to tissue source site (institution), after controlling
for genomic features?

**RQ2:** Does the estimated site variance component (ICC_site) vary
systematically across cancer types?

**RQ3:** Which genomic features lose prognostic signal when site is
explicitly modelled as a confounder — suggesting they may be proxies
for institutional rather than biological variation?

---

### 3.2 Data

**Cohorts:**
All TCGA cancer types with ≥ 50 death events in the OS endpoint.
Expected: ~18-22 cancer types, ~7,000-9,000 patients total.

**Survival labels:**
TCGA-CDR Table S1 (Liu et al. 2018). Endpoint: OS (overall survival).
`time` = OS.time (days). `event` = OS (1=dead, 0=censored).

**Genomic features:**
RNA-seq FPKM-UQ data from TCGA GDC portal.
- Download via GDC API (Python `requests` library or `gdc-client`)
- Filter to primary tumour samples only (barcode positions 14-15 = "01")
- Select cancer hallmark genes: 671 genes from the Hanahan framework
  (available via MSigDB hallmark gene sets)
- Preprocessing: log2(FPKM+1) → standardise per gene across patients
  within each cancer type

**Site labels:**
Extract TSS code (positions 2-3) from patient barcode. 
Join to GDC tissue source site lookup table.

---

### 3.3 Model Architecture

**Two models per cancer type:**

#### Model A — Naive (no site control)
```
log(t_i) = μ + β · x_i + ε_i

t_i     ~ Weibull(α, exp(-η_i))  [survival time]
η_i     = μ + β · x_i            [linear predictor]
ε_i     ~ Gumbel(0, 1/α)         [residual]

Priors:
μ       ~ Normal(6, 2)            [global intercept, ~log(400 days)]
α       ~ HalfNormal(1)           [Weibull shape]
τ       ~ HalfStudentT(3, 0, σ_τ) [horseshoe global shrinkage]
λ_j     ~ HalfCauchy(1)           [horseshoe local shrinkage, per gene]
β_j     ~ Normal(0, τ·λ_j)        [gene coefficients, horseshoe prior]
σ_τ     = p0/(p-p0) · σ/√n        [Piironen & Vehtari (2017) calibration]
                                   p0 = expected number of non-zero genes
                                   p  = total genes (671)
                                   n  = patients in cohort
```

#### Model B — Site-adjusted
```
log(t_i) = μ + β · x_i + u_{s(i)} + ε_i

Additional:
σ_site  ~ HalfNormal(1)           [site variance component]
u_s     ~ Normal(0, σ_site)        [site random intercepts, s=1..S]
         [non-centred: u_s = σ_site · z_s, z_s ~ Normal(0,1)]
```

**ICC computation (from posterior samples):**
```
σ²_residual = π² / (6 · α²)       [Weibull/Gumbel residual variance]
σ²_site     = σ_site²
σ²_genomic  = Var(β · X)           [empirical variance of linear predictor]

ICC_site    = σ²_site / (σ²_site + σ²_genomic + σ²_residual)
ICC_genomic = σ²_genomic / (σ²_site + σ²_genomic + σ²_residual)
ICC_residual = σ²_residual / (σ²_site + σ²_genomic + σ²_residual)
```

---

### 3.4 Gene Attenuation Analysis (RQ3)

For each gene j and each cancer type k, compute:

```
attenuation_jk = |β_jk(A)| - |β_jk(B)|
```

Where β_jk(A) is the posterior mean of gene j's coefficient in Model A,
and β_jk(B) is the posterior mean in Model B.

A positive attenuation means the gene's apparent effect shrinks when site
is controlled for — evidence that part of its signal was institutional
rather than biological.

**Reporting:**
- Rank genes by mean attenuation across cancer types
- Flag genes where attenuation > 1 posterior SD (credibly attenuated)
- Group attenuated genes by cancer hallmark category
- Key claim: "Genes in [hallmark X] show systematic attenuation under
  site adjustment across [N] cancer types"

---

### 3.5 Model Comparison

For each cancer type, compare Model A vs Model B using:
- `az.compare()` with LOO-CV (ArviZ)
- Report ELPD difference and standard error
- A positive ELPD difference for Model B = site adjustment improves
  predictive accuracy = site carries real information beyond genomics

---

### 3.6 Implementation Plan

#### Infrastructure
- **Local laptop:** Data download, preprocessing, single-cohort testing
- **Hetzner CCX53** (32 vCPU, 128GB RAM, ~€0.40/hr): Full model runs
  across all cancer types. Run in tmux, detach, check next morning.
- **Estimated compute cost:** ~€20-50 total for all runs

#### Tech stack
```
Python 3.11 (not 3.14 — PyMC incompatibility)
PyMC 5.x
ArviZ (model comparison, diagnostics)
pytensor
pandas, numpy
requests / gdc-client (RNA-seq download)
scikit-survival (C-index validation)
matplotlib / seaborn (figures)
```

#### Directory structure
```
tcga-site-genomics-survival/
├── data/
│   ├── raw/          # downloaded RNA-seq + CDR files
│   ├── processed/    # clean patient tables per cohort
├── models/
│   ├── model_a.py    # naive model
│   ├── model_b.py    # site-adjusted model
├── notebooks/
│   ├── 01_data_prep.ipynb
│   ├── 02_eda.ipynb
│   ├── 03_model_fit.ipynb
│   ├── 04_results.ipynb
├── results/
│   ├── traces/       # saved PyMC trace objects (.nc files)
│   ├── figures/
├── paper/
│   ├── main.tex
│   ├── figures/
└── README.md
```

---

### 3.7 Build Stages

#### Stage 0 — Data pipeline (1-2 weeks)
- [ ] Download TCGA RNA-seq for all viable cohorts via GDC API
- [ ] Download TCGA-CDR survival labels
- [ ] Extract TSS from barcodes, build site lookup
- [ ] Filter to tumour samples, select hallmark genes (671)
- [ ] Preprocess: log2(FPKM+1), standardise per gene
- [ ] Join RNA-seq + survival + site into patient tables
- [ ] Sanity checks: event rates, site distributions, gene variance
- **Deliverable:** One clean CSV per cancer type

#### Stage 1 — Single cohort pilot (1 week)
- [ ] Fit Model A (naive) on LUSC — check convergence
- [ ] Fit Model B (site-adjusted) on LUSC — check convergence
- [ ] Compute ICC, gene attenuation for LUSC
- [ ] Compare LOO-CV scores
- **Deliverable:** Working model on one cohort

#### Stage 2 — Full pan-cancer run (1-2 weeks)
- [ ] Port pipeline to Hetzner VPS
- [ ] Run Model A + Model B across all cohorts
- [ ] Save traces as NetCDF (.nc) files with ArviZ
- [ ] Collect convergence diagnostics — flag any problematic cohorts
- **Deliverable:** Complete trace files for all cohorts

#### Stage 3 — Analysis and figures (1 week)
- [ ] ICC per cancer type — ranked bar chart with 95% HDIs
- [ ] Three-component ICC decomposition (site / genomic / residual)
- [ ] Gene attenuation ranking — top 20 attenuated genes across cohorts
- [ ] LOO-CV comparison Model A vs B per cancer type
- [ ] Kaplan-Meier validation: do the high/low risk groups from Model B
      separate better than Model A?
- **Deliverable:** All paper figures

#### Stage 4 — Paper (1-2 weeks)
- [ ] Write in LaTeX using bioRxiv/PLOS template
- [ ] Sections: Abstract, Introduction, Methods, Results, Discussion,
      Limitations
- [ ] Post to arXiv (need endorsement first — approach LSE professor
      at end of Stage 3 with near-complete draft)
- **Deliverable:** arXiv preprint

---

### 3.8 Figures for the Paper

**Figure 1 — ICC decomposition per cancer type**
Stacked bar chart: for each cancer type, three bars showing
ICC_site / ICC_genomic / ICC_residual. Sorted by ICC_site descending.
Main finding: which cancer types are most site-confounded.

**Figure 2 — Gene attenuation under site adjustment**
Scatter plot: x = |β(A)| (naive), y = |β(B)| (site-adjusted).
Points below the diagonal = attenuated. Colour by hallmark category.
Named points = top 10 most attenuated genes.

**Figure 3 — LOO-CV comparison**
Forest plot: ELPD difference (Model B - Model A) per cancer type
with 95% CI. Positive = site adjustment improves prediction.

**Figure 4 — Survival discrimination**
C-index comparison: Model A vs Model B under site-preserved CV,
per cancer type. (Your existing pipeline from the first project.)

---

### 3.9 Limitations to acknowledge

- TCGA cohorts are small (200-500 patients per type). Wide credible
  intervals on gene coefficients are expected and honest.
- RNA-seq batch effects are not fully separable from site effects.
  ComBat correction as a sensitivity analysis.
- Hallmark gene selection introduces prior biological assumptions.
  Top-variance gene selection as a robustness check.
- Single-country, non-representative cohort (TCGA is predominantly
  US-based with known demographic skews).
- No external validation — a limitation for all TCGA-based analyses.

---

## PART 4: TIMELINE

| Week | Task |
|------|------|
| 1 | Read papers 1-5. Study horseshoe prior. Set up GDC API access. |
| 2 | Download RNA-seq data. Build preprocessing pipeline. |
| 3 | Fit pilot model on LUSC. Debug convergence. |
| 4 | Port to Hetzner. Run full pan-cancer Model A. |
| 5 | Run full pan-cancer Model B. Collect diagnostics. |
| 6 | Analysis: ICC, attenuation, LOO-CV. Build all figures. |
| 7 | Write paper draft. |
| 8 | Approach LSE professor for endorsement. Submit to arXiv. |

**Total: 8 weeks. Total Hetzner cost: ~€50.**
