# Novelty & Prior-Work Assessment: Bayesian Site-vs-Genomic Variance Decomposition of TCGA Survival

## TL;DR
- **The integrated project is novel as a whole**: no existing work computes a Bayesian, per-posterior-draw ICC variance decomposition of TCGA *survival outcomes* into tissue-source-site vs. genomic vs. residual components using a Weibull AFT model with site random intercepts and a regularized horseshoe prior over hallmark genes. Every individual building block has precedent, but their combination — and specifically the site-vs-genomic decomposition in *outcome space* — appears genuinely unoccupied.
- **The "published-signature-collapse-under-site-CV for RNA-seq" experiment has NOT been done**: site-preserved CV inflation has been demonstrated for histology images (Howard et al. 2021) and, modestly, for deep-learning omic features (Howard et al. 2022 Cancer Cell), but no one has taken a *named, previously-published bulk RNA-seq prognostic signature* and shown its C-index collapsing toward chance under leave-one-site-out CV within TCGA. This is a real, citable gap.
- **Closest works and how they differ**: Howard et al. 2021/2022 (establishes the site-leakage mechanism, but on images/DL-omics not expression signatures, and measures discrimination not variance); Samorodnitsky et al. 2019/2022 (Bayesian hierarchical pan-cancer survival, but on mutations/omic-factors with no site effect); Yu et al. 2020 "coxlmm" (partitions survival variance into clinical-vs-genomic but has no site term and is frequentist); and Ni/Qin "BatMan" (batch-induced spurious survival C-index, but miRNA, shipment-batch, simulation).

## Key Findings

**Verdict on each component:**

1. **Bayesian hierarchical survival with TSS random effects to decompose site vs genomic variance on TCGA RNA-seq — NO precedent found.** This is the novel core.
2. **ICC-style variance decomposition of survival OUTCOMES attributed to TCGA site — NO precedent found.** Survival-outcome variance decomposition exists in healthcare institution-profiling (Saarela; Chen et al.), and within TCGA it is done only for clinical-vs-genomic (coxlmm) or in *expression/embedding space* (tissue-of-origin), never for outcome-space site attribution.
3. **Horseshoe prior + site/batch random effects in a survival model on TCGA — NO combined precedent found.** Each ingredient exists separately; the combination is unoccupied.
4. **Leave-one-site-out CV for RNA-seq survival prognosis in TCGA — essentially NOT done.** Site-stratified CV is standard for *histology* (Howard 2021) and is sometimes applied to multimodal models that include expression, but a dedicated leave-one-site-out evaluation of *expression-only* survival prognosis is absent.
5. **Published expression signature collapsing under site-preserved CV — NOT done for RNA-seq.** Confirmed gap.
6. **Closest existing work**: Howard et al. 2022 Cancer Cell (closest on mechanism/molecular side) and Yu et al. 2020 coxlmm (closest on the variance-decomposition side).

## Details

### Components 1 & 2 — Bayesian hierarchical survival + ICC outcome-space decomposition

The strongest near-miss on the *variance-decomposition* axis is **Yu, Wang, Huang & Zeng 2020, "How Can Gene-Expression Information Improve Prognostic Prediction in TCGA Cancers: An Empirical Comparison Study on Regularization and Mixed Cox Models" (Frontiers in Genetics 11:920; PMC7472843)** — note the lead author is Xinghao Yu, not "Zhou." They fit a linear mixed-effects Cox model (coxlmm) across 32 TCGA cancers and explicitly "partition the survival variance into the relative contribution of clinical and transcriptomic components," defining PCE (proportion explained by clinical) and PGE (proportion explained by genomic). They report that "the average prediction gain was 2.4% for coxlasso, 4.2% for coxenet, and 7.2% for coxlmm across 16 low-censored cancers," and that the transcriptomic contribution "was higher for some cancers (e.g., LGG, CESC, PAAD, SKCM, and SARC) and lower for others (e.g., BRCA, COAD, KIRC, and STAD)." **Critically, there is no tissue-source-site term** — the partition is clinical-vs-genomic, it is frequentist (coxme/Laplace approximation, not Bayesian posterior draws), and it uses a Gaussian random effect over a genetic relationship matrix rather than a horseshoe over named hallmark genes. This is the paper whose framing is closest to RQ1/RQ2, but it answers a different question. (Their per-cancer PVE/PCE/PGE figure is the place to verify any exact decimals; the abstract reports prediction *gains* rather than a single PVE range.)

On the *Bayesian hierarchical pan-cancer survival* axis, the **Samorodnitsky et al.** line is closest:
- Samorodnitsky, Hoadley & Lock 2019/2020, "A Pan-Cancer and Polygenic Bayesian Hierarchical Model for the Effect of Somatic Mutations on Survival" (arXiv 1910.03447) — hierarchical Weibull/log-normal/exponential AFT across 27 cancer types, effects shrunk toward a common mean, **but on somatic mutations, not expression, and with no site random effect.**
- Samorodnitsky, Hoadley & Lock 2022, "A hierarchical spike-and-slab model for pan-cancer survival using pan-omic data" (BMC Bioinformatics 23:235; PMC9204947) — Bayesian hierarchical survival across 29 cancer types and 4 omics sources, variable selection via spike-and-slab over BIDIFAC+ factors, **but the hierarchical structure is over cancer types and omic factors, not over contributing institutions, and there is no ICC site decomposition.**

Other Bayesian-survival-on-TCGA precedents establish that the horseshoe/AFT pieces are individually standard: a Bayesian log-normal/AFT with horseshoe for circadian gene selection integrating CNV+RNA-seq (PMC8775911); Bayesian AFT + horseshoe for joint survival/binary RPPA protein modeling (Maity et al., Biometrics 2020 / PMC7729996); and rstanarm's documented multilevel survival model with hospital-site random intercepts and ICC interpretation on the AFT log-time scale (arXiv 2002.09633). Stata's `mestreg` documentation likewise treats ICCs after multilevel AFT survival models as a known construct. **None of these is applied to decompose TCGA survival-outcome variance by tissue source site.**

The healthcare institutional-comparison literature (Saarela; Chen, Lawson, Finelli & Saarela 2020, "Causal variance decompositions for institutional comparisons in healthcare," Statistical Methods in Medical Research) supplies the methodological template for hospital/site variance decomposition and its ICC connection — but in quality-of-care profiling, not cancer genomics, and with no genomic signal term.

### Component 3 — Horseshoe + site random effects in survival

Horseshoe-in-survival is well established (Piironen & Vehtari 2017 is the regularized-horseshoe method; multiple TCGA applications above). Site/batch random effects in survival are established (frailty/mixed AFT). I found **no paper combining a horseshoe (or regularized horseshoe) prior over genes with site/batch random intercepts in a survival model on TCGA.** This combination is novel.

### Component 4 — Leave-one-site-out CV for RNA-seq survival in TCGA

The site-leakage mechanism and "site-preserved cross-validation" originate with **Howard, Dolezal, Kochanny et al. 2021, Nature Communications 12:4423**, on *digital histology*: "we demonstrate that these features vary substantially across tissue submitting sites in TCGA for over 3,000 patients with six cancer subtypes," and that site-specific signatures "lead to an overestimation of model accuracy if multiple sites are included in both the training and validation datasets" for predictions including survival, gene expression, mutations, and ancestry. The follow-up **Howard, Kather & Pearson 2022/2023 Cancer Cell** letter applied site-preserved CV to the multimodal PORPOISE model (which includes RNA-seq among its omic inputs). Site-stratified CV "to mitigate batch effect (Howard et al., 2021)" is now cited as standard practice in multimodal survival papers (e.g., Multimodal Prototyping, arXiv 2407.00224). However, a dedicated **leave-one-site-out evaluation of expression-only survival prognosis across pan-cancer TCGA is essentially absent** — existing transcriptomic survival benchmarks (the coxlmm study; deep-learning RNA-seq survival benchmarks) use Monte Carlo, k-fold, or LOOCV that do not preserve sites. A related but distinct effort is the clustering-stratified CV framework for omics survival (Saadati et al., BMC Med Res Methodol 2025; PMC12709853), which stratifies by transcriptomic clusters and ComBat, not by tissue source site.

### Component 5 — Published RNA-seq signature collapsing under site-preserved CV

**This specific experiment has not been published.** After dedicated searching, no study takes a named, previously-published bulk RNA-seq prognostic signature (with a pre-reported C-index ~0.65) and shows its discrimination collapsing toward ~0.5 under leave-one-site-out/site-preserved CV within TCGA attributed to tissue-source-site leakage. The conceptual building blocks exist separately:

- **Howard, Kather & Pearson 2022 Cancer Cell + the PORPOISE_SITE GitHub repo (github.com/fmhoward/PORPOISE_SITE)** is the closest on the molecular side. It re-ran PORPOISE under site-preserved 5-fold CV. For the genomic/omic-only model (a self-normalizing network over RNA-seq + CNV + mutation signatures), the average C-index declined from ≈0.582 under standard CV to ≈0.554 under site-preserved CV (per the repo's per-cancer tables) — a modest decline, **not** collapse; KIRC and LUAD were unperturbed or improved, while HNSC fell below 0.5 (0.534→0.476). The letter notes that "performance in some cancer subtypes, such as clear cell renal carcinoma (KIRC), was unperturbed by site-preserved cross-validation and will likely remain accurate in external patient cohorts." **Crucially, these are PORPOISE's own deep-learning omic features, not a named, previously-published gene-expression signature.** (Exact decimals are repo-sourced; the journal full text was bot-blocked.)
- **Ni, Liu & Qin, "BatMan" (arXiv 2209.03902)** shows in simulation with TCGA ovarian data that batch effects confounded with survival can induce spurious C-index above 0.5 (~0.52–0.56), with null-model C-index ~0.5 — but on **miRNA microarray data, with shipment/array batch (not tissue source site), building its own models rather than re-evaluating a named signature.**
- **Caballé Mestres, Berenguer Llergo & Stephan-Otto Attolini 2018 (bioRxiv 360495)** tests a named signature (MammaPrint 70-gene) and shows technical bias distorts gene-signature survival associations — but the metric is hazard-ratio bias / Type-I-error inflation in GEO cohorts, not C-index collapse under TCGA site-stratified CV.
- **Venet, Dumont & Detours 2011 (PLoS Comput Biol 7(10):e1002240)** established signature skepticism: of 47 published breast-cancer signatures, "28 of them (60%) were not significantly better outcome predictors than random signatures of identical size and 11 (23%) were worse than the median random signature," and "more than 90% of random signatures >100 genes were significant outcome predictors." This is the canonical "random signatures are prognostic" result — but not site-leakage and not TCGA site-preserved CV.

So the proposed Experiment is novel; it would unify these strands for the first time on a named bulk-mRNA signature with TSS leave-one-site-out in TCGA.

### Supporting context

TCGA tissue-source-site batch structure is extensively documented (MD Anderson TCGA Batch Effects Viewer; Rasnic et al. 2019; Choi et al. 2017 on somatic-variant batch effects; the PRPS paper, Nature Biotechnology 2022, on RNA-seq TSS/plate confounding). The TCGA Pan-Cancer Clinical Data Resource (**Liu, Lichtenberg, Hoadley et al. 2018, Cell 173(2):400–416**) explicitly warns: "Because endpoint-confounding factors from different TSS populations might include patient age, tumor stage/grade, and treatment, TSSs might serve as a proxy for these as well as other unmeasured differences, including incomplete clinical annotation." This directly motivates RQ1. Hallmark-gene survival without site adjustment or shrinkage is **Nagy, Munkácsy & Győrffy 2021 (Sci Rep 11:6047)**: "RNA-seq HTSeq counts and survival data from 26 different tumor types were acquired from the TCGA repository," analyzed with univariate Cox per gene (no multivariable shrinkage, no site term); "renal clear cell cancer and low grade gliomas harbored the most prognostic changes … while thyroid and glioblastoma were largely independent of hallmark genes."

## Recommendations

**Positioning (do this):**
1. **Claim novelty at the level of the integrated framework and specifically the outcome-space site/genomic ICC decomposition.** Frame the contribution as the "first Bayesian variance decomposition of TCGA *survival outcomes* into tissue-source-site vs genomic vs residual components." Explicitly contrast with Yu et al. 2020 (clinical-vs-genomic, no site, frequentist) and Samorodnitsky 2022 (hierarchy over cancer types/omic factors, no site ICC).
2. **Do not claim novelty for any single ingredient.** Horseshoe-in-survival, Bayesian AFT on TCGA, hierarchical pan-cancer Bayesian survival, ICC-after-AFT, and site-preserved CV all have precedent; cite them as the components you combine.
3. **For the signature-collapse experiment, claim it as the first RNA-seq demonstration.** State precisely: site-preserved CV collapse has been shown for histology (Howard 2021) and partially for DL-omic features (Howard 2022), but never for a *named, published bulk-mRNA expression signature* in TCGA. Choose a signature with a clearly reported C-index in a cancer type with many tissue source sites (e.g., LUAD, BRCA, or LGG/KIRC) to maximize site multiplicity — but note that KIRC was robust in PORPOISE_SITE, so it may be a poor choice for demonstrating collapse and a good choice for demonstrating robustness.

**Methodological guardrails (benchmarks that would change the framing):**
4. If KIRC-type cancers show site ICC near zero and stable C-index under leave-one-site-out (as PORPOISE_SITE found), report this as a feature: the method should show that site effects are *cancer-type-specific* (RQ2), not universal. A finding that site ICC is uniformly negligible would weaken the headline and should be reported honestly.
5. Calibrate the signature-collapse claim against an explicit null: include random-signature and site-only (intercept) Cox baselines so that a drop "toward chance" is referenced against an actual ~0.5 baseline (mirroring BatMan's null C-index and Venet's random-signature framework).
6. If a pre-submission literature sweep surfaces a preprint doing the exact RNA-seq signature-collapse experiment, downgrade the Experiment to "independent confirmation/extension" and shift emphasis to the Bayesian outcome-space decomposition, which remains unoccupied.

## Caveats
- Absence of evidence is not proof of absence: the negative findings (RQs 1–5) rest on extensive but not exhaustive search; a very recent preprint could exist. The signature-collapse and site-ICC-outcome claims are the most confidently novel.
- The Howard 2022 Cancer Cell per-cancer molecular C-index numbers come from the authors' own GitHub repo tied to the peer-reviewed letter; the journal full text was bot-blocked, so exact decimals should be treated as repo-sourced and re-verified.
- "Variance decomposition" is used inconsistently across fields; some TCGA papers use the phrase for expression/embedding space (tissue-of-origin variance) rather than outcome space — be precise in your own framing to avoid a reviewer conflating them.
- Caballé Mestres et al. 2018 is a preprint version; verify peer-reviewed status before citing as established.
- The exact OV/LGG PVE percentages should be read off Yu et al. 2020's PVE/PCE/PGE figure rather than the abstract, which reports prediction gains rather than a single PVE range.

### Key sources (with URLs)
- Howard et al. 2021, Nat Commun 12:4423 — https://www.nature.com/articles/s41467-021-24698-1
- Howard, Kather & Pearson 2022/2023, Cancer Cell 41(1):5–6 — https://pubmed.ncbi.nlm.nih.gov/36368319/ ; repo: https://github.com/fmhoward/PORPOISE_SITE
- Samorodnitsky et al. 2019, arXiv 1910.03447 — https://arxiv.org/abs/1910.03447
- Samorodnitsky, Hoadley & Lock 2022, BMC Bioinformatics 23:235 — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9204947/
- Yu et al. 2020 (coxlmm), Front Genet 11:920 — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7472843/
- Nagy, Munkácsy & Győrffy 2021, Sci Rep 11:6047 — https://www.nature.com/articles/s41598-021-84787-5
- Venet, Dumont & Detours 2011, PLoS Comput Biol 7(10):e1002240 — https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002240
- Liu et al. 2018 (TCGA-CDR), Cell 173(2):400–416 — https://www.cell.com/cell/fulltext/S0092-8674(18)30229-0
- Ni/Liu/Qin "BatMan", arXiv 2209.03902 — https://arxiv.org/abs/2209.03902
- Caballé Mestres et al. 2018, bioRxiv 360495 — https://www.biorxiv.org/content/10.1101/360495
- Maity et al. 2020 / Bayesian AFT+horseshoe, PMC7729996 — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7729996/
- Bayesian survival w/ rstanarm (site random intercepts, ICC), arXiv 2002.09633 — https://arxiv.org/abs/2002.09633
- Chen, Lawson, Finelli & Saarela 2020, Stat Methods Med Res (causal variance decomposition, institutions) — https://journals.sagepub.com/doi/abs/10.1177/0962280219880571
- Saadati et al. 2025 (clustering-stratified CV), BMC Med Res Methodol — https://pmc.ncbi.nlm.nih.gov/articles/PMC12709853/