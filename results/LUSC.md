# LUSC — our result vs. the literature

## What we found (TCGA-LUSC, 494 patients, OS endpoint)

Bayesian Weibull AFT with a regularized horseshoe over the **MSigDB Hallmark**
panel (3,867 protein-coding genes expressed in LUSC), site (TSS) random
intercepts for Model B. nutpie NUTS, 0–9 divergences, r-hat ≈ 1.00.

| quantity | value (89% interval) |
|---|---|
| ICC_site | **0.099** (0.032 – 0.197) |
| ICC_genomic | 0.015 (0.000 – 0.049) |
| ICC_residual | 0.886 (0.786 – 0.960) |
| prognostic gene hits (ETI excl. 0) | **0 / 3,867** (both models) |
| LOSO C-index, Model A (naive) | **0.531** |
| LOSO C-index, Model B (site-adj) | **0.534** |
| LOSO C-index gap (B − A) | +0.003 (−0.006 – +0.011) |

Headline: site explains ~6× the *explainable* survival variance that the gene
panel does, and the genes have **no out-of-sample discriminative value**
(C-index ≈ 0.53, i.e. ≈ chance) under leave-one-site-out evaluation.

## What the literature reports (LUSC RNA-seq prognostic models)

Published, optimised LUSC expression signatures report **modest** discrimination:

- Aging-related gene signature: C-index **0.628** (95% CI 0.586–0.671). [PMC8921527]
- Five-lncRNA prognostic index: AUC **0.652**. [PMC5584202]
- Seven-lncRNA signature: AUC **0.66–0.67** (1/3/5-yr). [PMC11973350]

So the *published* discrimination ceiling for LUSC OS from expression is
~0.62–0.67 — already low — and our pre-specified-panel, site-honest estimate
is ~0.53.

## The difference

**~0.53 (ours) vs ~0.65 (literature)** — a gap of roughly 0.12 in C-index/AUC.
The question is whether the literature is *inflated*, our estimate is
*under-powered/mis-specified*, or both. Possibilities, most-to-least likely:

### A. The literature numbers are optimistic (inflation)
1. **Outcome-driven gene selection.** Published signatures are *built* by
   screening the same cohort for survival-associated genes, then reporting fit
   on (often) that cohort. This selection-on-the-outcome inflates apparent
   performance; our panel is fixed a priori (hallmark), so it cannot benefit.
2. **No site-preserved validation.** Their CV is random/internal, which lets
   institutional (TSS) signatures leak across the train/test split and inflates
   accuracy — the exact effect Howard et al. (2021) demonstrated. Our
   leave-one-site-out design removes this, which *necessarily* lowers the number.
3. **In-sample / time-point AUC vs out-of-sample Harrell C.** Several reports
   are time-dependent AUC on the development cohort, not pooled out-of-fold
   Harrell C. Different, generally more generous, metrics.

### B. Our estimate is conservative or mis-targeted (deflation)
4. **Different feature set.** The cited signatures lean on **lncRNAs** and
   bespoke gene lists that are **not in the Hallmark collection** (protein-coding
   only). We may simply be missing the features they used — this is a genuine
   non-leakage explanation and should not be dismissed.
5. **Heavy Bayesian shrinkage.** The regularized horseshoe pulls weak
   coefficients to ~0, so we forgo chance correlations a penalised/stepwise Cox
   signature would exploit. Conservative by design (the `p0` sweep shows the
   null is robust, but shrinkage still lowers absolute discrimination).
6. **AFT vs Cox / endpoint framing.** We model OS as Weibull-AFT; many
   signatures use Cox PH. Specification differences can shift C modestly.

### C. The signal really is weak (both are "right")
7. **LUSC OS is genuinely hard to predict from expression.** Even the optimistic
   literature tops out at ~0.65. Venet, Dumont & Detours (2011) showed that even
   biologically *unrelated* signatures (e.g. postprandial laughter, mouse social
   defeat, skin-fibroblast localization) are significantly associated with
   breast-cancer outcome — i.e. outcome association is largely non-specific and
   easy to obtain by chance. Our ~0.53 is what an honest, non-cherry-picked
   estimate looks like when the underlying signal is thin.

## Interpretation for the paper

The gap between ~0.65 (literature) and ~0.53 (ours) is **not a failure of our
model — it is the quantity of interest.** It is consistent with selection
optimism + site leakage inflating published estimates, exactly the bias this
study targets. The most important *alternative* to rule out is **(4)**: that the
shortfall is just our restricted hallmark panel, not bias. 

**Discriminating tests (future work):**
- Re-run with a **top-variance** panel (plan's robustness option) and with a
  superset that includes the genes/lncRNAs from the cited signatures — if
  C-index stays ~0.53, (4) is rejected and the inflation story strengthens.
- Reproduce a published signature **under our leave-one-site-out CV** — if its
  ~0.65 collapses toward ~0.53, that directly demonstrates the leakage/selection
  inflation (the cleanest possible result).
- The pan-cancer run (RQ2) shows whether LUSC is unusually null or typical.

## Sources
- Howard et al. 2021, *Nat Commun* — site leakage / preserved-site CV: https://www.nature.com/articles/s41467-021-24698-1
- Aging-related signature (C 0.628): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8921527/
- Five-lncRNA index (AUC 0.652): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5584202/
- Seven-lncRNA signature (AUC 0.66–0.67): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11973350/
- Venet, Dumont & Detours 2011, *PLoS Comput Biol* e1002240 — most random gene signatures associate with breast-cancer outcome: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002240
