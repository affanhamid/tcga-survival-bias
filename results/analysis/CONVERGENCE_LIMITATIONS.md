# Convergence limitation (stated)

Convergence was triaged for every saved trace. A pre-specified, mechanical rule was applied identically to all cohorts: a fit's **ICC-driver scalars** ({mu, alpha, tau, c2, sigma_site}) must satisfy max R-hat ≤ 1.05 and min bulk-ESS ≥ 100; otherwise the fit is flagged *convergence-severe* and the corresponding cohort×endpoint is **excluded** from the ICC cross-cohort analyses (Model B gates ICC; Model A is used only for gene-hit comparison).

The regularized-horseshoe Weibull-AFT has a residual funnel geometry that higher `target_accept` (0.999), longer tuning (2000), and more draws (2000) did not resolve for a minority of small / low-event cohorts — a limit of the model parameterization, not of compute. Rather than alter the prior post hoc, these fits are excluded under the rule above and reported transparently.

**Excluded (3 cohort×endpoints):** SARC/PFI, BLCA/PFI, PAAD/PFI.

Cross-cohort correlations are reported on both the full and the clean (post-exclusion) set (`results/analysis/crosscohort_r.csv`); see below. With ~20 cohorts and two patient-sharing endpoints these remain exploratory.


| relationship | r (full) | n | r (clean) | n |
|---|---|---|---|---|
| ICC_genomic vs out-of-sample C-index (Model B) | 0.767 | 39 | 0.875 | 36 |
| ICC_site vs out-of-sample C-index (Model B) | -0.142 | 39 | -0.155 | 36 |
| ICC_site OS vs PFI (paired) | 0.487 | 18 | 0.527 | 15 |
