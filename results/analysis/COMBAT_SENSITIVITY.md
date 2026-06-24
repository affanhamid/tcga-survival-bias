# ComBat sensitivity (matched-baseline)

The ComBat contrast is reported as **dropped-uncorrected (`hallmark+dropsing`) vs dropped-corrected (`hallmark+combat`)** — the same singleton-site-filtered patient set — so the change in ICC_site is attributable to ComBat, not to the dropped sites. The full-cohort `hallmark` ICC is shown for context only. The same mechanical convergence exclusion (scalar R-hat>1.05 or bulk-ESS<100) is applied to the drop/combat Model-B traces.


| cohort/endpoint | ICC_site full | drop | combat | Δcombat (drop→combat) | excluded |
|---|---|---|---|---|---|
| BLCA/OS | 0.0943 | nan | nan | nan |  |
| BRCA/OS | 0.1186 | 0.1197 | 0.1293 | 0.0096 | yes |
| CESC/OS | 0.1309 | nan | nan | nan |  |
| COAD/OS | 0.1334 | nan | nan | nan |  |
| ESCA/OS | 0.2058 | 0.2169 | 0.2175 | 0.0006 |  |
| GBM/OS | 0.0189 | nan | nan | nan |  |
| HNSC/OS | 0.0781 | nan | nan | nan |  |
| KIRC/OS | 0.1786 | 0.1942 | 0.1942 | 0.0 |  |
| LGG/OS | 0.0525 | nan | nan | nan |  |
| LIHC/OS | 0.2004 | nan | nan | nan |  |
| LUAD/OS | 0.0772 | nan | nan | nan |  |
| LUSC/OS | 0.099 | 0.0981 | 0.1 | 0.0019 |  |
| MESO/OS | 0.2739 | nan | nan | nan |  |
| OV/OS | 0.1707 | nan | nan | nan |  |
| PAAD/OS | 0.1645 | nan | nan | nan |  |
| SARC/OS | 0.0422 | nan | nan | nan |  |
| STAD/OS | 0.1481 | nan | nan | nan |  |
| UCEC/OS | 0.1123 | nan | nan | nan |  |
| BLCA/PFI | 0.0417 | nan | nan | nan |  |
| BRCA/PFI | 0.0129 | 0.0127 | 0.016 | 0.0033 |  |
| CESC/PFI | 0.2069 | nan | nan | nan |  |
| COAD/PFI | 0.0956 | nan | nan | nan |  |
| ESCA/PFI | 0.2973 | 0.294 | 0.2934 | -0.0006 |  |
| GBM/PFI | 0.0444 | nan | nan | nan |  |
| HNSC/PFI | 0.0552 | nan | nan | nan |  |
| KIRC/PFI | 0.0916 | 0.094 | 0.1237 | 0.0297 |  |
| KIRP/PFI | 0.0339 | nan | nan | nan |  |
| LGG/PFI | 0.0973 | nan | nan | nan |  |
| LIHC/PFI | 0.0838 | nan | nan | nan |  |
| LUAD/PFI | 0.0173 | nan | nan | nan |  |
| LUSC/PFI | 0.0468 | 0.0461 | 0.0553 | 0.0092 |  |
| MESO/PFI | 0.1565 | nan | nan | nan |  |
| OV/PFI | 0.0635 | nan | nan | nan |  |
| PAAD/PFI | 0.0706 | nan | nan | nan |  |
| PRAD/PFI | 0.0618 | nan | nan | nan |  |
| SARC/PFI | 0.099 | nan | nan | nan |  |
| STAD/PFI | 0.1632 | nan | nan | nan |  |
| THCA/PFI | 0.2635 | nan | nan | nan |  |
| UCEC/PFI | 0.02 | nan | nan | nan |  |

Headline (converged set):
  PFI: ComBat changes ICC_site by median +0.006 (range -0.001..+0.030) across 4 converged cohorts [same-patient drop->combat].
  OS: ComBat changes ICC_site by median +0.001 (range +0.000..+0.002) across 3 converged cohorts [same-patient drop->combat].
