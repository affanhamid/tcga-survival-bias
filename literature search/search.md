
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
