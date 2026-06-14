
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
