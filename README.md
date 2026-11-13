# H5N1-infected human astrocytes: microarray differential expression in R

Re-analysis of the public GEO dataset **[GSE66597](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66597)** with `limma`. The dataset covers human astrocytes infected with H5N1 avian influenza virus and uninfected controls, measured at **6 h, 12 h and 24 h** after infection (3 replicates each, 18 arrays). The arrays are Agilent one-colour 4x44K.

The goal was to find genes that change with infection, test how sensitive the results are to analysis choices, and check whether the original paper's findings could be reproduced.

## What the analysis does

All code is in [`H5N1_astrocytes_microarray_limma.R`](H5N1_astrocytes_microarray_limma.R). The first part (up to `papermethod`) is the step-by-step exercise. The second part is the cleaned-up full pipeline.

1. **Read raw data.** Unzips the GEO raw files, builds a targets table (condition × time) and reads the Agilent green-channel intensities with `read.maimages`.
2. **Explore raw data.** Boxplots and density plots of log2 intensities show why normalisation is needed.
3. **Normalise (two variants):**
   - **A:** `normexp` background correction (offset 16) + **quantile** normalisation (the paper's approach)
   - **B:** no background correction + **cyclic loess**
4. **Quality control.** An MDS plot coloured by time and shaped by condition shows that samples separate by infection mainly at 24 h.
5. **Filter.** Removes Agilent control probes and averages replicate probes (`avereps`).
6. **Pooled differential expression.** Model `~ Time + Condition` (infection effect, adjusted for time), fitted with `lmFit` + `eBayes`, Benjamini–Hochberg FDR. Three variants are compared:
   - A + ordinary least squares
   - B + ordinary least squares
   - A + **robust** regression

   A Venn diagram shows the overlap between them, and a volcano plot marks 2-fold and p = 0.01 cut-offs.
7. **GO enrichment.** `topGO` (Biological Process, Fisher's exact test, `org.Hs.eg.db`) on genes with p < 0.05.
8. **Reproduce the paper.** Separate Infected-vs-Control contrasts at 6 h, 12 h and 24 h, using the paper's rule of **FDR < 0.05 and |log2FC| ≥ 1** (≥ 2-fold). Includes up/down counts per time point, a Venn diagram across time points, and a check against the paper's three 24 h gene tables (antiviral/interferon-stimulated genes, pro-inflammatory genes, neuroactive ligand–receptor genes).

## Key results (24 h, FDR < 0.05 and ≥ 2-fold)

- **97 probes** were differentially expressed at 24 h. Almost all were **up-regulated** (96 up, 1 down).
- The strongest responders were classic **interferon-stimulated antiviral genes**:

| Gene | log2FC | adj. P |
|---|---|---|
| RSAD2 (viperin) | 4.16 | 7.4e-07 |
| IFIT2 | 4.06 | 9.5e-08 |
| IFIT3 | 3.95 | 2.2e-07 |
| IFIT1 | 3.31 | 7.4e-07 |
| OASL | 3.20 | 2.6e-07 |
| MX1 | 2.72 | 8.9e-06 |
| CXCL10 | 1.44 | 1.0e-04 |

This agrees with the paper: astrocytes mount a strong innate antiviral (type I interferon) response that peaks at 24 h.

## Packages

`limma`, `statmod`, `topGO`, `org.Hs.eg.db`, `R.utils` (Bioconductor: `BiocManager::install(c("limma","statmod","topGO","org.Hs.eg.db"))`)

## Data

The raw data are not included. Download `GSE66597_RAW.tar` from GEO, extract it into a `raw/` folder and update `setwd()` at the top of the script.
