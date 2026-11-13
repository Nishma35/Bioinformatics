setwd("C:/Users/nisma/OneDrive/Desktop/Bioinformatics")
file.info("GSE66597_RAW.tar")$size / 1e6   # size in MB
list.files("raw")                          # what's inside the raw folder
install.packages("R.utils")
gz <- list.files("raw", pattern = "\\.txt\\.gz$", full.names = TRUE)
for (f in gz) R.utils::gunzip(f, overwrite = TRUE)
list.files("raw")

targets <- data.frame(
          FileName  = paste0(c("GSM1625995_GE1","GSM1625996_GE2","GSM1625997_GE3",
                               "GSM1625998_GE4","GSM1625999_GE5","GSM1626000_GE6",
                               "GSM1626001_GE7","GSM1626002_GE8","GSM1626003_GE9",
                               "GSM1626004_GE10","GSM1626005_GE11","GSM1626006_GE12",
                               "GSM1626007_GE13","GSM1626008_GE14","GSM1626009_GE15",
                               "GSM1626010_GE16","GSM1626011_GE17","GSM1626012_GE18"), ".txt"),
          Name      = c("con6h_r1","con6h_r2","con6h_r3","con12h_r1","con12h_r2","con12h_r3",
                        "con24h_r1","con24h_r2","con24h_r3","hm6h_r1","hm6h_r2","hm6h_r3",
                        "hm12h_r1","hm12h_r2","hm12h_r3","hm24h_r1","hm24h_r2","hm24h_r3"),
          Condition = rep(c("Control","Infected"), each = 9),
          Time      = rep(rep(c("6h","12h","24h"), each = 3), times = 2),
          stringsAsFactors = FALSE
)

write.table(targets, "Targets.txt", sep = "\t", quote = FALSE, row.names = FALSE)
read.delim("Targets.txt")   # check it prints 18 rows correctly

library(limma)

targets <- readTargets("Targets.txt")

x <- read.maimages(targets$FileName,
                   source     = "agilent",
                   green.only = TRUE,
                   path       = "raw",
                   names      = targets$Name)


grp <- factor(targets$Condition)

boxplot(log2(x$E), las = 2, col = as.integer(grp),
        main = "Raw log2 intensity per array", ylab = "log2 intensity")

dev.copy2pdf(file = "01_raw_boxplot.PDF")


plotDensities(log2(x$E), group = grp, legend = "topright",
              main = "Raw intensity densities")

dev.copy2pdf(file = "02_raw_density.PDF")


#################################normalization##############3
yA <- backgroundCorrect(x, method = "normexp", offset = 16)
yA <- normalizeBetweenArrays(yA, method = "quantile")

boxplot(yA$E, las = 2, col = as.integer(grp),
        main = "Variant A: normexp + quantile", ylab = "log2 intensity")

dev.copy2pdf(file = "03_normA_boxplot.PDF")


########################3normalizing######################
yB <- backgroundCorrect(x, method = "none")
yB <- normalizeBetweenArrays(yB, method = "cyclicloess")

boxplot(yB$E, las = 2, col = as.integer(grp),
        main = "Variant B: cyclic loess", ylab = "log2 intensity")

dev.copy2pdf(file = "04_normB_boxplot.PDF")


plotMDS(yA, labels = targets$Name, col = as.integer(grp),
        main = "MDS - Variant A (colour = condition)")

dev.copy2pdf(file = "05_mds_A.PDF")



############task4#####################3
keep <- yA$genes$ControlType == 0
yA <- yA[keep, ]
yB <- yB[keep, ]

yA <- avereps(yA, ID = yA$genes$ProbeName)
yB <- avereps(yB, ID = yB$genes$ProbeName)

dim(yA)   # fewer rows now: real probes only, duplicates merged



##################design matrix#################3
Condition <- factor(targets$Condition, levels = c("Control", "Infected"))
Time      <- factor(targets$Time,      levels = c("6h", "12h", "24h"))

design <- model.matrix(~ Time + Condition)
design

###########################################
run_de <- function(y, fit.method = "ls") {
          fit <- lmFit(y, design, method = fit.method)
          fit <- eBayes(fit)
          topTable(fit, coef = "ConditionInfected", adjust.method = "BH", number = Inf)
}

ttA        <- run_de(yA)                # Variant A: normexp+quantile, ordinary lmFit
ttB        <- run_de(yB)                # Variant B: cyclic-loess, ordinary lmFit
ttA_robust <- run_de(yA, "robust")      # Variant A prep, but robust lmFit


cat("Significant genes (adj.P.Val < 0.05):\n",
    " A  normexp+quantile / ls    :", sum(ttA$adj.P.Val        < 0.05), "\n",
    " B  cyclicloess / ls         :", sum(ttB$adj.P.Val        < 0.05), "\n",
    " A  normexp+quantile / robust:", sum(ttA_robust$adj.P.Val < 0.05), "\n")


hitsA  <- rownames(ttA)[ttA$adj.P.Val < 0.05]
hitsB  <- rownames(ttB)[ttB$adj.P.Val < 0.05]
hitsAr <- rownames(ttA_robust)[ttA_robust$adj.P.Val < 0.05]

all.ids <- union(union(hitsA, hitsB), hitsAr)
vmat <- cbind(A_lsquant = all.ids %in% hitsA,
              B_loess   = all.ids %in% hitsB,
              A_robust  = all.ids %in% hitsAr)

vennDiagram(vennCounts(vmat), main = "DEG overlap across variants")
dev.copy2pdf(file = "06_venn_variants.PDF")



fitA <- eBayes(lmFit(yA, design))
volcanoplot(fitA, coef = "ConditionInfected", highlight = 10,
            names = fitA$genes$GeneName, main = "Volcano - Variant A")
dev.copy2pdf(file = "07_volcano_A.PDF")



write.csv(ttA[ttA$adj.P.Val < 0.05, ], "DEG_variantA.csv")
############cutoffline and volocanoplot######################
volcanoplot(fitA, coef = "ConditionInfected", main = "Volcano - Variant A")
abline(v = c(-1, 1), col = "red", lty = 2)   # |log2FC| = 1  (2-fold)
abline(h = -log10(0.01), col = "blue", lty = 2)  # p = 0.01
dev.copy2pdf(file = "07_volcano_A.PDF")

#################part5###########################

BiocManager::install(c("topGO", "org.Hs.eg.db"), update = FALSE)
library(topGO)
library(org.Hs.eg.db)

# ---- Data preparation ----
# Named vector of p-values, named by gene SYMBOL. Your results table is ttA,
# your symbol column is GeneName, your p-value column is P.Value.
geneList_full <- ttA$P.Value
names(geneList_full) <- ttA$GeneName
geneList_full <- geneList_full[!is.na(names(geneList_full)) & names(geneList_full) != ""]
geneList_full <- tapply(geneList_full, names(geneList_full), min)  # keep smallest p per symbol

# "interesting" genes = significant at p < 0.05
topDiffGenes <- function(allScore) { return(allScore < 0.05) }
sum(topDiffGenes(geneList_full))   # how many genes selected

# ---- Build the topGOdata object (human = org.Hs.eg.db) ----
GOdata <- new("topGOdata",
              description = "GO enrichment - H5N1 infected vs control astrocytes",
              ontology    = "BP",
              allGenes    = geneList_full,
              geneSel     = topDiffGenes,
              nodeSize    = 5,
              annot       = annFUN.org,
              mapping     = "org.Hs.eg.db",
              ID          = "symbol")
GOdata

# ---- Run Fisher's exact test (classic algorithm) ----
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher

# ---- Summarise top terms ----
goResults <- GenTable(GOdata, classicFisher = resultFisher,
                      orderBy = "classicFisher", topNodes = 15)
goResults

write.csv(goResults, file = "GO_enrichment_results.csv", row.names = FALSE)



##############papermethod#########################


## ===========================================================================
## GSE66597 — H5N1-infected vs control human astrocytes (Agilent 1-colour 4x44K)
## Full pipeline: read -> explore -> normalise -> DE (pooled) -> compare
## variants -> GO enrichment -> reproduce paper (per time point + 2-fold).
## Save any on-screen plot with: dev.copy2pdf(file = "name.PDF")
## ===========================================================================

setwd("C:/Users/nisma/OneDrive/Desktop/Bioinformatics")

## ---- packages (run once) --------------------------------------------------
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("limma","statmod","topGO","org.Hs.eg.db"), update=FALSE)
# install.packages("R.utils")

library(limma)

## ===========================================================================
## 0. UNZIP RAW FILES  (only needed once, if still .txt.gz)
## ===========================================================================
gz <- list.files("raw", pattern = "\\.txt\\.gz$", full.names = TRUE)
for (f in gz) R.utils::gunzip(f, overwrite = TRUE)
list.files("raw")

## ===========================================================================
## 1. TARGETS FILE  (condition/time from the GEO sample titles)
## ===========================================================================
targets <- data.frame(
          FileName  = paste0(c("GSM1625995_GE1","GSM1625996_GE2","GSM1625997_GE3",
                               "GSM1625998_GE4","GSM1625999_GE5","GSM1626000_GE6",
                               "GSM1626001_GE7","GSM1626002_GE8","GSM1626003_GE9",
                               "GSM1626004_GE10","GSM1626005_GE11","GSM1626006_GE12",
                               "GSM1626007_GE13","GSM1626008_GE14","GSM1626009_GE15",
                               "GSM1626010_GE16","GSM1626011_GE17","GSM1626012_GE18"), ".txt"),
          Name      = c("con6h_r1","con6h_r2","con6h_r3","con12h_r1","con12h_r2","con12h_r3",
                        "con24h_r1","con24h_r2","con24h_r3","hm6h_r1","hm6h_r2","hm6h_r3",
                        "hm12h_r1","hm12h_r2","hm12h_r3","hm24h_r1","hm24h_r2","hm24h_r3"),
          Condition = rep(c("Control","Infected"), each = 9),
          Time      = rep(rep(c("6h","12h","24h"), each = 3), times = 2),
          stringsAsFactors = FALSE
)
write.table(targets, "Targets.txt", sep = "\t", quote = FALSE, row.names = FALSE)
targets <- readTargets("Targets.txt")

## ===========================================================================
## 2. READ RAW DATA  (one-colour Agilent)
## ===========================================================================
x <- read.maimages(targets$FileName,
                   source     = "agilent",
                   green.only = TRUE,
                   path       = "raw",
                   names      = targets$Name)
dim(x)                                   # ~45015 probes x 18 arrays

grp <- factor(targets$Condition)

## ===========================================================================
## 3. EXPLORE RAW DATA  (justifies normalising)
## ===========================================================================


plotDensities(log2(x$E), group = grp, legend = "topright",
              main = "Raw intensity densities")
dev.copy2pdf(file = "02_raw_density.PDF")
grpcol <- c(Control = "steelblue", Infected = "firebrick")[targets$Condition]

boxplot(log2(x$E), las = 2, col = grpcol,
        main = "Raw log2 intensity per array", ylab = "log2 intensity")
dev.copy2pdf(file = "01_raw_boxplot.PDF")

## ===========================================================================
## 4. NORMALISATION  (two variants)
## ===========================================================================
## Variant A: normexp background + quantile  (matches the paper's quantile step)
yA <- backgroundCorrect(x, method = "normexp", offset = 16)
yA <- normalizeBetweenArrays(yA, method = "quantile")

grpcol <- c(Control = "steelblue", Infected = "firebrick")[targets$Condition]

boxplot(yA$E, las = 2, col = grpcol,
        main = "Variant A: normexp + quantile", ylab = "log2 intensity")
dev.copy2pdf(file = "03_normA_boxplot.PDF")
grpcol <- c(Control = "steelblue", Infected = "firebrick")[targets$Condition]

plotDensities(yA$E, col = grpcol, legend = "topright",
              main = "Densities after normalization (Variant A)")
dev.copy2pdf(file = "09_normA_density.PDF")


## Variant B: cyclic-loess, no separate background step
yB <- backgroundCorrect(x, method = "none")
yB <- normalizeBetweenArrays(yB, method = "cyclicloess")
grpcol <- c(Control = "steelblue", Infected = "firebrick")[targets$Condition]

boxplot(yB$E, las = 2, col = grpcol,
        main = "Variant B: cyclic loess", ylab = "log2 intensity")
dev.copy2pdf(file = "04_normB_boxplot.PDF")

plotDensities(yB$E, col = grpcol, legend = "topright",
              main = "Densities after normalization (Variant B)")
dev.copy2pdf(file = "10_normB_density.PDF")

## QC: MDS (samples separate by condition mostly at 24h -> time-dependent)
# colour by TIME, symbol by CONDITION
timecol <- c("6h" = "#1b9e77", "12h" = "#d95f02", "24h" = "#7570b3")[targets$Time]
condpch <- c(Control = 1, Infected = 19)[targets$Condition]   # open vs filled

plotMDS(yA,
        pch = condpch, col = timecol, cex = 2,
        main = "MDS - Variant A")

legend("topleft",
       legend = c("6h","12h","24h","Control (open)","Infected (filled)"),
       col = c("#1b9e77","#d95f02","#7570b3","black","black"),
       pch = c(15,15,15,1,19), pt.cex = 1.4, bty = "n")

dev.copy2pdf(file = "05_mds_A_better.PDF")
dev.copy2pdf(file = "05_mds_A.PDF")

## ===========================================================================
## 5. FILTER controls + collapse replicate probes
## ===========================================================================
keep <- yA$genes$ControlType == 0
yA <- yA[keep, ];  yB <- yB[keep, ]
yA <- avereps(yA, ID = yA$genes$ProbeName)
yB <- avereps(yB, ID = yB$genes$ProbeName)
dim(yA)                                   # ~41000 probes

## ===========================================================================
## 6. DIFFERENTIAL EXPRESSION — pooled (Infected vs Control, blocking on Time)
##    + compare 3 analysis variants
## ===========================================================================
Condition <- factor(targets$Condition, levels = c("Control","Infected"))
Time      <- factor(targets$Time,      levels = c("6h","12h","24h"))
design    <- model.matrix(~ Time + Condition)

run_de <- function(y, fit.method = "ls") {
          fit <- lmFit(y, design, method = fit.method)
          fit <- eBayes(fit)
          topTable(fit, coef = "ConditionInfected", adjust.method = "BH", number = Inf)
}
ttA        <- run_de(yA)                  # A: normexp+quantile, ls
ttB        <- run_de(yB)                  # B: cyclicloess, ls
ttA_robust <- run_de(yA, "robust")        # A prep, robust lmFit

cat("Significant genes (adj.P.Val < 0.05):\n",
    " A normexp+quantile / ls    :", sum(ttA$adj.P.Val        < 0.05), "\n",
    " B cyclicloess / ls         :", sum(ttB$adj.P.Val        < 0.05), "\n",
    " A normexp+quantile / robust:", sum(ttA_robust$adj.P.Val < 0.05), "\n")

## Venn of the three variants
hitsA <- rownames(ttA)[ttA$adj.P.Val<0.05]
hitsB <- rownames(ttB)[ttB$adj.P.Val<0.05]
hitsAr<- rownames(ttA_robust)[ttA_robust$adj.P.Val<0.05]
allid <- union(union(hitsA,hitsB),hitsAr)
vmat  <- cbind(A_lsquant=allid%in%hitsA, B_loess=allid%in%hitsB, A_robust=allid%in%hitsAr)
vennDiagram(vennCounts(vmat), main = "DEG overlap across variants")
dev.copy2pdf(file = "06_venn_variants.PDF")

## Volcano (Variant A) with 2-fold and p=0.01 cutoff lines
fitA <- eBayes(lmFit(yA, design))
volcanoplot(fitA, coef = "ConditionInfected", main = "Volcano - Variant A")
abline(v = c(-1, 1), col = "red",  lty = 2)      # |log2FC| = 1 (2-fold)
abline(h = -log10(0.01), col = "blue", lty = 2)  # p = 0.01
dev.copy2pdf(file = "07_volcano_A.PDF")

write.csv(ttA[ttA$adj.P.Val < 0.05, ], "DEG_variantA.csv")

## ===========================================================================
## 7. GO ENRICHMENT (bonus, topGO)
## ===========================================================================
library(topGO); library(org.Hs.eg.db)

geneList_full <- ttA$P.Value
names(geneList_full) <- ttA$GeneName
geneList_full <- geneList_full[!is.na(names(geneList_full)) & names(geneList_full) != ""]
geneList_full <- tapply(geneList_full, names(geneList_full), min)

topDiffGenes <- function(allScore) allScore < 0.05
GOdata <- new("topGOdata",
              description = "GO enrichment - H5N1 infected vs control astrocytes",
              ontology = "BP", allGenes = geneList_full, geneSel = topDiffGenes,
              nodeSize = 5, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
goResults <- GenTable(GOdata, classicFisher = resultFisher,
                      orderBy = "classicFisher", topNodes = 15)
goResults
write.csv(goResults, "GO_enrichment_results.csv", row.names = FALSE)

## ===========================================================================
## 8. REPRODUCE THE PAPER — per time point + 2-fold threshold
##    Paper rule: FDR < 0.05 AND |log2FC| >= 1 (fold change >= 2), each time
##    point analysed separately (Infected - Control within 6h/12h/24h).
## ===========================================================================
grp2 <- factor(paste(targets$Condition, targets$Time, sep = "_"),
               levels = c("Control_6h","Infected_6h",
                          "Control_12h","Infected_12h",
                          "Control_24h","Infected_24h"))
design2 <- model.matrix(~ 0 + grp2); colnames(design2) <- levels(grp2)

cm <- makeContrasts(
          h6  = Infected_6h  - Control_6h,
          h12 = Infected_12h - Control_12h,
          h24 = Infected_24h - Control_24h,
          levels = design2)

fit2 <- lmFit(yA, design2)
fit2 <- contrasts.fit(fit2, cm)
fit2 <- eBayes(fit2)

count_deg <- function(coef) {
          tt  <- topTable(fit2, coef = coef, number = Inf, adjust.method = "BH")
          sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1, ]
          c(Up = sum(sig$logFC > 0), Down = sum(sig$logFC < 0), Total = nrow(sig))
}
paper_table <- rbind(`6h`  = count_deg("h6"),
                     `12h` = count_deg("h12"),
                     `24h` = count_deg("h24"))
paper_table    # compare to paper Fig 2A: 192 / 233 / 721

## 24h list (their headline) + check the late chemokines return
tt24  <- topTable(fit2, coef = "h24", number = Inf, adjust.method = "BH")
deg24 <- tt24[tt24$adj.P.Val < 0.05 & abs(tt24$logFC) >= 1, ]
write.csv(deg24, "DEG_24h_2fold.csv")
deg24[deg24$GeneName %in% c("CXCL10","CXCL11","RSAD2","IFIT1","MX1"),
      c("GeneName","logFC","adj.P.Val")]

## Venn of the three time points (reproduces paper Fig 2B)
getSig <- function(coef){
          tt <- topTable(fit2, coef=coef, number=Inf)
          rownames(tt)[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1]
}
s6 <- getSig("h6"); s12 <- getSig("h12"); s24 <- getSig("h24")
ids <- union(union(s6,s12),s24)
vm  <- cbind(`6h`=ids%in%s6, `12h`=ids%in%s12, `24h`=ids%in%s24)
vennDiagram(vennCounts(vm), main = "DE genes per time point (2-fold)")
dev.copy2pdf(file = "08_venn_timepoints.PDF")






count_raw <- function(coef) {
          tt <- topTable(fit2, coef = coef, number = Inf)
          sig <- tt[tt$P.Value < 0.05 & abs(tt$logFC) >= 1, ]   # raw p, not adj
          c(Up = sum(sig$logFC>0), Down = sum(sig$logFC<0), Total = nrow(sig))
}
rbind(`6h`=count_raw("h6"), `12h`=count_raw("h12"), `24h`=count_raw("h24"))


tt24 <- topTable(fit2, coef="h24", number=Inf)
tt24[tt24$GeneName %in% c("CXCL10","CXCL11","RSAD2","IFIT1","MX1"),
     c("GeneName","logFC","adj.P.Val")]


## ---- Compare your 24h results to the paper's three gene tables ----
tt24 <- topTable(fit2, coef = "h24", number = Inf)

check_genes <- function(genes, label) {
          sub <- tt24[tt24$GeneName %in% genes,
                      c("GeneName","logFC","adj.P.Val")]
          # keep the strongest (smallest adj.P) row per gene
          sub <- sub[order(sub$adj.P.Val), ]
          sub <- sub[!duplicated(sub$GeneName), ]
          sub$hit <- ifelse(sub$adj.P.Val < 0.05 & abs(sub$logFC) >= 1, "YES", "no")
          sub$dir <- ifelse(sub$logFC > 0, "UP", "DOWN")
          cat("\n====", label, "====\n")
          cat("Found", nrow(sub), "of", length(genes),
              "| significant (FDR<0.05 & >=2fold):", sum(sub$hit=="YES"), "\n")
          print(sub, row.names = FALSE)
}

## Table 1 - top 20 up-regulated at 24h (antiviral / ISGs)
t1 <- c("RSAD2","IFIT2","IFIT3","CXCL10","OASL","MX2","BATF2","CXCL11","IFIT1",
        "DDX58","OAS2","IFI44L","RTP4","MX1","FOSB","OAS1","VPREB1","BST2",
        "SAMD9","USP18")

## Table 2 - pro-inflammatory (24h)
t2 <- c("IL17A","INHBC","IL12RB1","CXCL9","CCL5","CXCL10","CCL4","CCR4","IFNB1",
        "LEPR","TNFSF15","IL10RA","IL29","CD274","CTF1","CXCL11","TNFSF10",
        "CCL11","CSF3R")

## Table 3 - neuroactive ligand-receptor (24h)
t3 <- c("HTR1A","GHRHR","CHRNG","GRM7","SSTR2","CGA","P2RY13","OPRM1","GABRA1",
        "MC3R","CNR2","LEPR","CCKAR","HRH2","MCHR2","MC2R","HTR1B","NPB",
        "MLNR","TAAR8")

check_genes(t1, "Table 1: antiviral / ISGs")
check_genes(t2, "Table 2: pro-inflammatory")
check_genes(t3, "Table 3: neuroactive ligand-receptor")



