############################################################
# AML Transcriptomics Analysis Pipeline
# Reproducible R workflow for differential expression
# and visualization of public RNA-seq datasets
#
# Author: Daniel Muteb Muyey
# Repository: AML-transcriptomics-analysis
# Purpose: Minimal, reusable pipeline for GitHub + Zenodo
############################################################
AML-transcriptomics-integrative-analysis/
│
├── R/
│   ├── 01_load_data.R
│   ├── 02_preprocess.R
│   ├── 03_DE_analysis.R
│   ├── 04_ChIP_mapping.R
│   ├── 05_PCA_analysis.R
│   ├── 06_integration.R
│   ├── 07_visualization.R
│
├── data/
├── results/
├── figures/
├── README.md
└── run_pipeline.R

#1) MASTER PIPELINE (run_pipeline.R)
############################################################
# AML Multi-Omics Integration Pipeline (v1.0)
############################################################

source("R/01_load_data.R")
source("R/02_preprocess.R")
source("R/03_DE_analysis.R")
source("R/04_ChIP_mapping.R")
source("R/05_PCA_analysis.R")
source("R/06_integration.R")
source("R/07_visualization.R")

#2) DATA LOADING (01_load_data.R)
library(readr)
library(readxl)
library(dplyr)

# -------------------------
# RNA-seq datasets
# -------------------------
gse303795 <- read.table("data/GSE303795.txt", header=TRUE, sep="\t")
gse273642 <- read.table("data/GSE273642.txt", header=TRUE, sep="\t")

# -------------------------
# ChIP-seq dataset
# -------------------------
gse240526 <- read.table("data/GSE240526_peaks.txt", header=TRUE, sep="\t")

# -------------------------
# Ferroptosis gene set
# -------------------------
ferroptosis_genes <- read.csv("data/ferroptosis_genes.csv")

#3) PREPROCESSING (02_preprocess.R)
library(dplyr)

# Gene harmonization function
harmonize_genes <- function(df){
  df$gene <- toupper(df$gene)
  df <- df %>% filter(!is.na(gene))
  return(df)
}

gse303795 <- harmonize_genes(gse303795)
gse273642 <- harmonize_genes(gse273642)

# Log2 transformation (RNA-seq)
log_transform <- function(df){
  expr <- df[, -1]
  df[, -1] <- log2(expr + 0.5)
  return(df)
}

gse303795 <- log_transform(gse303795)
gse273642 <- log_transform(gse273642)

#4) DIFFERENTIAL EXPRESSION (03_DE_analysis.R)
safe_ttest <- function(x, group){
  x1 <- x[group == levels(group)[1]]
  x2 <- x[group == levels(group)[2]]

  if(sd(x1)==0 & sd(x2)==0) return(NA_real_)
  tryCatch(t.test(x1, x2)$p.value, error=function(e) NA_real_)
}

run_DE <- function(expr, group){
  pvals <- apply(expr, 1, safe_ttest, group=group)

  logFC <- rowMeans(expr[, group==levels(group)[1]]) -
           rowMeans(expr[, group==levels(group)[2]])

  FDR <- p.adjust(pvals, method="BH")

  data.frame(
    gene = rownames(expr),
    logFC = logFC,
    pval = pvals,
    FDR = FDR
  )
}

deg_303795 <- run_DE(gse303795[, -1], factor(c("Ctrl","Ctrl","HF","HF")))
deg_273642 <- run_DE(gse273642[, -1], factor(c("LSC","HSC")))

#5) CHROMATIN MAPPING (04_ChIP_mapping.R)
# Map peaks → nearest gene
map_peaks_to_genes <- function(peaks){
  peaks$gene <- toupper(peaks$nearest_gene)
  return(peaks)
}

chip_mapped <- map_peaks_to_genes(gse240526)

chip_genes <- unique(chip_mapped$gene)

#6) PCA ANALYSIS (05_PCA_analysis.R)
library(FactoMineR)

run_PCA <- function(expr){
  pca <- prcomp(t(expr), scale.=TRUE)
  return(pca)
}

pca_303795 <- run_PCA(gse303795[, -1])

#7) INTEGRATION CORE (06_integration.R)
# -------------------------
# Gene set intersections
# -------------------------

deg_genes_303 <- deg_303795$gene[deg_303795$FDR < 0.05]
deg_genes_273 <- deg_273642$gene[deg_273642$FDR < 0.05]

chip_genes <- chip_mapped$gene

ferroptosis <- ferroptosis_genes$gene

# -------------------------
# Cross-dataset overlap
# -------------------------

overlap_DE_chip <- intersect(deg_genes_303, chip_genes)
overlap_DE_ferro <- intersect(deg_genes_303, ferroptosis)
overlap_all <- Reduce(intersect, list(deg_genes_303,
                                      chip_genes,
                                      ferroptosis,
                                      deg_genes_273))

integration_table <- data.frame(
  gene = overlap_all
)

#8) VISUALIZATION (07_visualization.R)
library(ggplot2)

# Volcano plot example
ggplot(deg_303795, aes(logFC, -log10(pval))) +
  geom_point(alpha=0.5) +
  theme_minimal()

# Overlap summary plot
barplot(c(
  DE_Chip = length(overlap_DE_chip),
  DE_Ferro = length(overlap_DE_ferro),
  ALL = length(overlap_all)
))
#_gsva_pathway_analysis.R

############################################################
# AML Multi-Omics Integration Pipeline (v1.0)
# Reproducible transcriptomics + epigenomics + pathway + PPI
############################################################

#############################
# 1. LIBRARIES
#############################

library(readr)
library(readxl)
library(dplyr)
library(tidyverse)

library(VennDiagram)
library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(STRINGdb)

#############################
# 2. DATA LOADING
#############################

# RNA-seq datasets
gse303795 <- read.table("data/GSE303795.txt", header=TRUE, sep="\t")
gse273642 <- read.table("data/GSE273642.txt", header=TRUE, sep="\t")

# ChIP-seq dataset
gse240526 <- read.table("data/GSE240526_peaks.txt", header=TRUE, sep="\t")

# Ferroptosis gene set
ferroptosis_genes <- read.csv("data/ferroptosis_genes.csv")

#############################
# 3. DEG INPUTS (precomputed DESeq2 results)
#############################

deg1 <- read.csv("results/DEGs/DESeq2_GSE307476_A_vs_NC.csv")
deg2 <- read.csv("results/DEGs/DESeq2_GSE303795_HF50_vs_Control.csv")
deg3 <- read.csv("results/DEGs/DESeq2_GSE107490_AML_vs_Healthy.csv")

#############################
# 4. SIGNIFICANT GENE FILTER
#############################

get_sig_genes <- function(df){
  df$gene[df$padj < 0.05 & abs(df$log2FoldChange) > 1]
}

sig1 <- get_sig_genes(deg1)
sig2 <- get_sig_genes(deg2)
sig3 <- get_sig_genes(deg3)

#############################
# 5. VENN OVERLAP ANALYSIS
#############################

dir.create("figures/overlap", recursive = TRUE, showWarnings = FALSE)
dir.create("results/overlaps", recursive = TRUE, showWarnings = FALSE)

venn.diagram(
  x = list(
    GSE307476 = sig1,
    GSE303795 = sig2,
    GSE107490 = sig3
  ),
  filename = "figures/overlap/DEG_venn.png",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2
)

# Core intersection genes
core_genes <- Reduce(intersect, list(sig1, sig2, sig3))

write.csv(core_genes,
          "results/overlaps/core_AML_genes.csv",
          row.names = FALSE)

#############################
# 6. GSVA PATHWAY ACTIVITY
#############################

expr <- read.csv("normalized_expression_matrix.csv", row.names = 1)

gene_sets <- getGmt("h.all.v2024.symbols.gmt")

gsva_scores <- gsva(
  as.matrix(expr),
  gene_sets,
  method = "gsva"
)

dir.create("results/GSVA", recursive = TRUE, showWarnings = FALSE)
write.csv(gsva_scores, "results/GSVA/gsva_scores.csv")

#############################
# 7. GSEA ANALYSIS
#############################

geneList <- deg1$log2FoldChange
names(geneList) <- deg1$gene
geneList <- sort(geneList, decreasing = TRUE)

gsea_res <- GSEA(
  geneList = geneList,
  TERM2GENE = gene_sets
)

dir.create("results/GSEA", recursive = TRUE, showWarnings = FALSE)
saveRDS(gsea_res, "results/GSEA/GSE307476_GSEA.rds")

#############################
# 8. CHIP-seq → GENE MAPPING (simplified)
#############################

# If peaks already annotated:
gse240526$gene <- toupper(gse240526$nearest_gene)

chip_genes <- unique(na.omit(gse240526$gene))

#############################
# 9. FERROPTOSIS OVERLAP ANALYSIS
#############################

ferro_genes <- toupper(ferroptosis_genes$gene)

overlap_DE_chip   <- intersect(sig1, chip_genes)
overlap_DE_ferro  <- intersect(sig1, ferro_genes)
overlap_core_all  <- Reduce(intersect, list(sig1, chip_genes, ferro_genes))

integration_table <- data.frame(
  gene = overlap_core_all
)

write.csv(integration_table,
          "results/overlaps/integrated_core_genes.csv",
          row.names = FALSE)

#############################
# 10. PPI NETWORK (STRING)
#############################

string_db <- STRINGdb$new(version="11.5", species=9606)

core <- read.csv("results/overlaps/core_AML_genes.csv")

mapped <- string_db$map(core,
                        "gene",
                        removeUnmappedRows = TRUE)

ppi <- string_db$get_interactions(mapped$STRING_id)

dir.create("results/networks", recursive = TRUE, showWarnings = FALSE)
write.csv(ppi,
          "results/networks/PPI_network.csv",
          row.names = FALSE)

#############################
# 11. FINAL SUMMARY OUTPUT
#############################

cat("Pipeline completed successfully!\n")
cat("Core genes:", length(core_genes), "\n")
cat("Integrated genes:", length(overlap_core_all), "\n")

[![DOI](https://zenodo.org/badge/1261440380.svg)](https://doi.org/10.5281/zenodo.20573108)





