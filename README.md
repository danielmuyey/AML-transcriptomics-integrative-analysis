# AML-transcriptomics-analysis
# Create a minimal, reusable analysis workflow template (R)
# This writes a script you can upload to Zenodo as the reproducible pipeline entry point.

AML-transcriptomics-integrative-analysis/
│
├── scripts/
│   ├── 01_qc_preprocessing.R
│   ├── 02_deseq2_all_datasets.R
│   ├── 03_cross_dataset_DEG_overlap.R      NEW
│   ├── 04_gsva_pathway_analysis.R
│   ├── 05_gsea_enrichment.R
│   ├── 06_hub_gene_network.R                STRING/Cytoscape
│   ├── 07_final_figures.R
│
├── data/
├── results/
│   ├── DEGs/
│   ├── overlaps/
│   ├── GSVA/
│   ├── GSEA/
│
├── figures/
│   ├── volcano/
│   ├── heatmaps/
│   ├── pathway/
│   ├── networks/
│
├── README.md
└── renv.lock

library(DESeq2)
library(tidyverse)

run_deseq <- function(count_file, meta, contrast_name, group1, group2) {

  counts <- read.table(count_file, header = TRUE, row.names = 1)

  coldata <- data.frame(
    row.names = colnames(counts),
    condition = meta
  )

  dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)

  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", group1, group2))
  res <- lfcShrink(dds, coef = 2, res = res)

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)

  write.csv(res_df,
            paste0("results/DEGs/DESeq2_", contrast_name, ".csv"))

  return(res_df)
}


# Dataset 1
res1 <- run_deseq(
  "GSE307476_counts.txt",
  meta = c("A","A","A","NC","NC","NC"),
  contrast_name = "GSE307476_A_vs_NC",
  group1 = "A",
  group2 = "NC"
)

# Dataset 2
res2 <- run_deseq(
  "GSE303795_counts.txt",
  meta = c("Control","Control","HF50","HF50"),
  contrast_name = "GSE303795_HF50_vs_Control",
  group1 = "HF50",
  group2 = "Control"
)

# Dataset 3
res3 <- run_deseq(
  "GSE107490_counts.txt",
  meta = c("Healthy","MDS","AML"),
  contrast_name = "GSE107490_AML_vs_Healthy",
  group1 = "AML",
  group2 = "Healthy"
)

library(VennDiagram)
library(tidyverse)

deg1 <- read.csv("results/DEGs/DESeq2_GSE307476_A_vs_NC.csv")
deg2 <- read.csv("results/DEGs/DESeq2_GSE303795_HF50_vs_Control.csv")
deg3 <- read.csv("results/DEGs/DESeq2_GSE107490_AML_vs_Healthy.csv")

sig1 <- deg1$gene[deg1$padj < 0.05 & abs(deg1$log2FoldChange) > 1]
sig2 <- deg2$gene[deg2$padj < 0.05 & abs(deg2$log2FoldChange) > 1]
sig3 <- deg3$gene[deg3$padj < 0.05 & abs(deg3$log2FoldChange) > 1]

venn.diagram(
  x = list(GSE307476 = sig1,
           GSE303795 = sig2,
           GSE107490 = sig3),
  filename = "figures/overlap/DEG_venn.png",
  fill = c("red","blue","green"),
  alpha = 0.5
)

# Core shared genes
core_genes <- Reduce(intersect, list(sig1, sig2, sig3))
write.csv(core_genes, "results/overlaps/core_AML_genes.csv")

library(GSVA)
library(GSEABase)

expr <- read.csv("normalized_expression_matrix.csv", row.names = 1)

gene_sets <- getGmt("h.all.v2024.symbols.gmt")

gsva_scores <- gsva(as.matrix(expr),
                    gene_sets,
                    method = "gsva")

write.csv(gsva_scores, "results/GSVA/gsva_scores.csv")

library(clusterProfiler)

gsea_res <- GSEA(
  geneList = sort(res1$log2FoldChange, decreasing = TRUE),
  TERM2GENE = gene_sets
)

saveRDS(gsea_res, "results/GSEA/GSE307476_GSEA.rds")

library(STRINGdb)

string_db <- STRINGdb$new(version="11.5", species=9606)

core <- read.csv("results/overlaps/core_AML_genes.csv")

mapped <- string_db$map(core, "gene", removeUnmappedRows = TRUE)

ppi <- string_db$get_interactions(mapped$STRING_id)

write.csv(ppi, "results/networks/PPI_network.csv")




