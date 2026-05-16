# setwd("~/TNBC/scripts")

library(dplyr)
library(tidyverse)
library(GEOquery)
library(DESeq2)

dat <- read.csv(file = "../data/GSE309617_BulkRNASeq_Human_Gene_Counts.tsv")
dim(dat)

counts <- read.delim(
  "../data/GSE309617_BulkRNASeq_Human_Gene_Counts.tsv",
  row.names = 1
)

head(counts)

colnames(counts)

counts_matrix <- counts[, -c(1,2)]

rownames(counts_matrix) <- counts$geneSym

head(counts_matrix)

sample_names <- colnames(counts_matrix)

condition <- ifelse(
  grepl("CR", sample_names),
  "Resistant",
  "Sensitive"
)

metadata <- data.frame(
  row.names = sample_names,
  condition = condition
)

metadata

table(metadata$condition)

metadata$condition <- as.factor(metadata$condition)

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(dds)

head(res)

resOrdered <- res[order(res$padj), ]

head(resOrdered)


write.csv(
  as.data.frame(resOrdered),
  "../results/DEG_results.csv"
)

dir.create("../results", showWarnings = FALSE)

write.csv(
  as.data.frame(resOrdered),
  "../results/DEG_results.tsv"
)


res_df <- as.data.frame(resOrdered)

res_df$gene <- rownames(res_df)

head(res_df)


gene_symbols <- make.unique(as.character(counts$geneSym))

counts_matrix <- counts[, -c(1,2)]

rownames(counts_matrix) <- gene_symbols


dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(dds)

resOrdered <- res[order(res$padj), ]

res_df <- as.data.frame(resOrdered)

head(res_df)


res_df[grep("DNMT", rownames(res_df)), ]

sig_genes <- subset(
  res_df,
  padj < 0.05 & abs(log2FoldChange) > 1
)

head(sig_genes)

nrow(sig_genes)

write.csv(
  sig_genes,
  "../results/significant_genes.csv"
)

library(EnhancedVolcano)
EnhancedVolcano(
  res_df,
  lab = rownames(res_df),
  x = 'log2FoldChange',
  y = 'padj'
)

EnhancedVolcano(
  res_df,
  lab = rownames(res_df),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1,
  xlim = c(-8,8),
  ylim = c(0,10)
)
res_df[grep("DNMT", rownames(res_df)), ]

sig_genes <- subset(
  res_df,
  padj < 0.05 & abs(log2FoldChange) > 1
)

nrow(sig_genes)

head(sig_genes)

res_df[grep("DNMT", rownames(res_df)), ]

write.csv(
  sig_genes,
  "../results/significant_genes.csv"
)

gene_list <- rownames(sig_genes)

head(gene_list)

library(pheatmap)
vsd <- vst(dds, blind = FALSE)
top_genes <- head(rownames(sig_genes), 50)
mat <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat)))

annotation_col <- data.frame(
  Condition = metadata$condition
)

rownames(annotation_col) <- rownames(metadata)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6
)

top_genes <- head(rownames(sig_genes), 25)
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA
)

plotPCA(vsd, intgroup = "condition")


gene_list <- rownames(sig_genes)

head(gene_list)

write.table(
  gene_list,
  "../results/gene_list.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

res_df[grep("DNMT", rownames(res_df)), ]

res_df[grep("HIST|HDAC|KMT|EZH|SMAR|BRCA|ATM|ATR", rownames(res_df)), ]
