# ============================================================
# Script: 04_heatmap.R
# Purpose: Heatmap of immune checkpoint ligands and cytokines (Fig. 3G)
# Input:   gene_expression.xls (TPM matrix from BGI)
# ============================================================

library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)

# 1. Read expression matrix
expr <- read.table("data/gene_expression.xls",
                   header = TRUE, sep = "\t", quote = "",
                   comment.char = "", stringsAsFactors = FALSE)

# IMPORTANT: In this dataset, "N" samples are shCsf3 and "S" samples are NC
target_samples <- c("tpm_N1", "tpm_N2", "tpm_N3",   
                    "tpm_S1", "tpm_S2", "tpm_S3")    
expr_tpm <- expr %>%
  select(gene_symbol, all_of(target_samples)) %>%
  column_to_rownames("gene_symbol")

# Rename columns to real group labels
colnames(expr_tpm) <- c("shCsf3_1", "shCsf3_2", "shCsf3_3",
                        "NC_1", "NC_2", "NC_3")

sample_info <- data.frame(
  Group = factor(c(rep("shCsf3", 3), rep("NC", 3)),
                 levels = c("NC", "shCsf3")),
  row.names = colnames(expr_tpm))

# 2. Focus genes 
focus_genes <- c("Cd274", "Pdcd1lg2", "Icosl", "Vsir", "Cd80", "Cd86",
                 "Fgl1", "Tnfsf4", "Tnfsf9", "Cd70", "Hhla2", "B7h4",
                 "B7h3", "Cd160", "Nectin2", "Pvr", "Btn2a2",
                 "Cxcl12", "Ccl2", "Ccl5", "Cxcl1", "Cxcl2", "Cxcl10",
                 "Il6", "Il10", "Tgfb1", "Csf1", "Ccl22", "Ccl17", "Csf3")
genes_heat <- intersect(focus_genes, rownames(expr_tpm))

# 3. Log2 transform and row z-score 
mat <- as.matrix(expr_tpm[genes_heat, ])
mat_scaled <- t(scale(t(log2(mat + 1))))

# Remove genes with non-finite values 
mat_scaled <- mat_scaled[apply(mat_scaled, 1,
                               function(x) all(is.finite(x))), ]

# 4. Draw heatmap 
png("Fig3G_heatmap.png", width = 11, height = 6, units = "in", res = 300)
pheatmap(mat_scaled,
         cluster_rows     = TRUE,
         cluster_cols     = FALSE,
         annotation_col   = sample_info["Group"],
         annotation_colors = list(Group = c(NC = "#90bff9", shCsf3 = "#ff8080")),
         color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
         main  = "Cytokines & Immune checkpoint ligands",
         fontsize = 12, angle_col = 0, border_color = NA)
dev.off()
