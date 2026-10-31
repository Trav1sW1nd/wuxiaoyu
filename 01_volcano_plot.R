# ============================================================
# Script: 01_volcano_plot.R
# Purpose: Differential expression volcano plot (Fig. 3A)
# Input:   gene_diff.xls (DESeq2 output from BGI)
# Output:  volcano plot
# ============================================================

# Load packages 
library(dplyr)
library(ggplot2)
library(ggrepel)

# Set input directory 
data_dir <- "data"

# 1. Read differential expression table 
#    this .xls is actually a tab-delimited text file 
diff <- read.table(file.path(data_dir, "gene_diff.xls"),
                   header = TRUE, sep = "\t", quote = "",
                   comment.char = "")

# 2. Align log2FC direction
# The original BGI output reports log2FC as NC vs shCsf3 (NC as reference).
# We invert the sign so that positive values indicate up-regulation in shCsf3.
diff_corrected <- diff %>%
  select(gene      = gene_symbol,
         logFC_raw = diffexp_log2fc_NC.vs.SH,
         pvalue    = diffexp_deseq2_pvalue_NC.vs.SH,
         padj      = diffexp_deseq2_qvalue_NC.vs.SH) %>%
  filter(!is.na(logFC_raw)) %>%
  mutate(logFC = -logFC_raw) %>%    
  arrange(desc(logFC))

# 3. Focus genes: immune checkpoints, chemokines and cytokines
focus_genes <- c("Cd274", "Pdcd1lg2", "Icosl", "Vsir", "Cd80", "Cd86",
                 "Fgl1", "Tnfsf4", "Tnfsf9", "Cd70", "Hhla2", "B7h4",
                 "B7h3", "Cd160", "Nectin2", "Pvr", "Btn2a2",
                 "Cxcl12", "Ccl2", "Ccl5", "Cxcl1", "Cxcl2", "Cxcl10",
                 "Il6", "Il10", "Tgfb1", "Csf1", "Ccl22", "Ccl17", "Csf3")

# 4. Assign color groups and labels 
diff_plot <- diff_corrected %>%
  mutate(label = ifelse(gene %in% focus_genes, gene, NA),
         color_group = case_when(
           logFC < -1 & pvalue < 0.05 ~ "Down",
           logFC >  1 & pvalue < 0.05 ~ "Up",
           TRUE ~ "NS"))

color_vec <- c(Down = "#ff8080", Up = "#90bff9", NS = "grey70")

# 5. Volcano plot 
p_volcano <- ggplot(diff_plot, aes(x = logFC, y = -log10(pvalue))) +
  geom_point(aes(color = color_group), size = 1.5, alpha = 0.6) +
  scale_color_manual(values = color_vec) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = subset(diff_plot, !is.na(label) & pvalue < 0.05 & abs(logFC) > 1),
    aes(label = gene), size = 5.25, max.overlaps = 20, box.padding = 0.5) +
  labs(x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# 6. Save figure 
