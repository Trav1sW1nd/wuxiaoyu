# ============================================================
# Script: 05_TF_activity_fgsea_dorothea.R
# Purpose: Transcription factor activity inferred by fgsea on Dorothea
#          regulons; top 20 TFs with reduced activity (Fig. 3J)
# Requires: ranked_genes from 02_GSEA_Hallmark.R
# ============================================================

# Load packages 
library(dplyr)
library(dorothea)
library(fgsea)
library(ggplot2)

# 1. Load Dorothea mouse regulons 
data(dorothea_mm, package = "dorothea")

# 2. Keep confidence levels A-C 
regulons_df <- dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C"))

# 3. Build a named list: TF -> target genes
all_regulons <- split(regulons_df$target, regulons_df$tf)

# 4. Deduplicate and keep TFs with >= 5 targets
all_regulons <- lapply(all_regulons, unique)
all_regulons <- all_regulons[lengths(all_regulons) >= 5]

# 5. Run fgsea 
set.seed(123)
gsea_all_tf <- fgseaMultilevel(
  pathways  = all_regulons,
  stats     = ranked_genes,
  minSize   = 5,
  maxSize   = 500,
  scoreType = "std")

# 6. Extract top 20 TFs with most reduced activity
top_20_down <- gsea_all_tf %>%
  arrange(NES) %>%
  slice_head(n = 20) %>%
  select(tf = pathway, NES, pval, padj, size)

# 7. Bar plot 
p_all <- ggplot(top_20_down,
                aes(x = reorder(tf, -NES), y = NES,
                    fill = ifelse(padj < 0.05, "FDR < 0.05", "NS"))) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("FDR < 0.05" = "#ff8080", "NS" = "grey70")) +
  coord_flip() +
  labs(x = "",
       y = "Normalized Enrichment Score (NES)",
       title = "Top 20 transcription factors with reduced activity upon CSF3 knockdown",
       subtitle = "Inferred by fgsea on Dorothea regulons (confidence A/B/C)") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 20),
        axis.text  = element_text(size = 16),
        legend.position = "top")

# 8. Save figure and result table 

