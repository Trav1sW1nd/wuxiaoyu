# ============================================================
# Script: 02_GSEA_Hallmark.R
# Purpose: GSEA with Hallmark gene sets
#          - bubble plot (Fig. 3B)
#          - enrichment plots (Fig. 3C-E)
# Requires: diff_corrected from 01_volcano_plot.R
# ============================================================

library(dplyr)
library(ggplot2)
library(fgsea)

# 1. Build ranked gene list 
#    ranked by descending logFC
ranked_genes <- diff_corrected %>%
  filter(!is.na(logFC)) %>%
  arrange(desc(logFC)) %>%
  tibble::deframe()

# 2. Load Hallmark gene sets 
#    GMT file can be downloaded from MSigDB 
gmt_file <- "data/mh.all.v2023.2.Mm.symbols.gmt"
hallmark_list <- gmtPathways(gmt_file)

# 3. Run fgsea 
set.seed(123)
gsea_res <- fgsea(pathways = hallmark_list,
                  stats    = ranked_genes,
                  minSize  = 5,
                  maxSize  = 500,
                  nperm    = 10000)

# 4. Bubble plot of significant pathways 
bubble_df <- gsea_res %>%
  as_tibble() %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES)) %>%
  mutate(pathway_name   = gsub("HALLMARK_", "", pathway),
         neg_log10_padj = -log10(padj)) %>%
  slice_head(n = 20)          # Top 20 pathways 

p_bubble <- ggplot(bubble_df, aes(x = NES, y = reorder(pathway_name, NES))) +
  geom_point(aes(size = size, color = neg_log10_padj), alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
  scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_line(color = "grey90"))

ggsave("Fig3B_GSEA_bubble.pdf", plot = p_bubble, width = 10, height = 7)

# 5. Enrichment plots for selected pathways
selected_pathways <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                       "HALLMARK_INFLAMMATORY_RESPONSE",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE")

for (pw in selected_pathways) {
  p <- plotEnrichment(hallmark_list[[pw]], ranked_genes) +
    labs(title = pw, x = "Rank", y = "Enrichment Score") +
    theme_minimal(base_size = 12)

  short_name <- gsub("HALLMARK_", "", pw)
  ggsave(paste0("Fig3_", short_name, ".pdf"), plot = p, width = 8, height = 6)
}

# 6. Save full GSEA results 

