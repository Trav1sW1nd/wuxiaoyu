# ============================================================
# Script: 03_GSEA_custom_checkpoint.R
# Purpose: GSEA with a custom immune checkpoint ligand gene set (Fig. 3F)
# Requires: ranked_genes from 02_GSEA_Hallmark.R
# ============================================================

library(fgsea)
library(ggplot2)

# 1. Define custom immune checkpoint ligand gene set
# (B7 family, checkpoint ligands and co-stimulatory molecules)
ic_ligands <- c("Cd274", "Pdcd1lg2", "Icosl", "Vsir", "Cd80", "Cd86",
                "Fgl1", "Tnfsf4", "Tnfsf9", "Cd70", "Hhla2",
                "Cd160", "Nectin2", "Pvr", "Btn2a2")

# Keep only genes present in the ranked list 
ic_ligands <- intersect(ic_ligands, names(ranked_genes))

# 2. Run fgsea 
set.seed(123)
gsea_custom <- fgsea(pathways = custom_list,
                     stats    = ranked_genes,
                     minSize  = 5,
                     nperm    = 10000)

# 3. Enrichment plot 
p_custom <- plotEnrichment(custom_list[["Immune_Checkpoint_Ligands"]], ranked_genes) +
  labs(title = "Immune Checkpoint Ligands",
       x = "Rank", y = "Enrichment Score") +
  theme_minimal(base_size = 12)

