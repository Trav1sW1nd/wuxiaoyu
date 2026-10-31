# ============================================================
# Script: 07_TF_target_network.R
# Purpose: Build and visualize TF-target gene regulatory network (Fig. 4A)
# Input:   dorothea_mm (Dorothea mouse regulon)
# Output:  Fig4A_TF_target_network.png
# ============================================================

# Load packages 
library(dplyr)
library(dorothea)
library(igraph)

# 1. Core transcription factors 
core_tfs <- c("Stat1", "Irf1", "Nfkb1", "Rela", "Stat3", "Cebpb", "Jun", "Fos")

# 2. Focus genes: secreted factors and checkpoint ligands
focus_genes <- c("Cd274", "Cxcl12", "Ccl2", "Ccl5", "Csf1", "Il6",
                 "Cd80", "Btn2a2", "Tnfsf9", "Vsir", "Icosl", "Fgl1")

# 3. Extract TF-target regulatory edges (
# Use Dorothea regulon with confidence levels A-C
edge_list <- data.frame()
for (tf in core_tfs) {
  targets <- dorothea_mm %>%
    filter(tf == !!tf, confidence %in% c("A", "B", "C")) %>%
    pull(target) %>% unique()

  # Keep only targets present in the focus gene set
  targets_focus <- intersect(targets, focus_genes)

  if (length(targets_focus) > 0) {
    edge_list <- rbind(edge_list,
                       data.frame(from = rep(tf, length(targets_focus)),
                                  to   = targets_focus))
  }
}

# 4. Build directed network 
g <- graph_from_data_frame(edge_list, directed = TRUE)

# Node color and size: TFs in blue, targets in pink
V(g)$color <- ifelse(V(g)$name %in% core_tfs, "#90bff9", "#ff8080")
V(g)$size  <- ifelse(V(g)$name %in% core_tfs, 25, 35)

# 5. Plot and save
png("Fig4A_TF_target_network.png", width = 8, height = 8, units = "in", res = 300)

set.seed(123)
layout_kk <- layout_with_kk(g)

plot(g,
     layout = layout_kk,
     main = "Core TF - target gene regulatory network",
     cex.main = 30,
     margin = c(0.2, 0.1, 0.5, 0.1),
     vertex.label.dist = 0,
     vertex.label.cex = 1.3,
     vertex.label.font = 2,
     vertex.label.color = "black",
     vertex.frame.color = "grey30",
     vertex.frame.width = 1.5,
     edge.color = "grey30",
     edge.width = 2.5,
     edge.arrow.size = 0.6,
     edge.arrow.width = 0.6)

dev.off()
