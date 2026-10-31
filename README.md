# R code for RNA-seq analysis of G-CSF in lung cancer

## Study

Tumor-derived G-CSF promotes immune evasion and lung metastasis via PD-L1 
upregulation and immunosuppressive microenvironment remodeling.

## Data

RNA-seq data: GEO accession [GSE343634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE343634)

- Differential expression analysis was performed by BGI Genomics using DESeq2.
- Downstream analyses (volcano plot, GSEA, transcription factor activity, 
  heatmap, and regulatory network) were performed in R on the BGI-provided 
  differential expression table and expression matrix.
- Reads were aligned to the mouse reference genome **GRCm39** 
  (NCBI GCF_000001635.27) using HISAT2, and gene expression was quantified 
  as FPKM using Bowtie2.
- In the differential expression table, "N" and "S" denote the sh-NC and 
  sh-Csf3 groups, respectively. log2FC values are oriented so that positive 
  values indicate up-regulation in sh-Csf3 relative to sh-NC.

## Scripts

| File | Description | Figure |
|------|-------------|--------|
| `01_volcano_plot.R` | Differential expression volcano plot | Fig. 3A |
| `02_GSEA_Hallmark.R` | GSEA with Hallmark gene sets: bubble plot and enrichment plots | Fig. 3B–E |
| `03_GSEA_custom_checkpoint.R` | GSEA with a custom immune checkpoint ligand gene set | Fig. 3F |
| `04_heatmap.R` | Heatmap of immune checkpoint ligands and cytokines | Fig. 3G |
| `05_TF_activity_fgsea_dorothea.R` | Transcription factor activity inferred by fgsea on Dorothea regulons | Fig. 3J |
| `06_TF_target_network.R` | Transcriptional regulatory network of core TFs and target genes | Fig. 4A |

## Software and packages

- R (version 4.2)
- DESeq2 (v1.34.0)
- fgsea
- clusterProfiler (v4.2.0)
- dorothea
- igraph
- ggplot2, ggrepel, dplyr, tibble, pheatmap, RColorBrewer

## Note

Statistical analyses and graphs for the experimental data (Figs. 1C–O, 2, 
3H–I, 4B–J, 5, and 6) were performed using GraphPad Prism (v8.0). Flow 
cytometry data were analyzed with FlowJo (v10.8.1).

## Contact

Corresponding author: wangtong@sdsmu.edu.cn
