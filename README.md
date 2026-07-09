# Cross-Species Analysis of Alzheimer's Disease

This repository contains the computational pipeline and gene modules generated for a systems-level, cross-species investigation of Alzheimer's Disease (AD). 
By integrating transcriptomic data from the mouse dataset and human post-mortem brain datasets, we identified conserved molecular networks driving neurodegeneration.

# Repository Structure

The repository is organized to ensure reproducibility of the computational analysis:

# Scripts
- `mus_musculus.R`: R script for the analysis of the mouse AD model.
- `homo_sapiens.R`: R script for the human MTG transcriptomic analysis.
- `Cross-Species_Analysis.R`: Implementation of the ortholog mapping, hypergeometric overlap testing, and consensus signature extraction.

# Data & Modules
- `modules_homo_sapiens/`: Human WGCNA module gene lists (Ensembl IDs).
- `modules_mouse_original/`: Native mouse WGCNA module gene lists.
- `modules_mouse_orthologs_IDs/`: Mapped ortholog modules in mouse Ensembl IDs.
- `metadata(1).txt`: Sample experimental design and phenotypic traits.
# Network & Annotation Files
- `string_interactions.tsv`: Protein-protein interaction data derived from the STRING database.
- `protein_functional_annotations.tsv`: Functional enrichment annotations.
- `protein_node_degrees.tsv`: Network centrality metrics.
- `node_table.csv`: Processed node attributes for visualization.

# Usage Instructions
1. Clone this repository to your local machine.
2. Ensure you have the required R packages (WGCNA, edgeR, gprofiler2, org.Mm.eg.db, pheatmap, flashClust, igraph, ggraph, tidygraph, reshape2).
3. Update the working directory paths in the provided R scripts to match your local environment, as indicated by the `[USER ACTION REQUIRED]` comments in the code.
4. Run the scripts in the following order:
    - `mus_musculus.R` / `homo_sapiens.R`
    - `Cross-Species_Analysis.R`
