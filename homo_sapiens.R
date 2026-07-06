# ROUTINE 1: DIFFERENTIAL EXPRESSION ANALYSIS (edgeR) - Homo sapiens MTG
#SECTION 1: DATA LOADING AND DESIGN MATRIX
# Load human MTG raw RNA-Seq counts, define experimental groups 
# (young_control, old_control, AD), and set up the model design matrix.

library(edgeR)
rm(list=ls())

# [USER ACTION REQUIRED] ENVIRONMENT SETUP
# Before running the script, please update the path below to match your local
# working directory where the input data files are stored.
setwd("C:/ΒΙΟΛΟΓΙΚΟ/PTYXIAKI/analysis/")

# Import the raw expression counts (Ensure the file is in your working directory)
a = read.table("alzheimer_homo_sapiens.tabular", h=TRUE, sep="\t", row.names=1)
head(a)
group = factor(c(rep("old_control",1), rep("AD", 1), rep("old_control", 2), rep("young_control",1), rep("AD",5), rep("young_control",1), rep("old_control",5), rep("AD",3), rep("young_control",2), rep("old_control",2), rep("AD",3), rep("young_control",4)), levels= c("young_control", "old_control", "AD"))
sample = colnames(a)
metadata = data.frame(
  sample = sample,
  group = group
)
rownames(metadata) = sample
metadata
dim(a)
design = model.matrix(~0 + group)
colnames(design) = levels(group)
design
head(a)

# SECTION 2: DATA FILTERING AND NORMALIZATION
# Apply TMM normalization to scale library sizes, filter out low-expressed 
# genes using 'filterByExpr', and estimate negative binomial dispersions.

y = DGEList(counts = a, group=group)
tail(sort(a[,1]))
y = normLibSizes(y)
keep = filterByExpr(y, design)
y = y[keep, , keep.lib.sizes=FALSE]
y = estimateDisp(y)
logcpm = cpm(y, log=TRUE)

# SECTION 3: EXPLORATORY DATA ANALYSIS (PCA & HEATMAP)
# Perform Principal Component Analysis (PCA) and generate a global heatmap 
# on logCPM values to check sample clustering patterns and profiles.

col = rep("white", length(group))
col[group == "young_control"] <- "blue"
col[group == "old_control"]   <- "black"
col[group == "AD"]            <- "red"

# PCA plot
pdf("PCA_exploratory_homo_sapiens.pdf")
pca = prcomp(t(logcpm), scale. = TRUE, center = TRUE)
PC1 = pca$x[,1]
PC2 = pca$x[,2]
plot(
  PC1, PC2,
  xlab = "PC1",
  ylab = "PC2",
  pch = 19,
  cex = 1.5,
  col = col
)
dev.off()

# Heatmap
library(pheatmap)
pdf("heatmap_exploratory_homo_sapiens.pdf")
pheatmap(
  logcpm,
  color = colorRampPalette(c("blue", "black", "red"))(50),
  scale = "row",
  show_rownames = FALSE,
  show_colnames = FALSE
)
dev.off()

# SECTION 4: GLM TESTING AND CONTRASTS
# Fit a Quasi-Likelihood GLM and execute pairwise tests for two contrasts: 
# AD vs Old Control (pathology) and Old vs Young Control (healthy aging).

fit = glmQLFit(y, design)

# Contrast 1: AD vs Old Control
contrast_AD_vs_Old = makeContrasts(AD_vs_Old = AD - old_control,levels = design)
contrast_AD_vs_Old
qlf_AD = glmQLFTest(fit, contrast = contrast_AD_vs_Old)
res = topTags(qlf_AD)
res

# Extract and plot significant DEGs (FDR < 0.05)
genes = rownames(res$table[res$table['FDR']<0.05,])
genes
pheatmap(logcpm[genes,])
dev.off()

# Contrast 2: Old vs Young Control
contrast_aging = makeContrasts(Old_vs_Young = old_control - young_control,levels = design)
contrast_aging
qlf_age = glmQLFTest(fit, contrast = contrast_aging)
res = topTags(qlf_age)
res


# ROUTINE 2: WEIGHTED GENE CO-EXPRESSION NETWORK ANALYSIS (WGCNA) - Homo sapiensWEIGHTED GENE CO-EXPRESSION NETWORK ANALYSIS (WGCNA) - Homo sapiensWEIGHTED GENE CO-EXPRESSION NETWORK ANALYSIS (WGCNA) - Homo sapiens
# SECTION 1: ENVIRONMENT SETUP AND DATA QUALITY CONTROL
# Prepare expression matrix for WGCNA, enable multi-threading, and filter out
# low-quality samples or genes with zero variance using 'goodSamplesGenes'.

library(WGCNA)
library(flashClust)
allowWGCNAThreads()
datExpr = as.data.frame(t(logcpm))
gsg = goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr)

# SECTION 2: SOFT-THRESHOLDING POWER SELECTION
# Evaluate scale-free topology fit indices and mean connectivity across powers (1-15)
# to select the optimal soft-thresholding power for network construction.

powers = c(1:15)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, RsquaredCut=0.8)

# Plot Scale-Free Topology Fit
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (Power)", ylab = "Scale-Free Topology Model Fit", type = "n", main = "Scale-Free Topology")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab = "Soft Threshold (Power)", ylab = "Mean Connectivity", type = "n", main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "blue")

# SECTION 3: NETWORK CONSTRUCTION AND MODULE DETECTION
# Compute adjacency and Topographic Overlap Matrix (TOM), perform hierarchical
# clustering, and apply dynamic tree cut to define co-expression modules.

softPower = 9
adjacency = adjacency(datExpr, power = softPower)

# Load or compute TOM similarity matrix
load(file='tom.RData')

TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

# [USER ACTION REQUIRED] Save computed TOM to local working directory
save(TOM, file="tom.RData")

# Hierarchical Clustering and Dynamic Tree Cutting
geneTree = flashClust(as.dist(dissTOM), method = "average")
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 100)
moduleColors = labels2colors(dynamicMods)

# Plot Dendrogram with Module Colors
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03)

# SECTION 4: MODULE GENE EXPORT
# Map gene names to their designated co-expression modules, clean Ensembl IDs
# by removing version suffixes, and export gene lists to text files.

geneNames = colnames(datExpr)
length(geneNames)
length(dynamicMods)
names(moduleColors) = geneNames
moduleGeneList = split(geneNames, moduleColors)
moduleGeneList

# Save text modules to local directory.
for (mod in names(moduleGeneList)) {
  genes = moduleGeneList[[mod]]
  genes_clean = gsub("\\..*$", "", genes)  
  writeLines(genes_clean, paste0("module_homosapiens_", mod, "_genes.txt"))
}

# SECTION 5: TRAIT CORRELATION AND MODULE-TRAIT HEATMAP
# Calculate module eigengenes (MEs), compute Pearson correlations with clinical 
# phenotypes (Young, Old, AD), and visualize results via labeled heatmaps.

MEs = moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEX = moduleEigengenes(datExpr, colors = moduleColors, nPC=10)

# [USER ACTION REQUIRED] Ensure path targets your correct local metadata file.
metadata = read.table("metadata1.txt", header = TRUE, sep = "\t", row.names = 1)
str(metadata)
metadata$Young = as.numeric(metadata$Young)
metadata$Old = as.numeric(metadata$Old)
metadata$AD = as.numeric(metadata$AD)

# Compute Pearson correlation matrix
moduleTraitCor = cor(MEs, metadata, use = "p")

# Generate and save Module-Trait Correlation Heatmap
pdf("module_trait_correlations_Homo_sapiens.pdf", width = 12, height = 10)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metadata),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               main = "Module-Trait Correlations_Homo_sapiens")
dev.off()


# ROUTINE 3: CO-EXPRESSION NETWORK CONSTRUCTION AND VISUALIZATION (igraph/ggraph)
# SECTION 1: NETWORK DIMENSIONALITY REDUCTION AND THRESHOLDING
# Map gene names to the Topographic Overlap Matrix (TOM), apply an adjacency 
# threshold (TOM > 0.1) to retain only strong connections, and optionally subset 
# the top 1000 most variable genes to optimize computational efficiency.

library(igraph)
library(ggraph)
library(tidygraph)
rownames(TOM) = colnames(datExpr)
colnames(TOM) = colnames(datExpr)

# Apply topological overlap similarity threshold
tom_threshold = 0.1
TOM_sparse = TOM
TOM_sparse[TOM_sparse < tom_threshold] = 0
cat("Remaining edges after thresholding:", sum(TOM_sparse > 0) / 2, "\n")

# Subset top variable genes if network scale exceeds threshold limits
if (ncol(datExpr) > 2000) {
  gene_var = apply(datExpr, 2, var)
  top_genes = names(sort(gene_var, decreasing = TRUE))[1:1000]
  TOM_sparse = TOM_sparse[top_genes, top_genes]
  moduleColors = moduleColors[top_genes]
}

# SECTION 2: GRAPH INTERACTION GENERATION AND NODE ATTRIBUTE MAPPING
# Construct an undirected, weighted igraph network object from the sparse TOM 
# matrix, remove redundant loops or multiple edges, and assign WGCNA module 
# membership colors as discrete node attributes.

graph = graph_from_adjacency_matrix(
  TOM_sparse,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

# Map module colors to graph vertices
V(graph)$module = moduleColors[match(V(graph)$name, names(moduleColors))]

# Simplify network architecture
graph = simplify(graph, remove.loops = TRUE, remove.multiple = TRUE)

# SECTION 3: TOPOLOGY LAYOUT GRAPH PLOTTING AND EXPORT
# Generate a network plot utilizing the Fruchterman-Reingold (FR) force-directed 
# layout, scale edge widths by statistical weight, color code nodes by module 
# identity, and export the high-resolution visualization to a PDF.

set.seed(123)
network_plot = ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.4, color = "gray") +
  geom_node_point(aes(color = module), size = 3) +
  theme_void() +
  ggtitle("WGCNA Gene Co-expression Network-Homo sapiens") +
  theme(legend.position = "none")

# Export graph topology visualization
ggsave("WGCNA_network_Homo_sapiens.pdf", plot = network_plot, width = 12, height = 12)