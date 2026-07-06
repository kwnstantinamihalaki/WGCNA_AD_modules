# ROUTINE 1: DIFFERENTIAL EXPRESSION ANALYSIS (edgeR) - Mouse Hippocampus
# SECTION 1: DATA LOADING, METADATA AND EXPERIMENTAL DESIGN
# Import mouse raw RNA-Seq counts and establish factors for genotype (wt/ko),
# disease status (h/ad), and sex (f/m). Construct a detailed sample metadata
# matrix and set up the experimental design groups.

library(edgeR)
library(pheatmap)
rm(list=ls())

# [USER ACTION REQUIRED] ENVIRONMENT SETUP
# Before running the script, please update the path below to match your local
# working directory where the input data files are stored.
setwd("C:/ΒΙΟΛΟΓΙΚΟ/PTYXIAKI/analysis/")

# Import the raw expression counts (Ensure the file is in your working directory)
a = read.table("alzheimer_newdata.tabular", h=TRUE, sep="\t", row.names=1)
genotype = factor(c(rep("wt", 6), rep("ko", 6)), levels=c('wt', 'ko'))
disease = factor(c(rep(rep(c("h", "ad"), each=3), 2)), levels=c('h', 'ad'))
gender = factor(c('f', 'f', 'f', 'm', 'm', 'm', 'f', 'm', 'm', 'f', 'f', 'f'), levels = c('f', 'm'))
metadata = data.frame(sample=1:length(genotype) , batch=gender, factor1=genotype, factor2=disease)
metadata
group = factor(c(paste0(genotype, disease, gender)), levels=c("wthf", "wtadf", "kohf", "koadf", "wthm", "wtadm", "kohm", "koadm"))
group
index=1:ncol(a)
dim(a)
oldnames = names(a)
newnames = c(paste(genotype, disease, index, sep=""))
colnames(a) = newnames
genotype = factor(genotype, levels=c("wt", "ko"))
disease = factor(disease, levels = c("h", "ad"))
group
design = model.matrix(~0 + group)

# SECTION 2: DATA FILTERING AND TMM NORMALIZATION
# Create a DGEList object, apply TMM normalization, filter out low-expressed
# genes across the design matrix using 'filterByExpr', and calculate empirical
# Bayes negative binomial dispersions.

y = DGEList(counts = a, group=group)
y = normLibSizes(y)
keep = filterByExpr(y, design)## filter only based on genotype not disease
y <- y[keep, , keep.lib.sizes=FALSE]
y = estimateDisp(y)
logcpm = cpm(y, log=TRUE)

# SECTION 3: EXPLORATORY DATA ANALYSIS (PCA & GLOBAL HEATMAP)
# Vectorize metadata group colors to evaluate sample segregation. Execute PCA
# and plot sample profiles alongside a global row-scaled expression heatmap.

# Define display colors based on disease phenotype.
col = rep('blue', length(gender))
col[gender=='f'] = 'red'
col = rep('blue', length(genotype))
col[genotype=='ko'] = 'red'
col = rep('blue', length(disease))
col[disease=='ad'] = 'red'
pdf("exploratory.pdf")

# Principal Component Analysis Plot.
pca = prcomp(t(logcpm), scale.=TRUE, center=TRUE)
PC1 <- pca$x[,1]
PC2 <- pca$x[,2]
plot(PC1, PC2, 
     xlab="PC1", ylab="PC2", 
     pch=19, cex=1.5, col='white')
text(PC1, PC2, labels=newnames, col=col, cex=1.2)

# Global Exploratory Heatmap.
pheatmap(logcpm, color=colorRampPalette(c('blue', 'white', 'red'))(50), scale = 'row', show_rownames = FALSE)
dev.off()

# SECTION 4: ADDITIVE MODEL ANALYSIS (GENDER, GENOTYPE, AND DISEASE EFFECTS)
# Re-construct the design matrices using additive formulations to identify 
# independent DEGs for demographic and clinical factors via Quasi-Likelihood F-tests.

# Evaluation of Gender-Specific DEGs
design = model.matrix(~0+genotype+disease+gender)

colnames(design) = gsub(pattern="genotype", "", x=colnames(design))
colnames(design) = gsub(pattern="disease", "", x=colnames(design))
colnames(design) = gsub(pattern=":", "", x=colnames(design))
colnames(design) = gsub(pattern="gender", "", x=colnames(design))

# Re-initialize DGEList for additive model filtering and dispersion mapping.
rm(y)
y = DGEList(counts = a)
y = normLibSizes(y)
keep = filterByExpr(y, design) 
y = y[keep, , keep.lib.sizes=FALSE]
y = estimateDisp(y)
fit = glmQLFit(y, design)
mycontrasts1 = makeContrasts(kowt = ko-wt, ad = ad, gender=m, levels=design)
mycontrasts1
qlf = glmQLFTest(fit, contrast=mycontrasts1[,"gender"])
res = topTags(qlf)
res

genesGender = rownames(res$table[res$table['FDR']<0.05,])
genesGender

# Generate Gender Expression Heatmap.
pheatmap(logcpm[genesGender,])
dev.off()


# Evaluation of Genotype-Specific DEGs.

design = model.matrix(~0+disease+gender + genotype)
design
colnames(design) = gsub(pattern="genotype", "", x=colnames(design))
colnames(design) = gsub(pattern="disease", "", x=colnames(design))
colnames(design) = gsub(pattern=":", "", x=colnames(design))
colnames(design) = gsub(pattern="gender", "", x=colnames(design))

rm(y)
y = DGEList(counts = a)
y = normLibSizes(y)
keep = filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes=FALSE]
y = estimateDisp(y)

fit = glmQLFit(y, design)
mycontrasts2 = makeContrasts(had = ad-h, gender=m, genotype=ko, levels=design)

# Test and extract Genotype (KO vs WT) DEGs
qlf = glmQLFTest(fit, contrast=mycontrasts2[,"genotype"])
topTags(qlf)
res = topTags(qlf)
genesGenotype = rownames(res$table[res$table['FDR']<0.05,])

# Generate Genotype Expression Heatmap.
pheatmap(logcpm[genesGenotype,])
dev.off()

# Evaluation of Disease-Specific DEGs.

# Test and extract Disease (AD vs Healthy) DEGs
qlf = glmQLFTest(fit, contrast=mycontrasts2[,"had"])
topTags(qlf)
res = topTags(qlf)
genesDisease = rownames(res$table[res$table['FDR']<0.05,])

# Generate Disease Expression Heatmap.
pheatmap(logcpm[genesDisease,])
dev.off()

# SECTION 5: INTERACTION MODEL ANALYSIS (DISEASE × GENDER EFFECTS)
# Implement an interaction design framework to intercept and extract genes 
# whose transcriptomic response to AD is conditional upon biological sex.

design = model.matrix(~0+disease*gender + genotype)
design
colnames(design) = gsub(pattern="genotype", "", x=colnames(design))
colnames(design) = gsub(pattern="disease", "", x=colnames(design))
colnames(design) = gsub(pattern=":", "", x=colnames(design))
colnames(design) = gsub(pattern="gender", "", x=colnames(design))

fit = glmQLFit(y, design)
mycontrasts3 = makeContrasts(had = ad-h, gender=m, genotype=ko, adm=adm, levels=design)

qlf_int = glmQLFTest(fit, contrast=mycontrasts3[,"adm"])
topTags(qlf_int)
res_int = topTags(qlf_int)
genesInteraction = rownames(res_int$table[res_int$table['FDR']<0.05,])

# Conditional Heatmap Output based on statistical significance thresholds
pdf("degs_interaction_heatmap.pdf")
if (length(genesInteraction) > 0) {
  pheatmap(logcpm[genesInteraction,])
} else {
  message("No significant genes detected for the interaction effect (disease × gender).")
}
dev.off()


# ROUTINE 2: WEIGHTED GENE CO-EXPRESSION NETWORK ANALYSIS (WGCNA) - Mouse
# SECTION 1: ENVIRONMENT SETUP AND DATA QUALITY CONTROL
# Prepare the expression matrix for WGCNA, enable parallel processing threads, 
# and filter out low-variance or missing gene/sample entries.

library(WGCNA)
library(flashClust)
allowWGCNAThreads()
datExpr = as.data.frame(t(logcpm))

# Check data profiles for outlying genes or samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr)

# SECTION 2: SOFT-THRESHOLDING POWER EVALUATION
# Screen a series of power values (1-15) to evaluate the network scale-free 
# topology fit index and mean connectivity patterns.
powers <- c(1:15)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, RsquaredCut=0.8)

# Plot Scale-Free Topology Fit and Mean Connectivity
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (Power)", ylab = "Scale-Free Topology Model Fit", type = "n", main = "Scale-Free Topology")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab = "Soft Threshold (Power)", ylab = "Mean Connectivity", type = "n", main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "blue")

# SECTION 3: TOPOLOGICAL OVERLAP MATRIX (TOM) AND CLUSTERING
# Construct adjacency network and convert to Topological Overlap Matrix (TOM).
# Apply average linkage hierarchical clustering and dynamic tree cutting.
softPower = 5
adjacency = adjacency(datExpr, power = softPower)

# [NOTE] Load or create local structural TOM cache object
load(file='tom.RData')

TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

# [USER ACTION REQUIRED] Save computed TOM to local working directory
save(TOM, file="tom.RData")

# Hierarchical Clustering and Module Partitioning
geneTree = flashClust(as.dist(dissTOM), method = "average")
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 100)
moduleColors = labels2colors(dynamicMods)

# Plot Gene Dendrogram and Module Assignment Colors
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03)

# SECTION 4: EXPORTING MODULE GENE LISTS
# Extract gene identities belonging to each co-expression module, remove 
# Ensembl dot notation suffixes, and save lists to flat text files.
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
  writeLines(genes_clean, paste0("module_", mod, "_genes.txt"))
}

# SECTION 5: TRAIT CORRELATION ANALYSIS AND HEATMAP EXPORT
# Calculate module eigengenes (MEs) and compute Pearson correlation coefficients 
# against discrete clinical phenotypic traits loaded from the sample metadata.
moduleColors = moduleColors[colnames(datExpr)]
stopifnot(ncol(datExpr) == length(moduleColors))
MEs = moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEX = moduleEigengenes(datExpr, colors = moduleColors, nPC=10)

# [USER ACTION REQUIRED] Ensure path targets your correct local metadata file.
metadata = read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1)
metadata$Genotype = as.numeric(as.factor(metadata$Genotype))
metadata$Gender = as.numeric(as.factor(metadata$Gender))

# Compute Pearson correlation matrix
moduleTraitCor = cor(MEs, metadata, use = "p")

# Generate and save Module-Trait Correlation Heatmap
pdf("module_trait_correlations.pdf", width = 12, height = 10)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metadata),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               main = "Module-Trait Correlations")
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
  ggtitle("WGCNA Gene Co-expression Network") +
  theme(legend.position = "none")

# Export graph topology visualization
ggsave("WGCNA_network.pdf", plot = network_plot, width = 12, height = 12)