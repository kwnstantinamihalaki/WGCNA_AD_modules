# ROUTINE 1: CROSS-SPECIES GENE ORTHOLOG MAPPING (gprofiler2)
# SECTION 1: IMPORT HUMAN MODULE LISTS
# Scan the working directory for human WGCNA module gene files, read their 
# content into a structured list object, and clean file string names to retain 
# discrete module color identities.

# [NOTE] This step expects human module text files generated from previous steps 
# to be located inside the current working directory.
human_files = list.files(pattern = "^module_homosapiens_.*_genes.txt$")
human_files
human_modules = lapply(human_files, function(f){
  scan(f, what="character")
})

# Refine list names by stripping prefix and suffix text formats
names(human_modules) = gsub("module_homosapiens_|_genes.txt", "", human_files)

# SECTION 2: HOMOLOGY CONVERSION AND MOUSE ORTHOLOG EXPORT
# Utilize the 'gprofiler2' package to query human Ensembl symbols against the 
# mouse genome (mmusculus), isolate unique converted mouse ortholog names, 
# and export the mapped gene lists into species-specific text files.

library(gprofiler2)

# Create an empty list to store mouse orthologs.
mouse_modules = list()

# Execute loop structure to convert genes across species barriers
for (mod in names(human_modules)) {
  
  genes = human_modules[[mod]]  
  
  # Convert human genes to mouse ortholog symbols.
  orth = gorth(
    query = genes,
    source_organism = "hsapiens",
    target_organism = "mmusculus"
  )
  
  # Store and filter out unique mouse ortholog names for the active module.
  mouse_modules[[mod]] = unique(orth$ortholog_name)
}

# View converted mouse modules.
mouse_modules

# Save converted mouse modules to files.
for (mod in names(mouse_modules)) {
  writeLines(mouse_modules[[mod]],
             paste0("module_mouse_", mod, "_genes.txt"))
}


# ROUTINE 2: CROSS-SPECIES MODULE OVERLAP ANALYSIS (HYPERGEOMETRIC TEST)
# SECTION 1: DATA LOADING AND GENE ID CONVERSION
# Load mouse-mapped human orthologs and native mouse WGCNA modules. Convert 
# mouse gene Symbols to Ensembl IDs using 'org.Mm.eg.db' to standardize formats.

library(pheatmap)
library(reshape2)
library(org.Mm.eg.db)

# [USER ACTION REQUIRED] Update these folder paths to match the exact names/locations
# of your local directories containing the human-to-mouse ortholog files and 
# the native original mouse module text files.
human_ortho_path = "modules_mouse_orthologs"
mouse_study_path = "modules_mouse_ptyxiakis"

# Function to read all files in a directory into a list.
read_gene_lists = function(path) {
  files = list.files(path, pattern = "\\.txt$", full.names = TRUE)
  gene_lists = lapply(files, function(f) scan(f, what = character(), quiet = TRUE))
  names(gene_lists) = gsub(".txt", "", basename(files))
  return(gene_lists)
}

list_m_symbols = read_gene_lists(human_ortho_path) # human to mouse ortholog modules
list_n = read_gene_lists(mouse_study_path) # original mouse modules

# Convert gene Symbols to Ensembl IDs for consistency
list_m = lapply(list_m_symbols, function(syms) {
  if (length(syms) > 0) {
    id_map = mapIds(org.Mm.eg.db, 
                    keys = syms, 
                    column = "ENSEMBL", 
                    keytype = "SYMBOL", 
                    multiVals = "first")
    return(unique(as.character(id_map[!is.na(id_map)])))
  } else {
    return(character(0))
  }
})


# SECTION 2: STATISTICAL OVERLAP TESTING (HYPERGEOMETRIC ENRICHMENT)
# Establish a shared gene background universe and perform pairwise statistical 
# evaluation using the hypergeometric distribution (phyper) to assess overlap significance.

universe = unique(c(unlist(list_m), unlist(list_n)))
N = length(universe)

# Initialize the P-value Matrix
pval_matrix = matrix(NA, nrow = length(list_m), ncol = length(list_n))
rownames(pval_matrix) = names(list_m)
colnames(pval_matrix) = names(list_n)
pval_matrix_neg = pval_matrix

# Calculate dual-tail hypergeometric p-values across all module combinations
for (i in names(list_m)) {
  for (j in names(list_n)) {
    
    
    genes_i = list_m[[i]]
    genes_j = list_n[[j]]
    q = length(intersect(genes_i, genes_j))
    m = length(genes_i)
    n = N - m
    k = length(genes_j)
    pval_matrix_neg[i,j] = phyper(q, m, n, k, lower.tail = FALSE)
    pval_matrix[i, j] = phyper(q - 1, m, n, k, lower.tail = FALSE)
    if(pval_matrix[i,j]> pval_matrix_neg[i,j]){ pval_matrix[i,j] = -pval_matrix_neg[i,j]}
  }
}

# SECTION 3: SIGNIFICANCE VISUALIZATION (OVERLAP HEATMAP)
# Transform p-values using -log10 transformation, clean display labels, 
# and generate a clustered significance heatmap using pheatmap.

log_pval_matrix = -log10(abs(pval_matrix) + 1e-30)

# Clean module names
names(list_m) = gsub("module_mouse_|_genes", "", names(list_m))
names(list_n) = gsub("module_mouse_|_genes", "", names(list_n))
rownames(log_pval_matrix) = names(list_m)
colnames(log_pval_matrix) = names(list_n)

# Fix pheatmap error if values are identical.
log_pval_matrix = log_pval_matrix + runif(length(log_pval_matrix), 0, 1e-10)

# Export high-resolution Overlap Significance Heatmap
pdf("module_overlap_heatmap.pdf", width = 10, height = 8)
pheatmap(log_pval_matrix,
         main = "-log10(P-value) of Module Overlaps",
         color = colorRampPalette(c("white", "orange", "red"))(50),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         angle_col = 45,
         treeheight_row = 30,
         treeheight_col = 30)
dev.off() 


# ROUTINE 3: CONSENSUS SIGNATURE EXTRACTION AND ID FORMATTING
# SECTION 1: ENSEMBL MODULE REPOSITORY EXPORT
# Create a local output directory and write out the mapped human-to-mouse 
# module files containing standardized Ensembl Gene IDs.

dir.create("modules_mouse_orthologs_IDs", showWarnings = FALSE)
for (mod_name in names(list_m)) {
  file_name = paste0("modules_mouse_orthologs_IDs/", mod_name, "_IDs.txt")
  writeLines(list_m[[mod_name]], con = file_name)
}
cat("Storage pipeline complete. Exported to: modules_mouse_orthologs_IDs\n")

# SECTION 2: TRANS-SPECIES ROBUST SIGNATURE EXTRACTION
# Isolate intersection segments from targeted highly conserved pairs 
# (Yellow-Yellow, Pink-Magenta, Brown-Darkturquoise) to extract the unique 
# core conserved cross-species AD gene footprint.

# [NOTE] Ensure that the specific module target keys requested below match 
# the cleaned entry identities inside your active list objects

# Pair 1: Yellow (Mouse) - Yellow (Human Ortholog).
inter1 = intersect(list_m[["yellow"]], list_n[["module_yellow"]])

# Pair 2: Pink (Mouse) - Magenta (Human Ortholog).
inter2 = intersect(list_m[["magenta"]], list_n[["module_pink"]])

# Pair 3: Brown (Mouse) - Darkturquoise (Human Ortholog).
inter3 = intersect(list_m[["darkturquoise"]], list_n[["module_brown"]])

# Combine overlapping sets into a unified, non-redundant consensus pool.
combined_genes = unique(c(inter1, inter2, inter3))
length(combined_genes)

# Export the final conserved molecular signature text file.
writeLines(combined_genes, "combined_AD_ortholog_overlap.txt")                         

