# Find list of expressed TFs
# Use count data from CAGE to filter TFs and TF-complexes


# Load libs
library(biomaRt)

# Argument parsing
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 4) {
  args <- c("--help")
}
## Help section
if("-h" %in% args ||"--help" %in% args) {
  cat("
      Please pass 4 arguments to the script in the following order:
      <cage_expression> cage expression count table
      <row_sum_cutoff> the row sum cutoff for the analysis (500 recommended)
      <motif_dir> directory containing the JASPAR motif files of available TFs
      <output_file> name of the output file")
  stop()
}




cage_expression <- args[1]
row_sum_cutoff <- as.numeric(args[2])
motif_dir <- args[3]
output <- args[4]
# 
# row_sum_cutoff = 500
# motif_dir = "~/resources/TF_interactions/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar/"
# cage_expression = "~/rimod/CAGE/cage_analysis/CAGE_annotated_counts_from_ct_2018-03-06_08.36.39/merged_annot_gene_counts_2018-03-06_08.36.39.txt"
# output = "expressed_TFs.txt"


# Load expression data 
cage <- read.table(cage_expression, sep="\t", header = T, row.names = 1)

# Filter for lowly expressed genes with row_sum_cutoff
row_sums <- apply(cage, 1, sum)
cage <- cage[row_sums > row_sum_cutoff,]

# Get gene symbols with biomaRt
print("Retrieving gene symbols from biomaRt")
genes <- as.character(sapply(rownames(cage), function(x){strsplit(x, split="[.]")[[1]][[1]]}))
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = genes, mart = ensembl)
expr_genes <- bm$hgnc_symbol

# Get name of all TFs (in JASPAR core)
print("Getting available TFs")
motif_files <- list.files(motif_dir, full.names = T)
tfs <- c()
# Copy motifs of expressed TFs to output directory
for (m in motif_files){
  mot = read.table(m, fill = T)
  gene = as.character(mot[1,2])
  tfs <- c(tfs, gene)
}

print("Calculating overlap")
# Dividie in TFs and TF-complexes
tf_comp = tfs[grepl("::", tfs)]
tfs_sgl = tfs[!grepl("::", tfs)]

# Function testing if tf is in a list
# Compare at uppercase level to avoid matching problems 
tf_in_list = function(tf, tf_list){
  if (toupper(tf) %in% toupper(tf_list)){
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

# Get TFs that are expressed in count table
# Compare at uppercase level to avoid matching problems
tfs_sgl_expr <- c()
for (tf in tfs_sgl) {
  if (tf_in_list(tf, expr_genes)){
    tfs_sgl_expr <- c(tfs_sgl_expr, tf)
  }
}

# Get all complexes of which all components are expressed in count table
tf_comp_expr = c()
for (tf in tf_comp){
  comps <- strsplit(tf, split="::")[[1]]
  complex_expressed = TRUE
  for (comp in comps){
    if (!tf_in_list(comp, expr_genes)){
      complex_expressed = FALSE
    }
  }
  if (complex_expressed){
    tf_comp_expr <- c(tf_comp_expr, tf)
  }
}

# Save result
all_expr_tfs <- c(tfs_sgl_expr, tf_comp_expr)
print(paste("Found ", as.character(length(all_expr_tfs)), " expressed TFs", sep=""))
write.table(all_expr_tfs, output, row.names = FALSE, quote = F)
