# Script to collect motif files given as input the names of expressed TFs

## Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3) {
  args <- c("--help")
}
## Help section
if("-h" %in% args ||"--help" %in% args) {
  cat("
      Necessary arguments:
      <motif directory> directory containin all motifs from JASPAR CORE
      <expr_tfs> table containing expressed TFs
      <output_dir> director for ouptut
      ")
  stop()
}

motif_dir <- args[1]
expr_tfs <- args[2]
output_dir <- args[3]


# For testing
#motif_dir = "~/resources/TF_interactions/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar/"
#expr_tfs <- "~/rimod/CAGE/temporal/expressed_TFs_temporal_500rs.txt"
#output_dir <- "~/rimod/CAGE/temporal/test/"

# Get TFs
tfs <- read.table(expr_tfs, sep="\t")
tfs <- as.character(tfs$V1)

# Get list of JASPAR motif files
motif_files <- list.files(motif_dir, full.names = T)

# Copy motifs of expressed TFs to output directory
for (m in motif_files){
  mot = read.table(m, fill = T)
  gene = as.character(mot[1,2])
  if (gene %in% tfs){
    file.copy(m, output_dir)
  }
}