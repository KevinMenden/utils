#######################################################
## Script to extract gene names of regulons
## and save them to separate files in a directory
########################################################

args <- commandArgs(trailingOnly = TRUE)
reg_file <- args[1]

regs <- readRDS(reg_file)

dir.create("go_enrichment")

tfs <- names(regs)
for (i in 1:length(tfs)) {
  tf <- tfs[i]
  genes <- as.character(unlist(regs[i]))
  write.table(genes, paste("go_enrichment/",tf,"_regulon_genes.txt", sep = ""), quote=F, row.names=F, col.names=F)
}
