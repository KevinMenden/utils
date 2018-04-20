
args <- commandArgs(trailingOnly = TRUE)

deg_file <- as.character(args[1])
pval <- args[2]
lfc <- args[3]
print(deg_file)

deg <- read.table(deg_file, sep="\t", row.names = 1, header = T)

deg.sig <- deg[deg$padj <= pval,]
if (lfc > 0){
  deg.sig <- deg.sig[deg.sig$log2FoldChange >= lfc,]
} else{
  deg.sig <- deg.sig[deg.sig$log2FoldChange <= lfc,]
}
dim(deg.sig)

deg.bg <- deg[!rownames(deg) %in% rownames(deg.sig),]

createBedFile = function(df){
  genes <- rownames(df)
  chr <- as.character(sapply(genes, function(x){strsplit(x, split="_")[[1]][[1]]}))
  start <- as.character(sapply(genes, function(x){strsplit(x, split="_")[[1]][[2]]}))
  end <- as.character(sapply(genes, function(x){strsplit(x, split="_")[[1]][[3]]}))
  strand <- as.character(sapply(genes, function(x){strsplit(x, split="_")[[1]][[4]]}))
  new_df = data.frame(chr, start, end, strand)
  return(new_df)
}

sig.bed <- createBedFile(deg.sig)
bg.bed <- createBedFile(deg.bg)

# Make background twice as big as sig genes
cutoff = 2 * nrow(sig.bed)
bg.bed <- bg.bed[1:cutoff,]

write.table(sig.bed, "sig_genes.bed", sep="\t", quote=F, col.names = F, row.names = F)
write.table(bg.bed, "background.bed", sep="\t", quote=F, col.names = F, row.names = F)

