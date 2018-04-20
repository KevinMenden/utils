setwd("~/rimod/CAGE/")
library(ggplot2)
cage <- read.table("~/rimod/CAGE/analysis/CAGE_rLog_expression_values_080218.txt", sep="\t", header=T, row.names = 1)

groups <- as.character(sapply(colnames(cage), function(x){strsplit(x, split="_")[[1]][[6]]}))
groups[41] <- "C9orf72"

plotGene = function(x){
  gene <- as.numeric(cage[rownames(cage) == x,])
  df <- data.frame(exp = gene, groups = groups)
  p <- ggplot(df, aes(x = groups, y = exp, color = groups)) + geom_boxplot() + ggtitle(x)
  p
}
