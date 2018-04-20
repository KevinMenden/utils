###############################################################
## Collection of utility funtions that might be useful       ##
## in many other scripts                                     ##
###############################################################
library(limma)

# Filter features (e.g. genes) based on value (e.g. CPM)
# cutoff is the cutoff value below which features will be discarded
# numSampls is the number of samples that must be above the cutoff in 
# order to keep the feature
filterFeats <- function(mat, cutoff = 1, numSamps = 1){
  f = function(x){
    y <- x[x >= cutoff]
    return(length(y))
  }
  keep <- apply(mat, 1, f)
  keep <- keep >= numSamps
  mat.flt <- mat[keep,]
  return(mat.flt)
}

# Plot multiple MDS plots
plotMultipleMDS = function(x, labels, cols, gene.selection = "common", top = 1000){
  library(limma)
  library(RColorBrewer)
  pal <- brewer.pal(length(levels(factor(cols))), "Dark2")
  par(mfrow=c(2,2))
  plotMDS(x, top=top, gene.selection=gene.selection, labels=labels, col = pal[factor(cols)])
  plotMDS(x, top=top, gene.selection=gene.selection, labels=labels, col = pal[factor(cols)], dim=c(2,3))
  plotMDS(x, top=top, gene.selection=gene.selection, labels=labels, col = pal[factor(cols)], dim=c(1,3))
  plotMDS(x, top=top, gene.selection=gene.selection, labels=labels, col = pal[factor(cols)], dim=c(3,4))
  par(mfrow=c(1,1))
}


## Calculate the fold change
# Provide proper contrast names!
calcFoldChange <- function(x, groups, conts){
  groups <- make.names(groups)
  res.df <- data.frame(dummy = rep(1, nrow(x)))
  groupvec <- c(1:ncol(x))
  for (j in 1:length(conts)){
    cont <- conts[j]
    cont <- strsplit(cont, split="-")[[1]]
    g1 <- cont[1]
    g2 <- cont[2]
    g1.vec <- groupvec[groups == g1]
    g2.vec <- groupvec[groups == g2]
    g1.sub <- x[,g1.vec]
    g2.sub <- x[,g2.vec]
    g1.mean <- apply(g1.sub, 1, mean)
    g2.mean <- apply(g2.sub, 1, mean)
    mean.df <- data.frame(g1 = g1.mean, g2 = g2.mean)
    rownames(mean.df) <- names(g1.mean)
    cont.fcs <- c()
    f = function(r){log2(r[1]/r[2])}
    cont.fcs <- apply(mean.df, 1, f)
    res.df <- cbind(res.df, cont.fcs)
  }
  res.df <- res.df[,-1]
  if (length(conts) == 1){
    res.df <- data.frame(cont = res.df)
  }
  colnames(res.df) <- make.names(conts)
  rownames(res.df) <- rownames(x)
  return(res.df)
}

