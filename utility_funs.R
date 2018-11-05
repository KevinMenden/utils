###############################################################
## Collection of utility funtions that might be useful       ##
## in many other scripts                                     ##
###############################################################
library(limma)
library(viridis)

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

utils.ewce.plot <- function(total_res,mtc_method="bonferroni"){
  if(!mtc_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")){
    stop("ERROR: Invalid mtc_method argument. Please see '?p.adjust' for valid methods.")
  }
  multiList = TRUE
  if(is.null(total_res$list)){multiList = FALSE}
  
  # Multiple testing correction across all rows
  total_res$q = p.adjust(total_res$p,method=mtc_method)
  
  # Mark significant rows with asterixes
  ast_q = rep("",dim(total_res)[1])
  ast_q[total_res$q<0.05] = "*"
  total_res$ast_q = ast_q
  
  # GENERATE THE PLOT
  total_res$sd_from_mean[total_res$sd_from_mean<0]=0
  graph_theme = theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(panel.grid.major = element_line(size = .5, color = "grey"),
          axis.line = element_line(size=.7, color = "black"), text = element_text(size=14),
          axis.title.y = element_text(vjust = 0.6))# + theme(legend.position="none")
  
  #total_res$
  upperLim = max(abs(total_res$sd_from_mean))
  
  total_res$y_ast = total_res$sd_from_mean*1.05
  
  total_res$abs_sd = abs(total_res$sd_from_mean)
  
  #print(upperLim)
  if("Direction" %in% colnames(total_res)){
    #the_plot = ggplot(total_res) + geom_bar(aes(x=CellType,y=abs(sd_from_mean),fill=Direction),position="dodge",stat="identity") + graph_theme
    the_plot = ggplot(total_res) + geom_bar(aes_string(x='CellType',y='abs_sd',fill='Direction'),position="dodge",stat="identity") + graph_theme
  }else{
    #the_plot = ggplot(total_res) + geom_bar(aes(x=CellType,y=abs(sd_from_mean),fill="red"),stat="identity") + graph_theme +theme(legend.position="none")
    the_plot = ggplot(total_res) + geom_bar(aes_string(x='CellType',y='abs_sd'),fill=viridis(3)[2],stat="identity") + graph_theme +theme(legend.position="none")
  }
  
  the_plot = the_plot  +
    theme(plot.margin=unit(c(1,0,0,0),"mm"),axis.text.x = element_text(angle = 55, hjust = 1, size = 15))+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
    xlab("") +
    theme(strip.text.y = element_text(angle = 0)) +
    coord_cartesian(ylim = c(0,1.1*upperLim))+
    ylab("Std.Devs. from the mean") + theme(plot.margin = unit(c(0,0,0,1.5), "cm"))
  
  the_plot = the_plot + scale_y_continuous(breaks=c(0,ceiling(upperLim*0.66))) + geom_text(aes_string(label="ast_q",x="CellType",y="y_ast"),size=10)
  
  if(multiList){
    the_plot = the_plot + facet_grid("list ~ .",scales="free", space = "free_x")
    #the_plot = the_plot + facet_grid(facets="list ~ .",scale="free", space = "free_x")
  }
  
  return(the_plot)
}

