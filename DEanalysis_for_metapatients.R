###################################################################
# DIFFERENTIAL EXPRESSION (DE) ANALYSIS FOR METAPATIENTS  
#
#
# Generates and saves the table of DEGs for each metapatient defined 
# in a cluster .rds object for a given disease.
#
# Script parameters (adjust them for each run)
#
# dis               - Disease of the metapatients
# clusters_filename - Filename of the clusters .rds object
# results_dir       - Output directory
#
# Beatriz Urda García 2020
######################################################################

library(knitr)
library(RColorBrewer)
library(SummarizedExperiment)
library(geneplotter)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(edgeR)
library(corpcor) 
library(M3C)
library(calibrate)
library(EnhancedVolcano) 
library(sva) 
library(pvclust) # Hierarchical clustering with p-values
library(dendextend) # fancy plotting dendograms

#### RUN DEAnalysis for a set of diseases


# < Script parameters for the meta-patients obtained with all genes> (adjust them for each run) 
setwd("~/Desktop/ANALYSIS")
clusters_dir <- 'Metapatients/with_entire_count_matrix/final_clusters/'
clusters_list <- list.files(clusters_dir)
dis_list <- gsub("_normalized_counts.rds","",clusters_list)
dis_list <- gsub("_combat_counts.rds.*","",dis_list)
metap_method = 'allgenes'
outputdir <- "Metapatients/with_entire_count_matrix/DEAnalysis/"
output_rnk <- 'Metapatients/Ranked_metapatients/Ranked_metapatients_allgenes/'
# </>

# < Script parameters for sDVGs and increased variability> (adjust them for each run)
setwd("~/Desktop/ANALYSIS")
clusters_dir <- 'Metapatients/using_sDVGs/final_clusters_increased/'
clusters_list <- list.files(clusters_dir)
dis_list <- gsub("_clusters_.*","",dis_list)
metap_method = 'sDVGs_increased'
outputdir <- "Metapatients/using_sDVGs/DEAnalysis/DEA_metapatients_increased/"
output_rnk <- 'Metapatients/Ranked_metapatients/Ranked_metapatients_increased/'
# </>

run_DEAnalysis_for_dis_metapatients <- function(dis,clusters_filename,results_dir){
  # dev.off()
  # pdf(paste(results_dir,"/",dis,".pdf",sep=""))
  # Importing metadata, se and dge objects
  meta <- readRDS("new_meta_jon_grein_good_quality.rds")
  se <- readRDS(paste('DL_RESULTS_GREIN/normalized_counts/',dis,'_se_normalized.rds',sep=""))
  dge <- readRDS(paste('DL_RESULTS_GREIN/normalized_counts/',dis,'_dge_normalized.rds',sep=""))
  
  #### DEFINING METAPATIENTS SE & DGE
  # Selecting the controls of the disease
  controls <- meta[meta$disease == dis,]
  controls <- as.character(controls$sample) # 61 samples
  
  # Importing the clusters
  clus <- readRDS(clusters_filename) 
  
  # Number of clusters
  n_clusters = dim(as.matrix(table(clus)))[1]
  
  # Creating a list with n_clusters elements. Each element contains a list with its samples -always patients
  clus_samples <- c()
  se_list <- c()
  dge_list <- c()
  for(k in 1:n_clusters){
    cpatients <- names(clus[clus == k])
    clus_samples <- append(clus_samples, list(cpatients))
    se_list <- append(se_list,list(se[,append(cpatients,controls)]))
    dge_list <- append(dge_list,list(dge[,append(cpatients,controls)]))
  }
  
  colnames(se_list[[1]])
  
  for(k in 1:n_clusters){
    se.filt <- se_list[[k]]
    dge.filt.norm <- dge_list[[k]]
    
    # Without correcting
    se.filt$type <- relevel(se.filt$type, ref="control")
    mod <- model.matrix(~ se.filt$type, colData(se.filt))
    mod0 <- model.matrix(~ 1, colData(se.filt))
    pv <- f.pvalue(assays(se.filt)$logCPM.norm, mod, mod0)
    sum(p.adjust(pv, method="fdr") < 0.5)
    
    # Distribution of expected p-values
    par(mfrow=c(1, 1))
    try(hist(pv, main="Distribution of expected p-values", las=1, ylim = c(0, 8000)))
    
    # Correcting for the series
    if (length(unique(se$series)) > 1){
      mod <- model.matrix(~type + series, colData(se.filt))
      mod0 <- model.matrix(~series, colData(se.filt))
      colnames(mod)
      pValues <- f.pvalue(assays(se.filt)$logCPM.norm,mod,mod0)
      sum(p.adjust(pValues, method="fdr") < 0.05)
    } else {
      print("There is only one series for this disease")
    }
    
    # Distribution of expected p-values after correcting for the series
    if (length(unique(se$serie)) > 1){
      try(hist(pValues, main="Distribution of expected p-values after correcting for series", las=1))
    }
    
    # Mean-variance relationship
    FDRcutoff <- 0.05
    par(mar=c(1,1,1,1))
    v <- voom(dge.filt.norm, mod, plot = TRUE) # Voom is applied to the normalized and filtered DGEList object
    fit <- lmFit(v, mod)
    fit<- eBayes(fit)
    res <- decideTests(fit, p.value = FDRcutoff)
    summary(res)
    
    # Add gene metadata (gene symbol) and fetch the table of results. Print first 15 sDEgenes
    rowRanges(se.filt)
    #genesmd <- data.frame(chr = as.character(seqnames(rowRanges(se.filt))), symbol = rowData(se.filt)[, 1], stringsAsFactors = FALSE)
    genesmd <- data.frame(symbol = rownames(res), stringsAsFactors = FALSE)
    fit$genes <- genesmd
    tt <- topTable(fit, coef = 2, n = Inf)
    head(tt, 15)
    
    write.table(tt,file=paste(results_dir,"/",dis,"_DEGs_",metap_method,n_clusters,"_",k,".txt",sep=""),sep="\t",row.names=F, quote=FALSE)
    for_jon <- tt[,c(1,2)]
    for_jon <- for_jon[order(-for_jon$logFC),]
    # write.table(for_jon,file=paste(results_dir,"/",dis,"_DEGs_",metap_method,n_clusters,"_",k,".rnk",sep=""),sep="\t",col.names=F,row.names=F, quote=FALSE)
    write.table(for_jon,file=paste(output_rnk,"/",dis,"_DEGs_",metap_method,n_clusters,"_",k,".rnk",sep=""),sep="\t",col.names=F,row.names=F, quote=FALSE)
    
    # Total number of DE genes (considering p-value)
    sum(tt$P.Value < 0.05)
    
    # When considering adjusted p-value, the number os significantly DEgenes is:
    sum(tt$adj.P.Val < 0.05)
    
    # Q-Q plot
    par(mar=c(1,1,1,1))
    {qqt(fit$t[, 2], df = fit$df.prior + fit$df.residual, main = "", pch = ".", cex = 3, las = 1)
      qqline(fit$t[, 2], col = "red")}
    
    # Volcano plot
    par(mar=c(1,1,1,1))
    par(xpd=T, mar=par()$mar+c(0,0,0,5))
    with(tt, plot(logFC, -log10(adj.P.Val), pch = 20, main = "", xlim=c(-6,7), ylim=c(0,15))) 
    
    
    with(subset(tt, logFC > 1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="green"))
    with(subset(tt, logFC < -1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="orange"))
    with(subset(tt, -log10(adj.P.Val) < 1.3), points(logFC, -log10(adj.P.Val), pch = 20, col="grey"))
    
    legend("bottomright", cex = .75, inset = c(-0.32,0.75), xjust = 2, yjust = 1,pch = 20,c("sDE & overexpressed", "sDE & underexpressed", "sDE", "non-sDE"), col = c("green", "orange", "grey", "black"),bg="white", bty= "n")
    with(subset(tt, -log10(adj.P.Val)>1.3 & abs(logFC)>1), textxy(logFC, -log10(adj.P.Val), labs = symbol, cex=.7, col= "#35445b"))
    par(xpd=F)
    abline(h= 1.3, col = "blue", lty= 2, lwd = 1)
    abline(v= c(-1,1), col = "blue", lty= 2, lwd = 1)
    
    # Select genes with a |logFC| higher than 1
    sum((tt$adj.P.Val < 0.05) & (abs(tt$logFC)) >= 1)
  }
  
  # dev.off()
  
}

wo_metapatients <- c()
for (dis in dis_list){
  # dis <- 'Schizophrenia'
  print(dis)
  cclustes_name <- clusters_list[grep(paste(dis,"_.*_normalized",sep=""),clusters_list)]
  if(length(cclustes_name) == 0){
    cclustes_name <- clusters_list[grep(paste(dis,"_.*_combat",sep=""),clusters_list)]
  }
  print(cclustes_name)
  
  if(length(cclustes_name) == 0){
    wo_metapatients <- append(wo_metapatients, dis)
    next
  }
  clusters_filename = paste(clusters_dir,cclustes_name,sep="")

  # Preparing output directory
  results_dir = paste(outputdir,dis,sep="")
  if(!dir.exists(results_dir)){
     dir.create(results_dir)
  }
  try(dev.off())
  pdf(paste(results_dir,"/",dis,".pdf",sep=""),width = 8.267, height = 11.692)
  par(mar=c(1,1,1,1))
  run_DEAnalysis_for_dis_metapatients(dis,clusters_filename,results_dir)
  dev.off()
  
}

wo_metapatients
# Diseases wo clusters for sDVGs > 0 :  "Keratoconus"   "Schizophrenia" "Smoker" 



#### RUN DEAnalysis for a given disease

# < Script parameters > (adjust them for each run)
dis = "Schizophrenia"

clusters_filename = 'Metapatients/Schizophrenia_metapatients_RNAseq/Schizophrenia_clusters_k2_dm_transcripts.rds'
metap_method = 'dm'

clusters_filename = 'Metapatients/Schizophrenia_metapatients_RNAseq/Schizophrenia_clusters_k2_transcripts.rds'
metap_method = 'allgenes'

# Preparing output directory
results_dir = paste("Metapatients/DEAnalysis/",dis,sep="")
if(!dir.exists(results_dir)){
  dir.create(results_dir)
}
# </>

# Importing metadata, se and dge objects
meta <- readRDS("new_meta_jon_grein_good_quality.rds")
se <- readRDS(paste('DL_RESULTS_GREIN/normalized_counts/',dis,'_se_normalized.rds',sep=""))
dge <- readRDS(paste('DL_RESULTS_GREIN/normalized_counts/',dis,'_dge_normalized.rds',sep=""))

#### DEFINING METAPATIENTS SE & DGE
# Selecting the controls of the disease
controls <- meta[meta$disease == dis,]
controls <- as.character(controls$sample) # 61 samples

# Importing the clusters
clus <- readRDS(clusters_filename) 

# Number of clusters
n_clusters = dim(as.matrix(table(clus)))[1]

# Creating a list with n_clusters elements. Each element contains a list with its samples -always patients
clus_samples <- c()
se_list <- c()
dge_list <- c()
for(k in 1:n_clusters){
  cpatients <- names(clus[clus == k])
  clus_samples <- append(clus_samples, list(cpatients))
  se_list <- append(se_list,list(se[,append(cpatients,controls)]))
  dge_list <- append(dge_list,list(dge[,append(cpatients,controls)]))
}

colnames(se_list[[1]])

for(k in 1:n_clusters){
  se.filt <- se_list[[k]]
  dge.filt.norm <- dge_list[[k]]
  
  # Without correcting
  se.filt$type <- relevel(se.filt$type, ref="control")
  mod <- model.matrix(~ se.filt$type, colData(se.filt))
  mod0 <- model.matrix(~ 1, colData(se.filt))
  pv <- f.pvalue(assays(se.filt)$logCPM.norm, mod, mod0)
  sum(p.adjust(pv, method="fdr") < 0.5)
  
  # Distribution of expected p-values
  par(mfrow=c(1, 1))
  hist(pv, main="Distribution of expected p-values", las=1, ylim = c(0, 8000))
  
  # Correcting for the series
  if (length(unique(se$series)) > 1){
    mod <- model.matrix(~type + series, colData(se.filt))
    mod0 <- model.matrix(~series, colData(se.filt))
    colnames(mod)
    pValues <- f.pvalue(assays(se.filt)$logCPM.norm,mod,mod0)
    sum(p.adjust(pValues, method="fdr") < 0.05)
  } else {
    print("There is only one series for this disease")
  }
  
  # Distribution of expected p-values after correcting for the series
  if (length(unique(se$serie)) > 1){
    hist(pValues, main="Distribution of expected p-values after correcting for series", las=1)
  }
  
  # Mean-variance relationship
  FDRcutoff <- 0.05
  v <- voom(dge.filt.norm, mod, plot = TRUE) # Voom is applied to the normalized and filtered DGEList object
  fit <- lmFit(v, mod)
  fit<- eBayes(fit)
  res <- decideTests(fit, p.value = FDRcutoff)
  summary(res)
  
  # Add gene metadata (gene symbol) and fetch the table of results. Print first 15 sDEgenes
  rowRanges(se.filt)
  #genesmd <- data.frame(chr = as.character(seqnames(rowRanges(se.filt))), symbol = rowData(se.filt)[, 1], stringsAsFactors = FALSE)
  genesmd <- data.frame(symbol = rownames(res), stringsAsFactors = FALSE)
  fit$genes <- genesmd
  tt <- topTable(fit, coef = 2, n = Inf)
  head(tt, 15)
  
  write.table(tt,file=paste(results_dir,"/",dis,"_DEGs_",metap_method,n_clusters,"_",k,".txt",sep=""),sep="\t",row.names=F, quote=FALSE)
  for_jon <- tt[,c(1,2)]
  for_jon <- for_jon[order(-for_jon$logFC),]
  write.table(for_jon,file=paste(results_dir,"/",dis,"_DEGs_",metap_method,n_clusters,"_",k,".rnk",sep=""),sep="\t",col.names=F,row.names=F, quote=FALSE)
  
  # Total number of DE genes (considering p-value)
  sum(tt$P.Value < 0.05)
  
  # When considering adjusted p-value, the number os significantly DEgenes is:
  sum(tt$adj.P.Val < 0.05)
  
  # Q-Q plot
  {qqt(fit$t[, 2], df = fit$df.prior + fit$df.residual, main = "", pch = ".", cex = 3, las = 1)
    qqline(fit$t[, 2], col = "red")}
  
  # Volcano plot
  par(xpd=T, mar=par()$mar+c(0,0,0,5))
  with(tt, plot(logFC, -log10(adj.P.Val), pch = 20, main = "", xlim=c(-6,7), ylim=c(0,15))) 
  
  
  with(subset(tt, logFC > 1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="green"))
  with(subset(tt, logFC < -1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="orange"))
  with(subset(tt, -log10(adj.P.Val) < 1.3), points(logFC, -log10(adj.P.Val), pch = 20, col="grey"))
  
  legend("bottomright", cex = .75, inset = c(-0.32,0.75), xjust = 2, yjust = 1,pch = 20,c("sDE & overexpressed", "sDE & underexpressed", "sDE", "non-sDE"), col = c("green", "orange", "grey", "black"),bg="white", bty= "n")
  with(subset(tt, -log10(adj.P.Val)>1.3 & abs(logFC)>1), textxy(logFC, -log10(adj.P.Val), labs = symbol, cex=.7, col= "#35445b"))
  par(xpd=F)
  abline(h= 1.3, col = "blue", lty= 2, lwd = 1)
  abline(v= c(-1,1), col = "blue", lty= 2, lwd = 1)
  
  # Select genes with a |logFC| higher than 1
  sum((tt$adj.P.Val < 0.05) & (abs(tt$logFC)) >= 1)
}

