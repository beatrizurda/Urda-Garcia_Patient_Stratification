###################################################################
# DEA FOR RANDOM METAPATIENTS  
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
# Beatriz Urda García 2021
######################################################################



library(knitr)
library(SummarizedExperiment)
library(devtools)
library(ggbiplot)
library(edgeR)
library(sva)

start_time <- Sys.time()

# < Script parameters for the meta-patients obtained with all genes> (adjust them for each run) 
# To test
# setwd("~/Desktop/ANALYSIS")
# dis <- 'AdenomatousPolyps'
# nrandom = 10
# clusters_dir <- paste0('Metapatients/with_entire_count_matrix/Randomization/final_clusters/',dis,"/")
# clusters_list <- list.files(clusters_dir);  clusters_list <- clusters_list[order(clusters_list)];
# clusters_list <- clusters_list[1:nrandom];  head(clusters_list)
# outputdir <- paste0("Metapatients/with_entire_count_matrix/Randomization/DEA_results/",dis,"/")

args <- commandArgs(trailingOnly = TRUE)
nrandom = args[1]
dis = args[2]
print(dis)
clusters_dir <- paste0('Randomization/final_clusters/',dis,"/")
clusters_list <- list.files(clusters_dir);  clusters_list <- clusters_list[order(clusters_list)];
clusters_list <- clusters_list[1:nrandom];  head(clusters_list)
outputdir <- paste0("Randomization/DEA_results/",dis,"/")

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
    genesmd <- data.frame(symbol = rownames(res), stringsAsFactors = FALSE)
    fit$genes <- genesmd
    tt <- topTable(fit, coef = 2, n = Inf)
    tt <- tt[, c("symbol", "logFC", "adj.P.Val")]
    
    write.table(tt,file=paste(results_dir,"_DEGs_",k,".txt",sep=""),sep="\t",row.names=F, quote=FALSE)
  }
  
  # dev.off()
  
}

### RUN DEA for all the metapatients of a given disease
for (cluster in clusters_list){
  cclustes_name = cluster; print(cclustes_name)
  clusters_filename = paste(clusters_dir,cclustes_name,sep="")
  
  # Preparing output directory
  results_dir = paste(outputdir,cluster,sep="")
  
  # try(dev.off())
  # pdf(paste0(results_dir,".pdf"),width = 8.267, height = 11.692)
  # par(mar=c(1,1,1,1))
  run_DEAnalysis_for_dis_metapatients(dis,clusters_filename,results_dir)
  # dev.off()
  
}

end_time <- Sys.time()
print(end_time - start_time)



