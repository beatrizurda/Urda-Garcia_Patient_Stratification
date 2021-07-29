#####################################################################################
# DEFINING META-PATIENTS BASED ON COUNTSTHE NORMALIZED (AND BATCH EFFECT REMOVED) COUNTS
#
#   We use PAM or WARD algorithms to find clusters for each disease based on their
#   normalized (and if necessary) batch effect removal counts
# 
# Beatriz Urda García 2020
######################################################################################

setwd("~/Desktop/ANALYSIS")

library(cluster)
library(factoextra)
library(ggplot2)
library(SummarizedExperiment)
library(dplyr)
library(fitdistrplus)
library(data.table)


meta <- readRDS("new_meta_jon_grein_good_quality.rds")

dismeta <- read.csv("dl_general_info_umls_ncolors.csv",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
ch1 <- which(dismeta$tis_cat_colors == 'fuchsia')
dismeta$tis_cat_colors[ch1] <- 'magenta3'

ch1 <- which(dismeta$dis_cat_colors == 'grey')
dismeta$dis_cat_colors[ch1] <- 'black'

PAM_clutering <- function(distmatrix,inputdata,max_k = 15, filename, return_sil = FALSE, summ_filename = FALSE){
  
  pdf(paste(filename,"_clustering_results.pdf",sep=""))
  
  ave_sil_df <- data.frame(matrix(ncol=2,nrow=0))
  print('Maximum k used')
  cmax_k <- min(max_k,dim(inputdata)[1]-1)
  print(cmax_k)
  if(cmax_k == 1){
    return(NA)
  }
  for(k in 2:cmax_k){
    pamres <- pam(dist(distmatrix), k, diss=TRUE)
    str(si <- silhouette(pamres))
    (ssi <- summary(si))
    plot(si)
    pamres$data <- inputdata
    cp <- fviz_cluster(pamres) + ggtitle(paste("PAM ",filename," - k = ",k,sep=""))
    print(cp)
    avg_sil <- mean(si[,3]) # Average silhouette
    ave_sil_df <- rbind(ave_sil_df, c(k, avg_sil))
  }
  colnames(ave_sil_df) <- c('N_clusters', 'Silhouette')
  best <- ave_sil_df[which.max(ave_sil_df$Silhouette),1]
  p <- ggplot(ave_sil_df, aes(N_clusters,Silhouette)) + geom_line() +
    labs(title = paste("Average Silhouette per n clusters (k). Optimal k =  ",as.character(best),sep=""))
  # title(paste("Average Silhouette per nº clusters (k). Optimal k =  ",str(best),sep=""))
  print(p)
  if (summ_filename != FALSE){
    write(paste(dis,best,sep='\t'),file=summ_filename,append=TRUE)
  }
  p2 <- fviz_dist(distmatrix, gradient = list(low = "#FC4E07", mid = "white", high = "#00AFBB")) +
    theme(axis.text=element_text(size=8))
  print(p2)
  
  dev.off()
  if(return_sil == TRUE){
    return(ave_sil_df$Silhouette)
  }
  
}

PAM_clutering_given_k <- function(distmatrix,inputdata,k = 15, filename){
  
  pdf(paste(filename,"_bestk_clustering.pdf",sep=""))
  
  pamres <- pam(dist(distmatrix), k, diss=TRUE)
  str(si <- silhouette(pamres))
  (ssi <- summary(si))
  plot(si)
  pamres$data <- inputdata
  cp <- fviz_cluster(pamres) + ggtitle(paste("PAM ",filename," - k = ",k,sep=""))
  print(cp)
  avg_sil <- mean(si[,3]) # Average silhouette
  
  p2 <- fviz_dist(distmatrix, gradient = list(low = "#FC4E07", mid = "white", high = "#00AFBB")) +
    theme(axis.text=element_text(size=8))
  print(p2)
  
  dev.off()
  return(pamres$clustering)
}

set.seed(5)

combat_files <- list.files(path = "DL_RESULTS_GREIN/after_combat_counts/")
combat_diseases <- gsub("_combat_counts.rds","",combat_files)

normalized_files <- list.files(path = "DL_RESULTS_GREIN/normalized_counts/", pattern = "_se_normalized.rds$")
disease_list <- gsub("_se_normalized.rds","",normalized_files)
disease_list <- disease_list

### OBTAIN METAPATIENTS FOR ALL DISEASES USING LAST COUNTS AND SPEARMAN DISTANCE - OBTAIN BEST SILHOUETTE

k = 25
sildf <- data.frame(matrix(ncol=length(disease_list) + 1,nrow=k-1))
sildf[,1] = c(2:k)
colnames(sildf) <- append('N_clusters',disease_list)
count = 2

fileConn<-file("Metapatients/with_entire_count_matrix/Summary_best_silhouette_metapatients.txt")
writeLines(paste("Disease","Best_silhouette",sep='\t'), fileConn)
close(fileConn)
  
for(dis in disease_list){
  
  print(dis)
  
  cmeta <- meta[meta$disease == dis,]
  cpatients <- cmeta[cmeta$type == 'control',]
  ccontrols <- as.character(cpatients$sample)
  
  # Get counts
  if(dis %in% combat_diseases){
    counts <- readRDS(paste("DL_RESULTS_GREIN/after_combat_counts/",dis,"_combat_counts.rds",sep=""))
    used_genes='combat_counts'
    # Removing the controls
    to_remove <- which(colnames(counts) %in% ccontrols)
    counts <- counts[,-to_remove]
    # hey <- intersect(colnames(counts), ccontrols) # Checking
  }else{
    cobject <- readRDS(paste("DL_RESULTS_GREIN/normalized_counts/",dis,"_se_normalized.rds",sep=""))
    used_genes='normalized_counts'
    # Removing the controls
    counts <- assays(cobject)$logCPM.norm
    to_remove <- which(colnames(counts) %in% ccontrols)
    counts <- counts[,-to_remove]
    # hey <- intersect(colnames(counts), ccontrols) # Checking
  }
  counts <- t(counts)
  print(dim(counts))
  
  # Apply PAM algorithm to obtain metapatients
  method='spearman'

  # Obtaining the distance matrix
  distmatrix <- get_dist(counts, method=method, stand=TRUE)
  
  # Performing PAM clustering and plotting the results with PCA
  filen <- paste("Metapatients/with_entire_count_matrix/all_k/",dis,"_",used_genes,sep="")
  csil <- PAM_clutering(distmatrix,counts,k,filen,return_sil = TRUE, summ_filename = 'Metapatients/with_entire_count_matrix/Summary_best_silhouette_metapatients.txt')
  csil_len <- length(csil)
  if(csil_len < k-1){
    csil <- append(csil,rep(NA,(k-1)-csil_len))
  }
  sildf[,count] = csil
  count = count + 1
  
  # Plot - Obtaining a heatmap using the distance matrix
  # p <- fviz_dist(distmatrix, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) +
  #   theme(axis.text=element_text(size=8))
  # p
}

### Plotting all Silhouette distribution for each method in one plot
melted <- melt(sildf,id.var='N_clusters')
tgp <- ggplot(melted,aes(x=N_clusters,y=value, col=variable)) + geom_line() +
  xlab("Number of clusters") + ylab("Average Silhouette") +
  labs(title = 'Silhouette distribution in patient clustering', color="Distance and gene set")
tgp

# Hacer un plot del máximo

### OBTAIN CLUSTERING FOR ALL DISEASES AND BEST K
best_sil <- read.csv('Metapatients/with_entire_count_matrix/Summary_best_silhouette_metapatients.txt',header=T, sep="\t",stringsAsFactors = F)
dim(best_sil) # Falta "SessileSerratedPolyp" because it only has two patients
cdisease_list <- intersect(disease_list, best_sil$Disease)
k <- 1
for(dis in cdisease_list){
  # dis = disease_list[1]
  
  print(dis)
  cbestk <- best_sil$Best_silhouette[k]
    
  cmeta <- meta[meta$disease == dis,]
  cpatients <- cmeta[cmeta$type == 'control',]
  ccontrols <- as.character(cpatients$sample)
  
  # Get counts
  if(dis %in% combat_diseases){
    counts <- readRDS(paste("DL_RESULTS_GREIN/after_combat_counts/",dis,"_combat_counts.rds",sep=""))
    used_genes='combat_counts'
    # Removing the controls
    to_remove <- which(colnames(counts) %in% ccontrols)
    counts <- counts[,-to_remove]
    # hey <- intersect(colnames(counts), ccontrols) # Checking
  }else{
    cobject <- readRDS(paste("DL_RESULTS_GREIN/normalized_counts/",dis,"_se_normalized.rds",sep=""))
    used_genes='normalized_counts'
    # Removing the controls
    counts <- assays(cobject)$logCPM.norm
    to_remove <- which(colnames(counts) %in% ccontrols)
    counts <- counts[,-to_remove]
    # hey <- intersect(colnames(counts), ccontrols) # Checking
  }
  counts <- t(counts)
  print(dim(counts))
  
  # Apply PAM algorithm to obtain metapatients
  method='spearman'
  
  # Obtaining the distance matrix
  distmatrix <- get_dist(counts, method=method, stand=TRUE)
  
  # Performing PAM clustering and plotting the results with PCA
  filen <- paste("Metapatients/with_entire_count_matrix/best_k/",dis,"_",used_genes,sep="")
  clusters <- PAM_clutering_given_k(distmatrix,counts,k = cbestk, filen)
  saveRDS(clusters, file = paste("Metapatients/with_entire_count_matrix/final_clusters/",dis,"_",used_genes,".rds",sep=""))
  k <- k + 1 
  
}


### OBTAIN CLUSTERING FOR A GIVEN DISEASE AND K
dis = 'Schizophrenia'

print(dis)

cmeta <- meta[meta$disease == dis,]
cpatients <- cmeta[cmeta$type == 'control',]
ccontrols <- as.character(cpatients$sample)

# Get counts
if(dis %in% combat_diseases){
  counts <- readRDS(paste("DL_RESULTS_GREIN/after_combat_counts/",dis,"_combat_counts.rds",sep=""))
  used_genes='combat_counts'
  # Removing the controls
  to_remove <- which(colnames(counts) %in% ccontrols)
  counts <- counts[,-to_remove]
  # hey <- intersect(colnames(counts), ccontrols) # Checking
}else{
  cobject <- readRDS(paste("DL_RESULTS_GREIN/normalized_counts/",dis,"_se_normalized.rds",sep=""))
  used_genes='normalized_counts'
  # Removing the controls
  counts <- assays(cobject)$logCPM.norm
  to_remove <- which(colnames(counts) %in% ccontrols)
  counts <- counts[,-to_remove]
  # hey <- intersect(colnames(counts), ccontrols) # Checking
}
counts <- t(counts)
print(dim(counts))

# Apply PAM algorithm to obtain metapatients
method='spearman'

# Obtaining the distance matrix
distmatrix <- get_dist(counts, method=method, stand=TRUE)

# Performing PAM clustering and plotting the results with PCA
filen <- paste("Metapatients/single_k/",dis,"_",used_genes,sep="")
clusters <- PAM_clutering_given_k(distmatrix,counts,k = 2, filen)
saveRDS(clusters, file = "Schizophrenia_clusters_k2_transcripts.rds")


################# OBTAIN CLUSTERING FOR A GIVEN DISEASE USING WARD ALGORITHM #################
library(ggplot2)
library(pvclust)
library(factoextra)
library(gplots) 

library(cluster)
library(factoextra)
library(ggplot2)
library(SummarizedExperiment)
library(dplyr)
library(fitdistrplus)
library(data.table)


meta <- readRDS("new_meta_jon_grein_good_quality.rds")

dismeta <- read.csv("dl_general_info_umls_ncolors.csv",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
ch1 <- which(dismeta$tis_cat_colors == 'fuchsia')
dismeta$tis_cat_colors[ch1] <- 'magenta3'

ch1 <- which(dismeta$dis_cat_colors == 'grey')
dismeta$dis_cat_colors[ch1] <- 'black'

set.seed(5)

combat_files <- list.files(path = "DL_RESULTS_GREIN/after_combat_counts/")
combat_diseases <- gsub("_combat_counts.rds","",combat_files)

normalized_files <- list.files(path = "DL_RESULTS_GREIN/normalized_counts/", pattern = "_se_normalized.rds$")
disease_list <- gsub("_se_normalized.rds","",normalized_files)
disease_list <- disease_list

dis = 'Schizophrenia'
inputname <- paste(dis,"metapatients_ward")

print(dis)

cmeta <- meta[meta$disease == dis,]
cpatients <- cmeta[cmeta$type == 'control',]
ccontrols <- as.character(cpatients$sample)

# Get counts
if(dis %in% combat_diseases){
  counts <- readRDS(paste("DL_RESULTS_GREIN/after_combat_counts/",dis,"_combat_counts.rds",sep=""))
  used_genes='combat_counts'
  # Removing the controls
  to_remove <- which(colnames(counts) %in% ccontrols)
  counts <- counts[,-to_remove]
  # hey <- intersect(colnames(counts), ccontrols) # Checking
}else{
  cobject <- readRDS(paste("DL_RESULTS_GREIN/normalized_counts/",dis,"_se_normalized.rds",sep=""))
  used_genes='normalized_counts'
  # Removing the controls
  counts <- assays(cobject)$logCPM.norm
  to_remove <- which(colnames(counts) %in% ccontrols)
  counts <- counts[,-to_remove]
  # hey <- intersect(colnames(counts), ccontrols) # Checking
}
counts <- t(counts)
print(dim(counts))

# Apply PAM algorithm to obtain metapatients
method='spearman'

# Use ward
# ccounts <- counts[,1:100]
h1 <- heatmap.2(ccounts,scale="none",
                hclustfun = function(data=ccounts){
                  return(hclust(dist(data), method="ward.D2")) # default: euclidean
                },
                # col=color.palette,breaks=palette.breaks,
                dendrogram = "row",trace="none",
                Colv=FALSE,cexRow = 0.6,cexCol = 0.4,"key"=FALSE,margins = c(15,18),
                # RowSideColors = as.vector(parent_col$Color),
                xlab="Ward2 clustering | euclidean distance | all diseases",
                ylab="Reactome pathways")
h1

cluster_eucl <-pvclust(ccounts,method.dist="euclidean", method.hclust = "ward.D2", nboot=10000)   # 10000
pdf(file=paste("Metapatients/Ward/",inputname,".pdf", sep=""))
plot(cluster_eucl,hang = -1, cex = 0.5)
pvrect(cluster_eucl)
dev.off()


# ANALYZING OVERLAP BETWEEN WARD AND PAM SCHIZOPHRENIA PATIENTS
ward1 <- c("GSM2081503","GSM2081499","GSM2081504","GSM2081483","GSM1556374","GSM1556375","GSM2081501",
           "GSM2081502","GSM2081481","GSM2081482","GSM2081484","GSM2081480","GSM2081485")

ward2 <- c("GSM2439188","GSM1556376","GSM1556377","GSM2081465","GSM2081461","GSM2081462","GSM2081460",
           "GSM2081454","GSM2081455","GSM2081452","GSM2081457","GSM2081451","GSM2081463",
           "GSM1556379","GSM2439190","GSM1556378","GSM2439189")

length(ward1) # 13
length(ward2) # 17

length(clusters) # Resultado de ejectuat PAM algorithm para Schizophrenia k = 2 con counts
pam1 <- names(clusters[clusters == 1])
pam2 <- names(clusters[clusters == 2])

# Checking that everything is ok!
length(unique(intersect(union(ward1,ward2),union(pam1,pam2)))) # 30 --> Perfect!

library(VennDiagram)
ven <- venn.diagram(x=list(pam1,pam2,ward1,ward2), 
                    category.names = c('pam1','pam2','ward1','ward2'),
                    filename=NULL)
grid.draw(ven)








