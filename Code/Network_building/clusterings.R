#################################################################################
# EXPLORING DISEASE SIMILARITIES in the new diseases 
#
#
# - Performing PAM to obtain disease clusters --> Plot the with PCA
# - Obtain heatmaps of the disease correlations
#
####################################################################################

library(cluster)
library(factoextra)
library(ggplot2)
library(reshape2)

setwd("~/Desktop/ANALYSIS/Network_building/")
set.seed(5)

dismeta <- read.csv("../dl_general_info_umls_ncolors.csv",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
ch1 <- which(dismeta$tis_cat_colors == 'fuchsia')
dismeta$tis_cat_colors[ch1] <- 'magenta3'
ch1 <- which(dismeta$dis_cat_colors == 'grey')
dismeta$dis_cat_colors[ch1] <- 'black'
# Saving the table with the definitive colors
# write.table(dismeta,file="../dl_general_info_umls_new_ncolors.csv",sep='\t',row.names=FALSE, quote = FALSE)

# A) TO USE NORMALIZED DATA + ONLY sDEGs
# alldf <- read.csv("all_genes_all_diseases_logFC_table_normalized.csv",header=T, sep=",",stringsAsFactors = F, row.names = 1)
# sdf <- read.csv2("sDEGs_all_diseases_logFC_table_normalized.csv",header=T, sep=",",stringsAsFactors = F, row.names = 1)

# Preparing the normalized data
# alldf <- alldf[,-dim(alldf)[2]]
# sdf <- sdf[,-dim(sdf)[2]]
# rownames <- row.names(sdf)
# sdf <- sapply(sdf[,1:dim(sdf)[2]],as.numeric)
# row.names(sdf) <- rownames

# B) TO USE UNNORMALIZED DATA + UNION OF sDEGs
# Unnormalized and with the union of sDEGs for sdf
alldf <- read.csv("all_genes_all_diseases_logFC_table.csv",header=T, sep=",",stringsAsFactors = F, row.names = 1)
sdf <- read.csv2("union_sDEGs_all_diseases_logFC_table.csv",header=T, sep=",",stringsAsFactors = F, row.names = 1)
dmdf <- read.csv2("../Disease_gene_variability/all_diseases_dm_ratio.csv",header=T, sep="\t",stringsAsFactors = F, row.names = 1)
dmdf[1:5,1:5]
str(dmdf)
dmdf <- as.data.frame(t(dmdf))

# Preparing the unnormalized data
rownames <- row.names(sdf)
sdf <- sapply(sdf[,1:dim(sdf)[2]],as.numeric)
row.names(sdf) <- rownames

# Preparing the unnormalized data - for the dm difference
to_numeric <- function(x){
  return(as.numeric(as.character(x)))
}
rownames <- row.names(dmdf)
dmdf <- as.data.frame(sapply(dmdf[,1:dim(dmdf)[2]],to_numeric))
row.names(dmdf) <- rownames

PAM_clutering <- function(distmatrix,inputdata,max_k = 15, filename, return_sil = FALSE){
  
  pdf(paste(filename,"_clustering_results.pdf",sep=""))
  
  ave_sil_df <- data.frame(matrix(ncol=2,nrow=0))
  for(k in 2:max_k){
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
  which.max(ave_sil_df$Silhouette)
  best <- ave_sil_df[1,1]
  p <- ggplot(ave_sil_df, aes(N_clusters,Silhouette)) + geom_line() +
    labs(title = paste("Average Silhouette per n clusters (k). Optimal k =  ",as.character(best),sep=""))
    # title(paste("Average Silhouette per nº clusters (k). Optimal k =  ",str(best),sep=""))
  print(p)
  
  dev.off()
  if(return_sil == TRUE){
    return(ave_sil_df$Silhouette)
  }
  
}

#### Obtain Silhouette for all methods ####

# Preparing df to store Silhouette values for all methods
k = 25
sildf <- data.frame(matrix(ncol=10,nrow=k-1))
sildf[,1] = c(2:k)
colnames(sildf) <- c('N_clusters','Euclidean_sde','Euclidean_allg','Pearson_sde','Pearson_allg','Spearman_sde','Spearman_allg','Hamming_sde','Hamming_allg','Jaccard_sde')
jacdf <- read.csv('diseases_jaccard_distance_sDEGs.csv',header=T, sep=",",stringsAsFactors = F, row.names = 1)
ham_sde <- read.csv('diseases_hamming_distance_sDEGs.csv',header=T, sep=",",stringsAsFactors = F, row.names = 1)
ham_allg <- read.csv('diseases_hamming_distance_all_genes.csv',header=T, sep=",",stringsAsFactors = F, row.names = 1)

methods = c('euclidean','pearson','spearman')
used_genesl = c('sDEGs','allgenes')
count = 2
for (method in methods){
  for(used_genes in used_genesl){
    if(used_genes == 'sDEGs'){
      input_data = sdf
    }else{
      input_data = alldf
    }
    # Obtaining the distance matrix
    distmatrix <- get_dist(input_data, method=method, stand=TRUE)
    
    # Performing PAM clustering and plotting the results with PCA
    csil <- PAM_clutering(distmatrix,input_data,k,paste(method,"_",used_genes,sep=""), return_sil = TRUE)
    sildf[,count] = csil
    count = count + 1
    
  }
}

method = 'hamming'
used_genes = 'sDEGs'
distmatrix <- as.dist(t(ham_sde), diag=TRUE, upper=TRUE)
csil <- PAM_clutering(distmatrix,sdf,k,paste(method,"_",used_genes,sep=""), return_sil = TRUE)
sildf[,count] = csil
count = count + 1

used_genes = 'all_genes'
distmatrix <- as.dist(t(ham_allg), diag=TRUE, upper=TRUE)
csil <- PAM_clutering(distmatrix,alldf,k,paste(method,"_",used_genes,sep=""), return_sil = TRUE)
sildf[,count] = csil
count = count + 1

method = 'jaccard'
used_genes = 'sDEGs'
distmatrix <- as.dist(t(jacdf), diag=TRUE, upper=TRUE)
csil <- PAM_clutering(distmatrix,sdf,k,paste(method,"_",used_genes,sep=""), return_sil = TRUE)
sildf[,count] = csil
count = count + 1

### Plotting all Silhouette distribution for each method in one plot
colnames(sildf) <- c('N_clusters','Euclidean_sde','Euclidean_allg','Pearson_sde','Pearson_allg','Spearman_sde','Spearman_allg','Hamming_sde','Hamming_allg','Jaccard_sde')
melted <- melt(sildf,id.var='N_clusters')
tgp <- ggplot(melted,aes(x=N_clusters,y=value, col=variable)) + geom_line() +
  xlab("Number of clusters") + ylab("Average Silhouette") +
  labs(title = 'Silhouette distribution with different distances and gene sets', color="Distance and gene set") +
  theme(plot.title = element_text(hjust = 0.5))
tgp

# Same plot without the Jaccard distance
other <- sildf
other$Jaccard_sde <- NULL
other_melted <- melt(other,id.var='N_clusters')
tgp <- ggplot(other_melted,aes(x=N_clusters,y=value, col=variable)) + geom_line() +
  xlab("Number of clusters") + ylab("Average Silhouette") +
  labs(title = 'Silhouette distribution with different distances and gene sets', color="Distance and gene set") +
  theme(plot.title = element_text(hjust = 0.5))
tgp


#### EUCLIDEAN, PEARSON & SPEARMAN FOR sDEGs
# Obtaining distance matrix
# distmatrix <- dist(alldf)
# distmatrix <- get_dist(alldf, method='spearman')
method='euclidean'
method='pearson'
method='spearman'
used_genes = 'sDEGs'

# Obtaining the distance matrix
distmatrix <- get_dist(sdf, method=method, stand=TRUE)

# Performing PAM clustering and plotting the results with PCA
PAM_clutering(distmatrix,sdf,25,paste(method,"_",used_genes,sep=""))

# Plot - Obtaining a heatmap using the distance matrix
p <- fviz_dist(distmatrix, gradient = list(low = "#FC4E07", mid = "white", high = "#00AFBB")) +
  theme(axis.text=element_text(size=8))
p

cdis_order <- gsub('-',"",as.character(p$data$Var1[1:45]))
correct <- data.frame("disease_name" = cdis_order)
to_merge <- dismeta[,c('disease_name','dis_cat_colors','tis_cat_colors')]
correct <- merge(correct,to_merge,by='disease_name', sort=FALSE)


p + theme(axis.text.x = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Disease category", sep=""))
p + theme(axis.text.x = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Tissue category", sep=""))

# FILENAME: sDEGs_spearman_dis.png

#### EUCLIDEAN, PEARSON & SPEARMAN FOR ALL GENES
# distmatrix <- dist(alldf)
method='euclidean'
method='pearson'
method='spearman'
used_genes = 'all_genes'
distmatrix <- get_dist(alldf, method=method, stand=TRUE)
PAM_clutering(distmatrix,alldf,25,paste(method,"_",used_genes,sep=""))

# Plot
p <- fviz_dist(distmatrix, gradient = list(low = "#FC4E07", mid = "white", high = "#00AFBB")) +
  theme(axis.text=element_text(size=8))
p

cdis_order <- gsub('-',"",as.character(p$data$Var1[1:45]))
correct <- data.frame("disease_name" = cdis_order)
to_merge <- dismeta[,c('disease_name','dis_cat_colors','tis_cat_colors')]
correct <- merge(correct,to_merge,by='disease_name', sort=FALSE)

p + theme(axis.text.x = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Disease category", sep=""))
p + theme(axis.text.x = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Tissue category", sep=""))

# FILENAME: allgenes_spearman_dis.png

### Using other distances to perform the clusterings
jacdf <- read.csv('diseases_jaccard_distance_sDEGs.csv',header=T, sep=",",stringsAsFactors = F, row.names = 1)
ham_sde <- read.csv('diseases_hamming_distance_sDEGs.csv',header=T, sep=",",stringsAsFactors = F, row.names = 1)
ham_allg <- read.csv('diseases_hamming_distance_all_genes.csv',header=T, sep=",",stringsAsFactors = F, row.names = 1)

colnames(jacdf)
row.names(sdf)

#### JACARD DISTANCE sDEGs
method='jaccard'
used_genes = 'sDEGs'
distmatrix <- as.dist(t(jacdf), diag=TRUE, upper=TRUE)
PAM_clutering(distmatrix,sdf,25,paste(method,"_",used_genes,sep=""))

p <- fviz_dist(distmatrix, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) +
  theme(axis.text=element_text(size=8))
p

cdis_order <- gsub('-',"",as.character(p$data$Var1[1:52]))
cdis_order[which(cdis_order == 'HIVARTtreated')] <- 'HIVART-treated'
correct <- data.frame("disease_name" = cdis_order)
to_merge <- dismeta[,c('disease_name','dis_cat_colors','tis_cat_colors')]
correct <- merge(correct,to_merge,by='disease_name', sort=FALSE)

p + theme(axis.text.x = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Disease category", sep=""))
p + theme(axis.text.x = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Tissue category", sep=""))

#### HAMMING DISTANCE sDEGs
method='hamming'
used_genes = 'sDEGs'
distmatrix <- as.dist(t(ham_sde), diag=TRUE, upper=TRUE)
PAM_clutering(distmatrix,sdf,25,paste(method,"_",used_genes,sep=""))
p <- fviz_dist(distmatrix, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) +
  theme(axis.text=element_text(size=8)) 

cdis_order <- gsub('-',"",as.character(p$data$Var1[1:52]))
cdis_order[which(cdis_order == 'HIVARTtreated')] <- 'HIVART-treated'
correct <- data.frame("disease_name" = cdis_order)
to_merge <- dismeta[,c('disease_name','dis_cat_colors','tis_cat_colors')]
correct <- merge(correct,to_merge,by='disease_name', sort=FALSE)

p + theme(axis.text.x = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Disease category", sep=""))
p + theme(axis.text.x = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Tissue category", sep=""))

#### HAMMING DISTANCE ALL GENES
method='hamming'
used_genes = 'all_genes'
distmatrix <- as.dist(t(ham_allg), diag=TRUE, upper=TRUE)
PAM_clutering(distmatrix,alldf,25,paste(method,"_",used_genes,sep=""))
p <- fviz_dist(distmatrix, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) +
  theme(axis.text=element_text(size=8)) 

cdis_order <- gsub('-',"",as.character(p$data$Var1[1:52]))
cdis_order[which(cdis_order == 'HIVARTtreated')] <- 'HIVART-treated'
cdis_order[which(cdis_order == 'HIVART.treated')] <- 'HIVART-treated'
correct <- data.frame("disease_name" = cdis_order)
to_merge <- dismeta[,c('disease_name','dis_cat_colors','tis_cat_colors')]
correct <- merge(correct,to_merge,by='disease_name', sort=FALSE)

p + theme(axis.text.x = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Disease category", sep=""))
p + theme(axis.text.x = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Tissue category", sep=""))

#### EUCLIDEAN, PEARSON & SPEARMAN FOR DEVGs CONSIDERING THE INTERSECTION
# Obtaining distance matrix
# distmatrix <- dist(alldf)
# distmatrix <- get_dist(alldf, method='spearman')
method='euclidean'
method='pearson'
method='spearman'
used_genes = 'DVGs'

# keep only genes that are in the intersection
is_in_all_diseases <- function(x){
  v <- is.na(x)
  return(all.equal(v,rep(FALSE,length(v))))
  # if(all.equal(v,rep(FALSE,length(v)))){
  #   return(TRUE)
  # }else{
  #   return(FALSE)
  # }
}

to_keep <- as.data.frame(sapply(dmdf[,1:dim(dmdf)[2]],is_in_all_diseases))
colnames(to_keep) <- c("result")
to_keep <- which(to_keep$result == TRUE)
length(to_keep) # 2974 genes

intercdmdf <- dmdf[,to_keep]

# Obtaining the distance matrix
distmatrix <- get_dist(intercdmdf, method=method, stand=TRUE)

# Performing PAM clustering and plotting the results with PCA
PAM_clutering(distmatrix,intercdmdf,25,paste(method,"_",used_genes,sep=""))

# Plot - Obtaining a heatmap using the distance matrix
p <- fviz_dist(distmatrix, gradient = list(low = "#FC4E07", mid = "white", high = "#00AFBB")) +
  theme(axis.text=element_text(size=8))
p

cdis_order <- gsub('-',"",as.character(p$data$Var1[1:52]))
cdis_order[which(cdis_order == 'HIVARTtreated')] <- 'HIVART-treated'
correct <- data.frame("disease_name" = cdis_order)
to_merge <- dismeta[,c('disease_name','dis_cat_colors','tis_cat_colors')]
correct <- merge(correct,to_merge,by='disease_name', sort=FALSE)


p + theme(axis.text.x = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Disease category", sep=""))
p + theme(axis.text.x = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Tissue category", sep=""))

# FILENAME: allgenes_spearman_dis.png

#### EUCLIDEAN, PEARSON & SPEARMAN FOR DEVGs PUTTING ZEROS IF NEEDED
# Obtaining distance matrix
# distmatrix <- dist(alldf)
# distmatrix <- get_dist(alldf, method='spearman')
method='euclidean'
method='pearson'
method='spearman'
used_genes = 'DVGs'

dmdf[1:5,1:5]
dmdf[is.na(dmdf)] <- 0
dmdf[1:5,1:5]

# Obtaining the distance matrix
distmatrix <- get_dist(dmdf, method=method, stand=TRUE)

# Performing PAM clustering and plotting the results with PCA
PAM_clutering(distmatrix,dmdf,25,paste(method,"_",used_genes,sep=""))

# Plot - Obtaining a heatmap using the distance matrix
p <- fviz_dist(distmatrix, gradient = list(low = "#FC4E07", mid = "white", high = "#00AFBB")) +
  theme(axis.text=element_text(size=8))
p

cdis_order <- gsub('-',"",as.character(p$data$Var1[1:52]))
cdis_order[which(cdis_order == 'HIVARTtreated')] <- 'HIVART-treated'
cdis_order[which(cdis_order == 'HIVART.treated')] <- 'HIVART-treated'
correct <- data.frame("disease_name" = cdis_order)
to_merge <- dismeta[,c('disease_name','dis_cat_colors','tis_cat_colors')]
correct <- merge(correct,to_merge,by='disease_name', sort=FALSE)

p + theme(axis.text.x = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Disease category", sep=""))
p + theme(axis.text.x = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Tissue category", sep=""))

#### PAIRWISE SPEARMAN DISTANCE WITH ALLGENES, SDEGs and DM.
#### ALL GENES
method='pairwise_spearman'
used_genes = 'all_genes'
corr_table <- read.csv('diseases_pairwise_spearman_all_genes.csv',header=T, sep=",",stringsAsFactors = F, row.names = 1)
corr_table <- 1-corr_table

distmatrix <- as.dist(t(corr_table), diag=TRUE, upper=TRUE)
# PAM_clutering(distmatrix,alldf,25,paste(method,"_",used_genes,sep=""))
# Plot
p <- fviz_dist(distmatrix, gradient = list(low = "#FC4E07", mid = "white", high = "#00AFBB")) +
  theme(axis.text=element_text(size=8))
p

cdis_order <- gsub('-',"",as.character(p$data$Var1[1:52]))
cdis_order[which(cdis_order == 'HIVARTtreated')] <- 'HIVART-treated'
cdis_order[which(cdis_order == 'HIVART.treated')] <- 'HIVART-treated'
correct <- data.frame("disease_name" = cdis_order)
to_merge <- dismeta[,c('disease_name','dis_cat_colors','tis_cat_colors')]
correct <- merge(correct,to_merge,by='disease_name', sort=FALSE)

p + theme(axis.text.x = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Disease category", sep=""))
p + theme(axis.text.x = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Tissue category", sep=""))

#### sDEGs
method='pairwise_spearman'
used_genes = 'sDEGs'
corr_table <- read.csv('diseases_pairwise_spearman_sDEGs.csv',header=T, sep=",",stringsAsFactors = F, row.names = 1)
corr_table <- 1-corr_table

distmatrix <- as.dist(t(corr_table), diag=TRUE, upper=TRUE)
# PAM_clutering(distmatrix,alldf,25,paste(method,"_",used_genes,sep=""))
# Plot
p <- fviz_dist(distmatrix, gradient = list(low = "#FC4E07", mid = "white", high = "#00AFBB")) +
  theme(axis.text=element_text(size=8))
p

cdis_order <- gsub('-',"",as.character(p$data$Var1[1:52]))
cdis_order[which(cdis_order == 'HIVARTtreated')] <- 'HIVART-treated'
cdis_order[which(cdis_order == 'HIVART.treated')] <- 'HIVART-treated'
correct <- data.frame("disease_name" = cdis_order)
to_merge <- dismeta[,c('disease_name','dis_cat_colors','tis_cat_colors')]
correct <- merge(correct,to_merge,by='disease_name', sort=FALSE)

p + theme(axis.text.x = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Disease category", sep=""))
p + theme(axis.text.x = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Tissue category", sep=""))

#### DM
method='pairwise_spearman'
used_genes = 'DM'
corr_table <- read.csv('diseases_pairwise_spearman_DM.csv',header=T, sep=",",stringsAsFactors = F, row.names = 1)
corr_table <- 1-corr_table

distmatrix <- as.dist(t(corr_table), diag=TRUE, upper=TRUE)
# PAM_clutering(distmatrix,alldf,25,paste(method,"_",used_genes,sep=""))
# Plot
p <- fviz_dist(distmatrix, gradient = list(low = "#FC4E07", mid = "white", high = "#00AFBB")) +
  theme(axis.text=element_text(size=8))
p

cdis_order <- gsub('-',"",as.character(p$data$Var1[1:52]))
cdis_order[which(cdis_order == 'HIVARTtreated')] <- 'HIVART-treated'
cdis_order[which(cdis_order == 'HIVART.treated')] <- 'HIVART-treated'
correct <- data.frame("disease_name" = cdis_order)
to_merge <- dismeta[,c('disease_name','dis_cat_colors','tis_cat_colors')]
correct <- merge(correct,to_merge,by='disease_name', sort=FALSE)

p + theme(axis.text.x = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$dis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Disease category", sep=""))
p + theme(axis.text.x = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  theme(axis.text.y = element_text(colour = as.character(correct$tis_cat_colors), size = 7, face = 'bold')) +
  labs(title = paste("Disease similarity based on ",method," distance - ",used_genes," | Tissue category", sep=""))




