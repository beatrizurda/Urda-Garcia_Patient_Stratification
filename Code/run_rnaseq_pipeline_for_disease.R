#########################################################################################
# RUN RNA-SEQ PIPELINE FOR A GIVEN DISEASE FROM SE OBJECT (grein counts + metadata)
# 
# FINAL VERSION STORED IN ~/Desktop/ANALYSIS
#
# Beatriz Urda García 2020
#########################################################################################

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

# For Functional Enrichment Analysis
# library(dplyr)
# library(org.Hs.eg.db)
# library(GSEABase)
# library(GSVAdata)
# library(GSVA)
# library(Category)

# File to test in this computer
# input_filename <- "ParkinsonsDisease_grein_se.rds"

# File to test in Mare Nostrum
# input_filename <- "GREIN_SE_good_quality/ParkinsonsDisease_grein_se.rds"
# input_filename <- "SE_icd9_good_quality/153_se.rds"

# OPTION A) TO RUN THE PIPELINE AT THE DISEASE LEVEL
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 1){
  # OPTION A) TO RUN THE PIPELINE AT THE DISEASE LEVEL
  input_filename <- paste("GREIN_SE_good_quality/",args[1],sep="")
  
  disease_name <- gsub('_grein.*',"",input_filename)
  disease_name <- gsub('GREIN_SE_good_quality/',"",disease_name)
  
  results_dir <- paste("DL_RESULTS_GREIN/",disease_name,sep="")
  parent_res_dir <- "DL_RESULTS_GREIN/"
  
  # If needed, create a folder for the disease
  if(!disease_name %in% list.files("DL_RESULTS_GREIN/")){
    dir.create(results_dir)
  }
  skip_unnecessary <- TRUE
  do_hierarchichal <- FALSE
}else if (length(args) > 1){
  if(args[2] == 'icd9'){
    # OPTION B) TO RUN THE PIPELINE AT THE ICD9 LEVEL
    print("runnning at the icd level")
    input_filename <- paste("SE_icd9_good_quality/",args[1],sep="")
    
    disease_name <- gsub('_se.rds',"",input_filename)
    disease_name <- gsub('SE_icd9_good_quality/',"",disease_name)
    
    results_dir <- paste("ICD9_RESULTS_GREIN/",disease_name,sep="")
    parent_res_dir <- "ICD9_RESULTS_GREIN/"
    
    # If needed, create a folder for the disease
    if(!disease_name %in% list.files("ICD9_RESULTS_GREIN/")){
      dir.create(results_dir)
    }
    
    skip_unnecessary <- TRUE
    do_hierarchichal <- FALSE
  }
}

se <- readRDS(input_filename)
head(se)

# Exploring the metadata (colData)
dim(colData(se))
colData(se)[1:5,]
nsamples <- length(se$sample)

# We explore the feature (gene) data
rowData(se)
rowRanges(se)

# Data exploration
metadata <- as.data.frame(colData(se))
fplot <- ggplot(as.data.frame(metadata))


# Exploratory plots about the metadata
pdf(file=paste(results_dir,"/",disease_name,"_analysis.pdf",sep=""))
# pdf(file=paste(disease_name,"_anlaysis.pdf",sep="")) # in the mare
fplot + geom_bar(aes(type, fill = type)) + xlab("Type") + ylab("Number of Samples") + 
  ggtitle("Sample type distribution")

fplot + geom_bar(aes(series, fill = type)) + xlab("Series") + 
  ylab("Number of Samples") + ggtitle("Controls and patients across series")

fplot + geom_bar(aes(tissue, fill = series)) + xlab("Tissue") + 
  ylab("Number of Samples") + ggtitle("Distribution of samples across tissue and series") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

fplot + geom_bar(aes(type, fill = series)) + xlab("Disease") + 
  ylab("Number of Samples") + ggtitle("Distribution of samples across disease and series") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

fplot + geom_bar(aes(series, fill = tissue)) + xlab("Serie") + 
  ylab("Number of Samples") + ggtitle("Distribution of samples across tissue and series")

# Quality assessment and normalization

## Preprocessing data
dge <- DGEList(counts = assays(se)$counts, genes = as.data.frame(mcols(se)), group = se$type)
names(dge)

assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
assays(se)$logCPM[1:5, 1:5]

dge$counts[1:5, 1:5]
head(dge$samples$lib.size)

## Examining the sequencing depth
ord <- order(dge$samples$lib.size)
head(dge$samples$lib.size)

# Barplot patients/controls
barplot(dge$samples$lib.size[ord]/1e+06, las = 1, main="Barplot - library size", ylab = "Million of counts", xlab = "Samples", col = c("blue","red")[(se$type[ord] == "patient") +1], border = NA)
legend("topleft", c("patient","control"), fill = c("red","blue"), inset = 0.01)

# Barplot series
palette <- RColorBrewer::brewer.pal(length(unique(se$series)),name = 'Set3')
head(palette)
se$color <- palette[as.factor(se$series)]
barplot(dge$samples$lib.size[ord]/1e+06, las = 1, ylab = "Million of counts", xlab = "Samples", 
        col = se$color[ord], border = NA)
legend("topleft", legend=se$series[!duplicated(se$series)], fill = se$color[!duplicated(se$color)], inset = 0.01)

# Barplot tissue
se$tissue <- as.character(se$tissue)
palette <- RColorBrewer::brewer.pal(length(unique(se$tissue)),name = 'Set3')
head(palette)
se$color <- palette[as.factor(se$tissue)]
barplot(dge$samples$lib.size[ord]/1e+06, las = 1, ylab = "Million of counts", xlab = "Samples", 
        col = se$color[ord], border = NA)
legend("topleft", legend=se$tissue[!duplicated(se$tissue)], fill = se$color[!duplicated(se$color)], inset = 0.01)

# QQplot
sampledepth <- round(dge$sample$lib.size / 1e6, digits=1)

qqnorm(sampledepth)
qqline(sampledepth)

# Distribution of expression levels among genes
avgexp <- rowMeans(assays(se)$logCPM)
hist(avgexp, xlab="log2 CPM", main="Distribution of the genes average expression", las=1)
abline(v=1, col="red", lwd=2)

# Filtering lowly expressed genes
sample_cutoff <- 0.2
nsamples_cutoff <- sample_cutoff*nsamples
nsamples_cutoff

logcpm_cutoff <- 1
mask <- rowSums(assays(se)$logCPM <= logcpm_cutoff) >= nsamples_cutoff

dim(se)
se.filt <- se[!mask, ]
dim(se.filt)
dge.filt <- dge[!mask, ]
dim(dge.filt)
kept_genes2b <- dim(se.filt)[1]

# Between sample normalization
dge.filt.norm <- calcNormFactors(dge.filt)
assays(se.filt)$logCPM.norm <- cpm(dge.filt.norm, log = TRUE, prior.count = 3, normalized.lib.sizes = TRUE)

saveRDS(se.filt, file.path(paste(parent_res_dir,"normalized_counts",sep=""), paste(disease_name,"_se_normalized.rds",sep="")))
saveRDS(dge.filt.norm, file.path(paste(parent_res_dir,"normalized_counts",sep=""), paste(disease_name,"_dge_normalized.rds",sep="")))

# Plots: no normalized, filtered, filtered and normalized
par(mfrow=c(3, 2))
multidensity(as.list(as.data.frame(assays(se[, se$type == "patient"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)
multidensity(as.list(as.data.frame(assays(se[, se$type == "control"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Control samples", las=1)

multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$type == "patient"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)
multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$type == "control"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Control samples", las=1)

multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$type == "patient"])$logCPM.norm)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)
multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$type == "control"])$logCPM.norm)),
             xlab="log 2 CPM", legend=NULL, main="Control samples", las=1)

# PCA based on logCPM
pca_logcpm <- prcomp(as.data.frame(t(assays(se.filt)$logCPM)) , center=TRUE, scale=TRUE)
ggbiplot(pca_logcpm, var.axes=FALSE, groups=se.filt$type ) + ggtitle("PCA based on logCPM and colored by type")

# PCA based on normalized counts
pca_norm <- prcomp(as.data.frame(t(assays(se.filt)$logCPM.norm)) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm, var.axes=FALSE, groups=se.filt$type ) + ggtitle("PCA based on normalized data and colored by type")

# tsne based on logCPM
max_perplex <- trunc(((nsamples - 1)/3)-1)
optimal_perplex = round(nsamples^(1/2))
#max_perplex = nsamples/3
perplex1 = round(optimal_perplex/2)
perplex3 = mean(c(optimal_perplex,max_perplex))

se.filt$simplified_type <- se.filt$type
se.filt$simplified_type <- gsub("control","C",se.filt$simplified_type)
se.filt$simplified_type <- gsub("atient","",se.filt$simplified_type)

# tsnes on the logCPM counts
if(skip_unnecessary == FALSE){
  if(perplex1 <= max_perplex){
    print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$type, perplex=perplex1, seed=TRUE)  + ggtitle("tsne based on logCPM perplex1"))
    print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$series, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex1"))
    print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$tissue, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex1"))

    if(optimal_perplex <= max_perplex){
      print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$type, perplex=perplex1, seed=TRUE)  + ggtitle("tsne based on logCPM perplex2"))
      print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$series, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex2"))
      print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$tissue, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex2"))

      if(perplex3 <= max_perplex){
        print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$type, perplex=perplex3, seed=TRUE)  + ggtitle("tsne based on logCPM perplex3"))
        print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$series, perplex=perplex3, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
        print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$tissue, perplex=perplex3, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
      }
    }
  }

  if(max_perplex %in% c(perplex1,optimal_perplex,perplex3)){
    print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$type, perplex=max_perplex, seed=TRUE)  + ggtitle("tsne based on logCPM max_perplex"))
    print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$series, perplex=max_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM max_perplex"))
    print(tsne(as.data.frame(assays(se.filt)$logCPM), labels=se.filt$tissue, perplex=max_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM max_perplex"))
  }

}


# tsne based on normalized counts
if(perplex1 <= max_perplex){
  print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$type, perplex=perplex1, seed=TRUE) + ggtitle("tsne based on normalized counts perplex1"))
  print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$series, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on normalized counts perplex1"))
  print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$tissue, perplex=perplex1, seed=TRUE) + ggtitle("tsne based on normalized counts perplex1"))
  
  if(optimal_perplex <= max_perplex){
    print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$type, perplex=optimal_perplex, seed=TRUE) + ggtitle("tsne based on normalized counts perplex2"))
    print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$series, perplex=optimal_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on normalized counts perplex2"))
    print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$tissue, perplex=optimal_perplex, seed=TRUE) + ggtitle("tsne based on normalized counts perplex2"))
    
    if(perplex3 <= max_perplex){
      print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$type, perplex=perplex3, seed=TRUE) + ggtitle("tsne based on normalized counts perplex3"))
      print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$series, perplex=perplex3, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on normalized counts perplex3"))
      print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$tissue, perplex=perplex3, seed=TRUE) + ggtitle("tsne based on normalized counts perplex3"))
    }
  }
}

if(max_perplex %in% c(perplex1,optimal_perplex,perplex3)){
  print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$type, perplex=max_perplex, seed=TRUE)  + ggtitle("tsne based on normalized counts max_perplex"))
  print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$series, perplex=max_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
  print(tsne(as.data.frame(assays(se.filt)$logCPM.norm), labels=se.filt$tissue, perplex=max_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
}

# spectrum clustering based on logCPM


# spectrum clustering based on normalized counts

# Hierarchichal clustering based on filtered and normalized counts using spearman correlation

## Hierarchichal clustering

plot_hierarchichal_clustering <- function(data,info, save_clusters=FALSE, key_word="unkown", n_boot=100){
  # checkiing
  # data <- logCPM
  # info <- 'based on normalized counts'
  # n_boot <- 50
  
  # pvclust
  pvclust_result <- pvclust(as.matrix(data), method.hclust = "average", method.dist = "correlation", nboot=n_boot, parallel=TRUE, iseed = 1)
  
  # Basic
  par(mfrow=c(1, 1))
  plot(pvclust_result)
  clusters <- pvpick(pvclust_result) # Returns the clusters and GSMs in it.
  pvrect(pvclust_result)
  
  # Saving the obtained clusters
  # if(save_clusters == TRUE){
  #   saveRDS(clusters, file.path(results_dir, paste(disease_name,"_clusters_",key_word,".rds",sep="")))
  #   # write(clusters,file=paste(results_dir,"/",disease_name,"_clusters_",key_word,".rds",sep=""),append=TRUE)
  # }
  # 
  # Colouring by type and series
  
  # Colouring by type
  par(mfrow=c(1, 1))
  colors=c("orange","blue")
  colors = colors[se.filt$type]
  dend <- as.dendrogram(pvclust_result)
  colors = colors[order.dendrogram(dend)]
  labels_colors(dend)<-rep(0,nsamples)
  dend %>% set("leaves_pch", 19) %>%  # node point type
    set("leaves_cex", 0.6) %>%  # node point size
    set("leaves_col", colors) %>%
    # set("leaves(dend)", colors[order.dendrogram(dend)]) %>% ### NEW
    # set("leaves_col", colors[as.integer(factor(se.filt$type))]) %>% # node point color
    plot(main = paste("Hierarchichal clustering",info,sep=" "))
  legend("topright",legend = unique(se.filt$type[order.dendrogram(dend)]), fill = unique(colors), cex=0.6, box.lty=0)
  pvrect(pvclust_result)
  
  # Colouring by series
  par(mfrow=c(1, 1))
  colors = brewer.pal(length(unique(metadata$series)),"Set3")
  colors = colors[as.factor(se.filt$series)]
  colors = colors[order.dendrogram(dend)]
  labels_colors(dend)<-rep(0,nsamples) # Changing the color of the text labels to white
  dend %>% set("leaves_pch", 19) %>%  # node point type
    set("leaves_cex", 0.6) %>%  # node point size
    # set("leaves_col", colors[as.factor(se.filt$series)]) %>%
    set("leaves_col", colors) %>%
    # set("leaves_col", colors[as.integer(factor(se.filt$type))]) %>% # node point color
    plot(main = paste("Hierarchichal clustering",info,sep=" ")) 
  legend("topright",legend = unique(se.filt$series[order.dendrogram(dend)]), fill = unique(colors), cex=0.6, box.lty=0)
  pvrect(pvclust_result)
  
  return(clusters)
}

logCPM <- cpm(dge.filt.norm, log=TRUE, prior.count=3)
if(do_hierarchichal == TRUE){
  normalized_clusters <- plot_hierarchichal_clustering(logCPM,'based on normalized counts', n_boot=1000)
}


# Multidimensional plot

# A) Colour:type and text:tissue
par(mfrow=c(1, 1))
logCPM <- cpm(dge.filt.norm, log=TRUE, prior.count=3)
batch <- as.integer(factor(se.filt$type))
names(batch) <- colnames(se.filt)
outcome <- substr(colnames(se.filt), 9, 12)
names(outcome) <- colnames(se.filt)

plotMDS(dge.filt.norm, labels=se.filt$tissue, col=batch)
legend("topleft",  paste(levels(factor(se.filt$type))), fill=sort(unique(batch)), inset = 0.05, cex = 0.7)


###################### DEAnalysis ############################################
# Without correcting
se.filt$type <- relevel(se.filt$type, ref="control")
mod <- model.matrix(~ se.filt$type, colData(se.filt))
mod0 <- model.matrix(~ 1, colData(se.filt))
pv <- f.pvalue(assays(se.filt)$logCPM.norm, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.5)
# 7882 DE genes

# Distribution of expected p-values
par(mfrow=c(1, 1))
hist(pv, main="Distribution of expected p-values", las=1, ylim = c(0, 8000))

# Correcting for the series
if (length(unique(se$series)) > 1){
  mod <- model.matrix(~type + series, colData(se.filt))
  mod0 <- model.matrix(~series, colData(se.filt))
  colnames(mod)
  pValues <- f.pvalue(logCPM,mod,mod0)
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

write.table(tt,file=paste(results_dir,"/",disease_name,"_DEGs.txt",sep=""),sep="\t",row.names=F)
for_jon <- tt[,c(1,2)]
for_jon <- for_jon[order(-for_jon$logFC),]
write.table(for_jon,file=paste(parent_res_dir,"/","ranked_gene_lists","/" ,disease_name,"_DEGs.rnk",sep=""),sep="\t",col.names=F,row.names=F, quote=FALSE)

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

# general_info_header <- paste('disease_name','kept_genes','n_samples','Number_DEGs','Number_sDEGs','Percentage_aligned','Library_size','perplex1','optimal_perplex','perplex3',sep="\t")
average_perc_mapped <- mean(metadata$Percentage)
average_libsize <- mean(dge$samples$lib.size)
general_info <- paste(disease_name,kept_genes2b,nsamples,sum(tt$P.Value < 0.05),sum(tt$adj.P.Val < 0.05),average_perc_mapped,average_libsize,perplex1,optimal_perplex,perplex3,sep="\t")
write(general_info,file=paste(parent_res_dir,"dl_general_info.txt",sep=""),append=TRUE)

###################### Batch effect removal and PCA ############################################
# Combat
if (length(unique(se$series)) > 1){
  batch <- se.filt$series
  mod <- model.matrix(~type, colData(se.filt))
  combatexp <- ComBat(assays(se.filt)$logCPM.norm, batch, mod) # matrix with the batch effect corrected
  
  class(combatexp)
  
  dim(combatexp)
  
  saveRDS(combatexp, file.path(paste(parent_res_dir,"after_combat_counts",sep=""), paste(disease_name,"_combat_counts.rds",sep="")))
}


if (length(unique(se$series)) > 1){
  # PCA
  pca_logcpm <- prcomp(as.data.frame(t(combatexp)) , center=TRUE, scale=TRUE)
  ggbiplot(pca_logcpm, var.axes=FALSE, groups=se.filt$type ) + ggtitle("PCA after Combat batch effect removal")
  
  # tsne after Combat batch effect correction
  if(perplex1 <= max_perplex){
    print(tsne(as.data.frame(combatexp), labels=se.filt$type, perplex=perplex1, seed=TRUE) + ggtitle("tsne after Combat batch effect removal perplex1"))
    print(tsne(as.data.frame(combatexp), labels=se.filt$series, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after Combat batch effect removal perplex1"))
    print(tsne(as.data.frame(combatexp), labels=se.filt$tissue, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after Combat batch effect removal perplex1"))
    
    if(optimal_perplex <= max_perplex){
      print(tsne(as.data.frame(combatexp), labels=se.filt$type, perplex=optimal_perplex, seed=TRUE) + ggtitle("tsne after Combat batch effect removal perplex2"))
      print(tsne(as.data.frame(combatexp), labels=se.filt$series, perplex=optimal_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after Combat batch effect removal perplex2"))
      print(tsne(as.data.frame(combatexp), labels=se.filt$tissue, perplex=optimal_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after Combat batch effect removal perplex2"))
      
      if(perplex3 <= max_perplex){
        print(tsne(as.data.frame(combatexp), labels=se.filt$type, perplex=perplex3, seed=TRUE) + ggtitle("tsne after Combat batch effect removal perplex3"))
        print(tsne(as.data.frame(combatexp), labels=se.filt$series, perplex=perplex3, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after Combat batch effect removal perplex3"))
        print(tsne(as.data.frame(combatexp), labels=se.filt$tissue, perplex=perplex3, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after Combat batch effect removal perplex3"))
        
      }
    }
  }
  
  if(max_perplex %in% c(perplex1,optimal_perplex,perplex3)){
    print(tsne(as.data.frame(combatexp), labels=se.filt$type, perplex=max_perplex, seed=TRUE)  + ggtitle("tsne based on logCPM perplex3"))
    print(tsne(as.data.frame(combatexp), labels=se.filt$series, perplex=max_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
    print(tsne(as.data.frame(combatexp), labels=se.filt$tissue, perplex=max_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
  }
  
  if(do_hierarchichal == TRUE){
  combat_clusters <- plot_hierarchichal_clustering(combatexp, 'after Combat batch effect removal',save_clusters=TRUE, key_word="combat", n_boot=1000)
  }
}

# QR Decomposition
if (length(unique(se$series)) > 1){
  batch <- se.filt$series
  qrexp <- removeBatchEffect(assays(se.filt)$logCPM.norm, batch, design = mod)
  class(qrexp)
  
  saveRDS(qrexp, file.path(paste(parent_res_dir,"after_qrdecomp_counts",sep=""), paste(disease_name,"_qrdecomp_counts.rds",sep="")))
  
}


if (length(unique(se$series)) > 1){
  # PCA
  pca_logcpm <- prcomp(as.data.frame(t(qrexp)) , center=TRUE, scale=TRUE)
  ggbiplot(pca_logcpm, var.axes=FALSE, groups=se.filt$type ) + ggtitle("PCA after QR Decomposition batch effect removal")
  
  # tsne after QR Decomposition batch effect correction
  if(perplex1 <= max_perplex){
    print(tsne(as.data.frame(qrexp), labels=se.filt$type, perplex=perplex1, seed=TRUE) + ggtitle("tsne after QR Decomposition batch effect removal perplex1"))
    print(tsne(as.data.frame(qrexp), labels=se.filt$series, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after QR Decomposition batch effect removal perplex1"))
    print(tsne(as.data.frame(qrexp), labels=se.filt$tissue, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after QR Decomposition batch effect removal perplex1"))
    
    if(optimal_perplex <= max_perplex){
      print(tsne(as.data.frame(qrexp), labels=se.filt$type, perplex=optimal_perplex, seed=TRUE) + ggtitle("tsne after QR Decomposition batch effect removal perplex2"))
      print(tsne(as.data.frame(qrexp), labels=se.filt$series, perplex=optimal_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after QR Decomposition batch effect removal perplex2"))
      print(tsne(as.data.frame(qrexp), labels=se.filt$tissue, perplex=optimal_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after QR Decomposition batch effect removal perplex2"))
      
      if(perplex3 <= max_perplex){
        print(tsne(as.data.frame(qrexp), labels=se.filt$type, perplex=perplex3, seed=TRUE) + ggtitle("tsne after QR Decomposition batch effect removal perplex3"))
        print(tsne(as.data.frame(qrexp), labels=se.filt$series, perplex=perplex3, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after QR Decomposition batch effect removal perplex3"))
        print(tsne(as.data.frame(qrexp), labels=se.filt$tissue, perplex=perplex3, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after QR Decomposition batch effect removal perplex3"))
        
      }
    }
  }
  
  if(max_perplex %in% c(perplex1,optimal_perplex,perplex3)){
    print(tsne(as.data.frame(qrexp), labels=se.filt$type, perplex=max_perplex, seed=TRUE)  + ggtitle("tsne based on logCPM perplex3"))
    print(tsne(as.data.frame(qrexp), labels=se.filt$series, perplex=max_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
    print(tsne(as.data.frame(qrexp), labels=se.filt$tissue, perplex=max_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
  }
  
  if(do_hierarchichal == TRUE){
    qrdecomposition_clusters <- plot_hierarchichal_clustering(qrexp,'after QR Decomposition batch effect removal',save_clusters=TRUE, key_word="qr_decomposition",n_boot=1000)
    # qrdecomposition_clusters <- plot_hierarchichal_clustering(qrexp,'after QR Decomposition batch effect removal')
    print(qrdecomposition_clusters)
  }
}


# SVD - it removes unkown batch effect, it behaves the worst
s <- fast.svd(t(scale(t(assays(se.filt)$logCPM.norm), center = TRUE, scale = TRUE)))
pcSds <- s$d
pcSds[1] <- 0
svdexp <- s$u %*% diag(pcSds) %*% t(s$v)
colnames(svdexp) <- colnames(se.filt)
class(svdexp)
dim(svdexp)

d <- as.dist(1 - cor(svdexp, method = "spearman"))
sampleClustering <- hclust(d)

# tsne after SVD batch effect correction
if(perplex1 <= max_perplex){
  print(tsne(as.data.frame(svdexp), labels=se.filt$type, perplex=perplex1, seed=TRUE) + ggtitle("tsne after SVD batch effect removal perplex1"))
  print(tsne(as.data.frame(svdexp), labels=se.filt$series, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after SVD batch effect removal perplex1"))
  print(tsne(as.data.frame(svdexp), labels=se.filt$tissue, perplex=perplex1, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after SVD batch effect removal perplex1"))
  
  if(optimal_perplex <= max_perplex){
    print(tsne(as.data.frame(svdexp), labels=se.filt$type, perplex=optimal_perplex, seed=TRUE) + ggtitle("tsne after SVD batch effect removal perplex2"))
    print(tsne(as.data.frame(svdexp), labels=se.filt$series, perplex=optimal_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after SVD batch effect removal perplex2"))
    print(tsne(as.data.frame(svdexp), labels=se.filt$tissue, perplex=optimal_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after SVD batch effect removal perplex2"))
    
    if(perplex3 <= max_perplex){
      print(tsne(as.data.frame(svdexp), labels=se.filt$type, perplex=perplex3, seed=TRUE) + ggtitle("tsne after SVD batch effect removal perplex3"))
      print(tsne(as.data.frame(svdexp), labels=se.filt$series, perplex=perplex3, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after SVD batch effect removal perplex3"))
      print(tsne(as.data.frame(svdexp), labels=se.filt$tissue, perplex=perplex3, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne after SVD batch effect removal perplex3"))
      
    }
  }
}

if(max_perplex %in% c(perplex1,optimal_perplex,perplex3)){
  print(tsne(as.data.frame(svdexp), labels=se.filt$type, perplex=max_perplex, seed=TRUE)  + ggtitle("tsne based on logCPM perplex3"))
  print(tsne(as.data.frame(svdexp), labels=se.filt$series, perplex=max_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
  print(tsne(as.data.frame(svdexp), labels=se.filt$tissue, perplex=max_perplex, seed=TRUE, text=as.character(se.filt$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
}

if(do_hierarchichal == TRUE){
  svd_clusters <- plot_hierarchichal_clustering(svdexp,  'after SVD batch effect removal',n_boot=1000)
}

write(disease_name,file=paste(parent_res_dir,"is_finished.txt",sep=""),append=TRUE)


dev.off()
