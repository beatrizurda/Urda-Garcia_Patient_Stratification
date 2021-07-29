#####################################################################################
# PVCLUST TO UNRAVEL MOLECULAR MECHANISMS UNDER DISEASE (DIS)SIMILARITIES
#
#  NEW VERSION in ANALYSIS.
# 
#  Beatriz Urda García 2020
######################################################################################

# Importing libraries
library(ggplot2)
library(pvclust)
library(factoextra)
library(Hmisc) # function: capitalize

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

setwd("~/Desktop/ANALYSIS")

### ARGUMENTS
## Generate a table with the up- and down- regulated pathways
FEdir <- "FE_results/Pathways_DEGs_diseases"
# FEdir <- "FE_results/Pathways_DMs_diseases"

if(FEdir == "FE_results/Pathways_DEGs_diseases"){
  inputname = paste("logFC_8")
}else if(FEdir == "FE_results/Pathways_DMs_diseases"){
  inputname = paste("DM_8")
}

dist_method <- "spearman"
dist_method <- "euclidean"

### START of the script

# OLD COLORS
# dismeta <- read.csv2("dl_general_info_umls_ncolors_fixed.csv",stringsAsFactors = F,sep="\t",header=T)
# ch1 <- which(dismeta$tis_cat_colors == 'fuchsia')
# dismeta$tis_cat_colors[ch1] <- 'magenta3'
# ch1 <- which(dismeta$dis_cat_colors == 'grey')
# dismeta$dis_cat_colors[ch1] <- 'black'

# NEW COLORS
dismeta <- read.csv2("new_disease_metadata_final_names.txt",stringsAsFactors = F,sep="\t",header=T)
initial_dismeta <- dismeta
dismeta$dis_cat_colors <- dismeta$new_dis_cat_colors


# Adding number of patients and controls for each disease
newdismeta <- read.csv2("Disease_gene_variability/summary_disease_variability.csv",stringsAsFactors = F,sep="\t",header=T)
newdismeta <- newdismeta[, c(1,16,17)]
dismeta <- merge(dismeta, newdismeta, by.x="disease_name", by.y = "disease")
dismeta <- dismeta[,c("disease_name","n_samples","dis_cat_colors", "tis_cat_colors", "final_disease_name")]
# dismeta <- dismeta[,c(1,3,19,20)]

list_paths<-list()
pats2<-c()
wo_pathways <- c()

a <- 1
for(disdir in list.files(FEdir)){
  # disdir <- list.files(FEdir)[1]
  cdis <- gsub("\\..+","",disdir) # c disease name
  cdis <- gsub("_DEGs","",cdis) # c disease name    # Before we had _logFC
  cdis <- gsub("_delta_dm_values","",cdis) # c disease name    # Before we had _logFC
  cdir <- paste(FEdir,"/",disdir,"/",sep="")
  listf <- list.files(cdir)
  # Getting filenames for the FE results
  pos <- listf[intersect(grep("gsea_report_for_na_pos",listf),grep("xls",listf))]
  neg <- listf[intersect(grep("gsea_report_for_na_neg",listf),grep("xls",listf))]
  # Reading FE results
  posi <- read.csv2(paste(cdir,pos,sep=""),stringsAsFactors = F,sep="\t")
  negi <- read.csv2(paste(cdir,neg,sep=""),stringsAsFactors = F,sep="\t")
  posu<-posi[,c(1,6,8)]
  negu<-negi[,c(1,6,8)]
  posu <- posu[grep("REACTOME_", posu$NAME), ]
  negu <- negu[grep("REACTOME_", negu$NAME), ]
  # Keeping the names of significant pathways
  respos<-posu[which(as.numeric(as.character(posu$FDR.q.val))<=0.05),1]
  resneg<-negu[which(as.numeric(as.character(negu$FDR.q.val))<=0.05),1]
  list_paths[[cdis]]$pos<-respos
  list_paths[[cdis]]$neg<-resneg
  pats2 <- c(pats2,respos, resneg)
  if(length(respos)==0 && length(resneg)==0){
    wo_pathways<-c(wo_pathways,a)
  }
  a <- a + 1
}

pats<-unique(pats2)
length(pats) # 4173 con todo --> 800 solo con REACTOME (únicas) | 459 de REACTOME con DM
patst2<-sort(table(pats2),decreasing = T)

# Sort pats by parents
pcats<-read.csv2("Reactome/Reactome_parents.txt",stringsAsFactors = F,sep="\t",header=T)

# Add the ones that are missing
pcats <- rbind(pcats, c("REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PEROXISOME_PROLIFERATOR_ACTIVATED_RECEPTOR_ALPHA_PPARALPHA", "Metabolism"))
pcats <- rbind(pcats, c("REACTOME_MITOTIC_G1_G1_S_PHASES", "Cell Cycle"))
pcats <- rbind(pcats, c("REACTOME_MICROAUTOPHAGY", "Autophagy"))

pcats <- pcats[order(pcats$Parent),]
head(pcats)
dim(pcats) # 1501 2 --> after adding the parents: 1527

# Select from pcats the ones that are in pats
sorted_paths <- pcats[pcats$GSEA_pathway %in% pats,]
dim(sorted_paths) # 823

# Compare both
length(unique(pats))
length(unique(sorted_paths$GSEA_pathway))
setdiff(unique(pats), unique(sorted_paths$GSEA_pathway))

# Coloring by category
# OPTION 1
length(unique(sorted_paths$Parent))
ncols <- 27
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(ncols)
# OPTION 2 -- THE BEST I THINK
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(1,ncols), col=sample(color, ncols))
mycolors=sample(color, ncols)
pie(rep(1,ncols), col=mycolors)
pie(rep(1,ncols), labels=unique(as.character(sorted_paths$Parent)), col=mycolors)
# saveRDS(mycolors,file="mycolors1.rds")
mycolors <- readRDS("mycolors1.rds")
mycolors[16] <- "#cfc5e8"; mycolors[17] <- "#dfe66c"; mycolors[19] <- "#66d4cd";
mycolors[25] <- "#648f6e"; mycolors[27] <- "#daeddd"; 
if(length(unique(sorted_paths$Parent)) != 27){
  parent_colors <- readRDS("Reactome/parent_colors.rds")
  subset <- parent_colors[which(parent_colors$Parent  %in% sorted_paths$Parent),]
  mycolors <- subset$Color
}
length(mycolors)

# OPTION 3
# n <- ncols
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=sample(col_vector, n))
# mycolors=sample(col_vector, n)

cols <- data.frame(Parent = unique(sorted_paths$Parent), Color = mycolors)
with_cols <- merge(sorted_paths, cols, by='Parent', all.x = TRUE)

# Saving parent colors dataframe
parent_colors <- with_cols[,c(1,3)]
parent_colors <- parent_colors[!duplicated(parent_colors),]
parent_colors$Color <- as.character(parent_colors$Color)
# saveRDS(parent_colors, file="Reactome/parent_colors.rds")
# write.table(parent_colors,file="FE_results/parent_colors_DM.txt",sep='\t',row.names=FALSE, quote = FALSE)

# Saving ICD9 category colors dataframe
icd9_colors <- initial_dismeta[,c("disease_cat","new_dis_cat_colors")]
icd9_colors <- icd9_colors[!duplicated(icd9_colors),]

# Plotting the legend
pdf(file=paste("FE_results/Heatmaps/Legend_Pathway_category_",inputname,".pdf", sep=""))
legendcolor<-rep(0,length(parent_colors[,1]))
names(legendcolor)<-parent_colors[,1]
barplot(legendcolor,ylim=c(0,10),las=2,col=parent_colors[,2],
        legend=T,args.legend = list("cex"=0.5,"border"=NA))

legendicd9 <- rep(0,length(icd9_colors[,1]))
names(legendicd9) <- capitalize(icd9_colors[,1])
barplot(legendicd9,ylim=c(0,10),las=2,col=icd9_colors[,2],
        legend=T,args.legend = list("cex"=0.5,"border"=NA))
dev.off()

unik <- !duplicated(as.character(with_cols$Parent))
uniq_indexes <- which(unik == TRUE)
# hey <- I(nlevels(as.factor(sorted_paths$Parent)), mycolors)
# hey <- I(brewer.pal(nlevels(as.factor(sorted_paths$Parent)), name = 'Set2')(nb.cols))
# 
# levs <- levels(as.factor(sorted_paths$Parent))
# cols <- mycolors

# Using categories
matriz<-matrix(nrow=length(sorted_paths$GSEA_pathway),ncol=length(list_paths),0) ; rownames(matriz)<-sorted_paths$GSEA_pathway ; colnames(matriz)<-names(list_paths)
# Using pathways
# matriz<-matrix(nrow=length(pats),ncol=length(list_paths),0) ; rownames(matriz)<-pats ; colnames(matriz)<-names(list_paths)
for(a in 1:length(names(list_paths))){
  # a<-1
  if(length(list_paths[[a]]$pos)>0){
    matriz[list_paths[[a]]$pos,a]<-1
  }
  
  if(length(list_paths[[a]]$neg)>0){
    matriz[list_paths[[a]]$neg,a]<-(-1)
  }
}
matriz[1:5,1:5]

# Changing names to the new ones:
colnames(matriz) == dismeta$disease_name
colnames(matriz) <- dismeta$final_disease_name

tmatriz <- t(matriz)
tmatriz[1:5,1:5]

cenf<-matriz
palette.breaks=c(seq(-1,-0.01,length=10),seq(-0.001,0.001,length=10),seq(0.01,1,length=10))
color.palette<-colorRampPalette(c("#C81E17","#FFFFFF","#405191"))(length(palette.breaks)-1) # OLD COLORS: "#FF0000","#FFFFFF","#1022ea"
# "#C81E17","#FFFFFF","#3D4D8B"
cenf2<-cenf
enf<-cenf2

## PLOT SUMMARY HEATMAPS
parents <- c("REACTOME_AUTOPHAGY","REACTOME_CELL_CYCLE","REACTOME_CELL_CELL_COMMUNICATION","REACTOME_CELLULAR_RESPONSES_TO_EXTERNAL_STIMULI",
             "REACTOME_CHROMATIN_ORGANIZATION","REACTOME_CIRCADIAN_CLOCK","","REACTOME_DIGESTION_AND_ABSORPTION",
             "REACTOME_DNA_REPAIR","REACTOME_DNA_REPAIR","REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION","REACTOME_HEMOSTASIS",
             "REACTOME_NEURONAL_SYSTEM","REACTOME_PROGRAMMED_CELL_DEATH", "REACTOME_PROTEIN_LOCALIZATION","REACTOME_REPRODUCTION",
             "REACTOME_TRANSPORT_OF_SMALL_MOLECULES",
             "REACTOME_DEVELOPMENTAL_BIOLOGY","REACTOME_DISEASE","REACTOME_GENE_EXPRESSION_TRANSCRIPTION",
             "REACTOME_METABOLISM_OF_RNA","REACTOME_MUSCLE_CONTRACTION","REACTOME_ORGANELLE_BIOGENESIS_AND_MAINTENANCE",
             "REACTOME_VESICLE_MEDIATED_TRANSPORT")

####### IMPORTANT: Parents que no están en pats!!!: Immune System; Metabolism; Metabolism of proteins;

parent_index <- which(rownames(matriz) %in% parents)
sum_matrix <- matriz[parent_index, ]
parent_col <- with_cols[parent_index, ]

# Changing disease names to final names
# colnames(sum_matrix) == dismeta$disease_name
# colnames(sum_matrix) <- dismeta$final_disease_name

dev.off()
pdf(file=paste("FE_results/Heatmaps/Summary_Heatmap_wardD_pathways_",inputname,".pdf", sep=""))
h1 <- heatmap.2(sum_matrix,scale="none",
          hclustfun = function(data=sum_matrix){
            return(hclust(dist(data), method="ward.D2")) # default: euclidean
          },
          col=color.palette,breaks=palette.breaks,dendrogram = "col",trace="none",
          Rowv=FALSE,cexRow = 0.6,cexCol = 0.4,"key"=FALSE,margins = c(15,18),
          RowSideColors = as.vector(parent_col$Color),
          colCol = dismeta$dis_cat_colors,
          xlab="Ward2 clustering | euclidean distance | all diseases",
          ylab="Reactome pathways")
h1



# Adding the histogram
cmat <- sum_matrix
cdisorder <- cmat[rev(h1$rowInd), h1$colInd]
cdisorder <- colnames(cdisorder)
cdismeta <- dismeta[dismeta$final_disease_name %in% cdisorder, ]
cdismeta <- cdismeta[match(cdisorder, cdismeta$final_disease_name), ] # Sort the dataframe by external vector
cdismeta$order <- c(1:length(cdismeta$final_disease_name))

# Plotting the distribution of controls and patients in sample size
library(reshape2)
melted <- melt(cdismeta[,c(1,3,4,5)], id.vars= c("disease_name","order"))
shist <- ggplot(melted, aes(reorder(disease_name, order),value, fill=variable)) + 
  # geom_bar(stat="identity", position="dodge") +
  geom_bar(stat="identity") +
  theme(text = element_text(size=8), axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("") + scale_fill_discrete(name = "Sample size", labels = c("Patients", "Controls"))
shist

# Remove diseases with 0 enriched parents
is_empty <- apply(sum_matrix,2, function(x) all(x == 0))
empty_index <- which(is_empty == TRUE)
sum_matrix_wo <- sum_matrix[, -empty_index]

h1 <- heatmap.2(sum_matrix_wo,scale="none",
          hclustfun = function(data=sum_matrix_wo){
            return(hclust(dist(data), method="ward.D2")) # default: euclidean
          },
          col=color.palette,breaks=palette.breaks,dendrogram = "col",trace="none",
          Rowv=FALSE,cexRow = 0.6,cexCol = 0.4,"key"=FALSE,margins = c(15,18),
          RowSideColors = as.vector(parent_col$Color),
          # colCol = dismeta$dis_cat_colors,
          xlab="Ward2 clustering | euclidean distance | wo diseases with 0 pathways",
          ylab="Reactome pathways")
h1
# Adding the histogram
cmat <- sum_matrix_wo
cdisorder <- cmat[rev(h1$rowInd), h1$colInd]
cdisorder <- colnames(cdisorder)
cdismeta <- dismeta[dismeta$final_disease_name %in% cdisorder, ]
cdismeta <- cdismeta[match(cdisorder, cdismeta$final_disease_name), ] # Sort the dataframe by external vector
cdismeta$order <- c(1:length(cdismeta$final_disease_name))

# Plotting only sample size
# shist <- ggplot(cdismeta, aes(reorder(disease_name, order),n_samples)) + geom_bar(stat="identity") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #+ coord_flip()
# shist

# Plotting the distribution of controls and patients in sample size
library(reshape2)
melted <- melt(cdismeta[,c(1,3,4,5)], id.vars= c("disease_name","order"))
shist <- ggplot(melted, aes(reorder(disease_name, order),value, fill=variable)) + 
  # geom_bar(stat="identity", position="dodge") +
  geom_bar(stat="identity") +
  theme(text = element_text(size=12), axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("") + scale_fill_discrete(name = "Sample size", labels = c("Patients", "Controls"))
shist

# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")
# library(ggpubr)
# ggarrange(h1, shist, heights = c(2, 0.7),
#           ncol = 1, nrow = 2, align = "v")
# 
# library(gridExtra)
# grid.arrange(grobs=list(shist,shist), nrow=2, top="Title")
# grid.arrange(grobs=list(h1,shist), nrow=2, top="Title")

heatmap.2(t(sum_matrix_wo),scale="none",
          hclustfun = function(data=t(sum_matrix_wo)){
            return(hclust(dist(data), method="ward.D2")) # default: euclidean
          },
          col=color.palette,breaks=palette.breaks,dendrogram = "col",trace="none",
          Rowv=FALSE,cexRow = 0.6,cexCol = 0.6,"key"=FALSE,margins = c(18,18),
          ColSideColors = as.vector(parent_col$Color),
          # colCol = dismeta$dis_cat_colors,
          xlab="Ward2 clustering | euclidean distance | wo diseases with 0 pathways")


dev.off()

## PLOT HEATMAPS ##
dev.off()
pdf(file=paste("FE_results/Heatmaps/Heatmap_wardD_pathways_",inputname,".pdf", sep=""))
# pdf(file="PLOTS/Heatmaps/Based_on_pathways/Heatmap_wardD_pathways_logFC_2.pdf")
# pdf(file="PLOTS/Heatmaps/Based_on_pathways/Heatmap_wardD_pathways_DM_1.pdf")

# Changing disease names to final names
# colnames(cenf2) == dismeta$disease_name
# colnames(cenf2) <- dismeta$final_disease_name

h1 <- heatmap.2(cenf2,scale="none",
          hclustfun = function(data=cenf2){
            return(hclust(dist(data), method="ward.D2")) # default: euclidean
          },
          col=color.palette,breaks=palette.breaks,dendrogram = "col",trace="none",
          Rowv=FALSE,cexRow = 0.0001,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          RowSideColors = as.vector(with_cols$Color),
          colCol = dismeta$dis_cat_colors,
          sepcolor = 'black', rowsep = uniq_indexes, sepwidth = c(0.001,0.001),
          xlab="Ward2 clustering | euclidean distance | all diseases",
          ylab="Reactome pathways")

h1

# pie(rep(1,ncols), col=mycolors, clockwise = TRUE)
# pie(rep(1,ncols), labels=unique(as.character(sorted_paths$Parent)), col=mycolors, clockwise = TRUE)

# barplot(rep(1,n), names.arg=unique(as.character(sorted_paths$Parent)), col=mycolors, legend=TRUE)

# Adding the histogram
cmat <- cenf2
cdisorder <- cmat[rev(h1$rowInd), h1$colInd]
cdisorder <- colnames(cdisorder)
cdismeta <- dismeta[dismeta$final_disease_name %in% cdisorder, ]
cdismeta <- cdismeta[match(cdisorder, cdismeta$final_disease_name), ] # Sort the dataframe by external vector
cdismeta$order <- c(1:length(cdismeta$final_disease_name))

# Plotting the distribution of controls and patients in sample size
library(reshape2)
melted <- melt(cdismeta[,c(1,3,4,5)], id.vars= c("disease_name","order"))
shist <- ggplot(melted, aes(reorder(disease_name, order),value, fill=variable)) + 
  # geom_bar(stat="identity", position="dodge") +
  geom_bar(stat="identity") +
  theme(text = element_text(size=8), axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("") + scale_fill_discrete(name = "Sample size", labels = c("Patients", "Controls"))
shist

# Remove disease with 0 dysregulated pathways
colnames(cenf2)[wo_pathways]
# for logFC:  "HIV"
# for DM: 
# [1] "AlagilleSyndrome"                      "BorreliaBurdorferiInfection"          
# [3] "Cardiomyopathy"                        "CoeliacDisease"                       
# [5] "FamilialPulmonaryArterialHypertension" "HIVART"                               
# [7] "HyperplasticPolyps"                    "Ischemia"                             
# [9] "MyotonicDystrophy"                     "SessileSerratedPolyp"                 
# [11] "SLE" 
cenf2wo <- cenf2[,-wo_pathways]
cdiscols <- dismeta[dismeta$final_disease_name %in% colnames(cenf2wo), ]
h1 <- heatmap.2(cenf2wo,scale="none",
          hclustfun = function(data=cenf2wo){
            return(hclust(dist(data), method="ward.D2")) # default: euclidean
          },
          col=color.palette,breaks=palette.breaks,dendrogram = "col",trace="none",
          Rowv=FALSE,cexRow = 0.0001,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          RowSideColors = as.vector(with_cols$Color),
          sepcolor = 'black', rowsep = uniq_indexes, sepwidth = c(0.03,0.03),
          colCol = cdiscols$dis_cat_colors,
          xlab="Ward2 clustering | euclidean distance | wo diseases with 0 pathways")
h1

# Adding the histogram
cmat <- cenf2wo
cdisorder <- cmat[rev(h1$rowInd), h1$colInd]
cdisorder <- colnames(cdisorder)
cdismeta <- dismeta[dismeta$final_disease_name %in% cdisorder, ]
cdismeta <- cdismeta[match(cdisorder, cdismeta$final_disease_name), ] # Sort the dataframe by external vector
cdismeta$order <- c(1:length(cdismeta$final_disease_name))

# Plotting the distribution of controls and patients in sample size
library(reshape2)
melted <- melt(cdismeta[,c(1,3,4,5)], id.vars= c("disease_name","order"))
shist <- ggplot(melted, aes(reorder(disease_name, order),value, fill=variable)) + 
  # geom_bar(stat="identity", position="dodge") +
  geom_bar(stat="identity") +
  theme(text = element_text(size=8), axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("") + scale_fill_discrete(name = "Sample size", labels = c("Patients", "Controls"))
shist

# With Spearman
heatmap.2(cenf2wo,scale="none",
          distfun = function(dat=cenf2wo){
            return(get_dist(dat, method='spearman', stand=FALSE))
          },
          hclustfun = function(data=cenf2wo){
            return(hclust(dist(data), method="ward.D2"))
          },
          col=color.palette,breaks=palette.breaks,dendrogram = "col",trace="none",
          Rowv=FALSE,cexRow = 0.0001,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          RowSideColors = as.vector(with_cols$Color),
          sepcolor = 'black', rowsep = uniq_indexes, sepwidth = c(0.03,0.03),
          colCol = cdiscols$dis_cat_colors,
          xlab="Ward2 clustering | 1-spearman | wo diseases with 0 pathways")


### Pathways vs diseases
heatmap.2(tmatriz,scale="none",
          hclustfun = function(data=tmatriz){
            return(hclust(dist(data), method="ward.D2")) # default: euclidean
          },
          col=color.palette,breaks=palette.breaks,dendrogram = "col",trace="none",
          Rowv=FALSE,cexRow = 0.5,cexCol = 0.1,"key"=FALSE,margins = c(10,15),
          colCol = dismeta$dis_cat_colors,
          # RowSideColors = as.vector(with_cols$Color),
          # sepcolor = 'black', rowsep = uniq_indexes, sepwidth = c(0.03,0.03),
          xlab="t Ward2 clustering | euclidean distance | all diseases",
          ylab="Diseases")

rownames(tmatriz)[wo_pathways]
tmatrizwo <- tmatriz[-wo_pathways,]
heatmap.2(tmatrizwo,scale="none",
          hclustfun = function(data=tmatrizwo){
            return(hclust(dist(data), method="ward.D2")) # default: euclidean
          },
          col=color.palette,breaks=palette.breaks,dendrogram = "col",trace="none",
          Rowv=FALSE,cexRow = 0.5,cexCol = 0.1,"key"=FALSE,margins = c(10,15),
          # RowSideColors = as.vector(with_cols$Color),
          # sepcolor = 'black', rowsep = uniq_indexes, sepwidth = c(0.03,0.03),
          xlab="t Ward2 clustering | euclidean distance | wo diseases with 0 pathways",
          ylab="Diseases")

heatmap.2(tmatrizwo,scale="none",
          distfun = function(dat=tmatrizwo){
            return(get_dist(dat, method='spearman', stand=FALSE))
          },
          hclustfun = function(data=tmatrizwo){
            return(hclust(dist(data), method="ward.D2"))
          },
          col=color.palette,breaks=palette.breaks,dendrogram = "col",trace="none",
          Rowv=FALSE,cexRow = 0.0001,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          colCol = dismeta$dis_cat_colors,
          # RowSideColors = as.vector(with_cols$Color),
          # sepcolor = 'black', rowsep = uniq_indexes, sepwidth = c(0.03,0.03),
          xlab="Ward2 clustering | 1-spearman | wo diseases with 0 pathways")

dev.off()


### PLOT HEATMAPS BY PARENT CATEGORY
draw_heatmap <- function(matrix,parent,alldis=FALSE,cexrow=0.4,cexcol=0.5){
  if(alldis == F){
    alldis <- ""
  }else{
    alldis <- "| all diseases"
  }
  h1 <- heatmap.2(matrix,scale="none",
                  hclustfun = function(data=cenf2){
                    return(hclust(dist(data), method="ward.D2")) # default: euclidean
                  },
                  col=color.palette,
                  # breaks=palette.breaks,
                  dendrogram = "col",trace="none",
                  Rowv=FALSE,cexRow = cexrow,cexCol = cexcol, "key"=FALSE, margins = c(9,30),
                  # colCol = dismeta$dis_cat_colors,
                  # RowSideColors = as.vector(with_cols$Color),
                  # sepcolor = 'black', rowsep = uniq_indexes, sepwidth = c(0.03,0.03),
                  # xlab=paste("Ward2 euclidean clustering |  ",parent,alldis,sep=""),
                  ylab=parent) + coord_flip()
  return(h1)
}

remove_dis_with_zero_paths <- function(matrix){
  is_empty <- apply(matrix,2, function(x) all(x == 0))
  empty_index <- which(is_empty == TRUE)
  sum_matrix_wo <- matrix[, -empty_index]
  return(sum_matrix_wo)
}

pdf(file=paste("FE_results/Heatmaps/Parents_Heatmaps_wardD_pathways_",inputname,".pdf", sep=""),width=10)
parent_colors <- readRDS("Reactome/parent_colors.rds")
parent_colors <- parent_colors[parent_colors$Parent %in% unique(sorted_paths$Parent),]
rownames(cenf2) <- gsub("REACTOME_","",rownames(cenf2)) # Removing REACTOME_ del nombre de las pathways
length(uniq_indexes)
only_one <- c()
for(k in 1:length(uniq_indexes)){
  # k <- 1
  starti <- uniq_indexes[k]
  fini <- uniq_indexes[k+1]-1; if(is.na(fini) == T){fini <- dim(cenf2)[1]}
  cpaths <- cenf2[starti:fini,]
  if(class(cpaths) == 'matrix'){
    cpaths <- remove_dis_with_zero_paths(cpaths) # TO REMOVE DIS WITH 0 PATHS -- uncomment to undo
    if(class(cpaths) == 'matrix'){
      # Choose size of column labels
      if(dim(cpaths)[2] < 26){
        cexcol=0.5
      }else if((dim(cpaths)[2] >= 26) & (dim(cpaths)[2] < 40)){
        cexcol=0.4
      }else{
        cexcol=0.3
      }
      # Choose size of row labels
      if(dim(cpaths)[1] <= 30){
        cexrow=0.5
      }else if(dim(cpaths)[1] <= 60){
        cexrow=0.4
      }else{
        cexrow=0.3
      }
      h1 <- draw_heatmap(cpaths, parent=parent_colors$Parent[k], cexrow=cexrow, cexcol=cexcol, alldis=F) 
      h1
    }else{
      only_one <- append(only_one,parent_colors$Parent[k])
    }
  }else{
    only_one <- append(only_one,parent_colors$Parent[k])
  }
}
dev.off()

## DENDOGRAMS, removing those diseases with 0 pathways enriched ##

matriz2<-matriz[,-wo_pathways]

cluster_eucl <-pvclust(matriz2,method.dist="euclidean", method.hclust = "ward.D2", nboot=10000)
pdf(file=paste("FE_results/Heatmaps/Dendogram_wardD_pathways_",inputname,".pdf", sep=""))
plot(cluster_eucl,hang = -1, cex = 0.5)
pvrect(cluster_eucl)
dev.off()

pvclust_dis <- cluster_eucl$hclust$labels[cluster_eucl$hclust$order]
pvclust_order <- cluster_eucl$hclust$order
new_matrix <- tmatrizwo[match(pvclust_dis,rownames(tmatrizwo)), ]

pdf(file=paste("FE_results/Heatmaps/Heatmap_wardD_pathways_pvclust_",inputname,".pdf", sep=""))
# Without reordering
heatmap.2(t(new_matrix),Rowv=FALSE,Colv=FALSE,scale="none",
          col=color.palette,breaks=palette.breaks,dendrogram = "none",trace="none",
          cexRow = 0.0001,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          RowSideColors = as.vector(with_cols$Color),
          sepcolor = 'black', rowsep = uniq_indexes, sepwidth = c(0.03,0.03),
          colCol = cdiscols$dis_cat_colors[cluster_eucl$hclust$order],
          xlab="pvclust order | Ward2 clustering | euclidean distance | wo diseases with 0 pathways")

plot(cluster_eucl,hang = -1, cex = 0.5, print.num = FALSE, print.pv = FALSE)
dev.off()

# Reordering with pvclust -- doesn't work
# heatmap.2(t(new_matrix),Rowv=FALSE,scale="none",
#           hclustfun = function(data=t(new_matrix)){
#             return(pvclust(t(matriz2), method.dist="euclidean", method.hclust = "ward.D2", nboot=10000)) # default: euclidean
#           },
#           col=color.palette,breaks=palette.breaks,dendrogram = "none",trace="none",
#           cexRow = 0.0001,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
#           RowSideColors = as.vector(with_cols$Color),
#           sepcolor = 'black', rowsep = uniq_indexes, sepwidth = c(0.03,0.03),
#           colCol = cdiscols$dis_cat_colors,
#           xlab="pvclust order | Ward2 clustering | euclidean distance | wo diseases with 0 pathways")





# Clustering en Spearman -- no funciona parece
cluster<-pvclust(matriz2,
                 method.dist = function(dat=cenf2){
                               return(get_dist(dat, method='spearman', stand=TRUE))},
                 method.hclust = "ward.D2", nboot = 1000)
pdf(file="PLOTS/Heatmaps/Based_on_pathways/Dendrogram_wardD_pathways_logFC_spearman.pdf")
plot(cluster,hang = -1, cex = 0.5)
pvrect(cluster)
dev.off()


