#####################################################################################
# MERGING SE OBJECTS THAT CORRESPOND TO THE SAME DISEASE
#
#   
#   
# 
#
######################################################################################

library(SummarizedExperiment)

setwd("~/Desktop/ANALYSIS/")

# Breast Cancer
se1 <- readRDS("GREIN_SE_good_quality_initial/BreastCancer_EstrogenReceptorPositive_grein_se.rds") 
se2 <- readRDS("GREIN_SE_good_quality_initial/BreastCancer_TripleNegative_grein_se.rds")
se3 <- readRDS("GREIN_SE_good_quality_initial/DuctalBreastCancer_grein_se.rds")

se <- cbind(se1,se2)
se <- cbind(se,se3) # 164 samples

se <- se[,!duplicated(se$sample)]
length(se$sample)
length(unique(se$sample))

saveRDS(se,"GREIN_SE_good_quality/BreastCancer_grein_se.rds")

# HIV
se1 <- readRDS("GREIN_SE_good_quality_initial/HIV_grein_se.rds") 
se2 <- readRDS("GREIN_SE_good_quality_initial/HIVART-treated_grein_se.rds")

se <- cbind(se1,se2)

se <- se[,!duplicated(se$sample)]
length(se$sample)
length(unique(se$sample))

saveRDS(se,"GREIN_SE_good_quality/HIV_grein_se.rds")

# LUNG CANCER
se1 <- readRDS("GREIN_SE_good_quality_initial/LungCancer_grein_se.rds") 
se2 <- readRDS("GREIN_SE_good_quality_initial/NSCLC_grein_se.rds")

se <- cbind(se1,se2)

se <- se[,!duplicated(se$sample)]
length(se$sample)
length(unique(se$sample))

saveRDS(se,"GREIN_SE_good_quality/LungCancer_grein_se.rds")

# COLORRECTAL CANCER
se1 <- readRDS("GREIN_SE_good_quality_initial/ColorectalCancer_grein_se.rds") 
se2 <- readRDS("GREIN_SE_good_quality_initial/ColorectalCancer_Adenocarcinoma_grein_se.rds")

se <- cbind(se1,se2)

se <- se[,!duplicated(se$sample)]
length(se$sample)
length(unique(se$sample))

saveRDS(se,"GREIN_SE_good_quality/ColorectalCancer_grein_se.rds")

# LIVER CANCER, THAT HAS MATRICES WITH DIFFERENT NUMBER OF GENES

se1 <- readRDS("GREIN_SE_good_quality_initial/IntrahepaticCholangiocarcinoma_grein_se.rds")
se2 <- readRDS("GREIN_SE_good_quality_initial/LiverCancer_grein_se.rds")
se3 <- readRDS("GREIN_SE_good_quality_initial/HepatocellularCarcinoma_grein_se.rds")
length(se1$sample) # 14
length(se2$sample) # 64
length(se3$sample) # 99
str(se1)
class(assays(se1)$counts)
length(rownames(assays(se1)$counts)) # 28089
length(rownames(assays(se2)$counts)) # 27885
length(rownames(assays(se3)$counts)) # 27885
common <- intersect(rownames(assays(se2)$counts),rownames(assays(se2)$counts))
common <- intersect(common,rownames(assays(se3)$counts))
length(common) # 27885 common set of genes

counts1 <- assays(se1)$counts[which(rownames(assays(se1)$counts) %in% common),]
counts2 <- assays(se2)$counts[which(rownames(assays(se2)$counts) %in% common),]
counts3 <- assays(se3)$counts[which(rownames(assays(se3)$counts) %in% common),]
nse1 <- SummarizedExperiment(assays=list(counts=counts1), colData=colData(se1))
nse2 <- SummarizedExperiment(assays=list(counts=counts2), colData=colData(se2))
nse3 <- SummarizedExperiment(assays=list(counts=counts3), colData=colData(se3))
fse <- cbind(nse1,nse2); fse <- cbind(fse,nse3) # 177 samples
length(fse$sample)
fse <- fse[,!duplicated(fse$sample)]
length(fse$sample)
saveRDS(fse,file="GREIN_SE_good_quality/LiverCancer_grein_se.rds")


# REMOVING OLD PHENOTYPES! IMP: do not remove the ones that you have overwritten
to_remove <- c("BreastCancer_EstrogenReceptorPositive_grein_se.rds","BreastCancer_TripleNegative_grein_se.rds",
               "DuctalBreastCancer_grein_se.rds","HIVART-treated_grein_se.rds",
               "NSCLC_grein_se.rds","HepatocellularCarcinoma_grein_se.rds",
               "IntrahepaticCholangiocarcinoma_grein_se.rds","ColorectalCancer_Adenocarcinoma_grein_se.rds")
for (file in to_remove){
  cfile <- paste("GREIN_SE_good_quality/",file,sep="")
  if (file.exists(cfile)) 
    #Delete file if it exists
    file.remove(cfile)
}

