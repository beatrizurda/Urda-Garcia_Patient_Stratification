#####################################################################################
# COMPUTING OTHER MOLECULAR NETWORK OVERLAP AT THE ICD9 LEVEL
#
#   We compute the overlap of:
#     - microbiome network
#     - miRNA network
#     - PPI network
#     - plus the new molecular networks
# 
#   with the epidemiological network from Hidalgo et al. over the set of 
#   diseases contained in the DSN network.
#
#
# Beatriz Urda García 2021 Jan
######################################################################################

setwd("~/Desktop/ANALYSIS/Other_molecular_networks/")

library(igraph)

dsn_path <- '~/Desktop/Desktop_20200508/Network_building/Defined_networks/icd9_pairwise_union_spearman_distance_sDEGs_network.txt'
dsn <- read.csv(dsn_path,header=T, sep="\t",stringsAsFactors = F)
my_icd9 <- unique(union(dsn$Dis1, dsn$Dis2))

##### MICROBIOME NETWORK #####

cnetwork_filename <- "microbiome_network.txt"
icds_correspondance_filename <- "icds_traslation_microbiome.csv"
icds_correspondance <- read.csv(icds_correspondance_filename,header=T, sep="\t",stringsAsFactors = F)
icds_correspondance <- icds_correspondance[, c(1,2)]
cnetwork_format <- "pajek"

cnetwork <- read.graph(cnetwork_filename, format = cnetwork_format); cnetwork
plot(cnetwork)
cnetwork_df <- as_data_frame(cnetwork, what="edges"); head(cnetwork_df); dim(cnetwork_df) 

# Selecting the positive interactions (weight > 0) They are colored in red.
cnetwork_df_pos <- cnetwork_df[cnetwork_df$weight > 0, ]; dim(cnetwork_df_pos) # 112 positive interacions
cnetwork_df_pos <- cnetwork_df_pos[, c(1:4)]; head(cnetwork_df_pos)

# Number of nodes
length(unique(union(cnetwork_df_pos$from, cnetwork_df_pos$to))) 

# Translate to ICD9 and put in correct format
merged <- merge(cnetwork_df_pos, icds_correspondance, by.x = "from", by.y = "Disease.name", all.x = TRUE, all.y = FALSE)
colnames(merged) <- c("from", "to", "weight", "color", "Dis1")
dim(merged); dim(cnetwork_df_pos)

merged <- merge(merged, icds_correspondance, by.x = "to", by.y = "Disease.name", all.x = TRUE, all.y = FALSE)
dim(merged); dim(cnetwork_df_pos)
colnames(merged)[6] <- "Dis2"

microbiome_pos_network <- merged[, c("Dis1", "Dis2", "weight")]
head(microbiome_pos_network)
colnames(microbiome_pos_network) <- c("Dis1", "Dis2", "Distance")

# Remove duplicated interactions
microbiome_pos_network$dup <- paste(microbiome_pos_network$Dis1, microbiome_pos_network$Dis2)
if(length(unique(microbiome_pos_network$dup)) == length(microbiome_pos_network$dup)){
  print("No duplicated interactions") 
}else{
  print(length(microbiome_pos_network$dup))
  print(length(unique(microbiome_pos_network$dup)))
  microbiome_pos_network <- microbiome_pos_network[!duplicated(microbiome_pos_network$dup), ]
  microbiome_pos_network <- microbiome_pos_network[, c(1:3)]
  dim(microbiome_pos_network)
}
head(microbiome_pos_network); dim(microbiome_pos_network) # 92
length(unique(union(microbiome_pos_network$Dis1, microbiome_pos_network$Dis2)))

write.table(microbiome_pos_network, file= paste("Transformed_networks/","icd9_",cnetwork_filename,sep=""), sep="\t",row.names=F, quote = FALSE)

# Keep only interactions within MY ICD9s (in our study)
microbiome_pos_network <- microbiome_pos_network[microbiome_pos_network$Dis1 %in% my_icd9, ]
microbiome_pos_network <- microbiome_pos_network[microbiome_pos_network$Dis2 %in% my_icd9, ]
dim(microbiome_pos_network) 
length(unique(union(microbiome_pos_network$Dis1, microbiome_pos_network$Dis2)))

write.table(microbiome_pos_network, file= paste("Transformed_networks/","my_icd9_",cnetwork_filename,sep=""), sep="\t",row.names=F, quote = FALSE)



##### MICRO-RNA NETWORK #####
cnetwork_filename <- "Datos_microRNAs/disease2b.net"
icds_correspondance_filename <- "icds_traslation_microRNAs.csv"
icds_correspondance <- read.csv(icds_correspondance_filename,header=T, sep=",",stringsAsFactors = F)
icds_correspondance <- icds_correspondance[, c(1,2)]; head(icds_correspondance)
cnetwork_format <- "pajek"

cnetwork <- read.graph(cnetwork_filename, format = cnetwork_format); cnetwork
plot(cnetwork); dim(cnetwork); plot(cnetwork, layout=layout_in_circle(cnetwork))
cnetwork_df <- as_data_frame(cnetwork, what="edges"); head(cnetwork_df); dim(cnetwork_df) 
length(unique(union(cnetwork_df$from, cnetwork_df$to)))

# Selecting the positive interactions (weight > 0) They are colored in red.
cnetwork_df_pos <- cnetwork_df[cnetwork_df$weight > 0, ]; dim(cnetwork_df_pos) 
cnetwork_df_pos <- cnetwork_df_pos[, c(1:3)]; head(cnetwork_df_pos)

# Translate to ICD9 and put in correct format
merged <- merge(cnetwork_df_pos, icds_correspondance, by.x = "from", by.y = "Disease.name", all.x = TRUE, all.y = FALSE)
colnames(merged) <- c("from", "to", "weight", "Dis1")
dim(merged); dim(cnetwork_df_pos)

merged <- merge(merged, icds_correspondance, by.x = "to", by.y = "Disease.name", all.x = TRUE, all.y = FALSE)
dim(merged); dim(cnetwork_df_pos)
colnames(merged)[5] <- "Dis2"

microbiome_pos_network <- merged[, c("Dis1", "Dis2", "weight")]
head(microbiome_pos_network)
colnames(microbiome_pos_network) <- c("Dis1", "Dis2", "Distance")

# Remove duplicated interactions
microbiome_pos_network$dup <- paste(microbiome_pos_network$Dis1, microbiome_pos_network$Dis2)
if(length(unique(microbiome_pos_network$dup)) == length(microbiome_pos_network$dup)){
  print("No duplicated interactions") 
}else{
  print(length(microbiome_pos_network$dup))
  print(length(unique(microbiome_pos_network$dup)))
  microbiome_pos_network <- microbiome_pos_network[!duplicated(microbiome_pos_network$dup), ]
  microbiome_pos_network <- microbiome_pos_network[, c(1:3)]
  dim(microbiome_pos_network)
}
head(microbiome_pos_network); dim(microbiome_pos_network) 

# Remove interactions that have diseases without an ICD9 category
microbiome_pos_network <- microbiome_pos_network[is.na(microbiome_pos_network$Dis1) == FALSE, ]
microbiome_pos_network <- microbiome_pos_network[is.na(microbiome_pos_network$Dis2) == FALSE, ]
dim(microbiome_pos_network) 
length(unique(union(microbiome_pos_network$Dis1, microbiome_pos_network$Dis2))) 

write.table(microbiome_pos_network, file= paste("Transformed_networks/","icd9_microRNAs_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)  

# Keep only interactions within MY ICD9s (in our study)
microbiome_pos_network <- microbiome_pos_network[microbiome_pos_network$Dis1 %in% my_icd9, ]
microbiome_pos_network <- microbiome_pos_network[microbiome_pos_network$Dis2 %in% my_icd9, ]
dim(microbiome_pos_network) #46
length(unique(union(microbiome_pos_network$Dis1, microbiome_pos_network$Dis2))) 

write.table(microbiome_pos_network, file= paste("Transformed_networks/","my_icd9_microRNAs_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)  
  

##### PPI NETWORK ##### 
#Menche et al., Science 2015

cnetwork_filename <- "ppi_network.txt"
cnetwork_df <- read.csv(paste("Datos_PPI/",cnetwork_filename,sep=""),header=T, sep="\t",stringsAsFactors = F); head(cnetwork_df); dim(cnetwork_df) # 44551

# Select positive and significant interactions (Sab > 0 AND q-valye < 0.05)
cnetwork_df <- cnetwork_df[((cnetwork_df$q..full.rand. < 0.05) & (cnetwork_df$s_AB..observed. < 0)), ]; dim(cnetwork_df)
cnetwork_df <- cnetwork_df[, c("X..disease_A", "disease_B", "s_AB..observed.", "q..full.rand.")]; head(cnetwork_df)
colnames(cnetwork_df) <- c("DisA", "DisB", "sAB", "q-value")
length(unique(union(cnetwork_df$DisA, cnetwork_df$DisB)))

# Write unique disease names to create the translation table with icd9s
ppi_disease_names <- unique(union(cnetwork_df$DisA, cnetwork_df$DisB)); length(ppi_disease_names)
# write.table(c("Disease name", ppi_disease_names), file= "icds_traslation_ppi_names.csv", sep="\n",row.names=F, col.names=F, quote = FALSE)

icds_correspondance_filename <- "icds_traslation_ppi.csv"
icds_correspondance <- read.csv(icds_correspondance_filename,header=T, sep="\t",stringsAsFactors = F)
icds_correspondance <- icds_correspondance[, c(1,2)]

# Translate to ICD9 and put in correct format
merged <- merge(cnetwork_df, icds_correspondance, by.x = "DisA", by.y = "Disease.name", all.x = TRUE, all.y = FALSE)
colnames(merged) <- c("DisA", "DisB", "sAB", "q-value", "Dis1")
dim(merged); dim(cnetwork_df)

merged <- merge(merged, icds_correspondance, by.x = "DisB", by.y = "Disease.name", all.x = TRUE, all.y = FALSE)
dim(merged); dim(cnetwork_df)
colnames(merged)[6] <- "Dis2"

ppi_pos_network <- merged[, c("Dis1", "Dis2", "sAB")]
head(ppi_pos_network)
colnames(ppi_pos_network) <- c("Dis1", "Dis2", "Distance")

# Remove the interactions that do not match to an icd9 code
ppi_pos_network <- ppi_pos_network[((is.na(ppi_pos_network$Dis1) == F) & (is.na(ppi_pos_network$Dis2) == F)) , ]
dim(ppi_pos_network)

# Remove duplicated interactions
ppi_pos_network$dup <- paste(ppi_pos_network$Dis1, ppi_pos_network$Dis2)
if(length(unique(ppi_pos_network$dup)) == length(ppi_pos_network$dup)){
  print("No duplicated interactions") 
}else{
  print(length(ppi_pos_network$dup))
  print(length(unique(ppi_pos_network$dup)))
  ppi_pos_network <- ppi_pos_network[!duplicated(ppi_pos_network$dup), ]
  ppi_pos_network <- ppi_pos_network[, c(1:3)]
  dim(ppi_pos_network)
}
head(ppi_pos_network); dim(ppi_pos_network) # 816
length(unique(union(ppi_pos_network$Dis1, ppi_pos_network$Dis2))) # 138

write.table(ppi_pos_network, file= paste("Transformed_networks/","icd9_",cnetwork_filename,sep=""), sep="\t",row.names=F, quote = FALSE)

# Keep only interactions within MY ICD9s (in our study)
ppi_pos_network <- ppi_pos_network[ppi_pos_network$Dis1 %in% my_icd9, ]
ppi_pos_network <- ppi_pos_network[ppi_pos_network$Dis2 %in% my_icd9, ]
dim(ppi_pos_network) # 43
length(unique(union(ppi_pos_network$Dis1, ppi_pos_network$Dis2)))

write.table(ppi_pos_network, file= paste("Transformed_networks/","my_icd9_",cnetwork_filename,sep=""), sep="\t",row.names=F, quote = FALSE)

##### SUBCELLULAR LOCALIZATION NETWORK #####
cnetwork_filename <- "Rev1_other_molecular_networks/Subcelllar_localization_network.csv"
cnetwork <- read.csv(cnetwork_filename,header=T, sep="\t",stringsAsFactors = F); head(cnetwork); nrow(cnetwork) # 44551

colnames(cnetwork)[c(1,3)] <- c("Dis1", "Dis2")
cnetwork$Dis1 = gsub("\\..*","", cnetwork$Dis1)
cnetwork$Dis2 = gsub("\\..*","", cnetwork$Dis2)

cnetwork$Dis1 = as.numeric(cnetwork$Dis1)
cnetwork$Dis2 = as.numeric(cnetwork$Dis2)
cnetwork = sort_icd_interactions(cnetwork)

# Remove duplicated interactions
cnetwork$dup = paste(cnetwork$Dis1, cnetwork$Dis2, sep="_")
if(length(unique(cnetwork$dup)) == length(cnetwork$dup)){
  print("No duplicated interactions") 
}else{
  print(length(cnetwork$dup))
  print(length(unique(cnetwork$dup)))
  cnetwork <- cnetwork[!duplicated(cnetwork$dup), ]
  nrow(cnetwork)    
}

length(unique(union(cnetwork$Dis1, cnetwork$Dis2))) 
cnetwork$Distance = cnetwork$Loc_PCC

write.table(cnetwork, file= paste("Transformed_networks/","icd9_","subcellular_localization_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)

# Keep only interactions within MY ICD9s (in our study)
cnetwork <- cnetwork[cnetwork$Dis1 %in% my_icd9, ]
cnetwork <- cnetwork[cnetwork$Dis2 %in% my_icd9, ]
nrow(cnetwork) # 112
length(unique(union(cnetwork$Dis1, cnetwork$Dis2))) 

write.table(cnetwork, file= paste("Transformed_networks/","my_icd9_","subcellular_localization_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)


##### CELLULAR COMPONENT NETWORK #####
cnetwork_filename <- "Rev1_other_molecular_networks/Network_cellular.csv"
cnetwork <- read.csv(cnetwork_filename,header=T, sep="\t",stringsAsFactors = F); head(cnetwork); dim(cnetwork) 

# Formatting ICD9 names
colnames(cnetwork)[c(1,3)] <- c("Dis1", "Dis2")
cnetwork$Dis1 = gsub("\\..*","", cnetwork$Dis1)
cnetwork$Dis2 = gsub("\\..*","", cnetwork$Dis2)

cnetwork$Dis1 = gsub("\\[|\\]","",cnetwork$Dis1)
cnetwork$Dis2 = gsub("\\[|\\]","",cnetwork$Dis2)

# Changing to numeric and sorting the interactions
cnetwork$Dis1 = as.numeric(cnetwork$Dis1)
cnetwork$Dis2 = as.numeric(cnetwork$Dis2)
cnetwork = sort_icd_interactions(cnetwork)

# Remove duplicated interactions
cnetwork$dup = paste(cnetwork$Dis1, cnetwork$Dis2, sep="_")
if(length(unique(cnetwork$dup)) == length(cnetwork$dup)){
  print("No duplicated interactions") 
}else{
  print(length(cnetwork$dup))
  print(length(unique(cnetwork$dup)))
  cnetwork <- cnetwork[!duplicated(cnetwork$dup), ]
  nrow(cnetwork)    
}

length(unique(union(cnetwork$Dis1, cnetwork$Dis2))) 

cnetwork_ppis = cnetwork[cnetwork$NP >= 1, ]; cnetwork_ppis$Distance = cnetwork_ppis$NP; nrow(cnetwork_ppis)   
cnetwork_genes = cnetwork[cnetwork$NG >= 1, ]; cnetwork_genes$Distance = cnetwork_genes$NG; nrow(cnetwork_genes) 
cnetwork$Distance = cnetwork$NP + cnetwork$NG; nrow(cnetwork)

write.table(cnetwork, file= paste("Transformed_networks/","icd9_","cellular_components_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)
write.table(cnetwork_ppis, file= paste("Transformed_networks/","icd9_","cellular_components_network_ppis.txt",sep=""), sep="\t",row.names=F, quote = FALSE)
write.table(cnetwork_genes, file= paste("Transformed_networks/","icd9_","cellular_components_network_genes.txt",sep=""), sep="\t",row.names=F, quote = FALSE)

# Keep only interactions within MY ICD9s (in our study)
cnetwork <- cnetwork[cnetwork$Dis1 %in% my_icd9, ]
cnetwork <- cnetwork[cnetwork$Dis2 %in% my_icd9, ]
nrow(cnetwork) # 118
length(unique(union(cnetwork$Dis1, cnetwork$Dis2)))  

# same for ppis
cnetwork_ppis <- cnetwork_ppis[cnetwork_ppis$Dis1 %in% my_icd9, ]
cnetwork_ppis <- cnetwork_ppis[cnetwork_ppis$Dis2 %in% my_icd9, ]
nrow(cnetwork_ppis) # 104
length(unique(union(cnetwork_ppis$Dis1, cnetwork_ppis$Dis2))) 

# same for genes
cnetwork_genes <- cnetwork_genes[cnetwork_genes$Dis1 %in% my_icd9, ]
cnetwork_genes <- cnetwork_genes[cnetwork_genes$Dis2 %in% my_icd9, ]
nrow(cnetwork_genes) # 39
length(unique(union(cnetwork_genes$Dis1, cnetwork_genes$Dis2))) 


write.table(cnetwork, file= paste("Transformed_networks/","my_icd9_","cellular_components_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)
write.table(cnetwork_ppis, file= paste("Transformed_networks/","my_icd9_","cellular_components_network_ppis.txt",sep=""), sep="\t",row.names=F, quote = FALSE)
write.table(cnetwork_genes, file= paste("Transformed_networks/","my_icd9_","cellular_components_network_genes.txt",sep=""), sep="\t",row.names=F, quote = FALSE)

##### GEO Hu GSK Network #####
# Net1 based on correlations
cnetwork_filename <- "Rev1_other_molecular_networks/from_review/pone.0006536.s002.txt"
cnetwork <- read.csv(cnetwork_filename,header=T, sep="\t",stringsAsFactors = F); head(cnetwork); dim(cnetwork)
nrow(cnetwork) 

# Net2 based on enrichments
cnetwork_filename2 <- "Rev1_other_molecular_networks/from_review/pone.0006536.s004.txt"
cnetwork2 <- read.csv(cnetwork_filename2,header=T, sep="\t",stringsAsFactors = F); head(cnetwork); dim(cnetwork) 
nrow(cnetwork2) 

# Select interactions between diseases in both networks
cnetwork = cnetwork[(cnetwork$s_comparison == 'disease state') & (cnetwork$p_comparison == 'disease state'), ]
nrow(cnetwork) 
cnetwork2 = cnetwork2[(cnetwork2$s_comparison == 'disease state') & (cnetwork2$p_comparison == 'disease state'), ]
nrow(cnetwork2) 

# Get disease anmes
diseasenames_gsk = unique(union(cnetwork$s_sampleB, cnetwork$p_sampleB))
length(diseasenames_gsk)  

diseasenames_gsk2 = unique(union(cnetwork2$s_sampleB, cnetwork2$p_sampleB))
length(diseasenames_gsk2)   

setdiff(diseasenames_gsk2, diseasenames_gsk)
# write.table(diseasenames_gsk, file= "Rev1_other_molecular_networks/from_review/icd9_gsk_geo.txt", sep="\t",row.names=F, quote = FALSE)

dic_gsk = read.csv("Rev1_other_molecular_networks/from_review/icd9_gsk_geo.txt",header=T, sep=",",stringsAsFactors = F)
head(dic_gsk)

dic_gsk = dic_gsk[is.na(dic_gsk$icd9) == FALSE, ]; nrow(dic_gsk) 
head(cnetwork)

valid_disnames = dic_gsk$disease_name
print(length(unique(dic_gsk$icd9)))

#####  Transform to icd9 and format: cnetwork
cnetwork = cnetwork[((cnetwork$s_sampleB %in% valid_disnames) & (cnetwork$p_sampleB %in% valid_disnames)), ]; nrow(cnetwork) 

cnetwork = merge(cnetwork, dic_gsk, by.x = "s_sampleB", by.y = "disease_name", all.x = TRUE, all.y = FALSE); nrow(cnetwork) 
cnetwork = merge(cnetwork, dic_gsk, by.x = "p_sampleB", by.y = "disease_name", all.x = TRUE, all.y = FALSE); nrow(cnetwork) 

cnetwork = cnetwork[, c('icd9.x', 'icd9.y', 'corr')]
colnames(cnetwork) = c("Dis1", "Dis2", "Distance")

cnetwork = sort_icd_interactions(cnetwork); nrow(cnetwork) 

# Remove duplicated interactions
cnetwork$dup = paste(cnetwork$Dis1, cnetwork$Dis2, sep="_")
if(length(unique(cnetwork$dup)) == length(cnetwork$dup)){
  print("No duplicated interactions") 
}else{
  print(length(cnetwork$dup))
  print(length(unique(cnetwork$dup)))
  cnetwork <- cnetwork[!duplicated(cnetwork$dup), ]
  nrow(cnetwork)   
}

cnetwork = cnetwork[cnetwork$Dis1 != cnetwork$Dis2, ]; nrow(cnetwork)   

cnetworkpos = cnetwork[cnetwork$Distance >= 0, ]; nrow(cnetworkpos) 
cnetworkneg = cnetwork[cnetwork$Distance < 0, ]; nrow(cnetworkneg)   

##### Transform to icd9 and format: cnetwork2
cnetwork2 = cnetwork2[((cnetwork2$s_sampleB %in% valid_disnames) & (cnetwork2$p_sampleB %in% valid_disnames)), ]; nrow(cnetwork2) 

cnetwork2 = merge(cnetwork2, dic_gsk, by.x = "s_sampleB", by.y = "disease_name", all.x = TRUE, all.y = FALSE); nrow(cnetwork2) 
cnetwork2 = merge(cnetwork2, dic_gsk, by.x = "p_sampleB", by.y = "disease_name", all.x = TRUE, all.y = FALSE); nrow(cnetwork2)

cnetwork2 = cnetwork2[, c('icd9.x', 'icd9.y', 'score')]
colnames(cnetwork2) = c("Dis1", "Dis2", "Distance")

cnetwork2 = sort_icd_interactions(cnetwork2); nrow(cnetwork2) 

# Remove duplicated interactions
cnetwork2$dup = paste(cnetwork2$Dis1, cnetwork2$Dis2, sep="_")
if(length(unique(cnetwork2$dup)) == length(cnetwork2$dup)){
  print("No duplicated interactions") 
}else{
  print(length(cnetwork2$dup))
  print(length(unique(cnetwork2$dup)))
  cnetwork2 <- cnetwork2[!duplicated(cnetwork2$dup), ]
  nrow(cnetwork2)    
}

cnetwork2 = cnetwork2[cnetwork2$Dis1 != cnetwork2$Dis2, ]; nrow(cnetwork2)    

cnetwork2pos = cnetwork2[cnetwork2$Distance >= 0, ]; nrow(cnetwork2pos)  
cnetwork2neg = cnetwork2[cnetwork2$Distance < 0, ]; nrow(cnetwork2neg)   

allinter = rbind(cnetworkpos, cnetwork2pos) ; nrow(allinter)
allinter = allinter[!duplicated(allinter$dup), ]; nrow(allinter)

write.table(allinter, file= paste("Transformed_networks/","icd9_","geo_gsk_pos_union_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)
write.table(cnetworkpos, file= paste("Transformed_networks/","icd9_","geo_gsk_pos_1_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)
write.table(cnetwork2pos, file= paste("Transformed_networks/","icd9_","geo_gsk_pos_2_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)

# Keep only interactions within MY ICD9s (in our study)
allinter <- allinter[allinter$Dis1 %in% my_icd9, ]
allinter <- allinter[allinter$Dis2 %in% my_icd9, ]; nrow(allinter) 
allinter = allinter[allinter$Distance >= 0, ]; nrow(allinter) 
length(unique(union(allinter$Dis1, allinter$Dis2)))   

cnetworkpos <- cnetworkpos[cnetworkpos$Dis1 %in% my_icd9, ]
cnetworkpos <- cnetworkpos[cnetworkpos$Dis2 %in% my_icd9, ]; nrow(cnetworkpos)

cnetwork2pos <- cnetwork2pos[cnetwork2pos$Dis1 %in% my_icd9, ]
cnetwork2pos <- cnetwork2pos[cnetwork2pos$Dis2 %in% my_icd9, ]; nrow(cnetwork2pos)

write.table(allinter, file= paste("Transformed_networks/","my_icd9_","geo_gsk_pos_union_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)
write.table(cnetworkpos, file= paste("Transformed_networks/","my_icd9_","geo_gsk_pos_1_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)
write.table(cnetwork2pos, file= paste("Transformed_networks/","my_icd9_","geo_gsk_pos_2_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)