#####################################################################################
# COMPUTING OTHER MOLECULAR NETWORK OVERLAP AT THE ICD9 LEVEL
#
#   We compute the overlap of:
#     - microbiome network
#     - miRNA network
#     - PPI network
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
cnetwork_df <- as_data_frame(cnetwork, what="edges"); head(cnetwork_df); dim(cnetwork_df) # 180

# Selecting the positive interactions (weight > 0) They are colored in red.
cnetwork_df_pos <- cnetwork_df[cnetwork_df$weight > 0, ]; dim(cnetwork_df_pos) # 112 positive interacions
cnetwork_df_pos <- cnetwork_df_pos[, c(1:4)]; head(cnetwork_df_pos)

# Number of nodes
length(unique(union(cnetwork_df_pos$from, cnetwork_df_pos$to))) # 33

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
dim(microbiome_pos_network) # 9
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
cnetwork_df <- as_data_frame(cnetwork, what="edges"); head(cnetwork_df); dim(cnetwork_df) # 414
length(unique(union(cnetwork_df$from, cnetwork_df$to)))

# Selecting the positive interactions (weight > 0) They are colored in red.
cnetwork_df_pos <- cnetwork_df[cnetwork_df$weight > 0, ]; dim(cnetwork_df_pos) # 414 positive interacions - ALL
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
head(microbiome_pos_network); dim(microbiome_pos_network) # 363

# Remove interactions that have diseases without an ICD9 category
microbiome_pos_network <- microbiome_pos_network[is.na(microbiome_pos_network$Dis1) == FALSE, ]
microbiome_pos_network <- microbiome_pos_network[is.na(microbiome_pos_network$Dis2) == FALSE, ]
dim(microbiome_pos_network) #355
length(unique(union(microbiome_pos_network$Dis1, microbiome_pos_network$Dis2))) # 46

write.table(microbiome_pos_network, file= paste("Transformed_networks/","icd9_microRNAs_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)  

# Keep only interactions within MY ICD9s (in our study)
microbiome_pos_network <- microbiome_pos_network[microbiome_pos_network$Dis1 %in% my_icd9, ]
microbiome_pos_network <- microbiome_pos_network[microbiome_pos_network$Dis2 %in% my_icd9, ]
dim(microbiome_pos_network) #46
length(unique(union(microbiome_pos_network$Dis1, microbiome_pos_network$Dis2))) # 13

write.table(microbiome_pos_network, file= paste("Transformed_networks/","my_icd9_microRNAs_network.txt",sep=""), sep="\t",row.names=F, quote = FALSE)  
  

##### PPI NETWORK ##### 
#Menche et al., Science 2015

cnetwork_filename <- "ppi_network.txt"
cnetwork_df <- read.csv(paste("Datos_PPI/",cnetwork_filename,sep=""),header=T, sep="\t",stringsAsFactors = F); head(cnetwork_df); dim(cnetwork_df) # 44551

# Select positive and significant interactions (Sab > 0 AND q-valye < 0.05)
cnetwork_df <- cnetwork_df[((cnetwork_df$q..full.rand. < 0.05) & (cnetwork_df$s_AB..observed. < 0)), ]; dim(cnetwork_df) # 1383
cnetwork_df <- cnetwork_df[, c("X..disease_A", "disease_B", "s_AB..observed.", "q..full.rand.")]; head(cnetwork_df)
colnames(cnetwork_df) <- c("DisA", "DisB", "sAB", "q-value")
length(unique(union(cnetwork_df$DisA, cnetwork_df$DisB)))

# Write unique disease names to create the translation table with icd9s
ppi_disease_names <- unique(union(cnetwork_df$DisA, cnetwork_df$DisB)); length(ppi_disease_names) # Number of nodes: 289
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

