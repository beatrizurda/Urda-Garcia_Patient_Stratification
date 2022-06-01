###################################################################
#  Find newly detected comorbidities in SSN with respect to DSN
#
#
#   Beatriz Urda 2022
###################################################################

library(stringr)
source("analysis_library.R")

# ANOTHER APPROACH:
# DSN
dsn = read.csv("shiny_network_app/final_pairwise_union_spearman_distance_sDEGs_network.txt", sep="\t", stringsAsFactors = FALSE)
head(dsn)
nrow(dsn) # 658

# SSN
ssn = read.csv("shiny_network_app/final_metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt", sep="\t", stringsAsFactors = FALSE)
head(ssn)
nrow(ssn) # 9017
ssn$Dis1 = str_trim(ssn$Dis1)
ssn$Dis2 = str_trim(ssn$Dis2)

# Select only positive interactions in both networks
dsn = dsn[dsn$Distance >= 0, ]; nrow(dsn)  
ssn = ssn[ssn$Distance >= 0, ]; nrow(ssn)  

# Select D-M or M-D interactions
ssn = ssn[(((ssn$Dis1 == ssn$Corr_Dis1) & (ssn$Dis2 != ssn$Corr_Dis2)) | ((ssn$Dis1 != ssn$Corr_Dis1) & (ssn$Dis2 == ssn$Corr_Dis2))), ]

# Select epidemiological interactions
ssn = ssn[ssn$in_epidem == TRUE, ]; nrow(ssn) # 1006

# NEW interactions at the icd9 level
# Sort the interactions by icd9
dsn = sort_icd_interactions(dsn, Dis1 = "Dis1_icd9", Dis2 = "Dis2_icd9")
ssn = sort_icd_interactions(ssn, Dis1 = "Dis1_icd9", Dis2 = "Dis2_icd9")

# NEW interactions at the disease level
dsn = sort_icd_interactions(dsn)
ssn = sort_icd_interactions(ssn, Dis1 = "Corr_Dis1", Dis2 = "Corr_Dis2")

dsn$pairs = paste(dsn$Dis1, dsn$Dis2, sep="_")
ssn$pairs = paste(ssn$Corr_Dis1, ssn$Corr_Dis2, sep="_")

only_ssn = ssn[!(ssn$pairs %in% dsn$pairs), ]; nrow(only_ssn) 
dsn_ssn = ssn[ssn$pairs %in% dsn$pairs, ]; nrow(dsn_ssn) 

length(unique(only_ssn$pairs))  

new_unique_ssn = only_ssn[!duplicated(only_ssn$pairs), ]; nrow(new_unique_ssn) 
new_unique_ssn$icd9pairs = paste(new_unique_ssn$Dis1_icd9, new_unique_ssn$Dis2_icd9, sep="_")

df = as.data.frame(table(append(new_unique_ssn$Corr_Dis1, new_unique_ssn$Corr_Dis2)))

to_save = new_unique_ssn[, c("Corr_Dis1", "Corr_Dis2")]
write.table(to_save, "newly_comorbidities_SSN.txt", quote=FALSE, sep = "\t", row.names = FALSE)



