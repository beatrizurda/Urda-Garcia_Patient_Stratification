#####################################################################################
# COMPUTING NETWORK OVERLAP AT THE DISEASE LEVEL
#
#   
#  
#  
# 
# Beatriz Urda García 2020
######################################################################################

library(VennDiagram)

# PARAMETERS
setwd("~/Desktop/ANALYSIS")

de_network_filename <- "comparable_spearman_distance_sDEGs_pos_network.txt"
dv_network_filename <- "comparable_spearman_distance_sDVGs_pos_network.txt"

### FUNCTIONS
sort_icd_interactions <- function(df){
  # Sort interactions such that icd of Dis1 < icd of Dis2
  
  for (k in 1:length(df$Dis1)){
    if (df$Dis1[k] > df$Dis2[k]){
      dis1 <- df$Dis1[k]
      df$Dis1[k] <- df$Dis2[k]
      df$Dis2[k] <- dis1
    }
  }
  return(df)
}

filedir <- 'Network_building/Defined_networks/'
de_network <- read.csv(paste(filedir, de_network_filename, sep=""),header=T, sep="\t",stringsAsFactors = F)
dv_network <- read.csv(paste(filedir, dv_network_filename, sep=""),header=T, sep="\t",stringsAsFactors = F)

head(de_network); dim(de_network) # 505
head(dv_network); dim(dv_network)

# SORT INTERACTIONS
de_network <- sort_icd_interactions(de_network)
dv_network <- sort_icd_interactions(dv_network)

# Venn diagram of new vs unique interactions
DE_interactions <- paste(de_network$Dis1, de_network$Dis2, sep="_")
DV_interactions <- paste(dv_network$Dis1, dv_network$Dis2, sep="_")

venn.diagram(x=list(DE_interactions, DV_interactions), category.names = c("DE Interactions","DV Interactions"), filename="Network_building/Overall_overlap/venn_diagram", output=TRUE)

############# Networks transformed into ICD9:
secondfiledir_de <- 'Network_building/Overlapping_results/Shuffling_labels/DE_overlap/networks_DE_overlap/comparable_spearman_distance_sDEGs_pos_consistent_C_TRUE_final_network.txt'
secondfiledir_dv <- 'Network_building/Overlapping_results/Shuffling_labels/DV_overlap/networks_DV_overlap/comparable_spearman_distance_sDVGs_pos_consistent_C_TRUE_final_network.txt'

secondfiledir_de <- 'Network_building/Overlapping_results/Shuffling_labels/DE_overlap/networks_DE_overlap/comparable_spearman_distance_sDEGs_pos_allpos_C_TRUE_final_network.txt'
secondfiledir_dv <- 'Network_building/Overlapping_results/Shuffling_labels/DV_overlap/networks_DV_overlap/comparable_spearman_distance_sDVGs_pos_allpos_C_TRUE_final_network.txt'

overlap_condition <- gsub('Network_building/Overlapping_results/Shuffling_labels/DE_overlap/networks_DE_overlap/comparable_spearman_distance_sDEGs_pos_',"",secondfiledir_de)
overlap_condition <- gsub('_TRUE_final_network.txt',"",overlap_condition)

fbarabasi1 <- gsub("network.txt", "barabasi.txt", secondfiledir_de,)
fbarabasi2 <- gsub("network.txt", "barabasi.txt", secondfiledir_dv)
  
de_network <- read.csv(secondfiledir_de,header=T, sep="\t",stringsAsFactors = F); dim(de_network)
dv_network <- read.csv(secondfiledir_dv ,header=T, sep="\t",stringsAsFactors = F); dim(dv_network)

n_icds_de <-  length(unique(union(de_network$Dis1, de_network$Dis2))); n_icds_de # Option B: 41   | OPTION C: 30
n_icds_dv <- length(unique(union(dv_network$Dis1, dv_network$Dis2))); n_icds_dv # Option B: 41   | OPTION C: 37

common_icds <- intersect(unique(union(de_network$Dis1, de_network$Dis2)), unique(union(dv_network$Dis1, dv_network$Dis2)))
length(common_icds) # 29


barabasi_de <- read.csv(fbarabasi1 ,header=T, sep="\t",stringsAsFactors = F); dim(barabasi_de)
barabasi_dv <- read.csv(fbarabasi2 ,header=T, sep="\t",stringsAsFactors = F); dim(barabasi_dv)

# barabasi_de == barabasi_dv # THEY ARE EQUAL

de_network <- sort_icd_interactions(de_network); dim(de_network)
dv_network <- sort_icd_interactions(dv_network); dim(dv_network )
barabasi_de <- sort_icd_interactions(barabasi_de); dim(barabasi_de)
barabasi_dv <- sort_icd_interactions(barabasi_dv); dim(barabasi_dv)

common_de_network <- de_network[(de_network$Dis1 %in% common_icds) & (de_network$Dis2 %in% common_icds), ]; dim(common_de_network)
common_dv_network <- dv_network[(dv_network$Dis1 %in% common_icds) & (dv_network$Dis2 %in% common_icds), ]; dim(common_dv_network)
common_barabasi_de <- barabasi_de[(barabasi_de$Dis1 %in% common_icds) & (barabasi_de$Dis2 %in% common_icds), ]; dim(common_barabasi_de); rownames(common_barabasi_de) <- NULL
common_barabasi_dv <- barabasi_dv[(barabasi_dv$Dis1 %in% common_icds) & (barabasi_dv$Dis2 %in% common_icds), ]; dim(common_barabasi_dv); rownames(common_barabasi_dv) <- NULL

if(all.equal(common_barabasi_de,common_barabasi_dv)){
  common_barabasi <- barabasi_de
  print("Barabasi dv and de are equal")
}else{
  print("Barabasi dv and de are not equal")
}

DE_interactions <- paste(de_network$Dis1, de_network$Dis2, sep="_"); length(DE_interactions)
DV_interactions <- paste(dv_network$Dis1, dv_network$Dis2, sep="_"); length(DV_interactions)
Barabasi_interactions <- paste(barabasi_de$Dis1, barabasi_de$Dis2, sep="_"); length(Barabasi_interactions)
Barabasi_interactions_dv <- paste(barabasi_dv$Dis1, barabasi_dv$Dis2, sep="_"); length(Barabasi_interactions_dv)
Barabasi_total <- append(Barabasi_interactions, setdiff(Barabasi_interactions_dv, Barabasi_interactions)); length(setdiff(Barabasi_interactions_dv, Barabasi_interactions))
length(intersect(Barabasi_interactions, Barabasi_interactions_dv));length(Barabasi_total)

venn.diagram(x=list(DE_interactions, DV_interactions, Barabasi_total), 
             category.names = c("DE Interactions","DV Interactions", "Barabasi interactions"), 
             filename=paste("Network_building/Overall_overlap/venn_diagram_overall_barabasi_",overlap_condition), output=TRUE)


total_barabasi <- length(Barabasi_total); total_barabasi
total_de <- length(intersect(DE_interactions,Barabasi_total)); total_de; (total_de/total_barabasi)*100
total_dv <- length(intersect(DV_interactions,Barabasi_total)); total_dv; (total_dv/total_barabasi)*100
union_dv_de <- union(DE_interactions, DV_interactions); length(union_dv_de)
total_coverage <- length(intersect(union_dv_de,Barabasi_total)); total_coverage; (total_coverage/total_barabasi)*100

(total_de/length(DE_interactions))*100   # 38.84298
(total_dv/length(DV_interactions))*100   # 45.1613
(total_coverage/length(union_dv_de))*100 # 40%

# With the common ICDs
common_de_interactions <- paste(common_de_network$Dis1, common_de_network$Dis2, sep="_"); length(DE_interactions)
common_dv_interactions <- paste(common_dv_network$Dis1, common_dv_network$Dis2, sep="_"); length(DV_interactions)
common_barabasi_interactions <- paste(common_barabasi$Dis1, common_barabasi$Dis2, sep="_"); length(Barabasi_interactions)

venn.diagram(x=list(common_de_interactions, common_dv_interactions, common_barabasi_interactions), 
             category.names = c("DE Interactions","DV Interactions", "Barabasi interactions"), 
             filename=paste("Network_building/Overall_overlap/venn_diagram_overall_barabasi_common_icds_",overlap_condition), output=TRUE)

total_barabasi <- length(common_barabasi_interactions); total_barabasi
total_de <- length(intersect(common_de_interactions,common_barabasi_interactions)); total_de; (total_de/total_barabasi)*100
total_dv <- length(intersect(common_dv_interactions,common_barabasi_interactions)); total_dv; (total_dv/total_barabasi)*100
union_dv_de <- union(common_de_interactions, common_dv_interactions); length(union_dv_de)
total_coverage <- length(intersect(union_dv_de,common_barabasi_interactions)); total_coverage; (total_coverage/total_barabasi)*100

(total_de/length(common_de_interactions))*100   # 38.84298
(total_dv/length(common_dv_interactions))*100   # 45.1613
(total_coverage/length(union_dv_de))*100 # 40%


