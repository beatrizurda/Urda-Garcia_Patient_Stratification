#########################################################################################
# COMPARISON OF THE DSN WITH THE PPI from Menche et al. 2015
# 
# IN ~/Desktop/ANALYSIS/Other_molecular_networks
#
# Beatriz Urda García 2021
#########################################################################################

setwd("~/Desktop/ANALYSIS/Other_molecular_networks/")

library(igraph)
library(VennDiagram)
library(ggplot2)

dsn_path <- '~/Desktop/Desktop_20200508/Network_building/Defined_networks/icd9_pairwise_union_spearman_distance_sDEGs_network.txt'
dsn <- read.csv(dsn_path,header=T, sep="\t",stringsAsFactors = F)
my_icd9 <- unique(union(dsn$Dis1, dsn$Dis2)); length(my_icd9) # 41
dsn_pos <- dsn[dsn$Distance >= 0, ]

# Read and filter barabasi network
barabasi <- read.csv("~/Desktop/Desktop_20200508/Network_building/Overlapping_results/Shuffling_labels/icd9_pairwise_union_spearman_distance_sDEGs_pos_B_TRUE_final_barabasi.txt",
                     header=T, sep="\t",stringsAsFactors = F)
dim(barabasi) # 331
barabasi_graph <- graph_from_data_frame(barabasi, directed=FALSE)
plot(barabasi_graph)

ppinetwork_path <- 'Transformed_networks/Nodups_my_icd9_ppi_network.txt'
ppinetwork <- read.csv(ppinetwork_path,header=T, sep="\t",stringsAsFactors = F)
ppi_icd9 <- unique(union(ppinetwork$Dis1, ppinetwork$Dis2)); length(ppi_icd9) # 19
dim(ppinetwork)

big_ppinetwork_path <- 'Transformed_networks/Nodups_icd9_ppi_network.txt'
big_ppinetwork <- read.csv(big_ppinetwork_path,header=T, sep="\t",stringsAsFactors = F)
big_ppi_icd9 <- unique(union(big_ppinetwork$Dis1, big_ppinetwork$Dis2)); length(big_ppi_icd9) # 136

only_dsn <- setdiff(my_icd9, ppi_icd9); length(setdiff(my_icd9, ppi_icd9)) # 22
only_ppi <- setdiff(big_ppi_icd9, my_icd9); length(setdiff(big_ppi_icd9, my_icd9)) # 111

# venn.diagram(x = list(my_icd9, big_ppi_icd9), category.names = c("DSN", "PPI"), filename = "ICD9_venndiagram_big.pdf", output=TRUE)
# venn.diagram(x = list(my_icd9, ppi_icd9), category.names = c("DSN", "PPI"), filename = "ICD9_venndiagram.pdf_comparabe.pdf", output=TRUE)

# There is a discrepancy of two diseases:
# The next 6 diseases only have interactions with diseases that do not belong to the DSN
setdiff(intersect(my_icd9,big_ppi_icd9), ppi_icd9)

common_icds <- unique(intersect(my_icd9, ppi_icd9)); length(common_icds) # 19

dsn_common <- dsn_pos[((dsn_pos$Dis1 %in% common_icds) & (dsn_pos$Dis2 %in% common_icds)), ]
dim(dsn_common)

dsn_graph <- graph_from_data_frame(dsn_pos, directed=FALSE)
dsn_common_graph <- graph_from_data_frame(dsn_common, directed=FALSE)
ppi_graph <- graph_from_data_frame(ppinetwork, directed=FALSE)
plot(dsn_graph); gsize(dsn_graph) # 20
plot(dsn_common_graph); gsize(dsn_common_graph) # 73
plot(ppi_graph); gsize(ppi_graph) # 20

intersect <- intersection(dsn_common_graph,ppi_graph)
plot(intersect)
print(gsize(intersect)) # 6
intersect_in_epidem <- intersection(intersect, barabasi_graph); gsize(intersect_in_epidem) # 5
plot(intersect_in_epidem)
intersect_df <- as_data_frame(intersect_in_epidem)
perc_overlap_dsn <- (gsize(intersect)/dim(dsn_common)[1])*100 # 8.22%
print(paste('% of overlap from the DSN: ',perc_overlap_dsn, sep=''))
perc_overlap_ppi <- (gsize(intersect)/dim(ppinetwork)[1])*100 # 30%
print(paste('% of network from the PPI: ',perc_overlap_ppi, sep=''))

only_dsn_int <- difference(dsn_common_graph, ppi_graph); plot(only_dsn_int)
gsize(only_dsn_int) # 67
only_dsn_int_epidem <- intersection(difference(dsn_common_graph, ppi_graph), barabasi_graph); plot(only_dsn_int_epidem); 
gsize(only_dsn_int_epidem) # 29
only_dsn_int_epidem_df <- as_data_frame(only_dsn_int_epidem)

only_ppi_int <- difference(ppi_graph, dsn_common_graph); plot(only_ppi_int)
gsize(only_ppi_int) # 14
only_ppi_int_epidem <- intersection(difference(ppi_graph, dsn_common_graph), barabasi_graph); plot(only_ppi_int_epidem)
gsize(only_ppi_int_epidem) # 10
only_ppi_int_epidem_df <- as_data_frame(only_ppi_int_epidem)

only_dsn_int_all <- difference(dsn_graph, ppi_graph); plot(only_dsn_int_all)
gsize(only_dsn_int_all) # 341
only_dsn_int_epidem_all <- intersection(difference(dsn_graph, ppi_graph), barabasi_graph); plot(only_dsn_int_epidem_all); 
gsize(only_dsn_int_epidem_all) # 149
only_dsn_int_epidem_all_df <- as_data_frame(only_dsn_int_epidem_all)

# TO PLOT

### Epidemiologically known links in COMMON, only DNS, only PPI between the common ICD9s.
common_icds_intersect <- intersect_df[, c("from", "to")]; common_icds_intersect$network <- "both"; dim(common_icds_intersect) # 5
common_icds_only_ppi <- only_ppi_int_epidem_df[, c("from", "to")]; common_icds_only_ppi$network <- "PPI"; dim(common_icds_only_ppi) # 10
common_icds_only_dsn <- only_dsn_int_epidem_df[, c("from", "to")]; common_icds_only_dsn$network <- "DSN"; dim(common_icds_only_dsn) # 29

common_icds_networks <- rbind(common_icds_intersect, common_icds_only_dsn)
common_icds_networks <- rbind(common_icds_networks, common_icds_only_ppi)
dim(common_icds_networks) # 44
write.table(common_icds_networks,file="figures_comparison/comparison_with_ppi_network.txt",sep="\t",row.names=FALSE, quote=FALSE)

dismeta <- read.csv("../dl_general_info_umls_ncolors.csv",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)

# Add the icds that are missing
dismeta[dismeta$disease_name == 'ThyroidCancer_Papillary', ]$icd9 <- '193' 
dismeta[dismeta$disease_name == 'ThyroidCancer_Papillary', ]$icd10 <- 'C73'

# Put three letters where there are only 2 before the point.
dismeta[dismeta$disease_name == 'HIV', ]$icd9 <- '042' 
dismeta[dismeta$disease_name == 'BorreliaBurdorferiInfection', ]$icd9 <- '088' 

dismeta$icd9 <- substr(dismeta$icd9, start = 1, stop = 3)
dismeta$icd10 <- substr(dismeta$icd10, start = 1, stop = 3)

dismeta$icd9 <- as.numeric(dismeta$icd9)

colores<-c("#DD3232","#EC00BC","#FAC6CC","#FFB6AD","#F68B69","#E43E22","#FBAF5F","#FACF63","#FFEF6C","#CED75C",
           "#D2EFDB","#A1CE5E","#969A52","#DBDFC3","#AEC6CC","#A2D0CF","#1A86A8","#00619C","#0065A9","#002B54",
           "#985396","#805462","#391242","#CCC1DB","#E3DFD6","#B2B1A5","#58574B","#C5AB89","#572C29","#000000")
length(colores) # 30

length(unique(dismeta$disease_cat))
df_colors <- data.frame(disease_cat = unique(dismeta$disease_cat), num=c(1:10), colors = colores[1:10])
pp <- ggplot(df_colors, aes(reorder(disease_cat,num), num, fill=disease_cat))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=as.character(df_colors$colors)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pp

df_colors <- data.frame(cat = as.character(c(1:30)), num=c(1:30), colors = as.character(colores))
pp <- ggplot(df_colors, aes(cat, num, fill=cat))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=as.character(df_colors$colors)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pp


selected_colors <- as.character(colores[c(1,4,5,7,12,13,15,18,21,22,26)])
n_selected_colors <- length(selected_colors); n_selected_colors

df_colors <- data.frame(cat = as.character(c(1:n_selected_colors)), num=c(1:n_selected_colors), colors = as.character(selected_colors))
df_colors$colors <- as.character(df_colors$colors)
pp <- ggplot(df_colors, aes(cat, num, fill=cat))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=df_colors$colors) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pp

sorted_disease_cats <- c("diseases of the circulatory system", "neoplasms", 
                         "diseases of the skin and subcutaneous tissue","diseases of the digestive system",
                         "infectious and parasitic diseases", "diseases of the musculoskeletal system and connective tissue",
                         "congenital anomalies", "diseases of the respiratory system", 
                         "diseases of the nervous system and sense organs", "mental disorders")
new_colors <- data.frame(disease_cat = sorted_disease_cats, 
                         new_dis_cat_colors=colores[c(1,4,5,7,12,13,15,18,21,22)])
str(new_colors)
new_colors$disease_cat <- as.character(new_colors$disease_cat)
new_colors$new_dis_cat_colors <- as.character(new_colors$new_dis_cat_colors)

new_dismeta <- merge(dismeta, new_colors, by.x = "disease_cat", by.y = "disease_cat", all.x = T, all.y = F)
write.table(new_dismeta,file="../new_disease_metadata.txt",sep="\t",row.names=FALSE, quote=FALSE)

node_colors_to_save <- new_dismeta[new_dismeta$icd9 %in% common_icds, ]
node_colors_to_save <- node_colors_to_save[, c("icd9", "disease_cat", "disease_name", "new_dis_cat_colors")]
dim(node_colors_to_save)
length(unique(node_colors_to_save$icd9))

node_colors_to_save <- node_colors_to_save[!duplicated(node_colors_to_save$icd9), ]
dim(node_colors_to_save)
write.table(node_colors_to_save,file="figures_comparison/node_colors_ppi.txt",sep="\t",row.names=FALSE, quote=FALSE)

# pie(df_colors$num, labels=df_colors$disease_cat, col = df_colors$colors)

# Generate 10.000 iterations of: select random network with my nº of interactions and compute overlap with Barabasi
n_iters <- 10000 
randomizing_labels <- TRUE
random_overlaps <- c()
random_overlaps_second <- c()
set.seed(5) ### ESTABLECER SEMILLA

# Randomizing labels
if(randomizing_labels == TRUE){
  for(k in 1:n_iters){
    cgraph <- rewire(dsn_common_graph, with=keeping_degseq(loops=FALSE, niter = ecount(dsn_common_graph)*10))
    # str(cgraph)
    cintersect <- intersection(cgraph,ppi_graph)
    # cintersect
    gsize(cintersect)
    random_overlaps <- append(random_overlaps,gsize(cintersect))
    # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
  }
  for(k in 1:n_iters){ # it has to go into another loop to keep consistent p-values (reproducibility)
    # p-value for the comparison with the my network
    cgraph_second <- rewire(ppi_graph, with=keeping_degseq(loops=FALSE, niter = ecount(ppi_graph)*10))
    # str(cgraph)
    cintersect_second <- intersection(cgraph_second,dsn_common_graph)
    # cintersect
    gsize(cintersect_second)
    random_overlaps_second <- append(random_overlaps_second,gsize(cintersect_second))
    # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
  }
}

# Comparison with the ppi | P-value:
pval <- length(which(random_overlaps>=gsize(intersect)))/length(random_overlaps)
pval
# 0.9422625

# Comparison with our network | P-value:
pval_second <- length(which(random_overlaps_second>=gsize(intersect)))/length(random_overlaps_second)
pval_second
# 0.9448

hist(random_overlaps, main = "Random overlaps - pairwise spearman sDEGs")
abline(v = gsize(intersect), col="red")

hist(random_overlaps_second, main = "Random overlaps - pairwise spearman sDEGs")
abline(v = gsize(intersect), col="red")


