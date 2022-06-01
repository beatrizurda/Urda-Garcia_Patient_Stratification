#########################################################################################
# CHECK THE NUMBER OF INTERACTIONS THAT ARE ALSO CAPTURED WITH OTHER MOLECULAR LAYERS
#    - EI AND NEIS
#
# Beatriz Urda García 2022
#########################################################################################

setwd("~/Desktop/ANALYSIS/Other_molecular_networks/")

library(igraph)
library(VennDiagram)
library(ggplot2)
source("../analysis_library.R")


dsn_path <- '~/Desktop/Desktop_20200508/Network_building/Defined_networks/icd9_pairwise_union_spearman_distance_sDEGs_pos_network.txt'
dsn <- read.csv(dsn_path,header=T, sep="\t",stringsAsFactors = F)
my_icd9 <- unique(union(dsn$Dis1, dsn$Dis2)); length(my_icd9) # 41

# Read and filter barabasi network
barabasi <- read.csv("~/Desktop/Desktop_20200508/Network_building/Overlapping_results/Shuffling_labels/icd9_pairwise_union_spearman_distance_sDEGs_pos_B_TRUE_final_barabasi.txt",
                     header=T, sep="\t",stringsAsFactors = F)
dim(barabasi) # 331
barabasi_graph <- graph_from_data_frame(barabasi, directed=FALSE)
plot(barabasi_graph)

other_molecular_networks = c('Transformed_networks/Nodups_icd9_ppi_network.txt', 'Transformed_networks/Nodups_icd9_cellular_components_network_genes.txt',
                             'Transformed_networks/Nodups_icd9_cellular_components_network_ppis.txt', 'Transformed_networks/Nodups_icd9_microbiome_network.txt',
                             'Transformed_networks/Nodups_icd9_microRNAs_network.txt',
                             'Transformed_networks/Nodups_icd9_geo_gsk_pos_union_network.txt', 
                             "Transformed_networks/Nodups_icd9_microarrays_network_pos.txt",
                             'Transformed_networks/Nodups_icd9_subcellular_localization_network.txt')
names(other_molecular_networks) = c("PPIs_Menche", "Cellular_componets_genes", "Cellular_components_PPIs", "Microbiome", "MicroRNAs", 
                                    "Union_gene_expression", "Microarrays", "Subcellular_localization")


robustness_other_networks = function(other_path, only_epidemiology = FALSE){
  # Read the other molecular network
  # other_path = 'Transformed_networks/Nodups_icd9_ppi_network.txt'
  othernetwork <- read.csv(other_path,header=T, sep="\t",stringsAsFactors = F)
  other_icd9 <- unique(union(othernetwork$Dis1, othernetwork$Dis2)); length(other_icd9) # 19
  dim(othernetwork)
  
  # Take the molecular networks with common icd9s
  common_icds <- unique(intersect(my_icd9, other_icd9)); length(common_icds) # 19
  
  dsn_common <- dsn[((dsn$Dis1 %in% common_icds) & (dsn$Dis2 %in% common_icds)), ]
  nrow(dsn_common)
  othernetwork_common = othernetwork[((othernetwork$Dis1 %in% common_icds) & (othernetwork$Dis2 %in% common_icds)), ] # just in case
  nrow(othernetwork_common)
  
  dsn_graph <- graph_from_data_frame(dsn_common, directed=FALSE)
  other_graph <- graph_from_data_frame(othernetwork, directed=FALSE)
  plot(dsn_graph); gsize(dsn_graph) # 20
  plot(other_graph); gsize(other_graph) # 20
  
  if(only_epidemiology == TRUE){
    dsn_graph = intersection(dsn_graph, barabasi_graph)
    other_graph = intersection(other_graph, barabasi_graph)
    size_filtered_dsn = gsize(dsn_graph)
    size_filtered_other_network = gsize(other_graph)
  }
  
  intersect <- intersection(dsn_graph,other_graph)
  plot(intersect)
  print(gsize(intersect)) # 12
  intersect_in_epidem <- intersection(intersect, barabasi_graph); gsize(intersect_in_epidem) # 5
  plot(intersect_in_epidem)
  intersect_df <- as_data_frame(intersect_in_epidem)
  recall <- (gsize(intersect)/dim(dsn_common)[1])*100 # 8.22%
  print(paste('RECALL (% of interactions from the DSN): ',recall, sep=''))
  precision <- (gsize(intersect)/dim(othernetwork_common)[1])*100 # 30%
  print(paste('PRECISION (% of network from the other molecular network): ',precision, sep=''))
  
  # Generate 10.000 iterations of: select random network with my nº of interactions and compute overlap with Barabasi
  n_iters <- 10000 
  randomizing_labels <- TRUE
  random_overlaps <- c()
  random_overlaps_second <- c()
  set.seed(5) ### ESTABLECER SEMILLA
  
  # Randomizing labels
  if(randomizing_labels == TRUE){
    for(k in 1:n_iters){
      cgraph <- rewire(dsn_graph, with=keeping_degseq(loops=FALSE, niter = ecount(dsn_graph)*10))
      # str(cgraph)
      cintersect <- intersection(cgraph,other_graph)
      # cintersect
      gsize(cintersect)
      random_overlaps <- append(random_overlaps,gsize(cintersect))
      # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
    }
    for(k in 1:n_iters){ # it has to go into another loop to keep consistent p-values (reproducibility)
      # p-value for the comparison with the my network
      cgraph_second <- rewire(other_graph, with=keeping_degseq(loops=FALSE, niter = ecount(other_graph)*10))
      # str(cgraph)
      cintersect_second <- intersection(cgraph_second,dsn_graph)
      # cintersect
      gsize(cintersect_second)
      random_overlaps_second <- append(random_overlaps_second,gsize(cintersect_second))
      # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
    }
  }
  
  # PRECISION | P-value:
  pval <- length(which(random_overlaps>=gsize(intersect)))/length(random_overlaps)
  pval
  
  # RECALL | P-value:
  pval_second <- length(which(random_overlaps_second>=gsize(intersect)))/length(random_overlaps_second)
  pval_second
  
  if(only_epidemiology == TRUE){
    output = c(nrow(othernetwork), 
               size_filtered_dsn, size_filtered_other_network,
               length(common_icds), nrow(dsn_common), nrow(othernetwork_common), gsize(intersect),
               recall, precision, pval_second, pval)
  }else{
    output = c(nrow(othernetwork), length(common_icds), nrow(dsn_common), nrow(othernetwork_common), gsize(intersect),
               recall, precision, pval_second, pval)
  }
  
}

# ALL
outputdf = data.frame(matrix(NA, nrow=0, ncol=9))

for(k in 1:length(other_molecular_networks)){
  # k = 1
  cname = names(other_molecular_networks)[k]
  cnetwork_path = other_molecular_networks[k]
  crow = robustness_other_networks(cnetwork_path)
  outputdf = rbind(outputdf, crow)
}

outputdf$network = names(other_molecular_networks)  
colnames(outputdf) = c("Dim_network", "Common_icds", "Dim_dsn_common_icds", "Dim_net_common_icds", "Overlap", "Recall", "Precision", "Pval_recall", "Pval_precision", "Mol_network")
write.table(outputdf, "robustness_other_molecular_layers.txt", quote = FALSE, sep="\t", row.names = FALSE)

# ONLY IN THE EPIDEMIOLOGY
outputdf2 = data.frame(matrix(NA, nrow=0, ncol=11))

for(k in 1:length(other_molecular_networks)){
  # k = 1 
  cname = names(other_molecular_networks)[k]
  cnetwork_path = other_molecular_networks[k]
  crow = robustness_other_networks(cnetwork_path, only_epidemiology = TRUE)
  outputdf2 = rbind(outputdf2, crow)
}

outputdf2$network = names(other_molecular_networks)  
colnames(outputdf2) = c("Dim_network", "Dim_epidem_DSN","Dim_epidem_network", "Common_icds", "Dim_dsn_common_icds", "Dim_net_common_icds", 
                       "Overlap", "Recall", "Precision", "Pval_recall", "Pval_precision", "Mol_network")
write.table(outputdf2, "robustness_other_molecular_layers_in_epidem.txt", quote = FALSE, sep="\t", row.names = FALSE)

