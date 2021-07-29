#########################################################################################
# OVERLAP OF THE MICROARRAYS NETWORK WITH THE DSN BASED ON RNA-SEQ
# 
# FINAL VERSION STORED IN ~/Desktop/ANALYSIS
#
# Beatriz Urda García 2020
#########################################################################################

library(igraph)
library(ggplot2)
library(visNetwork)
library(stringr)

# COMMAND-LINE ARGUMENTS
# args = commandArgs(trailingOnly=TRUE)
# input_filename <- args[1]

setwd("~/Desktop/ANALYSIS")


# input_filename <- 'pairwise_union_spearman_distance_sDEGs_network.txt'
input_filename <- 'icd9_pairwise_union_spearman_distance_sDEGs_network.txt'
analyzing_variability <- FALSE
dismeta <- read.csv("dl_general_info_umls_ncolors.csv",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
dismeta$icd9 <- as.character(dismeta$icd9)
# Add the icds that are missing
dismeta[dismeta$disease_name == 'ThyroidCancer_Papillary', ]$icd9 <- '193' 
dismeta[dismeta$disease_name == 'ThyroidCancer_Papillary', ]$icd10 <- 'C73'

# Put three letters where there are only 2 before the point.
dismeta[dismeta$disease_name == 'HIV', ]$icd9 <- '042' 
dismeta[dismeta$disease_name == 'BorreliaBurdorferiInfection', ]$icd9 <- '088' 

dismeta$icd9 <- substr(dismeta$icd9, start = 1, stop = 3)
dismeta$icd10 <- substr(dismeta$icd10, start = 1, stop = 3)

icd_names <- data.frame('disease_name' = dismeta$disease_name,'icd9' = dismeta$icd9)
dups_icds <- unique(as.character(icd_names[which(duplicated(icd_names$icd9)), ]$icd9))
merged_icds <- c()
for(k in dups_icds){
  cdf <- icd_names[icd_names$icd9 == k, ]
  merged_icds <- rbind(merged_icds, c(paste(cdf$disease_name,collapse = "\n"),k))
}
colnames(merged_icds) <- c("disease_name", 'icd9')
icd_names <- icd_names[!(icd_names$icd9 %in% dups_icds), ]
icd_names <- rbind(icd_names, data.frame(merged_icds))

network <- read.csv(paste('Network_building/Defined_networks/',input_filename,sep=""),header=T, sep="\t",stringsAsFactors = F)
if(analyzing_variability == TRUE){
  icds_dis_one_one <- dismeta[, c('disease_name','icd9')]
  dim(icds_dis_one_one)
  new_name <- merge(network, icds_dis_one_one,  by.x = 'Dis1' , by.y= 'disease_name', all.x = TRUE, all.y = FALSE, sort=FALSE); dim(new_name);
  new_name <- merge(new_name, icds_dis_one_one,  by.x = 'Dis2' , by.y= 'disease_name', all.x = TRUE, all.y = FALSE, sort=FALSE); dim(new_name);
  network$Dis1 <- new_name$icd9.x
  network$Dis2 <- new_name$icd9.y
}

# Put three letters where there are only 2 before the point.
network[network$Dis1 == '42',]$Dis1 <- '042' 
network[network$Dis2 == '42',]$Dis2 <- '042' 
network[network$Dis1 == '88',]$Dis1 <- '088' 
network[network$Dis2 == '88',]$Dis2 <- '088'
# Remove interactions where icd1 = icd2 ¡¡NEW!! No hay pero por si!
dim(network) # 545
network <- network[!(network$Dis1 == network$Dis2), ]
dim(network)
head(network)
network$Dis1 <- as.numeric(network$Dis1)
network$Dis2 <- as.numeric(network$Dis2)

microarrays <- read.csv('Network_building/ICD9_RMS_Sanchez_et_al_2020.txt',header=T, sep="\t",stringsAsFactors = F)
microarrays$Disease_1 <- as.numeric(gsub('ICD9-','',microarrays$Disease_1))
microarrays$Disease_2 <- as.numeric(gsub('ICD9-','',microarrays$Disease_2))
colnames(microarrays) <- c('Dis1','Dis2','Distance')
head(microarrays)
dim(microarrays)


icds_rnaseq <- unique(union(network$Dis1,network$Dis2))
icds_microarrays <- unique(union(microarrays$Dis1,microarrays$Dis2))
common_icds <- intersect(icds_rnaseq, icds_microarrays)
length(common_icds) # 27

# Unique in RNA-seq analysisg
only_rnaseq <- setdiff(icds_rnaseq, icds_microarrays)
length(unique(only_rnaseq)) # 14
only_rnaseq
only_microarray <- setdiff(icds_microarrays, icds_rnaseq)
length(unique(only_microarray)) # 65

##### Select network based on icd9
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

network <- sort_icd_interactions(network)
microarrays <- sort_icd_interactions(microarrays)

# Remove duplicated interactions
remove_duplicated_interactions <- function(df, directed=FALSE){
  # df <- microarrays
  concat <- paste(df$Dis1,df$Dis2,sep="_")
  length(concat)             # 241
  length(unique(concat))     # 167
  df$concat <- concat
  if(directed==FALSE){
    df <- df[!duplicated(df$concat),]
    dim(df)
  }else{
    multiple <- df[duplicated(df$concat),]
    combs <- unique(multiple$concat) # combinations
    signs <- rep(NA,length(combs))
    for (k in 1:length(combs)){
      signs[k] <- prod(df[df$concat == combs[k],]$Distance)
    }
    to_remove <- data.frame(combs,signs)
    to_remove <- as.character(to_remove[to_remove$signs < 0, ]$combs)
    length(to_remove) # 19 disease pairs with opposite patterns
    head(to_remove)
    df <- df[!duplicated(df$concat),]
    df <- df[!(df$concat %in% to_remove),]
  }
  df$concat <- NULL
  return(df)
}

dim(network)
network <- remove_duplicated_interactions(network, directed=FALSE)
dim(network)
dim(microarrays)
microarrays <- remove_duplicated_interactions(microarrays, directed=TRUE)
dim(microarrays)

# Saving the entire DSN and the DSN with positive interactions
entire_dsn <- network; dim(network) # 545
entire_dsn_pos <- entire_dsn[entire_dsn$Distance > 0, ]; dim(entire_dsn_pos) # 347

# Filtering both networks with common icds
network <- network[(network$Dis1 %in% common_icds) & (network$Dis2 %in% common_icds), ] # Filtering my network 
microarrays <- microarrays[(microarrays$Dis1 %in% common_icds) & (microarrays$Dis2 %in% common_icds), ] # Filtering my network
dim(network)
dim(microarrays)

# mygraph <- graph_from_data_frame(network, directed=FALSE)
# microarraygraph <- graph_from_data_frame(microarrays, directed=FALSE)
# intersect <- intersection(mygraph,microarraygraph)
# gsize(intersect) # Number of edges: 101
# perc_overlap <- (gsize(intersect)/dim(network)[1])*100 # 40.23904
# perc_overlap
# perc_overlap <- (gsize(intersect)/dim(microarrays)[1])*100 # 74.81481
# perc_overlap

# Remove duplicated interactions
compute_overlap <- function(df1,df2){
  df1graph <- graph_from_data_frame(df1, directed=FALSE)
  df2graph <- graph_from_data_frame(df2, directed=FALSE)
  intersect <- intersection(df1graph,df2graph)
  print(gsize(intersect)) # Number of edges: 101
  perc_overlap1 <- (gsize(intersect)/dim(df1)[1])*100 # 40.23904
  print(paste('% of network 1 interactions: ',perc_overlap1, sep=''))
  perc_overlap2 <- (gsize(intersect)/dim(df2)[1])*100 # 74.81481
  print(paste('% of network 2 interactions: ',perc_overlap2, sep=''))
  return(c(gsize(intersect),perc_overlap1,perc_overlap2))
}

allinteractions <- compute_overlap(network,microarrays)
net_pos <- network[network$Distance >= 0, ] ; 
net_neg <- network[network$Distance < 0, ]   ; 
microarrays_pos <- microarrays[microarrays$Distance >= 0, ] ; 
microarrays_neg <- microarrays[microarrays$Distance < 0, ]
positive_interactions <- compute_overlap(net_pos,microarrays_pos) # 27.439% & 65.217%
negative_interactions <- compute_overlap(net_neg,microarrays_neg) # 20.69% & 27.69%
all_with_direction <- ((positive_interactions[1]+negative_interactions[1])/dim(microarrays)[1])*100
positive_interactions[1]+negative_interactions[1]
all_with_direction

# DEFINE GRAPHS
mygraph_entire_dsn_pos <- graph_from_data_frame(entire_dsn_pos, directed=FALSE)
mygraph_pos <- graph_from_data_frame(net_pos, directed=FALSE)
mygraph_neg <- graph_from_data_frame(net_neg, directed=FALSE)
graph_microarrays_pos <- graph_from_data_frame(microarrays_pos, directed=FALSE)
graph_microarrays_neg <- graph_from_data_frame(microarrays_neg, directed=FALSE)

#### SIGNIFICANCIA DEL SOLAPAMIENTO
# Generate 10.000 iterations of: select random network with my nº of interactions and compute overlap with Barabasi
n_iters <- 10000
random_overlaps_pos <- c()
random_overlaps_neg <- c()
random_overlaps_all <- c()
set.seed(5) ### ESTABLECER SEMILLA
for(k in 1:n_iters){
  cgraph_pos <- rewire(mygraph_pos, with=keeping_degseq(loops=FALSE, niter = ecount(mygraph_pos)*10))
  cgraph_neg <- rewire(mygraph_neg, with=keeping_degseq(loops=FALSE, niter = ecount(mygraph_neg)*10))
  # str(cgraph)
  cintersect_pos <- intersection(cgraph_pos,graph_microarrays_pos)
  cintersect_neg <- intersection(cgraph_neg,graph_microarrays_neg)
  # cintersect
  # gsize(cintersect)
  random_overlaps_pos <- append(random_overlaps_pos,gsize(cintersect_pos))
  random_overlaps_neg <- append(random_overlaps_neg,gsize(cintersect_neg))
  # random_overlaps_all <- append(random_overlaps_all,gsize(cintersect_pos)+gsize(cintersect_neg))
  # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
}

pval_pos <- length(which(random_overlaps_pos>=positive_interactions[1]))/length(random_overlaps_pos); pval_pos  # 0.002
pval_neg <- length(which(random_overlaps_neg>=negative_interactions[1]))/length(random_overlaps_neg); pval_neg  # 0.6243

random_overlaps_all <- apply(data.frame(random_overlaps_pos,random_overlaps_neg), 1, sum)
pval_overall <- length(which(random_overlaps_all>=(positive_interactions[1]+negative_interactions[1])))/length(random_overlaps_all); pval_overall # 0.0272


plot_network <- function(graph, df=NULL){
  # Visualize the communities
  nodes <- data.frame(id = V(graph)$name, title = V(graph)$name)
  nodes <- nodes[order(nodes$id, decreasing = F),]
  edges <- get.data.frame(graph, what="edges")[1:2]
  # edges$color <- rep("grey",length(edges$from))
  if (!is.null(df)){
    edges$value <- df$Distance
  }
  visNetwork(nodes, edges) %>%
    visExport(type = "pdf") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%   # list(enabled=TRUE, selected="BreastCancer")
    visIgraphLayout() %>%
    visInteraction(multiselect = T) %>%
    visEvents(select = "function(nodes) {
            Shiny.onInputChange('current_node_id', nodes.nodes);
            ;}")# edges$value <- df$Distance
}

network_in_disease_names <- function(graph){
  # graph <- common_transcriptomics
  # graph <- solo_epidem
  cvnames <- V(graph)$name; cvnames <- data.frame('icd9'=cvnames); dim(cvnames)
  # to_rename_v <- icd_names[icd_names$icd9 %in% cvnames, ]; 
  to_rename_v <- merge(cvnames,icd_names, by='icd9', all.x = TRUE, all.y = FALSE, sort=FALSE); dim(to_rename_v);
  # to_rename_v <- to_rename_v[!duplicated(to_rename_v$icd9), ]; dim(to_rename_v)
  head(cvnames); head(to_rename_v); 
  # V(common_transcriptomics)$name <- to_rename_v$disease_name
  graph <- set.vertex.attribute(graph, "name", value=as.character(to_rename_v$disease_name))
  return(graph)
}

### PLOTTING GRAPHS
plot_network(mygraph_pos,net_pos)

## FILTRAR RED DE BARABASI PARA PODER COMPARARLA CON LA MIA
# Read and filter barabasi network
barabasi <- read.csv("Network_building/PDN_3_digits.net",header=F, sep="\t",stringsAsFactors = F)
colnames(barabasi) <- c('Dis1','Dis2','Prev1','Prev2','Co-ocurrence','RR','RR_99perc_left_bound','RR_99perc_right_bound','phi_corr','p_value')
barabasi <- barabasi[,c('Dis1','Dis2','RR_99perc_left_bound')]
barabasi <- barabasi[which(barabasi$RR_99perc_left_bound > 1 ), ]
dim(barabasi)
barabasi <- sort_icd_interactions(barabasi)
barabasi  <- remove_duplicated_interactions(barabasi , directed=FALSE)
dim(barabasi)
# barabasi$Dis1 <- str_pad(as.character(barabasi$Dis1), 3, pad = "0")
# barabasi$Dis2 <- str_pad(as.character(barabasi$Dis2), 3, pad = "0")
barabasi_graph <- graph_from_data_frame(barabasi, directed=FALSE)

#### DE LO COMÚN EN AMBAS REDES, QUÉ SOLAPA Y QUÉ NO CON LA EPIDEMIOLOGÍA
common_transcriptomics_pos <- intersection(mygraph_pos, graph_microarrays_pos, byname = TRUE, keep.all.vertices = FALSE)
# E(common_transcriptomics_pos)$color <- 'red'
common_transcriptomics_neg <- intersection(mygraph_neg, graph_microarrays_neg, byname = TRUE, keep.all.vertices = FALSE)
# E(common_transcriptomics_neg)$color <- 'blue'
common_transcriptomics <- union(common_transcriptomics_pos, common_transcriptomics_neg)
gsize(common_transcriptomics) # 63
plot(common_transcriptomics) # intersection
plot_network(common_transcriptomics)
plot_network(network_in_disease_names(common_transcriptomics))
plot_network(network_in_disease_names(common_transcriptomics_neg))
write_graph(common_transcriptomics_neg, file=paste('Network_building/Microarrays_comparison/','common_transcriptomics_neg',".txt",sep=""), format = "ncol")
# adjacent_vertices(common_transcriptomics,'332')

# Solapa
solapa <- intersection(common_transcriptomics_pos, barabasi_graph, byname = TRUE, keep.all.vertices = FALSE)
gsize(solapa) # 13
plot(solapa)
plot_network(solapa)
plot_network(network_in_disease_names(solapa))
write_graph(common_transcriptomics_neg, file=paste('Network_building/Microarrays_comparison/','common_transcriptomics_pos_in_epidemiology',".txt",sep=""), format = "ncol")

# No Solapa
solo_common <- difference(common_transcriptomics_pos,barabasi_graph)
plot(solo_common)
plot_network(network_in_disease_names(solo_common))
solo_epidem <- difference(barabasi_graph,common_transcriptomics_pos)
plot(solo_epidem)
plot_network(solo_epidem)
plot_network(network_in_disease_names(solo_epidem))
write_graph(solo_common, file=paste('Network_building/Microarrays_comparison/','common_transcriptomics_not_in_epidemiology',".txt",sep=""), format = "ncol")

### DE LO NUEVO, QUÉ SOLAPA EN LA EPIDEMIOLOGÍA
solo_rnaseq_pos <- difference(mygraph_pos, graph_microarrays_pos)
solo_rnaseq_neg <- difference(mygraph_neg, graph_microarrays_neg)

solo_rnaseq_epidem <- intersection(solo_rnaseq_pos, barabasi_graph, byname = TRUE, keep.all.vertices = FALSE)
gsize(solo_rnaseq_epidem) # 42
plot_network(solo_rnaseq_epidem)
plot_network(network_in_disease_names(solo_rnaseq_epidem))
plot_network(network_in_disease_names(solo_rnaseq_neg))
write_graph(solo_rnaseq_epidem, file=paste('Network_building/Microarrays_comparison/','solo_rnaseq_in_epidemiology',".txt",sep=""), format = "ncol")
write_graph(solo_rnaseq_neg, file=paste('Network_building/Microarrays_comparison/','solo_rnaseq_neg',".txt",sep=""), format = "ncol")

visExport(solo_rnaseq_epidem,type = "pdf",name = "solo_rna_seq_pos_known",
  label = paste0("Export as ", type),
  background = "#fff",
  float = "right"
)

export_pdf_network(solo_rnaseq_epidem, df=NULL, network_name="solo_rnaseq_epidem")

# COMPARANDO LA DSN ENTERA CON LA DE MICROARRAYS: TOTAL NEW KNOWN INTERACTIONS 
entire_dsn_minus_microarray <- difference(mygraph_entire_dsn_pos, graph_microarrays_pos)
entire_dsn_minus_microarray_in_epidem <- intersection(entire_dsn_minus_microarray, barabasi_graph, byname = TRUE, keep.all.vertices = FALSE)
plot_network(entire_dsn_minus_microarray_in_epidem)
gsize(entire_dsn_minus_microarray_in_epidem) # 141

solo_microarrays <- difference(graph_microarrays_pos, mygraph_entire_dsn_pos)
solo_microarrays_in_epidem <- intersection(solo_microarrays, barabasi_graph, byname = TRUE, keep.all.vertices = FALSE)
gsize(solo_microarrays_in_epidem) # 7

# TO PLOT

### Epidemiologically known links in COMMON, only DNS, only PPI between the common ICD9s.
solo_microarray_epidem <- intersection(difference(graph_microarrays_pos, mygraph_pos), barabasi_graph,byname = TRUE, keep.all.vertices = FALSE)
plot(solo_microarray_epidem)
intersect_df <- as_data_frame(solapa)
only_microarray_int_epidem_df <- as_data_frame(solo_microarray_epidem)
only_dsn_int_epidem_df <- as_data_frame(solo_rnaseq_epidem)

common_icds_intersect <- intersect_df[, c("from", "to")]; common_icds_intersect$network <- "both"; dim(common_icds_intersect) # 13
only_dsn_int_epidem_df <- only_dsn_int_epidem_df[, c("from", "to")]; only_dsn_int_epidem_df$network <- "DSN"; dim(only_dsn_int_epidem_df) # 42
only_microarray_int_epidem_df <- only_microarray_int_epidem_df[, c("from", "to")]; only_microarray_int_epidem_df$network <- "Microarrays"; dim(only_microarray_int_epidem_df) # 7

common_icds_networks <- rbind(common_icds_intersect, only_dsn_int_epidem_df)
common_icds_networks <- rbind(common_icds_networks, only_microarray_int_epidem_df)
dim(common_icds_networks) # 62
length(unique(union(common_icds_networks$from, common_icds_networks$to))) # 25
length(intersect(unique(union(common_icds_networks$from, common_icds_networks$to)), common_icds)) # 25
setdiff(common_icds, unique(union(common_icds_networks$from, common_icds_networks$to)))

write.table(common_icds_networks,file="Other_molecular_networks/figures_comparison/comparison_with_microarray_network.txt",sep="\t",row.names=FALSE, quote=FALSE)

#  Create table with node colors
new_dismeta <- read.csv("new_disease_metadata.txt",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)

node_colors_to_save <- new_dismeta[new_dismeta$icd9 %in% common_icds, ]
node_colors_to_save <- node_colors_to_save[, c("icd9", "disease_cat", "disease_name", "new_dis_cat_colors")]
dim(node_colors_to_save)
length(unique(node_colors_to_save$icd9))

node_colors_to_save <- node_colors_to_save[!duplicated(node_colors_to_save$icd9), ]
dim(node_colors_to_save)
write.table(node_colors_to_save,file="Other_molecular_networks/figures_comparison/node_colors_microarrays.txt",sep="\t",row.names=FALSE, quote=FALSE)



### PLOT DSN with variability new interactions
new_in_variability <- read.csv('Disease_gene_variability/Newly_detected_epidemiological_interactions_with_sdvgs.txt',header=T, sep="\t",stringsAsFactors = F)
variab_graph <- graph_from_data_frame(new_in_variability, directed=FALSE)
plot(variab_graph)
plot_network(variab_graph)
