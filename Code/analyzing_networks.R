#########################################################################################
# ANALYZE NETWORK'S TOPOLOGY
# 
# 
#
# Beatriz Urda García 2020
#########################################################################################

library(igraph)
library(ggplot2)
library(gplots)
library(ggpubr)

setwd("~/Desktop/ANALYSIS/")
metadata <- read.csv("dl_general_info_umls_new_ncolors.csv",header=T, sep=",",stringsAsFactors = F); dim(metadata)
new_meta <- read.csv("new_disease_metadata_final_names.txt",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
names_correspondences <- new_meta[, c("disease_name","final_disease_name")]

# RUN ANALYSIS FOR SSN
network_file <- "metapatients_and_disease/metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt"
output_filename <- "SSN"
obtain_cliques <- FALSE; cex_axis <- 4

# RUN ANALYSIS FOR DSN
network_file <- "pairwise_union_spearman_distance_sDEGs_network.txt"
output_filename <- "DSN"
obtain_cliques <- TRUE; cex_axis <- 6

# RUN ANALYSIS FOR MICROARRAY'S DATA
barabasi <- read.csv("Network_building/PDN_3_digits.net",header=TRUE, sep="\t",stringsAsFactors = F)
colnames(barabasi) <- c('Dis1','Dis2','Prev1','Prev2','Co-ocurrence','RR','RR_99perc_left_bound','RR_99perc_right_bound','phi_corr','p_value')
barabasi <- barabasi[,c('Dis1','Dis2','RR_99perc_left_bound')]
barabasi <- barabasi[which(barabasi$RR_99perc_left_bound > 1 ), ]
colnames(barabasi) <- c('Dis1','Dis2','Distance')
barabasi$Dis1 <- as.character(barabasi$Dis1); barabasi$Dis2 <- as.character(barabasi$Dis2)
df <- barabasi
output_filename <- "Barabasi"
obtain_cliques <- FALSE; cex_axis <- 4

# RUN ANALYSIS FOR MICROARRAY'S DATA
# network_file <- "../ICD9_RMS_Sanchez_et_al_2020.txt"
# output_filename <- "Sanchez_et_al"
# obtain_cliques <- FALSE; cex_axis <- 4



net_filepath <- paste("Network_building/Defined_networks/",network_file, sep="")
df <- read.csv(net_filepath,header=T, sep="\t",stringsAsFactors = F); dim(df)

## If network is the DSN, change names to final names
if(output_filename == "DSN"){
  newdf <- merge(df, names_correspondences, by.x = "Dis1", by.y = "disease_name", all.x = TRUE, all.y = FALSE)
  newdf <- merge(newdf, names_correspondences, by.x = "Dis2", by.y = "disease_name", all.x = TRUE, all.y = FALSE)
  df <- newdf[, c("final_disease_name.x", "final_disease_name.y", "Distance", "pvalue", "adj_pvalue")]
  colnames(df)[1:2] <- c("Dis1", "Dis2")
}


dfpos <- df[df$Distance >= 0, ]
dfneg <- df[df$Distance < 0, ]; dfneg$Distance <- abs(dfneg$Distance)
graph_list <- list(df,dfpos,dfneg)


# Does the degree correlate with the number of sDEGs or Sample size
if(output_filename == "DSN"){
  pdf(paste("Network_building/Network_analysis/Correlations/",output_filename,"_correlation_plots.pdf",sep=""),width = 6, height = 5) # 9 & 7
  graph <- graph_from_data_frame(df, directed=FALSE)
  deg <- degree(graph, mode="all")
  df_i_corr <- metadata[, c("disease_name", "n_samples", "Number_sDEGs")]
  deg_df <- data.frame("disease_name"=names(deg), "deg"=deg); rownames(deg_df)<- NULL; 
  df_i_corr <- merge(deg_df, df_i_corr, by= "disease_name", all.x=TRUE, all.y=FALSE)
  trans_corp <- cor.test(df_i_corr$deg, df_i_corr$n_samples, method="pearson"); trans_corp # YES
  trans_corp <- cor.test(df_i_corr$deg, df_i_corr$Number_sDEGs, method="pearson"); trans_corp # YES
  
  ## Correlation DEGREE & Number of samples
  p1 <- ggplot(df_i_corr,aes(n_samples,deg))
  p1 <- p1 + geom_point() + geom_smooth(method=lm, se=FALSE) +
    theme_classic()+
    xlab("Sample size") + ylab("Node degree in the DSN") +
    theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
  print(p1)
  
  # remove outliers
  c_df_i_corr <- df_i_corr[df_i_corr$n_samples < 300, ]
  p1 <- ggplot(c_df_i_corr,aes(n_samples,deg))
  p1 <- p1 + geom_point() + geom_smooth(method=lm, se=FALSE) +
    theme_classic()+
    xlab("Sample size") + ylab("Node degree in the DSN") +
    theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
  print(p1)
  
  ## Correlation DEGREE & #sDEGs
  p1 <- ggplot(df_i_corr,aes(Number_sDEGs,deg))
  p1 <- p1 + geom_point() + geom_smooth(method=lm, se=FALSE) +
    theme_classic()+
    xlab("Number of sDEGs") + ylab("Node degree in the DSN") +
    theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
  print(p1)
  
  # remove outliers
  c_df_i_corr <- df_i_corr[df_i_corr$Number_sDEGs > 0, ]
  p1 <- ggplot(c_df_i_corr,aes(Number_sDEGs,deg))
  p1 <- p1 + geom_point() + geom_smooth(method=lm, se=FALSE) +
    theme_classic()+
    xlab("Number of sDEGs") + ylab("Node degree in the DSN") +
    theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
  print(p1)
  
  dev.off()
  
}



# RUN THE TOPOLOGICAL ANALYSIS FOR THE NETWORK (k=1), POS (k=2) AND NEG SUBNETWORK (k=3)
k <- 1
for(df in graph_list){
  # head(df)
  # class(df)
  pdf(paste("Network_building/Network_analysis/Topology/",output_filename,k,"_topology_plots.pdf",sep=""))
  graph <- graph_from_data_frame(df, directed=FALSE)
  E(graph)$weight <- df$Distance
  
  connected_components <- clusters(graph)
  n_connected_components <- connected_components$no
  size_connected_components <- connected_components$csize
  if(length(size_connected_components) > 1){
    size_connected_components <- paste(size_connected_components,collapse="_")
  }
  
  deg <- degree(graph, mode="all")
  mean_deg <- mean(deg) 
  
  #Density of a graph: The proportion of present edges from all possible edges in the network
  edge_density <- edge_density(graph, loops=F)
  edge_density
  
  # Reciprocity: it only makes sense for directed graphs
  reciprocity(graph)
  
  # dyad_census(graph) # Mutual, asymmetric, and nyll node pairs
  # 2*dyad_census(graph)$mut/ecount(graph) # Calculating reciprocity
  
  # Transitivity: Transitivity measures the probability that the adjacent vertices of a vertex are connected
  transitivity <- transitivity(graph, type="global")
  local_trans <- transitivity(graph, vids=V(graph), type="local")
  transitivity_table <- data.frame(nodes=names(V(graph)), transitivity=local_trans)
  print(ggplot(transitivity_table, aes(reorder(nodes,transitivity), transitivity)) + geom_bar(stat="identity") +
          theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis),plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
    # ggtitle("Transitivity") + 
    # xlab("Diseases") + 
    ylab("Transitivity"))
  
  ### Diameter: longest geodesic distance (length of the shortest path between two nodes) in the network
  # With weights
  diameter <- diameter(graph, directed=F, weights=abs(E(graph)$weight))
  get_diameter(graph, directed=F, weights=abs(E(graph)$weight))
  # Without weights
  diameter_wo_weights <- diameter(graph, directed=F, weights=NA)
  get_diameter(graph, directed=F, weights=NA)
  
  ###Degree
  deg <- degree(graph, mode="all")
  plot(graph, vertex.size=deg*0.4,layout=layout.circle)
  hist(deg, breaks=1:vcount(graph)-1, main="Histogram of node degree")
  
  # Degree distribution
  deg.dist <- degree_distribution(graph, cumulative=T, mode="all")
  
  plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
        xlab="Degree", ylab="Cumulative Frequency")
  
  ### CENTRALITY AND CENTRALIZATION
  # Degree (number of ties) ### hacer un plot de esto con ggplot quizá!
  degree_table <- data.frame('nodes'=names(degree(graph, mode="in")), 'feature'=degree(graph, mode="in"))
  centr_degree(graph, mode="in", normalized=T)
  print(ggplot(degree_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
    # ggtitle("Transitivity") + 
    # xlab("Diseases") + 
    xlab("Diseases") + ylab("Node degree"))
  
  # Closeness (centrality based on distance to others in the graph)
  #Inverse of the node’s average geodesic distance to others in the network.
  closeness(graph, mode="all", weights=NA) 
  centr_clo(graph, mode="all", normalized=T)
  degree_table <- data.frame('nodes'=names(closeness(graph, mode="all", weights=NA)), 'feature'=closeness(graph, mode="all", weights=NA))
  print(ggplot(degree_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
    # ggtitle("Transitivity") + 
    # xlab("Diseases") + 
      xlab("Diseases") + ylab("Closeness"))
  
  # Eigenvector (centrality proportional to the sum of connection centralities)
  # Values of the first eigenvector of the graph matrix.
  
  eigen_centrality(graph, directed=T, weights=NA)
  centr_eigen(graph, directed=T, normalized=T) 
  
  # Betweenness (centrality based on a broker position connecting others)
  # Number of geodesics that pass through the node or the edge.
  betweenness <- betweenness(graph, directed=FALSE, weights=abs(E(graph)$weight))
  betweenness_wo_weights <- betweenness(graph, directed=FALSE, weights=NA)
  edge_betweenness(graph, directed=T, weights=NA)
  centr_betw(graph, directed=T, normalized=T)
  betweenness <- data.frame('nodes'=names(betweenness), 'feature'=betweenness)
  print(ggplot(betweenness, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
    ggtitle("Betweenness with weights") +
    # xlab("Diseases") + 
      xlab("Diseases") + ylab("Betweenness"))
  
  betweenness <- data.frame('nodes'=names(betweenness_wo_weights), 'feature'=betweenness_wo_weights)
  print(ggplot(betweenness, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
    ggtitle("Betweenness wo weights") +
    # xlab("Diseases") + 
      xlab("Diseases") + ylab("Betweenness"))
  
  ### HUBS AND AUTHORITIES
  # HUBS with weights
  hs <- hub_score(graph, weights=abs(E(graph)$weight))$vector
  hubs_table <- data.frame('nodes'=names(hs), 'feature'=hs)
  print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
    ggtitle("Hubs with weights") +
    # xlab("Diseases") + 
      xlab("Diseases") + ylab("Hubs"))
  # wo weights
  hs <- hub_score(graph, weights=NA)$vector
  hubs_table <- data.frame('nodes'=names(hs), 'feature'=hs)
  print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
    ggtitle("Hubs wo weights") +
    # xlab("Diseases") + 
      xlab("Diseases") + ylab("Hubs"))
  
  # AUTHORITIES with weights
  as <- authority_score(graph, weights=abs(E(graph)$weight))$vector
  auth_table <- data.frame('nodes'=names(as), 'feature'=as)
  print(ggplot(auth_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
    ggtitle("Authorities with weights") +
    # xlab("Diseases") + 
      xlab("Diseases") + ylab("Authority"))
  
  # wo weights
  as <- authority_score(graph, weights=NA)$vector
  auth_table <- data.frame('nodes'=names(as), 'feature'=as)
  print(ggplot(auth_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
    ggtitle("Authorities wo weights") +
    # xlab("Diseases") + 
      xlab("Diseases") + ylab("Authority"))
  
  par(mfrow=c(1,2))
  plot(graph, vertex.size=hs*20, main="Hubs")
  plot(graph, vertex.size=as*10, main="Authorities")
  
  ### DISTANCES AND PATHS
  # Average path length: the mean of the shortest distance between 
  #each pair of nodes in the network (in both directions for directed graphs).
  mean_distance <- mean_distance(graph, directed=F)
  # distances(graph) # with edge weights # YIELDS AND ERROR
  distances(graph, weights=NA) # ignore weights
  
  ### FINDING COMMUNITIES
  # cliques: complete subgraphs of an undirected graph
  if(obtain_cliques == TRUE){
    graph.sym <- as.undirected(graph, mode= "collapse",
                               edge.attr.comb=list(weight="sum", "ignore"))
    cliques(graph) # list of cliques       
    sapply(cliques(graph), length) # clique sizes
    saveRDS(largest_cliques(graph), file=paste("Network_building/Network_analysis/Topology/Cliques/",output_filename,k,'.rds', sep=""))
  }
  
  
  # # Community detection algorithms
  # if (input$comm_algorithm == 'greedy'){
  #   # Using greedy optimization of modularity
  #   fc <- fastgreedy.community(graph)
  #   V(graph)$community <- fc$membership
  # }else if(input$comm_algorithm == 'rand_walks'){
  #   # Using random walks
  #   fc <- cluster_walktrap(graph)
  #   V(graph)$community <- fc$membership #membership(fc)
  # }
  
  # K-core decomposition
  # The k-core is the maximal subgraph in which every node has degree of 
  # at least k. The result here gives the coreness of each vertex in the 
  # network. A node has coreness D if it belongs to a D-core but not 
  # to (D+1)-core.
  kc <- coreness(graph, mode="all")
  # plot(graph, vertex.size=kc*6, vertex.label=kc, vertex.color=colrs[kc]) # FIXX
  kc_table <- data.frame('nodes'=names(kc), 'feature'=kc)
  print(ggplot(kc_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
          theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))+
    # ggtitle("Authorities with weights") +
    # xlab("Diseases") + 
      xlab("Diseases") + ylab("K-core decomposition"))
  
  ### ASSORTATIVITY AND HOMOPHILY
  # assortativity_nominal(graph, V(graph)$media.type, directed=F)
  assortativity_degree <- assortativity_degree(graph, directed=F)
  
  cfeat <- c('N conn comp','Size conn comp','Mean deg',
             'edge_density','transitivity','diameter',
             'diameter_wo_weights','mean_distance','assortativity_degree') 
  cvals <- c(n_connected_components,size_connected_components,mean_deg,
             edge_density,transitivity,diameter,
             diameter_wo_weights,mean_distance,assortativity_degree)
  cdf <- data.frame(Feature=cfeat,Value = cvals)
  
  dev.off()
  write.table(cdf,file=paste("Network_building/Network_analysis/Topology/",output_filename,k,"_topology.txt",sep=""),sep="\t",row.names=F, quote=FALSE)
  k <- k+1
}

##### ANALYZE CLIQUES #####
cliq1 <- readRDS(file="Network_building/Network_analysis/Topology/Cliques/DSN1.rds")
cliq1

cliq2 <- readRDS(file="Network_building/Network_analysis/Topology/Cliques/DSN2.rds")
cliq2

cliq3 <- readRDS(file="Network_building/Network_analysis/Topology/Cliques/DSN3.rds")
cliq3

##### CORRELACIONES ENTRE LAS REDES COMPARABLES ######
filename_ours <- "../Desktop_20200508/Network_building/Overlapping_results/Shuffling_labels/icd9_pairwise_union_spearman_distance_sDEGs_pos_B_TRUE_final_network.txt"
filename_barabasi <- "../Desktop_20200508/Network_building/Overlapping_results/Shuffling_labels/icd9_pairwise_union_spearman_distance_sDEGs_pos_B_TRUE_final_barabasi.txt"
output_filename <- "B_all_diseases"
# output_filename <- "C_filtering_diseases"
cex_axis <- 6

our_network <- read.csv(filename_ours,header=T, sep="\t",stringsAsFactors = F); dim(our_network) # 347
barabasi_network <- read.csv(filename_barabasi,header=T, sep="\t",stringsAsFactors = F); dim(barabasi_network) # 331

str(our_network)
str(barabasi_network)

my_icds <- unique(union(our_network$Dis1, our_network$Dis2)); length(my_icds) # 41
barab_icds <- unique(union(barabasi_network$Dis1, barabasi_network$Dis2)); length(barab_icds) # 40
common_icds <- intersect(my_icds, barab_icds)
length(common_icds); setdiff(my_icds, barab_icds) # 174 Cancer de Mama

our_network <- our_network[(our_network$Dis1 %in% common_icds) & (our_network$Dis2 %in% common_icds), ]; dim(our_network) # 327
barabasi_network <- barabasi_network[(barabasi_network$Dis1 %in% common_icds) & (barabasi_network$Dis2 %in% common_icds), ]; dim(barabasi_network) # 331 --> 329

our_graph <- graph_from_data_frame(our_network, directed=FALSE)
E(our_graph)$weight <- our_network$Distance

barabasi_graph <- graph_from_data_frame(barabasi_network, directed=FALSE)
E(barabasi_graph)$weight <- barabasi_network$RR_99perc_left_bound

# SIN WEIGHT
# our_graph <- graph_from_data_frame(our_network, directed=FALSE)
# E(our_graph)$weight <- NA
# 
# barabasi_graph <- graph_from_data_frame(barabasi_network, directed=FALSE)
# E(barabasi_graph)$weight <- NA

# Degree
deg1 <- degree(our_graph, mode="all")
mean_deg1 <- mean(deg1); mean_deg1
deg1 <- deg1[order(names(deg1))]

deg2 <- degree(barabasi_graph, mode="all")
mean_deg2 <- mean(deg2); mean_deg2
deg2 <- deg2[order(names(deg2))]

dev.off()
pdf(paste("Network_building/Network_analysis/Correlations/",output_filename,"_epidem_correlations_plots.pdf",sep=""),width = 9, height = 7)

to_plot <- data.frame(nodes=as.numeric(names(deg1)), our_degree = deg1, barb_degree = deg2)
pp <- ggplot(to_plot, aes(deg1, deg2)) + geom_point() + 
  geom_text(label=to_plot$nodes, nudge_x = 0.3, nudge_y = 0.3, check_overlap = T)
pp


deg_corp <- cor.test(deg1, deg2, method="pearson"); deg_corp
deg_corp <- cor.test(deg1/sum(deg1), deg2/sum(deg2), method="pearson")
# deg_cors <- cor.test(deg1, deg2, method="spearman")
t.test(deg1, deg2, paired=TRUE)

p1 <- ggplot(to_plot,aes(our_degree,barb_degree))
p1 <- p1 + geom_point() + geom_smooth(method=lm, se=TRUE) +
  xlab("our_degree") + ylab("barb_degree") +
  theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
print(p1)

t.test(to_plot$our_degree, to_plot$barb_degree, paired=TRUE)

#Density of a graph: The proportion of present edges from all possible edges in the network
edge_density1 <- edge_density(our_graph, loops=F)
edge_density1

edge_density2 <- edge_density(barabasi_graph, loops=F)
edge_density2

# Transitivity: Transitivity measures the probability that the adjacent vertices of a vertex are connected
transitivity <- transitivity(our_graph, type="global"); transitivity
local_trans <- transitivity(our_graph, vids=V(our_graph), type="local")
transitivity_table <- data.frame(nodes=names(V(our_graph)), transitivity=local_trans)
print(ggplot(transitivity_table, aes(reorder(nodes,transitivity), transitivity)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis),plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        # ggtitle("Transitivity") + 
        # xlab("Diseases") + 
        ylab("Transitivity"))

transitivity2 <- transitivity(barabasi_graph, type="global"); transitivity2
local_trans2 <- transitivity(barabasi_graph, vids=V(barabasi_graph), type="local")
transitivity_table2 <- data.frame(nodes=names(V(barabasi_graph)), transitivity=local_trans2)

transitivity_table <- transitivity_table[order(transitivity_table$nodes), ]
transitivity_table2 <- transitivity_table2[order(transitivity_table2$nodes), ]
trans_corp <- cor.test(transitivity_table$transitivity, transitivity_table2$transitivity, method="pearson"); trans_corp
#trans_corp <- cor.test(deg1/sum(deg1), deg2/sum(deg2), method="pearson")
trans_cors <- cor.test(transitivity_table$transitivity, transitivity_table2$transitivity, method="spearman"); trans_cors
t.test(transitivity_table$transitivity, transitivity_table2$transitivity)

to_plot <- data.frame(nodes=transitivity_table$nodes, our = transitivity_table$transitivity, barb = transitivity_table2$transitivity)
pp <- ggplot(to_plot, aes(our, barb)) + geom_point() + 
  geom_text(label=to_plot$nodes, nudge_x = 0.01, nudge_y = 0.01, check_overlap = T)
pp

t.test(to_plot$our, to_plot$barb, paired=TRUE)

p1 <- ggplot(to_plot,aes(our,barb))
p1 <- p1 + geom_point() + geom_smooth(method=lm, se=TRUE) +
  xlab("our_transitivity") + ylab("barb_transitivity") +
  theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
print(p1)


### Diameter: longest geodesic distance (length of the shortest path between two nodes) in the network
# With weights
diameter <- diameter(our_graph, directed=F, weights=abs(E(our_graph)$weight)); diameter
get_diameter(our_graph, directed=F, weights=abs(E(our_graph)$weight))
diameter2 <- diameter(barabasi_graph, directed=F, weights=abs(E(barabasi_graph)$weight)); diameter2
# Without weights
diameter_wo_weights <- diameter(our_graph, directed=F, weights=NA); diameter_wo_weights 
hey <- get_diameter(our_graph, directed=F, weights=NA)

diameter_wo_weights2 <- diameter(barabasi_graph, directed=F, weights=NA); diameter_wo_weights2
get_diameter(barabasi_graph, directed=F, weights=NA)


###Degree
deg <- degree(our_graph, mode="all")
plot(graph, vertex.size=deg*0.4,layout=layout.circle)
hist(deg, breaks=1:vcount(our_graph)-1, main="Histogram of node degree")

# Degree distribution
deg.dist <- degree_distribution(our_graph, cumulative=T, mode="all")

plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")

### CENTRALITY AND CENTRALIZATION
# Degree (number of ties) ### hacer un plot de esto con ggplot quizá!
degree_table <- data.frame('nodes'=names(degree(our_graph, mode="in")), 'feature'=degree(our_graph, mode="total"))
centr_degree(our_graph, mode="in", normalized=T)
print(ggplot(degree_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        # ggtitle("Transitivity") + 
        # xlab("Diseases") + 
        ylab("Degree"))

# Closeness (centrality based on distance to others in the graph)
#Inverse of the node’s average geodesic distance to others in the network.
mean(closeness(our_graph, mode="all", weights=NA))
mean(closeness(barabasi_graph, mode="all", weights=NA))

closeness(our_graph, mode="all", weights=NA) 
centr_clo(our_graph, mode="all", normalized=T)
degree_table <- data.frame('nodes'=names(closeness(our_graph, mode="all", weights=NA)), 'feature'=closeness(our_graph, mode="all", weights=NA))
print(ggplot(degree_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        # ggtitle("Transitivity") + 
        # xlab("Diseases") + 
        ylab("Closeness"))

degree_table1 <- data.frame('nodes'=names(closeness(our_graph, mode="all", weights=NA)), 'feature'=closeness(our_graph, mode="all", weights=NA))
degree_table2 <- data.frame('nodes'=names(closeness(barabasi_graph, mode="all", weights=NA)), 'feature'=closeness(barabasi_graph, mode="all", weights=NA))
degree_table1 <- degree_table1[order(degree_table1$nodes), ]
degree_table2 <- degree_table2[order(degree_table2$nodes), ]
trans_corp <- cor.test(degree_table1$feature, degree_table2$feature, method="pearson"); trans_corp
#trans_corp <- cor.test(deg1/sum(deg1), deg2/sum(deg2), method="pearson")
trans_cors <- cor.test(degree_table1$feature, degree_table2$feature, method="spearman"); trans_cors
t.test(degree_table1$feature, degree_table2$feature, paired=TRUE)
t.test(degree_table1$feature, degree_table2$feature, paired=FALSE)

to_plot <- data.frame(nodes=degree_table1$nodes, our = degree_table1$feature, barb = degree_table2$feature)
pp <- ggplot(to_plot, aes(our, barb)) + geom_point() + 
  geom_text(label=to_plot$nodes, 
            nudge_x = 0.0001, nudge_y = 0.0001,
            check_overlap = T)
pp

p1 <- ggplot(to_plot,aes(our,barb))
p1 <- p1 + geom_point() + geom_smooth(method=lm, se=TRUE) +
  xlab("our_closeness") + ylab("barb_closeness") +
  theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
print(p1)


# Eigenvector (centrality proportional to the sum of connection centralities)
# Values of the first eigenvector of the graph matrix.

eigen_centrality(our_graph, directed=F, weights=NA)
centr_eigen(graph, directed=F, normalized=T) 

# Betweenness (centrality based on a broker position connecting others)
# Number of geodesics that pass through the node or the edge.
betweenness <- betweenness(our_graph, directed=FALSE, weights=abs(E(our_graph)$weight))
betweenness_wo_weights <- betweenness(our_graph, directed=FALSE, weights=NA)
edge_betweenness(our_graph, directed=T, weights=NA)
centr_betw(our_graph, directed=T, normalized=T)
betweenness <- data.frame('nodes'=names(betweenness), 'feature'=betweenness)
print(ggplot(betweenness, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        ggtitle("Betweenness with weights") +
        # xlab("Diseases") + 
        ylab("Betweenness"))

betweenness <- data.frame('nodes'=names(betweenness_wo_weights), 'feature'=betweenness_wo_weights)
print(ggplot(betweenness, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        ggtitle("Betweenness wo weights") +
        # xlab("Diseases") + 
        ylab("Betweenness"))

betweenness_wo_weights1 <- betweenness(our_graph, directed=FALSE, weights=NA)
betweenness_wo_weights2 <- betweenness(barabasi_graph, directed=FALSE, weights=NA)
degree_table1 <- data.frame('nodes'=names(betweenness_wo_weights1), 'feature'=betweenness_wo_weights1)
degree_table2 <- data.frame('nodes'=names(betweenness_wo_weights2), 'feature'=betweenness_wo_weights2)
degree_table1 <- degree_table1[order(degree_table1$nodes), ]
degree_table2 <- degree_table2[order(degree_table2$nodes), ]
trans_corp <- cor.test(degree_table1$feature, degree_table2$feature, method="pearson"); trans_corp
#trans_corp <- cor.test(deg1/sum(deg1), deg2/sum(deg2), method="pearson")
trans_cors <- cor.test(degree_table1$feature, degree_table2$feature, method="spearman"); trans_cors
t.test(degree_table1$feature, degree_table2$feature, paired=TRUE)

to_plot <- data.frame(nodes=degree_table1$nodes, our = degree_table1$feature, barb = degree_table2$feature)
pp <- ggplot(to_plot, aes(our, barb)) + geom_point() + 
  geom_text(label=to_plot$nodes, 
            nudge_x = 0.35, nudge_y = 0.35,
            check_overlap = T)
pp

mean(to_plot$our)
mean(to_plot$barb)
t.test(to_plot$our, to_plot$barb, paired=TRUE)

p1 <- ggplot(to_plot,aes(our,barb))
p1 <- p1 + geom_point() + geom_smooth(method=lm, se=TRUE) +
  xlab("our_betweenness") + ylab("barb_betweenness") +
  theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
print(p1)

p1 <- ggplot(to_plot,aes(our,barb))
p1 <- p1 + geom_point() + geom_smooth(method="loess", se=TRUE) +
  xlab("our_betweenness") + ylab("barb_betweenness") +
  theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
print(p1)



### HUBS AND AUTHORITIES
# HUBS with weights
hs <- hub_score(our_graph, weights=abs(E(our_graph)$weight))$vector
hubs_table <- data.frame('nodes'=names(hs), 'feature'=hs)
print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        ggtitle("Hubs with weights") +
        # xlab("Diseases") + 
        ylab("Hubs"))
# wo weights
hs <- hub_score(our_graph, weights=NA)$vector
hubs_table <- data.frame('nodes'=names(hs), 'feature'=hs)
print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        ggtitle("Hubs wo weights") +
        # xlab("Diseases") + 
        ylab("Hubs"))

hs2 <- hub_score(barabasi_graph, weights=NA)$vector
hubs_table2 <- data.frame('nodes'=names(hs2), 'feature'=hs2)
print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        ggtitle("Hubs wo weights") +
        # xlab("Diseases") + 
        ylab("Hubs"))

degree_table1 <- data.frame('nodes'=hubs_table$nodes, 'feature'=hubs_table$feature)
degree_table2 <- data.frame('nodes'=hubs_table2$nodes, 'feature'=hubs_table2$feature)
degree_table1 <- degree_table1[order(degree_table1$nodes), ]
degree_table2 <- degree_table2[order(degree_table2$nodes), ]

to_plot <- data.frame(nodes=degree_table1$nodes, our = degree_table1$feature, barb = degree_table2$feature)
pp <- ggplot(to_plot, aes(our, barb)) + geom_point() + 
  geom_text(label=to_plot$nodes, 
            nudge_x = 0.01, nudge_y = 0.01,
            check_overlap = T)
pp

p1 <- ggplot(to_plot,aes(our,barb))
p1 <- p1 + geom_point() + geom_smooth(method=lm, se=TRUE) +
  xlab("our_hubs") + ylab("barb_hubs") +
  theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
print(p1)

# # lm <- lm(formula = y ~ x) 
# p1 <- ggplot(to_plot,aes(our,barb))
# p1 <- p1 + geom_point() + geom_smooth(method="loess", se=FALSE) +
#   xlab("our_hubs") + ylab("barb_hubs") +
#   theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
#   theme(plot.title = element_text(hjust = 0.5))+
#   stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
# print(p1)


# AUTHORITIES with weights
as <- authority_score(graph, weights=abs(E(graph)$weight))$vector
auth_table <- data.frame('nodes'=names(as), 'feature'=as)
print(ggplot(auth_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        ggtitle("Authorities with weights") +
        # xlab("Diseases") + 
        ylab("Authority"))

# wo weights
as <- authority_score(graph, weights=NA)$vector
auth_table <- data.frame('nodes'=names(as), 'feature'=as)
print(ggplot(auth_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        ggtitle("Authorities wo weights") +
        # xlab("Diseases") + 
        ylab("Authority"))

par(mfrow=c(1,2))
plot(graph, vertex.size=hs*20, main="Hubs")
plot(graph, vertex.size=as*10, main="Authorities")

### DISTANCES AND PATHS
# Average path length: the mean of the shortest distance between 
#each pair of nodes in the network (in both directions for directed graphs).
mean_distance <- mean_distance(our_graph, directed=F); mean_distance
mean_distance2 <- mean_distance(barabasi_graph, directed=F); mean_distance2
# distances(graph) # with edge weights # YIELDS AND ERROR
distances(graph, weights=NA) # ignore weights

# K-core decomposition (degeneracy)
# The k-core is the maximal subgraph in which every node has degree of 
# at least k. The result here gives the coreness of each vertex in the 
# network. A node has coreness D if it belongs to a D-core but not 
# to (D+1)-core.
kc1 <- coreness(our_graph, mode="all")
# plot(graph, vertex.size=kc*6, vertex.label=kc, vertex.color=colrs[kc]) # FIXX
kc_table1 <- data.frame('nodes'=names(kc1), 'feature'=kc1)
print(ggplot(kc_table1, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        # ggtitle("Authorities with weights") +
        # xlab("Diseases") + 
        ylab("K-core decomposition"))

kc2 <- coreness(barabasi_graph, mode="all")
kc_table2 <- data.frame('nodes'=names(kc2), 'feature'=kc2)
print(ggplot(kc_table2, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
        # ggtitle("Authorities with weights") +
        # xlab("Diseases") + 
        ylab("K-core decomposition"))

mean(kc1); mean(kc2); median(kc1); median(kc2)

kc_table1 <- kc_table1[order(kc_table1$nodes), ]
kc_table2 <- kc_table2[order(kc_table2$nodes), ]
trans_corp <- cor.test(kc_table1$feature, kc_table2$feature, method="pearson"); trans_corp
#trans_corp <- cor.test(deg1/sum(deg1), deg2/sum(deg2), method="pearson")
trans_cors <- cor.test(kc_table1$feature, kc_table2$feature, method="spearman"); trans_cors
t.test(kc_table1$feature, kc_table2$feature, paired=TRUE)

to_plot <- data.frame(nodes=kc_table1$nodes, our = kc_table1$feature, barb = kc_table2$feature)
pp <- ggplot(to_plot, aes(our, barb)) + geom_point() + 
  geom_text(label=to_plot$nodes, 
            nudge_x = 0.05, nudge_y = 0.05,
            check_overlap = T)
pp

p1 <- ggplot(to_plot,aes(our,barb))
p1 <- p1 + geom_point() + geom_smooth(method=lm, se=TRUE) +
  xlab("our_kcore") + ylab("barb_kcore") +
  theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
print(p1)


### ASSORTATIVITY AND HOMOPHILY
# assortativity_nominal(graph, V(graph)$media.type, directed=F)
metadata <- read.csv("dl_general_info_umls_new_ncolors.csv",header=T, sep=",",stringsAsFactors = F); dim(metadata)
metadata <- metadata[, c("icd9", "disease_cat","disease_cat_number", "dis_cat_colors")]; head(metadata)
metadata$icd9 <- as.integer(metadata$icd9); head(metadata)
metadata <- metadata[!duplicated(metadata$icd9), ]; dim(metadata)
metadata <- metadata[metadata$icd9 %in% common_icds, ]; dim(metadata)

### FOR OUR GRAPH
V_order = as.integer(names(V(our_graph)))

# Sorting the dataframe in the same order as the nodes in the graph
match(V_order, metadata$icd9)
metadata2 <- metadata[match(V_order, metadata$icd9), ]; head(metadata2)

# This gives an error
# our_graph <- set_vertex_attr(our_graph, name="dis_category", value = metadata2$disease_cat)
# assortativity_degree <- assortativity_nominal(our_graph,types = V(our_graph)$dis_category, directed=F)

our_graph <- set_vertex_attr(our_graph, name="dis_category", value = metadata2$disease_cat_number)
assortativity_degree <- assortativity(our_graph, V(our_graph)$dis_category, directed=F) # 0.0230881
assortativity_degree

# ceb <- cluster_edge_betweenness(our_graph, weights = NA) 
# 
# dendPlot(ceb, mode="hclust")
# plot(ceb, net) 

### FOR BARABASI
V_order = as.integer(names(V(barabasi_graph)))

# Sorting the dataframe in the same order as the nodes in the graph
match(V_order, metadata$icd9)
metadata2 <- metadata[match(V_order, metadata$icd9), ]; head(metadata2)

# This gives an error
# our_graph <- set_vertex_attr(our_graph, name="dis_category", value = metadata2$disease_cat)
# assortativity_degree <- assortativity_nominal(our_graph,types = V(our_graph)$dis_category, directed=F)

barabasi_graph <- set_vertex_attr(barabasi_graph, name="dis_category", value = metadata2$disease_cat_number)
assortativity_degree <- assortativity(barabasi_graph, V(barabasi_graph)$dis_category, directed=F) # -0.03614256
assortativity_degree


cfeat <- c('N conn comp','Size conn comp','Mean deg',
           'edge_density','transitivity','diameter',
           'diameter_wo_weights','mean_distance','assortativity_degree') 
cvals <- c(n_connected_components,size_connected_components,mean_deg,
           edge_density,transitivity,diameter,
           diameter_wo_weights,mean_distance,assortativity_degree)
cdf <- data.frame(Feature=cfeat,Value = cvals)

dev.off()