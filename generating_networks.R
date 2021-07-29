#########################################################################################
# GENERATING AND ANALYZING NETWORKS FROM DISTANCE FILES
# 
# FINAL VERSION STORES IN ~/Desktop/ANALYSIS
#
# Beatriz Urda García 2020
#########################################################################################


library(igraph)
library(ggplot2)
library(gplots)


setwd("~/Desktop/ANALYSIS")
dismeta <- read.csv("dl_general_info_umls_new_ncolors.csv",header=T, sep=",",stringsAsFactors = F, row.names = NULL)
dismeta[dismeta$disease_name == 'ThyroidCancer_Papillary',]$icd9 <- 193
# dismeta$icd9 <- as.character(dismeta$icd9)

analyze_network <- function(dist_filename, output_filename, significant, posi, negat, draw_nets = FALSE){
  # Building the network
  graph <- graph_from_data_frame(significant, directed=FALSE)
  netpos <- graph_from_data_frame(posi, directed=FALSE)
  netneg <- graph_from_data_frame(negat, directed=FALSE)
  
  pdf(file=paste("Network_building/Network_analysis/",output_filename,"_networks.pdf",sep=""))
  par(mfrow=c(1,1))
  plot(graph, vertex.size = 5, main="Entire network")
  # plot(graph, vertex.size = 5, layout=layout.circle, main="Entire network")
  
  plot(netpos, vertex.size = 5, main="Positive interactions")
  # plot(netpos, vertex.size = 5, layout=layout.circle)
  
  plot(netneg, vertex.size = 5, main="Negative network")
  # plot(netneg, vertex.size = 5, layout=layout.circle)
  dev.off()
  
  ### NETWORK ANALYSIS ###
  # plot(graph, vertex.size=deg/5)
  # 
  # # Histogram of node degree
  # hist(deg, breaks=1:vcount(graph)-1, main="Histogram of node degree")
  # 
  # Degree distribution
  # deg <- degree(graph, mode="all")
  # deg.dist <- degree_distribution(graph, cumulative=T, mode="all")
  # png(file=paste("Network_building/Network_analysis/",output_filename,"_degdistribution.png",sep=""))
  # plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange",
  #       xlab="Degree", ylab="Cumulative Frequency")
  # dev.off()
  # 
  # Non-cumulative
  deg <- degree(graph, mode="all")
  df <- as.data.frame(table(deg))
  df$deg <- as.numeric(df$deg)
  p <- ggplot(df, aes(deg, Freq)) + geom_point() +
    xlab("Degree") + ylab("Number of nodes") + ggtitle(paste("Degree distribution | ",gsub("_distance.*","",output_filename),sep=""))
  p
  ggsave(paste("Network_building/Network_analysis/",output_filename,"_degdistribution.png",sep=""), width = 5, height = 5)
  # ggsave(paste("Network_building/Defined_networks/",specific_filename,"_pos_network.png",sep=""), width = 5, height = 5)
  # ggsave(paste("Network_building/Defined_networks/",specific_filename,"_neg_network.png",sep=""), width = 5, height = 5)
  # 
  # Número de componentes conexas
  # graph
  connected_components <- clusters(graph)
  n_connected_components <- connected_components$no
  size_connected_components <- connected_components$csize
  
  deg <- degree(graph, mode="all")
  mean_deg <- mean(deg) 
  
  # posi
  connected_components_pos <- clusters(netpos)
  n_connected_components_pos <- connected_components_pos$no
  size_connected_components_pos <- connected_components_pos$csize
  if(length(size_connected_components_pos) > 1){
    size_connected_components_pos <- paste(size_connected_components_pos[1],"_",size_connected_components_pos[2])
  }
  
  deg_pos <- degree(netpos, mode="all")
  mean_deg_pos <- mean(deg_pos) 
  
  # negat
  connected_components_neg <- clusters(netneg)
  n_connected_components_neg <- connected_components_neg$no
  size_connected_components_neg <- connected_components_neg$csize
  
  deg_neg <- degree(netneg, mode="all")
  mean_deg_neg <- mean(deg_neg) 
  
  if(length(size_connected_components)>1){
    size_connected_components <- paste(size_connected_components,collapse="_")
  }
  if(length(size_connected_components_pos)>1){
    size_connected_components_pos <- paste(size_connected_components_pos,collapse="_")
  }
  if(length(size_connected_components_neg)>1){
    size_connected_components_neg <- paste(size_connected_components_neg,collapse ="_")
  }
  
  cfeat <- c('N conn comp','Size conn comp','Mean deg',
             'N conn comp pos','Size conn comp pos','Mean deg pos',
             'N conn comp neg','Size conn comp neg','Mean deg neg') 
  cvals <- c(n_connected_components,size_connected_components,mean_deg,
             n_connected_components_pos,size_connected_components_pos,mean_deg_pos,
             n_connected_components_neg,size_connected_components_neg,mean_deg_neg)
  cdf <- data.frame(Feature=cfeat,Value = cvals)
   
  
  # Network representation
  netm <- get.adjacency(graph, attr="Distance", sparse=F)
  png(file=paste("Network_building/Network_analysis/",output_filename,"_clustering_heatmap.png",sep=""),width = 800, height = 800)
  heatmap.2(netm, col = bluered(100), 
            scale="none", margins=c(16,16), trace="none" )
  dev.off()
  
  netm <- get.adjacency(graph, attr="Distance", sparse=F)
  png(file=paste("Network_building/Network_analysis/",output_filename,"_heatmap.png",sep=""),width = 800, height = 800)
  
  heatmap.2(netm, col = bluered(100), 
            dendrogram = "none", Rowv=FALSE, Colv=FALSE, # Not reordering
            scale="none", margins=c(16,16), trace="none" )
  dev.off()
  # ggsave(paste("Network_building/Network_analysis/",specific_filename,"_heatmap.png",sep=""), width = 5, height = 5)
  return(cdf)
}

generate_network <- function(dist_filename, output_filename, analysis = TRUE){
  # dist_filename <- dist_filename_list[3]
  # output_filename <- output_filename_list[3]
  # print(dist_filename)
  # print(output_filename)
  
  dist <- read.csv(dist_filename,header=T, sep=",",stringsAsFactors = F)
  dim(dist)
  if(analysis == TRUE){
    cdf <- data.frame(Feature = "N posible int", Value = dim(dist)[1])
  }
  dist <- dist[is.na(dist$pvalue) == FALSE,]
  dim(dist)
  
  # Adding the adjusted p-value
  dist$adj_pvalue <- p.adjust(dist$pvalue, method="fdr")
  
  # Filtering: keeping the significant interactions
  # wo FDR
  significant <- dist[dist$pvalue <= 0.05,]
  dim(significant) # 1149    5
  # with FDR
  significant <- dist[dist$adj_pvalue <= 0.05,]
  nint_total <- dim(significant)[1] # All genes: 1141  # sDEGs:  328
  nint_total
  
  posi <- significant[significant$Distance > 0,]   
  negat <- significant[significant$Distance < 0,]   
  nint_posi <-dim(posi)[1]        # All genes: -  # sDEGs:  241
  nint_negat <- dim(negat)[1]     # All genes: -  # sDEGs:  87
  nint_posi
  nint_negat
  
  (nint_posi / nint_total)*100
  (nint_negat / nint_total)*100
  
  write.table(significant,file=paste("Network_building/Defined_networks/",output_filename,"_network.txt",sep=""),sep="\t",row.names=F)
  write.table(posi,file=paste("Network_building/Defined_networks/",output_filename,"_pos_network.txt",sep=""),sep="\t",row.names=F)
  write.table(negat,file=paste("Network_building/Defined_networks/",output_filename,"_neg_network.txt",sep=""),sep="\t",row.names=F)
  
  if(analysis == TRUE){
    cfeat <- c('N sign int','N pos int','N neg int','Perc pos int','Perc neg int' ) 
    cvals <- c(dim(significant)[1], nint_posi, nint_negat, (nint_posi / nint_total)*100, (nint_negat / nint_total)*100 )
    cdf <- rbind(cdf, data.frame(Feature=cfeat,Value = cvals))
    adf <- analyze_network(dist_filename, output_filename, significant, posi, negat)
    cdf <- rbind(cdf, adf)
    
  }  
  cdf$Feature <- as.character(cdf$Feature)
  return(cdf)
}


# A) GENERATE AND ANALYZE ICD9 LEVEL NETWORKS USING SDEGs 
dist_filename_list <- c("Network_building/icd9_distances/pairwise_spearman_distance_sDEGs.csv",
                        "Network_building/icd9_distances/pairwise_union_spearman_distance_sDEGs.csv",
                        "Network_building/icd9_distances/comparable_spearman_distance_sDEGs.csv")
output_filename_list <- c("icd9_pairwise_spearman_distance_sDEGs",
                     "icd9_pairwise_union_spearman_distance_sDEGs",
                     "icd9_comparable_spearman_distance_sDEGs")
specific_filename <- "icd9_sDEGs_level"

# B) GENERATE AND ANALYZE DISEASE LEVEL NETWORKS USING SDEGs 
dist_filename_list <- c("Network_building/distances/pairwise_spearman_distance_sDEGs.csv",
                        "Network_building/distances/pairwise_union_spearman_distance_sDEGs.csv",
                        "Network_building/distances/comparable_spearman_distance_sDEGs.csv")
output_filename_list <- c("pairwise_spearman_distance_sDEGs",
                          "pairwise_union_spearman_distance_sDEGs",
                          "comparable_spearman_distance_sDEGs")
specific_filename <- "disease_sDEGs_level"

# C) GENERATE AND ANALYZE DISEASE LEVEL AND ICD9 LEVEL NETWORKS USING SDEGs 
dist_filename_list <- c("Network_building/distances/pairwise_spearman_distance_sDEGs.csv",
                        "Network_building/distances/pairwise_union_spearman_distance_sDEGs.csv",
                        "Network_building/distances/comparable_spearman_distance_sDEGs.csv",
                        "Network_building/icd9_distances/pairwise_spearman_distance_sDEGs.csv",
                        "Network_building/icd9_distances/pairwise_union_spearman_distance_sDEGs.csv",
                        "Network_building/icd9_distances/comparable_spearman_distance_sDEGs.csv")

output_filename_list <- c("pairwise_spearman_distance_sDEGs",
                          "pairwise_union_spearman_distance_sDEGs",
                          "comparable_spearman_distance_sDEGs",
                          "icd9_pairwise_spearman_distance_sDEGs",
                          "icd9_pairwise_union_spearman_distance_sDEGs",
                          "icd9_comparable_spearman_distance_sDEGs")

specific_filename <- "disease_and_icd9_sDEGs_level"

# D) GENERATE AND ANALYZE DISEASE AND METAPATIENT LEVEL NETWORKS USING SDEGs 
dist_filename_list <- c("Network_building/distances/metapatients_and_disease/pairwise_spearman_distance_sDEGs.csv",
                        "Network_building/distances/metapatients_and_disease/pairwise_union_spearman_distance_sDEGs.csv",
                        "Network_building/distances/metapatients_and_disease/comparable_spearman_distance_sDEGs.csv")
output_filename_list <- c("metapatients_and_disease/metap_dis_pairwise_spearman_distance_sDEGs",
                          "metapatients_and_disease/metap_dis_pairwise_union_spearman_distance_sDEGs",
                          "metapatients_and_disease/metap_dis_comparable_spearman_distance_sDEGs")
specific_filename <- "metap_dis_sDEGs_level"

# E) GENERATE AND ANALYZE DISEASE LEVEL NETWORKS USING SDVGs 
dist_filename_list <- c("Network_building/distances/pairwise_spearman_distance_sDVGs.csv",
                        "Network_building/distances/pairwise_union_spearman_distance_sDVGs.csv",
                        "Network_building/distances/comparable_spearman_distance_sDVGs.csv")
output_filename_list <- c("pairwise_spearman_distance_sDVGs",
                          "pairwise_union_spearman_distance_sDVGs",
                          "comparable_spearman_distance_sDVGs")
specific_filename <- "disease_sDVGs_level"

# F) GENERATE AND ANALYZE ICD9 LEVEL NETWORKS USING SDVGs 
dist_filename_list <- c("Network_building/icd9_distances/pairwise_spearman_distance_sDVGs.csv",
                        "Network_building/icd9_distances/pairwise_union_spearman_distance_sDVGs.csv",
                        "Network_building/icd9_distances/comparable_spearman_distance_sDVGs.csv")
output_filename_list <- c("icd9_pairwise_spearman_distance_sDVGs",
                          "icd9_pairwise_union_spearman_distance_sDVGs",
                          "icd9_comparable_spearman_distance_sDVGs")
specific_filename <- "icd9_disease_sDVGs_level"
#------------------------------------------------------

rdf <- generate_network(dist_filename_list[1], output_filename_list[1], analysis = TRUE)
for(k in 2:length(dist_filename_list)){
  ccol <- generate_network(dist_filename_list[k], output_filename_list[k], analysis = TRUE)
  rdf <- merge(rdf, ccol, by="Feature", sort=FALSE)
}
colnames(rdf)[2:(length(output_filename_list)+1)] <- output_filename_list
rdf
write.table(rdf,file=paste("Network_building/Network_analysis/",specific_filename,"_net_analysis.txt",sep=""),sep="\t",row.names=F,quote=FALSE)

####################### ANALYZING BARABASI's NETWORK #####################---SA
dist_filename <- "Network_building/final_icd9_level_networks/icd9_pairwise_union_spearman_distance_sDEGs_pos_B_TRUE_final_barabasi.txt"
output_filename <- "analysis_icd9_pairwise_union_spearman_final_barabasi.txt"

significant <- read.csv(dist_filename,header=T, sep="\t",stringsAsFactors = F)
# hey <- analyze_network(dist_filename, output_filename, significant, posi, negat, draw_nets = FALSE)

nint_total <- dim(significant)[1] # All genes: 1141  # sDEGs:  328
nint_total

# Building the network
graph <- graph_from_data_frame(significant, directed=FALSE)

pdf(file=paste("Network_building/Network_analysis/",output_filename,"_networks.pdf",sep=""))
par(mfrow=c(1,1))
plot(graph, vertex.size = 5, main="Entire network")
# plot(graph, vertex.size = 5, layout=layout.circle, main="Entire network")
dev.off()

# Non-cumulative
deg <- degree(graph, mode="all")
df <- as.data.frame(table(deg))
df$deg <- as.numeric(df$deg)
p <- ggplot(df, aes(deg, Freq)) + geom_point() +
  xlab("Degree") + ylab("Number of nodes") + ggtitle(paste("Degree distribution | ",gsub("_distance.*","",output_filename),sep=""))
p
ggsave(paste("Network_building/Network_analysis/",output_filename,"_degdistribution.png",sep=""), width = 5, height = 5)
# ggsave(paste("Network_building/Defined_networks/",specific_filename,"_pos_network.png",sep=""), width = 5, height = 5)
# ggsave(paste("Network_building/Defined_networks/",specific_filename,"_neg_network.png",sep=""), width = 5, height = 5)
# 
# Número de componentes conexas
# graph
connected_components <- clusters(graph)
n_connected_components <- connected_components$no
size_connected_components <- connected_components$csize
n_connected_components
size_connected_components

deg <- degree(graph, mode="all")
mean_deg <- mean(deg) 
mean_deg

########### OBTAINING THE COLOURED HEATMAP FOR MY NETWORK ######################
dist <- read.csv("Network_building/distances/pairwise_union_spearman_distance_sDEGs.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "pairwise_union_spearman_sDEGs"
dim(dist)
dist <- dist[is.na(dist$pvalue) == FALSE,]
dist$adj_pvalue <- p.adjust(dist$pvalue, method="fdr")

# Filtering: keeping the significant interactions
# wo FDR
significant <- dist[dist$pvalue <= 0.05,]

# Sorting significant network
to_merge <- dismeta[,c("disease_name","icd9","dis_cat_colors")]
dim(significant)
dim(to_merge)
head(significant)
significant <- merge(significant,to_merge, by.x = 'Dis1', by.y = 'disease_name')
dim(significant)
head(significant)

# Sort the df by icd
significant <- significant[order(significant$icd9),]

# Put digestive system disorders in the end to color them properly
digestive_index <- which(significant$dis_cat_colors == 'springgreen1')
digestive <- significant[significant$dis_cat_colors == 'springgreen1', ]
dim(significant)
dim(digestive)
significant <- significant[-digestive_index,]
significant <- rbind(significant, digestive)

cdiseasecolors <- significant[!duplicated(significant$Dis1),]
cdiseasecolors <- cdiseasecolors$dis_cat_colors
cdiseasecolors <- append(cdiseasecolors, 'springgreen1')

significant$icd9 <- NULL ; significant$dis_cat_colors <- NULL

graph <- graph_from_data_frame(significant, directed=FALSE)
netm <- get.adjacency(graph, attr="Distance", sparse=F)
n_diseases <- dim(netm)[1]

# for(k in 1:n_diseases){
#   netm[k,k] <- 1
# }

png(file=paste("Network_building/Network_analysis/",specific_filename,"_coloured_heatmap.png",sep=""),width = 800, height = 800)
heatmap.2(netm, col = bluered(100), 
          dendrogram = "none", Rowv=FALSE, Colv=FALSE, # Not reordering,
          colCol = cdiseasecolors,
          colRow = cdiseasecolors,
          cexRow = 0.8,
          # key=FALSE,
          scale="none", margins=c(13,15), trace="none" )
dev.off()

pdf(file=paste("Network_building/Network_analysis/",specific_filename,"_coloured_heatmap_diagonal_fixed_pdf.pdf",sep=""))
heatmap.2(netm, col = bluered(100), 
          dendrogram = "none", Rowv=FALSE, Colv=FALSE, # Not reordering,
          colCol = cdiseasecolors,
          colRow = cdiseasecolors,
          cexRow = 0.5,cexCol = 0.5,
          key=FALSE,
          scale="none", margins=c(10,10), trace="none")
dev.off()

heatmap.2(cenf2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.5,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)


  
########### DOING THE ANALYSIS CONDITION BY CONDITION

# A) Spearman Pairwise All genes
dist <- read.csv("Network_building/distances/pairwise_spearman_distance_allgenes.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "pairwise_spearman_allgenes"
dim(dist)
# 1326

# B) Spearman Pairwise sDEGs
dist <- read.csv("Network_building/distances/pairwise_spearman_distance_sDEGs.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "pairwise_spearman_sDEGs"
dim(dist)
# 1326

# C) Spearman Union All genes
dist <- read.csv("Network_building/distances/comparable_spearman_distance_allgenes.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "union_spearman_allgenes"
dim(dist)
# 1326

# D) Spearman Union sDEGs
dist <- read.csv("Network_building/distances/comparable_spearman_distance_sDEGs.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "union_spearman_sDEGs"
dim(dist)
# 1326

# E) Spearman Pairwise Union sDEGs
dist <- read.csv("Network_building/distances/pairwise_union_spearman_distance_sDEGs.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "pairwise_union_spearman_sDEGs"
dim(dist)
# 1326

# F) Spearman Pairwise Union All genes
dist <- read.csv("Network_building/distances/pairwise_union_spearman_distance_allgenes.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "pairwise_union_spearman_sDEGs"
dim(dist)
# 1326

# E) Spearman Pairwise DM
dist <- read.csv("Network_building/distances/pairwise_spearman_distance_DM.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "pairwise_spearman_DM"
dim(dist)
# 1326

### WITH ICD9

dist <- read.csv("Network_building/icd9_distances/pairwise_spearman_distance_sDEGs.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "icd9_pairwise_spearman_distance_sDEGs"
dim(dist)
# 1326

dist <- read.csv("Network_building/icd9_distances/pairwise_union_spearman_distance_sDEGs.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "icd9_pairwise_spearman_distance_sDEGs"
dim(dist)
# 1326

dist <- read.csv("Network_building/icd9_distances/comparable_spearman_distance_sDEGs.csv",header=T, sep=",",stringsAsFactors = F)
specific_filename <- "icd9_pairwise_spearman_distance_sDEGs"
dim(dist)
# 1326

dist <- dist[is.na(dist$pvalue) == FALSE,]
dim(dist)
# 546






write.table(cdf,file=paste("Network_building/Network_analysis/",specific_filename,"_net_analysis.txt",sep=""),sep="\t",row.names=F)



#### ANALYZING NETWORKS ONE AT A TIME IN DETAIL

# Adding the adjusted p-value
dist$adj_pvalue <- p.adjust(dist$pvalue, method="fdr")

# Filtering: keeping the significant interactions
# wo FDR
significant <- dist[dist$pvalue <= 0.05,]
dim(significant) # 1149    5
# with FDR
significant <- dist[dist$adj_pvalue <= 0.05,]
nint_total <- dim(significant)[1] # All genes: 1141  # sDEGs:  328
nint_total

posi <- significant[significant$Distance > 0,]   
negat <- significant[significant$Distance < 0,]   
nint_posi <-dim(posi)[1]        # All genes: -  # sDEGs:  241
nint_negat <- dim(negat)[1]     # All genes: -  # sDEGs:  87
nint_posi
nint_negat

(nint_posi / nint_total)*100
(nint_negat / nint_total)*100


write.table(significant,file=paste("Network_building/Defined_networks/",specific_filename,"_network.txt",sep=""),sep="\t",row.names=F)
write.table(posi,file=paste("Network_building/Defined_networks/",specific_filename,"_pos_network.txt",sep=""),sep="\t",row.names=F)
write.table(negat,file=paste("Network_building/Defined_networks/",specific_filename,"_neg_network.txt",sep=""),sep="\t",row.names=F)

# Building the network
graph <- graph_from_data_frame(significant, directed=FALSE)
netpos <- graph_from_data_frame(posi, directed=FALSE)
netneg <- graph_from_data_frame(negat, directed=FALSE)

par(mfrow=c(1,1))
plot(graph, vertex.size = 5)
plot(graph, vertex.size = 5, layout=layout.circle)

plot(netpos, vertex.size = 5)
# plot(netpos, vertex.size = 5, layout=layout.circle)

plot(netneg, vertex.size = 5)
plot(netneg, vertex.size = 5, layout=layout.circle)

### NETWORK ANALYSIS ###
# graph <- netpos
# Node degrees
deg <- degree(graph, mode="all")
# Mean degree
mean_deg <- mean(deg) 
mean_deg

plot(graph, vertex.size=deg/5)

# Histogram of node degree
hist(deg, breaks=1:vcount(graph)-1, main="Histogram of node degree")

# Degree distribution
deg.dist <- degree_distribution(graph, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")

# Non-cumulative
df <- as.data.frame(table(deg))
df$deg <- as.numeric(df$deg)
p <- ggplot(df, aes(deg, Freq)) + geom_point() +
  xlab("Degree") + ylab("Number of nodes") + ggtitle("Degree distribution")
p
ggsave(paste("Network_building/Defined_networks/",specific_filename,"_network.png",sep=""), width = 5, height = 5)
ggsave(paste("Network_building/Defined_networks/",specific_filename,"_pos_network.png",sep=""), width = 5, height = 5)
ggsave(paste("Network_building/Defined_networks/",specific_filename,"_neg_network.png",sep=""), width = 5, height = 5)

# Número de componentes conexas
# graph
connected_components <- clusters(graph)
n_connected_components <- connected_components$no
n_connected_components
size_connected_components <- connected_components$csize
size_connected_components

deg <- degree(graph, mode="all")
mean_deg <- mean(deg) 
mean_deg

# posi
connected_components <- clusters(netpos)
n_connected_components <- connected_components$no
n_connected_components
size_connected_components <- connected_components$csize
size_connected_components

deg <- degree(netpos, mode="all")
mean_deg <- mean(deg) 
mean_deg

# negat
connected_components <- clusters(netneg)
n_connected_components <- connected_components$no
n_connected_components
size_connected_components <- connected_components$csize
size_connected_components

deg <- degree(netneg, mode="all")
mean_deg <- mean(deg) 
mean_deg

# Centrality - DO NOT INCLUDE IT IN THE TABLE
degree(graph, mode="in")

centr_degree(graph, mode="in", normalized=T)

closeness(graph, mode="all", weights=NA) 

centr_clo(graph, mode="all", normalized=T) 

# ...

# HUBS
hs <- hub_score(graph, weights=NA)$vector
as <- authority_score(graph, weights=NA)$vector

par(mfrow=c(1,2))
plot(graph, vertex.size=hs*15, main="Hubs")
plot(graph, vertex.size=as*15, main="Authorities")

# DISTANCES AND PATHS

# SUBGROUPS AND COMMUNITIES
# Cliques
cliques(graph)


# Network representation
netm <- get.adjacency(graph, attr="Distance", sparse=F)
heatmap(netm, Rowv = NA, Colv = NA, 
        scale="none", margins=c(10,10) )


heatmap.2(netm, col = bluered(100), 
        scale="none", trace="none" )




