#####################################################################################
# INTRA-INTERCATEGORY DENSITY OF THE DSN
#
#  
# 
#  Beatriz Urda García 2021    dsn_intra_inter_links.R
######################################################################################


library(igraph)


setwd("~/Desktop/ANALYSIS/Network_building/")

metadata =  read.csv2('../new_disease_metadata_final_names.txt',stringsAsFactors = F,sep="\t",header=T)
dsn_filename = "Defined_networks/pairwise_union_spearman_distance_sDEGs_pos_network.txt"
dsn = read.csv2(dsn_filename, header=TRUE, sep="\t",stringsAsFactors = F)

to_merge <- metadata[, c("disease_name", "disease_cat")]
dim(dsn)
dsn <- merge(dsn, to_merge, by.x = "Dis1", by.y = "disease_name", all.x=TRUE, all.y = FALSE ); dim(dsn)
dsn <- merge(dsn, to_merge, by.x = "Dis2", by.y = "disease_name", all.x=TRUE, all.y = FALSE ); dim(dsn)
colnames(dsn)[6:7] <- c("Dis1_cat", "Dis2_cat")

dsn <- dsn[, c("Dis1", "Dis2", "Distance", "adj_pvalue", "Dis1_cat", "Dis2_cat")]

discats <- unique(to_merge$disease_cat)
ndiscats <- length(discats); ndiscats

m <- matrix(NA, nrow=ndiscats, ncol=ndiscats)
for (k1 in 1:length(discats)){
  # k1 = 1; k2 = 2      # to test
  cat1 = discats[k1]
  for (k2 in 1:length(discats)){
    cat2 = discats[k2]
    subnet = dsn[((dsn$Dis1_cat == cat1) & (dsn$Dis2_cat == cat2)) | ((dsn$Dis1_cat == cat2) & (dsn$Dis2_cat == cat1)), ]
    if(dim(subnet)[1] != 0){
      subgraph = graph_from_data_frame(subnet, directed = FALSE)
      m[k1, k2] = edge_density(subgraph)
    }
  }
}

df <- as.data.frame(m)
colnames(df) <- discats; rownames(df) <- discats
write.table(df,file="Network_analysis/density_by_category_pair.txt", quote =FALSE, sep="\t", col.names = NA)

library(gplots)

m <- t(as.matrix(df))
g <- heatmap.2(m, scale="none",dendrogram = "none",
               cellnote = round(m,2), notecol="black", notecex = 0.7,
               trace="none",margins=c(18,18),
               cexRow = 0.8,cexCol = 0.8,
               main = "Precision of the DSN")




