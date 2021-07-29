#########################################################################################
# GENERATING AND ANALYZING NETWORKS FROM DISTANCE FILES
# 
# FINAL VERSION STORES IN ~/Desktop/ANALYSIS
#
# Beatriz Urda García 2020
########################################################################################

library(igraph)
library(ggplot2)
library(Hmisc) # function: capitalize

setwd("~/Desktop/ANALYSIS")

dismeta <- read.csv("new_disease_metadata_final_names.txt",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)

## Create Supl. Figures 2A y B. Correlation between #sDEGs with Sample size and library size
pdf(file=paste("Final_figures/Correlation_sDEGs_N_LibrarySize.pdf", sep=""), height = 5 , width = 6)

## Correlation sDEGs & Number of samples
p1 <- ggplot(dismeta ,aes(n_samples,Number_sDEGs))
p1 <- p1 + geom_point() + geom_smooth(method=lm, se=FALSE) +
  theme_classic()+
  xlab("Sample size") + ylab("Number of sDEGs") +
  theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme()+
  stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
print(p1)

wo_outlier <- dismeta[dismeta$n_samples < 300, ]
p1 <- ggplot(wo_outlier ,aes(n_samples,Number_sDEGs))
p1 <- p1 + geom_point() + geom_smooth(method=lm, se=FALSE) +
  theme_classic()+
  xlab("Sample size") + ylab("Number of sDEGs") +
  theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
print(p1)

## Correlation sDEGs & Average Library Size
p1 <- ggplot(dismeta ,aes(Library_size,Number_sDEGs))
p1 <- p1 + geom_point() + geom_smooth(method=lm, se=FALSE) +
  theme_classic()+
  xlab("Library Size") + ylab("Number of sDEGs") +
  theme(plot.margin = unit(c(1.3,1.1,1,1), "cm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_cor(method = "pearson", label.x.npc = 'center', label.y.npc = 0.4)
print(p1)

dev.off()

# Sort dismeta
dismeta <- dismeta[order(dismeta$disease_cat), ]
disease_cat_order <- c("congenital anomalies", "diseases of the circulatory system", "diseases of the digestive system",
                       "diseases of the musculoskeletal system and connective tissue",
                       "diseases of the respiratory system", "diseases of the nervous system and sense organs",
                       "mental disorders", "infectious and parasitic diseases", "neoplasms")
desired_order = data.frame("disease_cat"=disease_cat_order, "order"=c(1:length(disease_cat_order)))

dismeta <- merge(dismeta, desired_order, by="disease_cat", all.x=TRUE, all.y=FALSE)
dismeta <- dismeta[order(dismeta$order, dismeta$final_disease_name), ]

newnamesdf <- dismeta[, c("disease_cat", "disease_name","final_disease_name", "new_dis_cat_colors")]
names_correspondences <- dismeta[, c("disease_name","final_disease_name")]

network_path <- "Network_building/Defined_networks/pairwise_union_spearman_distance_sDEGs_network.txt"
output_filename <- "pairwise_union_spearman_distance_sDEGs"

network <- read.csv(network_path,header=T, sep="\t",stringsAsFactors = F)
dim(network)
network <- network[, c("Dis1", "Dis2", "Distance")]; head(network)

# Transform network into correct names
network <- merge(network, names_correspondences, by.x="Dis1", by.y="disease_name", all.x = TRUE, all.y = FALSE)
network <- merge(network, names_correspondences, by.x="Dis2", by.y="disease_name", all.x = TRUE, all.y = FALSE)
network <- network[, c("final_disease_name.x", "final_disease_name.y", "Distance")]
colnames(network) <- c("Dis1", "Dis2", "Distance"); head(network); dim(network)

graph <- graph_from_data_frame(network, directed=FALSE)
netm <- get.adjacency(graph, attr="Distance", sparse=F)

# Sort by ICD9 category
netmdf <- as.data.frame(netm)
netmdf <- netmdf[, newnamesdf$final_disease_name]
netmdf <- netmdf[newnamesdf$final_disease_name, ]

# palette.breaks=c(seq(-1,-0.01,length=10),seq(-0.001,0.001,length=10),seq(0.01,1,length=10))
# color.palette<-colorRampPalette(c("#C81E17","#FFFFFF","#405191"))(length(palette.breaks)-1)
pdf(file=paste("Final_figures/Heatmap_",output_filename,"_.pdf", sep=""), height = 10 , width = 10)
heatmap.2(as.matrix(netmdf),col=colorpanel(100, "#405191","#FFFFFF","#C81E17"), 
          dendrogram = "none", Rowv=FALSE, Colv=FALSE, # Not reordering
          scale="none", margins=c(15,15), trace="none",
          colCol = newnamesdf$new_dis_cat_colors,
          colRow = newnamesdf$new_dis_cat_colors)

heatmap.2(as.matrix(netmdf),col = bluered(100), 
          dendrogram = "none", Rowv=FALSE, Colv=FALSE, # Not reordering
          scale="none", margins=c(15,15), trace="none",
          colCol = newnamesdf$new_dis_cat_colors,
          colRow = newnamesdf$new_dis_cat_colors)

# Create legend
legendcolor <- dismeta[, c("disease_cat", "new_dis_cat_colors")]
legendcolor$disease_cat <- as.character(legendcolor$disease_cat)
legendcolor <- legendcolor[!duplicated(legendcolor$disease_cat), ]
legendcolor$disease_cat <- capitalize(legendcolor$disease_cat)
to_plot_legend <- rep(0,length(legendcolor[,1])); names(to_plot_legend) <- legendcolor$disease_cat
barplot(to_plot_legend,ylim=c(0,10),las=2,col=legendcolor$new_dis_cat_colors,
        legend=T,args.legend = list("cex"=0.5,"border"=NA))

dev.off()
# col = bluered(100)

