#####################################################################################
# COMPUTING NETWORK OVERLAP AT THE META-PATIENT LEVEL
#
#   
#  
#  
# 
# Beatriz Urda García 2021
######################################################################################

library(igraph)
library(ggplot2)

setwd("~/Desktop/ANALYSIS")

filedir <- 'Network_building/Defined_networks/metapatients_and_disease/'
filename <- 'metap_dis_pairwise_union_spearman_distance_sDEGs_pos_network.txt' # SSN with interactions: M-M, M-D and D-D
filepath <- paste(filedir,filename,sep="")
randomizing_labels = TRUE
interactions <- 'metap_metap'; interactions <- "all_metap"; interactions <- 'metap_dis'
# interactions <- 'all'
filtering <- 'wo_filtering'; filtering <- 'filtering'

outpath <- paste('Network_building/Overlapping_results/Shuffling_labels/SSN/',filtering,"/",gsub('_network.txt',"",filename),'_',interactions,'_',randomizing_labels,'_overlap.txt',sep="")

network <- read.csv(filepath,header=T, sep="\t",stringsAsFactors = F)
filecon <- file(outpath)
writeLines(paste('OVERLAP ANALYSIS OF THE NETWORK:\t',filename,sep=""),filecon)
close(filecon)

metadata =  read.csv2('new_disease_metadata_final_names.txt',stringsAsFactors = F,sep="\t",header=T)
dismeta <- read.csv("dl_general_info_umls_new_ncolors.csv",header=T, sep=",",stringsAsFactors = F, row.names = NULL)

extended_node_names =  read.csv2('Metapatients/with_entire_count_matrix/DEAnalysis/extended_node_names.txt',stringsAsFactors = F,sep="\t",header=T)


# Create a df with diseases and their icds, keeping only first 3 letters for the icd codes.
icds <- metadata[,c('disease_name','icd10','icd9')]

# Select interactions - it could also be done with metaps- and diseases
if(interactions == 'metap_metap'){       # M-M. Only interactions between meta-patients are considered
  fnetwork <- network[grepl("\\d", network$Dis1) & grepl("\\d", network$Dis2), ]  # meta_patients contain a digit. CUIDADO CON LA DISEASE QUE NO TIENE METAPS
}else if(interactions == 'all_metap'){   # M-M and M-D. Only interactions between meta-patients and metap-diseases are considered
  fnetwork <- network[grepl("\\d", network$Dis1) | grepl("\\d", network$Dis2), ]
}else if(interactions == 'metap_dis'){  ### M-D. 
  fnetwork1 <- network[!grepl("\\d", network$Dis1) & grepl("\\d", network$Dis2), ]      # Dis1 D and Dis2 M
  fnetwork2 <- network[grepl("\\d", network$Dis1) & (!grepl("\\d", network$Dis2)), ]    # Dis1 M and Dis2 D   --- No hay!
  # fnetwork <- network[grepl("\\d", network$Dis1), ]; unique(fnetwork$Dis2)           #                     comprobación en esta línea
  fnetwork <- rbind(fnetwork1, fnetwork2)
}else if(interactions == 'all'){
  fnetwork <- network
}

# Filtering interactions based on sDEGs
if(filtering == 'filtering'){
  # wo_sdegs <- sort(unique(extended_node_names[extended_node_names$sDEGs == 0, ]$old_names))
  with_sdegs <- unique(extended_node_names[extended_node_names$sDEGs > 0, ]$old_names)
  fnetwork <- fnetwork[(fnetwork$Dis1 %in% with_sdegs) & ((fnetwork$Dis2 %in% with_sdegs)), ]
}

## To compare it with the epidemiological network, use the ICD9 codes. 

# First, remove the <_numbers> from the meta-patients name
fnetwork$Dis1 <- gsub("_\\d+","",fnetwork$Dis1)
fnetwork$Dis2 <- gsub("_\\d+","",fnetwork$Dis2)
dim(fnetwork) # 3124    5

# Second, Obtain the networks in icds
merged <- merge(fnetwork, icds, by.x = "Dis1", by.y="disease_name", all.x = TRUE, all.y=FALSE)
dim(merged)
colnames(merged)[6:7] <- c("Dis1_icd10","Dis1_icd9")

merged <- merge(merged, icds, by.x = "Dis2", by.y="disease_name", all.x = TRUE, all.y=FALSE)
dim(merged)
colnames(merged)[8:9] <- c("Dis2_icd10","Dis2_icd9")
head(merged)
dim(merged)
dim(merged) # 3124    9

write(paste('\nN interactions in the network:\t',dim(merged)[1],sep=""),file=outpath, append=TRUE)


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

# Read and filter barabasi network
barabasi <- read.csv("Network_building/PDN_3_digits.net",header=F, sep="\t",stringsAsFactors = F)
colnames(barabasi) <- c('Dis1','Dis2','Prev1','Prev2','Co-ocurrence','RR','RR_99perc_left_bound','RR_99perc_right_bound','phi_corr','p_value')
barabasi <- barabasi[,c('Dis1','Dis2','RR_99perc_left_bound')]
barabasi <- barabasi[which(barabasi$RR_99perc_left_bound > 1 ), ]
dim(barabasi)
write(paste('Barabasi network initial interactions:\t',dim(barabasi)[1],sep=""),file=outpath, append=TRUE)
# Sort barabasi network
barabasi <- sort_icd_interactions(barabasi)

# Obtain icd9-based network with sorted interactions. That is: icd9 Dis1 < icd9 Dis2
cnetwork <- merged[,c('Dis1_icd9','Dis2_icd9','Distance')]
names(cnetwork) <- c('Dis1','Dis2','Distance')
cnetwork <- sort_icd_interactions(cnetwork)
length(cnetwork$Distance) # 3124

# Remove interactions where icd1 = icd2 ¡¡NEW!!
cnetwork <- cnetwork[!(cnetwork$Dis1 == cnetwork$Dis2), ]
dim(cnetwork)   # 2990
write(paste('Number of interactions between different diseases in the SSN:\t',dim(cnetwork)[1],sep=""),file=outpath, append=TRUE)

# Keep all interactions that are unique (Not consistent, but all interactions)
concat <- paste(cnetwork$Dis1,cnetwork$Dis2,sep="_")
length(concat)             # 2990
length(unique(concat))     # 517
cnetwork$concat <- concat
cnetwork <- cnetwork[!duplicated(cnetwork$concat),]
dim(cnetwork)

# Compute number of interactions
cnetwork <- cnetwork[,c("Dis1","Dis2","Distance")]
dim(cnetwork)
My_interactions <- length(cnetwork$Dis1)  # 517
write(paste('Unique icd9 interactions in the SSN:\t',My_interactions,sep=""),file=outpath, append=TRUE)

#### Obtain overlap with Barabasi

# Quitar las enfermedades que no existen en la otra red
barabasi_icds <- unique(union(barabasi$Dis1,barabasi$Dis2)); length(barabasi_icds)   # 995
cnetwork_icds <- unique(union(cnetwork$Dis1, cnetwork$Dis2)); length(cnetwork_icds)  # 41

common_icds <- intersect(barabasi_icds, cnetwork_icds); length(common_icds)  # 41
cnetwork <- cnetwork[(cnetwork$Dis1 %in% common_icds) & (cnetwork$Dis2 %in% common_icds), ]; dim(cnetwork)  # 517
barabasi <- barabasi[(barabasi$Dis1 %in% common_icds) & (barabasi$Dis2 %in% common_icds), ]; dim(barabasi)  # 331

# Saving the networks
finalnet_filename <- gsub('_overlap.txt',"_final_network.txt",outpath)
finalbarabasi_filename <- gsub('_overlap.txt',"_final_barabasi.txt",outpath)
write.table(cnetwork,file=finalnet_filename,sep="\t",row.names=F, quote = FALSE)
write.table(barabasi,file=finalbarabasi_filename,sep="\t",row.names=F, quote = FALSE)

mygraph <- graph_from_data_frame(cnetwork, directed=FALSE)
barabasigraph <- graph_from_data_frame(barabasi, directed=FALSE)
intersect <- intersection(mygraph,barabasigraph)
gsize(intersect) # Number of edges: 211

write(paste('Interactions in the final network:\t',dim(cnetwork)[1],sep=""),file=outpath, append=TRUE)
write(paste('Interactions in the final barabasi network:\t',dim(barabasi)[1],sep=""),file=outpath, append=TRUE)
write(paste('OVERLAP - Interactions in both networks:\t',gsize(intersect),sep=""),file=outpath, append=TRUE)
precision <- (gsize(intersect)/dim(cnetwork)[1])*100; precision # 43.06122
write(paste('PRECISION: Percentage of overlap of our network:\t',precision,sep=""),file=outpath, append=TRUE)
recall <- (gsize(intersect)/dim(barabasi)[1])*100; recall # 64.13374
write(paste('RECALL: Percentage of overlap of the epidemiology:\t',recall,sep=""),file=outpath, append=TRUE)

### OBTAIN THE P-VALUE
# Generate a network with all my possible interactions
N <- length(common_icds)
connected_net <- data.frame(Dis1=c(),Dis2=c(),Distance=c())
set.seed(5)
for (k1 in 1:(N-1)){
  for (k2 in (k1+1):N){
    if(is.na(common_icds[k2])== T){
      print(k1)
      print(k2)
      break
    }
    connected_net <- rbind(connected_net,data.frame(common_icds[k1],common_icds[k2],1))
  }
}
# Sort the connected network
colnames(connected_net) <- c('Dis1','Dis2','Distance')
connected_net <- sort_icd_interactions(connected_net)
dim(connected_net) # ((41*41)-41)/2
# 820 possible interactions
write(paste('Dimension of the connected network:\t',dim(connected_net)[1],sep=""),file=outpath, append=TRUE)

# Generate 10.000 iterations of: select random network with my nº of interactions and compute overlap with Barabasi
n_iters <- 10000 
random_overlaps <- c()
random_overlaps_second <- c()
set.seed(5) ### ESTABLECER SEMILLA

for(k in 1:n_iters){
  cgraph <- rewire(mygraph, with=keeping_degseq(loops=FALSE, niter = ecount(mygraph)*10))
  # str(cgraph)
  cintersect <- intersection(cgraph,barabasigraph)
  # cintersect
  gsize(cintersect)
  random_overlaps <- append(random_overlaps,gsize(cintersect))
  # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
}
for(k in 1:n_iters){ # it has to go into another loop to keep consistent p-values (reproducibility)
  # p-value for the comparison with the my network
  cgraph_second <- rewire(barabasigraph, with=keeping_degseq(loops=FALSE, niter = ecount(barabasigraph)*10))
  # str(cgraph)
  cintersect_second <- intersection(cgraph_second,mygraph)
  # cintersect
  gsize(cintersect_second)
  random_overlaps_second <- append(random_overlaps_second,gsize(cintersect_second))
  # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
}

pval <- length(which(random_overlaps>=gsize(intersect)))/length(random_overlaps)
pval
write(paste('Comparison with the epidemiology | P-value:\t',pval,sep=""),file=outpath, append=TRUE)
# 0.0234 --> with correct ulcer : 0.0187

pval_second <- length(which(random_overlaps_second>=gsize(intersect)))/length(random_overlaps_second)
pval_second
write(paste('Comparison with our network | P-value:\t',pval_second,sep=""),file=outpath, append=TRUE)
# 0.0264 --> with correct ulcer : 0.0253

hist(random_overlaps, main = "Random overlaps - pairwise spearman sDEGs")
abline(v = gsize(intersect), col="red")


# CONTINUE TO COMPUTE THE OVERLAP BY DISEASE CATEGORY PAIR
# Adding the disease category to the networks over the common set of icd9s
discatsdf <- metadata[, c("icd9","disease_cat")]; 
discatsdf <- discatsdf[!duplicated(discatsdf$icd9), ]
dim(cnetwork)[1]; dim(barabasi)[1]

cnetwork <- merge(cnetwork, discatsdf, by.x = "Dis1", by.y = "icd9", all.x = TRUE, all.y = FALSE)
cnetwork <- merge(cnetwork, discatsdf, by.x = "Dis2", by.y = "icd9", all.x = TRUE, all.y = FALSE)
barabasi <- merge(barabasi, discatsdf, by.x = "Dis1", by.y = "icd9", all.x = TRUE, all.y = FALSE)
barabasi <- merge(barabasi, discatsdf, by.x = "Dis2", by.y = "icd9", all.x = TRUE, all.y = FALSE)
colnames(cnetwork)[4:5] <- c("Dis1_cat","Dis2_cat")
colnames(barabasi)[4:5] <- c("Dis1_cat","Dis2_cat")
dim(cnetwork)[1]; dim(barabasi)[1] # checking that the dimension of the networks does not change. 

discats <- unique(discatsdf$disease_cat) 
ncats <- length(discats)  # 9

# head(icd9_discat)
overlap_df = data.frame()
for (k1 in 1:ncats){
  cat1 = discats[k1]; # cat1 = discats[5]; cat2 = discats[10]; cat1 = discats[2]; cat2 = discats[4]
  for(k2 in k1:ncats){
    cat2 = discats[k2]
    cbarabasi = barabasi[((barabasi$Dis1_cat == cat1 & barabasi$Dis2_cat == cat2) | (barabasi$Dis1_cat == cat2 & barabasi$Dis2_cat == cat1)),  ]
    cssn = cnetwork[((cnetwork$Dis1_cat == cat1 & cnetwork$Dis2_cat == cat2) | (cnetwork$Dis1_cat == cat2 & cnetwork$Dis2_cat == cat1)),  ]
    barabasi_g = graph_from_data_frame(cbarabasi, directed = FALSE)
    ssn_g = graph_from_data_frame(cssn, directed = FALSE)
    intersect = gsize(intersection(barabasi_g, ssn_g))
    overlap_df = rbind(overlap_df, data.frame(cat1, cat2, intersect, gsize(ssn_g), gsize(barabasi_g)))
  }
}

colnames(overlap_df)[4:5] = c('size_dsn', 'size_barabasi')

overlap_df$perc_de_barabasi = (overlap_df$intersect / overlap_df$size_barabasi)*100
overlap_df$perc_de_dsn = (overlap_df$intersect / overlap_df$size_dsn)*100
# head(overlap_df)

filt_overlap_df = overlap_df
abbreviatures = data.frame(dis_cat = discats, abbrev= c("Congenital","Circulatory","Digestive",
                                                        "Musculoskeletal","Nervous","Respiratory",
                                                        "Infectious and parasitic",
                                                        "Mental disorders","Neoplasms"))

graph = graph_from_data_frame(filt_overlap_df, directed=FALSE)
matrix1 = t(as.matrix(get.adjacency(graph,attr = "perc_de_barabasi")))
matrix2 = t(as.matrix(get.adjacency(graph,attr = "perc_de_dsn")))

# Adding column abrreviatures
colnames(matrix1) = abbreviatures$abbrev; rownames(matrix1) = abbreviatures$abbrev
colnames(matrix2) = abbreviatures$abbrev; rownames(matrix2) = abbreviatures$abbrev

m1 = apply(matrix1, 2, round, 1)
m2 = apply(matrix2, 2, round, 1)
matrix1[1:5, 1:2]
m1[1:5, 1:2]

mean1 = apply(matrix1, 1, mean, na.rm=TRUE)
mean2 = apply(matrix2, 1, mean, na.rm=TRUE)
meandf = data.frame(mean1, mean2)
meandf$discat = names(mean1); rownames(meandf) = NULL

meandf_dl <- read.csv("PLOTS/Pathway_counts/mean_precision_recall_dsn.txt",header=T, sep="\t",stringsAsFactors = F)

meandf <- merge(meandf, meandf_dl, by="discat")
colnames(meandf)[2:5] <- c("mean1","mean2","mean1_dl","mean2_dl")
meandf$increment1 <- meandf$mean1 - meandf$mean1_dl
meandf$increment2 <- meandf$mean2 - meandf$mean2_dl

# Getting triangular matrices
library(reshape2)

# matrix1[lower.tri(matrix1)] = NA; m1[lower.tri(m1)] = NA
# matrix2[lower.tri(matrix2)] = NA; m2[lower.tri(m2)] = NA
df1 = melt(matrix1); rounded1 = melt(m1)
df2 = melt(matrix2); rounded2 = melt(m2)



overlap_path <- gsub('_overlap.txt',"_overlap_per_discat.pdf",outpath)
pdf(file=overlap_path, width =6.5, height = 6)

orange = "#C9A624"
green = "#28A894"

ssn_recall = 46.20061
ssn_precision = 43.80403

precision - ssn_precision   # Ha perdido 0.7428055% de precisión
recall - ssn_recall         # Ha ganado 17.93313% de recall

# Respecto a Barabasi: RECALL (df1)
h1 <- ggplot(data=df1, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color='white', na.rm = TRUE)+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white", 
                       midpoint = ssn_recall, space = "Lab",limit=c(0,100),
                       name="Recall")+
  geom_text(data=rounded1, aes(x=Var1, y=Var2, label=value), size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  ggtitle(paste0('SSN Recall by disease category pairs (',round(ssn_recall,2), " %)"), subtitle = "(Percentage of the epimiological network explained by the SSN)") +
  xlab('ICD9 disease categories')+ylab("")
print(h1)

# Respecto a la DSN: PRECISION (df2)
h2 <- ggplot(data=df2, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color='white', na.rm = TRUE)+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white", 
                       midpoint = ssn_precision, space = "Lab",limit=c(0,100),
                       name="Precision")+
  geom_text(data=rounded2, aes(x=Var1, y=Var2, label=value), size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  ggtitle(paste0('SSN Precision by disease category pairs (',round(ssn_precision, 2)," %)"), subtitle = "(Percentage of the SSN described in the epimiological network)") +
  xlab('ICD9 disease categories')+ylab("")
print(h2)


# Respecto a Barabasi + puntos de DSN AL REVÉS
h1 <- ggplot(data=df2, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color='white', na.rm = TRUE)+ # df2 precision COLOR
  theme_minimal()+ 
  geom_point(data=df1, aes(Var1, Var2), size=df1$value/8, color='black', alpha=0.1)+    # df1 recall size
  scale_fill_gradient2(low = orange, high = green, mid = "white", 
                       midpoint = ssn_precision, space = "Lab",limit=c(0,100),
                       name="Precision")+
  # geom_text(data=rounded1, aes(x=Var1, y=Var2, label=value), size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  ggtitle('SSN Precision and Recall by disease category pairs') +
  xlab('ICD9 disease categories')+ylab("")
print(h1)

# All together 
df1$Var1 == df2$Var1; df1$Var2 == df2$Var2
all_df = df1; colnames(all_df)[3] = "recall"; all_df$precision = df2$value

h1 <- ggplot(data=all_df, aes(x=Var1, y=Var2)) + geom_tile(color='lightgrey',fill="white", na.rm = TRUE)+   # df2 precision COLOR
  theme_minimal()+ 
  geom_point(aes(Var1, Var2, size=recall), colour="grey", stroke=0.6)+  # 
  # scale_size_area()+
  geom_point(aes(Var1, Var2, colour=precision, size=recall))+   # , colour=df1$value
  scale_size_area(max_size=14.4)+
  # scale_size_continuous(range = c(0.5, 10))+
  geom_point(aes(Var1, Var2), size=0.05)+
  scale_colour_gradient2(low = orange, high = green, mid = "white", 
                         midpoint = ssn_precision, space = "Lab",limit=c(0,100),
                         name="Precision")+
  # geom_text(data=rounded1, aes(x=Var1, y=Var2, label=value), size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle('SSN Precision and Recall by disease category pairs') +
  xlab('ICD9 disease categories')+ylab("")+ labs(size= "Recall")
print(h1)
overlap_path_triangular <- gsub('_overlap.txt',"_overlap_by_category_pair_triangular_combined.pdf",outpath)
ggsave(overlap_path_triangular, width =7, height = 6.1)



meang1 <- ggplot(meandf, aes(reorder(discat, mean1), mean1, fill=mean1)) + 
  geom_bar(stat = "identity", color="lightgrey", size=0.05)+ylim(0,85)+
  geom_point(aes(reorder(discat,mean1), mean1_dl), shape=18, size=1.3,color="grey30")+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white",
                       midpoint = ssn_recall, space = "Lab",limit=c(0,100),
                       name="Mean recall")+
  geom_text(aes(x=reorder(discat, mean1), y=mean1+2, label=round(mean1,2)), size=2.5)+
  geom_text(aes(x=reorder(discat, increment1), y=2, label=round(increment1,1)), size=2.5, color="grey30")+
  geom_abline(intercept=ssn_recall, slope = 0, colour='darkgrey', size=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))+
  ggtitle('Mean recall') +
  xlab('ICD9 disease categories')+ylab("Mean recall")
meang1

meang2 <- ggplot(meandf, aes(reorder(discat, mean2), mean2, fill=mean2)) + 
  # geom_bar(stat="identity", color="lightgrey", aes(reorder(discat,mean2), mean2_dl), alpha=0.2)+
  geom_bar(stat = "identity", color="lightgrey", size=0.05)+ylim(0,85)+
  geom_point(aes(reorder(discat,mean2), mean2_dl), shape=18, size=1.3,color="grey30")+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white",
                       midpoint = ssn_precision, space = "Lab",limit=c(0,100),
                       name="Mean precision")+
  geom_text(aes(x=reorder(discat, mean2), y=mean2+2, label=round(mean2,2)), size=2.5)+
  geom_text(aes(x=reorder(discat, increment2), y=2, label=round(increment2,1)), size=2.5, color="grey30")+
  geom_abline(intercept=ssn_precision, slope = 0, colour='grey', size=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))+
  ggtitle('Mean precision') +
  xlab('ICD9 disease categories')+ylab("Mean precision")
meang2

meang1 <- ggplot(meandf, aes(reorder(discat, mean1), mean1, fill=mean1)) + 
  geom_bar(stat = "identity", color="lightgrey", size=0.05)+ylim(0,85)+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white",
                       midpoint = recall, space = "Lab",limit=c(0,100),
                       name="Mean recall")+
  geom_text(aes(x=reorder(discat, mean1), y=mean1+2, label=round(mean1,2)), size=2.5)+
  geom_abline(intercept=recall, slope = 0, colour='grey', size=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))+
  ggtitle('Mean recall') +
  xlab('ICD9 disease categories')+ylab("Mean recall")
meang1

meang2 <- ggplot(meandf, aes(reorder(discat, mean2), mean2, fill=mean2)) + 
  geom_bar(stat = "identity", color="lightgrey", size=0.05)+ylim(0,85)+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white",
                       midpoint = precision, space = "Lab",limit=c(0,100),
                       name="Mean precision")+
  geom_text(aes(x=reorder(discat, mean2), y=mean2+2, label=round(mean2,2)), size=2.5)+
  geom_abline(intercept=precision, slope = 0, colour='grey', size=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))+
  ggtitle('Mean precision') +
  xlab('ICD9 disease categories')+ylab("Mean precision")
meang2

dev.off()



