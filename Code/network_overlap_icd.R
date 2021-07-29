#####################################################################################
# COMPUTING NETWORK OVERLAP AT THE ICD9 LEVEL
#
#   
#  
#  
# 
# Beatriz Urda García 2020
######################################################################################


library(igraph)
library(ggplot2)

# COMMAND-LINE ARGUMENTS
args = commandArgs(trailingOnly=TRUE)

# Network filename
filename <- args[1] 

# Universe option
universe_option <- args[2]
randomizing_labels <-  args[3]


setwd("~/Desktop/Desktop_20200508/")

# Read the networks
# MANUAL ARGUMENTS
# filename <- 'pairwise_spearman_sDEGs_pos_network.txt'        # PEARWISE SPEARMAN sDEGs
# filename <- 'icd9_pairwise_spearman_distance_sDEGs_pos_network.txt'        # PEARWISE SPEARMAN sDEGs
# # # # filename <- 'union_spearman_sDEGs_pos_network.txt'           # SPEARMAN UNION of sDEGs for ALL DISEASES
# # # # filename <- 'pairwise_union_spearman_sDEGs_pos_network.txt'    # SPEARMAN UNION of sDEGs for ALL DISEASES
# filename <- 'icd9_pairwise_union_spearman_distance_sDEGs_pos_network.txt'
# universe_option <- "C"
# randomizing_labels <-  TRUE

# SCRIPT BEGINS
filedir <- 'Network_building/Defined_networks/'
filepath <- paste(filedir,filename,sep="")
if(randomizing_labels == TRUE){
  outpath <- paste('Network_building/Overlapping_results/Shuffling_labels/',gsub('_network.txt',"",filename),'_',universe_option,'_',randomizing_labels,'_overlap.txt',sep="")
}else{
  outpath <- paste('Network_building/Overlapping_results/Sampling_randomly/',gsub('_network.txt',"",filename),'_',universe_option,'_',randomizing_labels,'_overlap.txt',sep="")
}
# outpath <- paste('Network_building/Overlapping_results/',gsub('_network.txt',"",filename),'_',universe_option,'_overlap.txt',sep="")
network <- read.csv(filepath,header=T, sep="\t",stringsAsFactors = F)

# Save the network with the correct ICD9 for ULCER
if((707 %in% network$Dis1) | (707 %in% network$Dis2)){
  if(dim(network[network$Dis1 == 707, ])[1] > 0){
    network[network$Dis1 == 707, ]$Dis1 <- 569
  }
  if(dim(network[network$Dis2 == 707, ])[1] > 0){
    network[network$Dis2 == 707, ]$Dis2 <- 569
  }
  
  write.table(network,file=filepath,sep="\t",row.names=F, quote = FALSE)
  
}

# Change the ICD9 of Ulcer: from 707 to 569
dismeta <- read.csv("ICD9_RESULTS_GREIN/dl_general_info.txt",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
dismeta[dismeta$disease_name == 707, ]$disease_name <- '569'


filecon <- file(outpath)
writeLines(paste('OVERLAP ANALYSIS OF THE NETWORK:\t',filename,sep=""),filecon)
close(filecon)
write(paste('\n','Universe option: ',universe_option,sep=""),file=outpath, append=TRUE)

write(paste('\nN interactions in the network:\t',dim(network)[1],sep=""),file=outpath, append=TRUE)


if(length(grep("sDEGs",filename)) == 1){  # If the network is based on DExpression
  dis_wo_sDEGs <- dismeta[dismeta$Number_sDEGs == 0,]
  dis_wo_sDEGs <- dis_wo_sDEGs$disease_name
}else{     # If the network is based on DExpression
  # UNCOMMENT AFTER HAVING P-VALUES AND CREATE THE FILE
  # number_sDVGs <- read.csv("Number_over_under_sDVGs.txt",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
  # number_sDVGs <- number_sDVGs[, c(1,2)]; colnames(number_sDVGs) <- c("disease_name", "Number_sDVGs")
  # dis_wo_sDVGs <- number_sDVGs[number_sDVGs$Number_sDVGs == 0,]$disease_name
  # dismeta <- merge(dismeta, number_sDVGs, by = "disease_name", all.x = TRUE, all.y = FALSE)
}


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
cnetwork <- network[,c('Dis1','Dis2','Distance')]
cnetwork <- sort_icd_interactions(cnetwork)
length(cnetwork$Distance) # 241

# Remove interactions where icd1 = icd2 ¡¡NEW!! No hay pero por si!
dim(cnetwork)
cnetwork <- cnetwork[!(cnetwork$Dis1 == cnetwork$Dis2), ]
dim(cnetwork)
head(cnetwork)

# Compute number of interactions
length(cnetwork$Dis1) # 167
My_interactions <- length(cnetwork$Dis1)
write(paste('Unique icd9 interactions:\t',My_interactions,sep=""),file=outpath, append=TRUE)

# Obtain overlap with Barabasi
# my_icds <- union()
barabasi_icds <- unique(union(barabasi$Dis1,barabasi$Dis2))
length(barabasi_icds)

# Quitar las enfermedades que no existen en la otra red -- Remove from Barabasi the interactions that contain a disease that I don't have
icds_mine <- unique(union(cnetwork$Dis1,cnetwork$Dis2)); length(icds_mine)
icds_barabasi <- unique(union(barabasi$Dis1,barabasi$Dis2)); length(icds_barabasi)
common_icds <- intersect(icds_mine,icds_barabasi); length(common_icds)
barabasi <- barabasi[(barabasi$Dis1 %in% common_icds) & (barabasi$Dis2 %in% common_icds), ] # Filtering Barabasi's network
cnetwork <- cnetwork[(cnetwork$Dis1 %in% common_icds) & (cnetwork$Dis2 %in% common_icds), ] # Filtering my network 
dim(barabasi)
dim(cnetwork)

print("Before universe")
if(universe_option == "B"){
  # OPTION B. Considering the diseases in my network currently - The common set of icds
  write(paste('icd9 codes in our network (option B):\t',length(icds_mine),sep=""),file=outpath, append=TRUE)
}else if(universe_option == "C"){
  # OPTION C. Considering the diseases in my network with at least 3sDEGs
  min_sDEGs <- 1
  
  if(length(grep("sDEGs",filename)) == 1){    # If the network is based on DExpression
    power_dis <- dismeta[dismeta$Number_sDEGs >= min_sDEGs, ] ; power_dis <- unique(power_dis$disease_name)
  }else{    # If the network is based on DExpression
    power_dis <- dismeta[dismeta$Number_sDVGs >= min_sDEGs, ] ; power_dis <- unique(power_dis$disease_name)
  }
  
  length(power_dis) # 36 always --> 28 sin duplicates
  all_my_icds <-power_dis
  write(paste('icd9 codes in our network (option C):\t',length(power_dis),sep=""),file=outpath, append=TRUE)
  num_all_my_icds <- as.numeric(all_my_icds) ; num_all_my_icds <- unique(num_all_my_icds)
  cnetwork <- cnetwork[(cnetwork$Dis1 %in% num_all_my_icds) & (cnetwork$Dis2 %in% num_all_my_icds), ] # Filtering my network
  barabasi <- barabasi[(barabasi$Dis1 %in% num_all_my_icds) & (barabasi$Dis2 %in% num_all_my_icds), ] # Filtering Barabasi's network
  My_interactions <- length(cnetwork$Dis1)
}
print("After universe")

icds_mine <- unique(union(cnetwork$Dis1,cnetwork$Dis2)); length(icds_mine)
icds_barabasi <- unique(union(barabasi$Dis1,barabasi$Dis2)); length(icds_barabasi)
common_icds <- intersect(icds_mine,icds_barabasi); length(common_icds)
write(paste('common icd9 codes:\t',length(common_icds),sep=""),file=outpath, append=TRUE)
write(paste('Expected number of interactions in the connected network:\t',((length(common_icds)*length(common_icds))-length(common_icds))/2,sep=""),file=outpath, append=TRUE)

dim(barabasi) # 160
write(paste('Barabasi interaction within our icd9 codes:\t',dim(barabasi)[1],sep=""),file=outpath, append=TRUE)

cnetwork$Dis1 <- as.numeric(cnetwork$Dis1)
cnetwork$Dis2 <- as.numeric(cnetwork$Dis2)

# Saving the networks
finalnet_filename <- gsub('_overlap.txt',"_final_network.txt",outpath)
finalbarabasi_filename <- gsub('_overlap.txt',"_final_barabasi.txt",outpath)
write.table(cnetwork,file=finalnet_filename,sep="\t",row.names=F, quote = FALSE)
write.table(barabasi,file=finalbarabasi_filename,sep="\t",row.names=F, quote = FALSE)

mygraph <- graph_from_data_frame(cnetwork, directed=FALSE)
barabasigraph <- graph_from_data_frame(barabasi, directed=FALSE)
intersect <- intersection(mygraph,barabasigraph)
gsize(intersect) # Number of edges
write(paste('OVERLAP - Interactions in both networks:\t',gsize(intersect),sep=""),file=outpath, append=TRUE)
perc_overlap <- (gsize(intersect)/dim(cnetwork)[1])*100
write(paste('Percentage of overlap:\t',perc_overlap,sep=""),file=outpath, append=TRUE)
# 64 interactions in the olverlap
perc_overlap_from_epidemiology <- (gsize(intersect)/dim(barabasi)[1])*100
write(paste('Percentage of overlap of the epidemiology:\t',perc_overlap_from_epidemiology,sep=""),file=outpath, append=TRUE)

# Generate a network with all my possible interactions (52 diseases)
N <- length(common_icds) # 31
# n_possible_interactions <- sum(range(N-1)) # 60

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
dim(connected_net) # ((31*31) - 31)/2
# 465 possible interactions
write(paste('Dimension of the connected network:\t',dim(connected_net)[1],sep=""),file=outpath, append=TRUE)

# Generate 10.000 iterations of: select random network with my nº of interactions and compute overlap with Barabasi
n_iters <- 10000 
random_overlaps <- c()
random_overlaps_second <- c()
set.seed(5) ### ESTABLECER SEMILLA

if(randomizing_labels == FALSE){
  for(k in 1:n_iters){
    crandom_indexes <- sample(c(1:length(connected_net$Dis1)),My_interactions)
    crandom_net <- connected_net[crandom_indexes,]
    cgraph <- graph_from_data_frame(crandom_net, directed=FALSE)
    cintersect <- intersection(cgraph,barabasigraph)
    gsize(cintersect) # Number of edges
    random_overlaps <- append(random_overlaps,gsize(cintersect))
  }
}else{ # Randomizing labels
  print("Shuffling labels")
  for(k in 1:n_iters){
    # p-value for the comparison with the epidemiology
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
  
  
}

pval <- length(which(random_overlaps>=gsize(intersect)))/length(random_overlaps)
pval
write(paste('Comparison with epidemiology | P-value:\t',pval,sep=""),file=outpath, append=TRUE)
# 0.1076

pval_second <- length(which(random_overlaps_second>=gsize(intersect)))/length(random_overlaps_second)
pval_second
write(paste('Comparison with our network | P-value:\t',pval_second,sep=""),file=outpath, append=TRUE)

# TP, FP, FN, TN
write(paste('Total Number of interactions:\t',dim(connected_net)[1],sep=""),file=outpath, append=TRUE)

TP <- gsize(intersect)
FP <- gsize(mygraph) - gsize(intersect)
FP_2 <- gsize(difference(mygraph, barabasigraph))
FN <- gsize(barabasigraph) - gsize(intersect)
FN_2 <- gsize(difference(barabasigraph, mygraph))

connected_net_graph <- graph_from_data_frame(connected_net, directed=FALSE)
cnetwork_union_barabasi <- union(mygraph, barabasigraph, byname = TRUE)
print(cnetwork_union_barabasi); print(gsize(cnetwork_union_barabasi)); print(is_simple(cnetwork_union_barabasi)) 
TN_graph <- difference(connected_net_graph, cnetwork_union_barabasi, byname= TRUE)
TN <- gsize(TN_graph)


if(randomizing_labels == TRUE){
  hist_outpath <- paste('Network_building/Overlapping_results/Shuffling_labels/histograms/',gsub('_network.txt',"",filename),'_',universe_option,'_',randomizing_labels,'_hist.pdf',sep="")
}else{
  hist_outpath <- paste('Network_building/Overlapping_results/Sampling_randomly/histograms/',gsub('_network.txt',"",filename),'_',universe_option,'_',randomizing_labels,'_hist.pdf',sep="")
}
pdf(file=hist_outpath)
hist(random_overlaps, main = "Random overlaps - pairwise spearman sDEGs")
abline(v = gsize(intersect), col="red")
dev.off()

print("which_loop(mygraph)")
print(is_simple(mygraph)) # No loops or multiple edges between nodes
# print(which_loop(mygraph))
# print(count_multiple(mygraph))
print("which_loop(barabasigraph)")
print(is_simple(barabasigraph)) # No loops or multiple edges between nodes
# print(which_loop(barabasigraph))
# print(count_multiple(barabasigraph))
print("which_loop(connected_net_graph)")
print(is_simple(connected_net_graph)) # No loops or multiple edges between nodes


