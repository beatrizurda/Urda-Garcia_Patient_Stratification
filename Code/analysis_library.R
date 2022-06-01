###################################################################
#  Functions for the ANALYSIS folder: RNA-seq comorbidities
#
#
#   Beatriz Urda 2022
###################################################################

library(igraph)
# sort_icd_interactions <- function(df){
#   # Sort interactions such that icd of Dis1 < icd of Dis2
#   
#   for (k in 1:length(df$Dis1)){
#     if (df$Dis1[k] > df$Dis2[k]){
#       dis1 <- df$Dis1[k]
#       df$Dis1[k] <- df$Dis2[k]
#       df$Dis2[k] <- dis1
#     }
#   }
#   return(df)
# }

sort_icd_interactions <- function(df, Dis1 = "Dis1", Dis2 = "Dis2"){
  # Sort interactions such that icd of Dis1 < icd of Dis2
  
  for (k in 1:nrow(df)){
    if (df[k, Dis1] > df[k, Dis2]){
      disease1 <- df[k, Dis1]
      df[k, Dis1] <- df[k, Dis2]
      df[k, Dis2] <- disease1
    }
  }
  return(df)
}

get_sorted_significant_barabasi <- function(){
  barabasi <- read.csv("Network_building/PDN_3_digits.net",header=F, sep="\t",stringsAsFactors = F)
  colnames(barabasi) <- c('Dis1','Dis2','Prev1','Prev2','Co-ocurrence','RR','RR_99perc_left_bound','RR_99perc_right_bound','phi_corr','p_value')
  barabasi <- barabasi[,c('Dis1','Dis2','RR_99perc_left_bound')]
  barabasi <- barabasi[which(barabasi$RR_99perc_left_bound > 1 ), ]
  barabasi <- sort_icd_interactions(barabasi)
  return(barabasi)
}

compute_overlap <- function(net1, net2){
  # net1 = net; net2= west # to test
  
  # Sorting the interactions
  net1 = sort_icd_interactions(net1); nrow(net1) # 345
  net2 = sort_icd_interactions(net2); nrow(net2) # 25707
  
  # Getting common names
  names1 = unique(union(net1$Dis1, net1$Dis2)); length(names1)
  names2 = unique(union(net2$Dis1, net2$Dis2)); length(names2)
  
  common_names = intersect(names1, names2); length(common_names) # 40
  
  # Keeping the networks with common names
  net1 = net1[(net1$Dis1 %in% common_names) & (net1$Dis2 %in% common_names), ]; nrow(net1) # 332
  net2 = net2[(net2$Dis1 %in% common_names) & (net2$Dis2 %in% common_names), ]; nrow(net2) # 68
  
  net1$pairs = paste(net1$Dis1, net1$Dis2, sep="_"); length(unique(net1$pairs))
  net2$pairs = paste(net2$Dis1, net2$Dis2, sep="_"); length(unique(net2$pairs))
  net1 = net1[!duplicated(net1$pairs), ]; nrow(net1)
  net2 = net2[!duplicated(net2$pairs), ]; nrow(net2)
  
  length(intersect(net1$pairs, net2$pairs))
  
  # Getting the graphs
  g1 = graph_from_data_frame(net1, directed=FALSE)
  g2 = graph_from_data_frame(net2, directed=FALSE)
  
  # Computing overlapping interactions, recall and precision
  intersect <- intersection(g1,g2)
  overlap = gsize(intersect); print(overlap) # 24
  plot(intersect)
  
  recall <- (gsize(intersect)/gsize(g2))*100 # 53.33%
  print(paste('RECALL (% of interactions from the DSN): ',recall, sep=''))
  precision <- (gsize(intersect)/gsize(g1))*100 # 7.23%
  print(paste('PRECISION (% of network from the other molecular network): ',precision, sep=''))
  
  # Getting the p-value of the Recall and precision
  n_iters <- 10000 
  randomizing_labels <- TRUE
  random_overlaps <- c()
  random_overlaps_second <- c()
  set.seed(5) ### ESTABLECER SEMILLA
  
  # Randomizing labels
  if(randomizing_labels == TRUE){
    for(k in 1:n_iters){
      cgraph <- rewire(g1, with=keeping_degseq(loops=FALSE, niter = ecount(g1)*10))
      # str(cgraph)
      cintersect <- intersection(cgraph,g2)
      # cintersect
      gsize(cintersect)
      random_overlaps <- append(random_overlaps,gsize(cintersect))
      # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
    }
    for(k in 1:n_iters){ # it has to go into another loop to keep consistent p-values (reproducibility)
      # p-value for the comparison with the my network
      cgraph_second <- rewire(g2, with=keeping_degseq(loops=FALSE, niter = ecount(g2)*10))
      # str(cgraph)
      cintersect_second <- intersection(cgraph_second,g1)
      # cintersect
      gsize(cintersect_second)
      random_overlaps_second <- append(random_overlaps_second,gsize(cintersect_second))
      # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
    }
  }
  
  # PRECISION | P-value:
  pval <- length(which(random_overlaps>=gsize(intersect)))/length(random_overlaps)
  print(pval)
  #  0.2544
  
  # RECALL | P-value: 
  pval_second <- length(which(random_overlaps_second>=gsize(intersect)))/length(random_overlaps_second)
  print(pval_second)
  # 0.2644
  return(common_names)
}




