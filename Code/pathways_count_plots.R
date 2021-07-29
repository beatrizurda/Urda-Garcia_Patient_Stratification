#####################################################################################
# PLOT THE PATHWAY COUNTS BEHIND THE DSN: EPIDEMIOLOGICAL AND NOT EPIDEM. INTERACTIONS
#
#  Finding the molecular mechanisms behind comorbidities. WHY these occur. 
# 
#  Beatriz Urda García 2021
######################################################################################

# Importing libraries
library(ggplot2)

setwd("~/Desktop/ANALYSIS")
try(dev.off())

set.seed(5)

metadata =  read.csv2('new_disease_metadata_final_names.txt',stringsAsFactors = F,sep="\t",header=T)
discats = unique(metadata$disease_cat)

parent_level_filename = "PLOTS/Pathway_counts/parent_level_counts.csv"
parent_level = read.csv2(parent_level_filename,stringsAsFactors = F,sep="\t",header=T)

parent_level_unique = "PLOTS/Pathway_counts/parent_level_counts_explained_interactions.csv"
parent_level_unique = read.csv2(parent_level_unique,stringsAsFactors = F,sep="\t",header=T)

# SORT DFs by Pathway Parent
parent_level = parent_level[order(parent_level$Parent), ]
parent_level_unique = parent_level_unique[order(parent_level_unique$Parent), ]

# NORMALIZE INTERACTIONS BY NUMBER OF EPIDEMI. AND NOT EPIDEMIOM. INTERACTIONS
epidem_int = 154 #145
notepidem_int = 197 #185

# Parent level - Average number of pathways UP / DOWN
parent_level$epidem_up_mean = parent_level$epidem_up / epidem_int
parent_level$epidem_down_mean = parent_level$epidem_down / epidem_int

parent_level$not_epidem_up_mean = parent_level$not_epidem_up / notepidem_int
parent_level$not_epidem_down_mean = parent_level$not_epidem_down / notepidem_int

parent_level$ratio_up = parent_level$epidem_up_mean / parent_level$not_epidem_up_mean
parent_level$ratio_down = parent_level$epidem_down_mean / parent_level$not_epidem_down_mean

# Parent unique
parent_level_unique$epidem_up_perc = (parent_level_unique$epidem_up / epidem_int)*100
parent_level_unique$epidem_down_perc = (parent_level_unique$epidem_down / epidem_int)*100

parent_level_unique$not_epidem_up_perc = (parent_level_unique$not_epidem_up / notepidem_int)*100
parent_level_unique$not_epidem_down_perc = (parent_level_unique$not_epidem_down / notepidem_int)*100

library(viridis)
library(ggrepel)

pdf(file="PLOTS/Pathway_counts/pathways_behind_comorbidities_summary_plots.pdf", width = 8, height = 8)

p <- ggplot(parent_level_unique, aes(not_epidem_up_perc, epidem_up_perc)) + 
  theme_minimal()+
  # geom_point() +
  geom_point(aes(size=parent_level$epidem_up_mean, 
                 color=(parent_level$epidem_up_mean / parent_level$not_epidem_up_mean)), 
             # position=position_jitter(h=0.1, w=0.1), 
             alpha=0.7) +
  scale_size_area(max_size=7)+
  # geom_label(label=parent_level_unique$Parent, nudge_x = 0.02, nudge_y = 0.02, 
  #            check_overlap = F) +
  # geom_label(label=parent_level_unique$Parent, nudge_x = 0.1, nudge_y = 0.1) +
  geom_text_repel(aes(label=Parent), size=2.5) +
  # scale_color_viridis(option="D")+
  scale_color_gradient2(midpoint = 1, low = "blue", mid = "lightgrey",
                        high = "red") + # , space = "Lab"
  geom_abline(intercept=0, slope = 1, colour='grey')+
  ggtitle("Common overexpressed pathways") +
  xlab("% interactions not in epidemiology") + ylab("% interactions in epidemiology")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.box = "vertical") + 
  labs(size = "Mean overexpressed pathways in epidem. interactions", color="Ratio of mean overexpressed pathways epidem vs not epidem.")
p

p + theme(legend.position = "right", legend.box = "vertical") 


# In the underexpressed plot, the outlier is: 'Digestion and absorption' with 5.55 Ratio
tail(sort((parent_level$epidem_down_mean / parent_level$not_epidem_down_mean))) # el outlier es 3.84
outlier = subset(parent_level, (epidem_down_mean / not_epidem_down_mean) > 2) # 3.84

p <- ggplot(parent_level_unique, aes(not_epidem_down_perc, epidem_down_perc)) + 
  theme_minimal()+
  # geom_point() +
  geom_point(aes(size=parent_level$epidem_down_mean,
                 color=(parent_level$epidem_down_mean / parent_level$not_epidem_down_mean)), 
             # position=position_jitter(h=0.1, w=0.1),
             alpha=0.7) +
  scale_size_area(max_size=7)+
  geom_text_repel(aes(label=Parent), size=2.5) +
  geom_point(data=subset(parent_level_unique, Parent == 'Digestion and absorption'), color="darkred")+
  geom_text(data=subset(parent_level_unique, Parent == 'Digestion and absorption'),
            label="3.8",vjust=0.54,hjust=1.3, color="darkred", size=2)+  # hjust, distancia hacia la izq; vjust = distancia hacia arriba
  scale_color_gradient2(midpoint = 1, low = "blue", mid = "lightgrey",
                        high = "red", limits=c(0.4585886,1.619154)) + # , space = "Lab" 
  geom_abline(intercept=0, slope = 1, colour='grey')+
  ggtitle("Common underexpressed pathways") +
  xlab("% interactions not in epidemiology") + ylab("% interactions in epidemiology")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.box = "vertical") + 
  labs(size = "Mean underexpressed pathways in epidem. interactions", color="Ratio of mean underexpressed pathways epidem vs not epidem.")
p

p + theme(legend.position = "right", legend.box = "vertical") 



dev.off()

sum(parent_level$epidem_up_mean); sum(parent_level$not_epidem_up_mean)
sum(parent_level$epidem_down_mean); sum(parent_level$not_epidem_down_mean)

# ggplot(parent_level, aes(epidem_up_mean))+geom_histogram(aes(fill="pink", alpha=0.5))+
#   geom_histogram(aes(not_epidem_up_mean, fill="purple", alpha=0.5))
# 
# ggplot(parent_level, aes(epidem_down_mean))+geom_histogram(aes(fill="pink", alpha=0.5))+
#   geom_histogram(aes(not_epidem_down_mean, fill="purple", alpha=0.5))

ggplot(parent_level, aes(epidem_up_mean))+geom_density(aes(fill="pink", alpha=0.5))+
  geom_density(aes(not_epidem_up_mean, fill="purple", alpha=0.5))

ggplot(parent_level, aes(epidem_down_mean))+geom_density(aes(fill="pink", alpha=0.5))+
  geom_density(aes(not_epidem_down_mean, fill="purple", alpha=0.5))

##### Do epidemiological interactions have more enriched common pathways?
pathways_in_epidem = read.csv2("PLOTS/Pathway_counts/pathways_dsn_in_epidemiology.csv",stringsAsFactors = F,sep="\t",header=T)
pathways_not_epidem = read.csv2("PLOTS/Pathway_counts/pathways_dsn_not_in_epidemiology.csv",stringsAsFactors = F,sep="\t",header=T)

mean_up_epidem = sum(pathways_in_epidem$n_common_up)/epidem_int; mean_up_epidem
mean_down_epidem = sum(pathways_in_epidem$n_common_down)/epidem_int; mean_down_epidem

mean_up_not_epidem = sum(pathways_not_epidem$n_common_up)/notepidem_int; mean_up_not_epidem
mean_down_not_epidem = sum(pathways_not_epidem$n_common_down)/notepidem_int; mean_down_not_epidem

# And if we remove neoplasms?
pathways_in_epidem_wo_neoplasms = pathways_in_epidem[(pathways_in_epidem$dis1_category != 'neoplasms') & (pathways_in_epidem$dis2_category != 'neoplasms'), ]
pathways_not_epidem_wo_neoplasms = pathways_not_epidem[(pathways_not_epidem$dis1_category != 'neoplasms') & (pathways_not_epidem$dis2_category != 'neoplasms'), ]

mean_up_epidem = sum(pathways_in_epidem_wo_neoplasms$n_common_up)/epidem_int; mean_up_epidem
mean_down_epidem = sum(pathways_in_epidem_wo_neoplasms$n_common_down)/epidem_int; mean_down_epidem

mean_up_not_epidem = sum(pathways_not_epidem_wo_neoplasms$n_common_up)/notepidem_int; mean_up_not_epidem
mean_down_not_epidem = sum(pathways_not_epidem_wo_neoplasms$n_common_down)/notepidem_int; mean_down_not_epidem

((mean_up_not_epidem + mean_down_not_epidem)/(mean_up_epidem + mean_down_epidem))*100   # 51.75328
(mean_up_not_epidem/mean_up_epidem)*100       # 52.43032 --> 51.40498 (ulcer)  -->  53.14655 (infectious dis)
(mean_down_not_epidem/mean_down_epidem)*100   # 54.17311 --> 53.1418 (ulcer)  --> 56.75026 (infectious dis)


#### FOR EACH PAIR OF ICD9 CATEGORIES
# - Obtain the df with the interactions
# - Plot the #over/underexpressed common pathways

library(stringr)
ncats = length(discats)
nreactome_parents = length(unique(parent_level_unique$Parent))

interact_epidem =  read.csv2("PLOTS/Pathway_counts/pathways_dsn_in_epidemiology.csv",stringsAsFactors = F,sep="\t",header=T)
interact_not_epidem =  read.csv2("PLOTS/Pathway_counts/pathways_dsn_not_in_epidemiology.csv",stringsAsFactors = F,sep="\t",header=T)

pdf(file="PLOTS/Pathway_counts/pathways_behind_comorbidities_per_category_plots.pdf", width = 8.4, height = 8)

discat_pairsdf = data.frame()

discats_fnames = c("Congenital","Circulatory","Digestive","Musculoskeletal","Nervous",
                   "Respiratory","Infectious and parasitic",
                   "Mental disorders","Neoplasms")

for (k1 in 1:ncats){
  cat1 = discats[k1]; # cat1 = discats[5]; cat2 = discats[10]; cat1 = discats[1]; cat2 = discats[6]; cat1 = discats[3]; cat2 = discats[7]
  for(k2 in k1:ncats){
    # print(paste(k1,k2, sep="_"))
    cat2 = discats[k2]
    # print(k2)
    
    # Obtain #interactions in epidem and not epidem in the DSN.
    ninter_epidem = interact_epidem[((interact_epidem$dis1_category == cat1 & interact_epidem$dis2_category == cat2) | (interact_epidem$dis1_category == cat2 & interact_epidem$dis2_category == cat1)),  ]
    ninter_not_epidem = interact_not_epidem[((interact_not_epidem$dis1_category == cat1 & interact_not_epidem$dis2_category == cat2) | (interact_not_epidem$dis1_category == cat2 & interact_not_epidem$dis2_category == cat1)),  ]
    ninter_epidem = dim(ninter_epidem)[1]; ninter_not_epidem = dim(ninter_not_epidem)[1]
    print(data.frame(cat1, cat2, ninter_epidem, ninter_not_epidem))
    discat_pairsdf = rbind(discat_pairsdf, data.frame(cat1, cat2, ninter_epidem, ninter_not_epidem))
    
    ### Obtain df with the #interactions for the current disease category pair (cat1 & cat2) 
    cpattern = paste0(paste(cat1, cat2, sep="_"),"|",paste(cat2, cat1, sep="_"))
    print(cpattern)
    cdf = data.frame()
    for (i in 1:(dim(parent_level_unique)[1])){ # For each Parent Category
      # Counting interactions
      cup_epidem = str_count(parent_level_unique[i,"epidem_up_discat"], pattern=cpattern)
      cdown_epidem = str_count(parent_level_unique[i,"epidem_down_discat"], pattern=cpattern)
      cup_not_epidem = str_count(parent_level_unique[i,"not_epidem_up_discat"], pattern=cpattern)
      cdown_not_epidem = str_count(parent_level_unique[i,"not_epidem_down_discat"], pattern=cpattern)
      
      # Counting number of pathways
      cup_epidem_quant = str_count(parent_level[i,"epidem_up_discat"], pattern=cpattern)
      cdown_epidem_quant = str_count(parent_level[i,"epidem_down_discat"], pattern=cpattern)
      cup_not_epidem_quant = str_count(parent_level[i,"not_epidem_up_discat"], pattern=cpattern)
      cdown_not_epidem_quant = str_count(parent_level[i,"not_epidem_down_discat"], pattern=cpattern)
      
      cdf = rbind(cdf, data.frame(parent_level_unique[i, "Parent"],
                                  cup_epidem, cdown_epidem, cup_not_epidem, cdown_not_epidem,
                                  cup_epidem_quant, cdown_epidem_quant, cup_not_epidem_quant, cdown_not_epidem_quant))
      
    }
    colnames(cdf)[1] <- "Parent"
    head(cdf)
      
      
    ### Normalize the number of interactions
    if(ninter_epidem != 0){     # In Epidemiology
      # unique interactions
      cdf$cup_epidem_perc = (cdf$cup_epidem/ninter_epidem)*100
      cdf$cdown_epidem_perc = (cdf$cdown_epidem/ninter_epidem)*100
      
      # average number of pathways
      cdf$cup_epidem_mean = cdf$cup_epidem_quant/ninter_epidem
      cdf$cdown_epidem_mean = cdf$cdown_epidem_quant/ninter_epidem
    }else{
      # unique interactions
      cdf$cup_epidem_perc = rep(0,nreactome_parents)
      cdf$cdown_epidem_perc = rep(0,nreactome_parents)
      
      cdf$cup_epidem_mean = rep(0,nreactome_parents)
      cdf$cdown_epidem_mean = rep(0,nreactome_parents)
    }
    
    if(ninter_not_epidem != 0){     # Not in Epidemiology
      # unique interactions
      cdf$cup_not_epidem_perc = (cdf$cup_not_epidem/ninter_not_epidem)*100
      cdf$cdown_not_epidem_perc = (cdf$cdown_not_epidem/ninter_not_epidem)*100
      
      # average number of pathways
      cdf$cup_not_epidem_mean = cdf$cup_not_epidem_quant/ninter_not_epidem
      cdf$cdown_not_epidem_mean = cdf$cdown_not_epidem_quant/ninter_not_epidem
    }else{
      # unique interactions
      cdf$cup_not_epidem_perc = rep(0,nreactome_parents)
      cdf$cdown_not_epidem_perc = rep(0,nreactome_parents)
      
      # average number of pathways
      cdf$cup_not_epidem_mean = rep(0,nreactome_parents)
      cdf$cdown_not_epidem_mean = rep(0,nreactome_parents)
    }
    
    cdfup = cdf[!((cdf$cup_not_epidem_perc == 0) & (cdf$cup_epidem_perc == 0)), ]
    cdfdown = cdf[!((cdf$cdown_not_epidem_perc == 0) & (cdf$cdown_epidem_perc == 0)), ]
    
    ### PLOT the df
    # Overexpressed common pathways
    p <- ggplot(cdfup, aes(cup_not_epidem_perc, cup_epidem_perc)) +  # aes(not_epidem_up_perc, epidem_up_perc)
      theme_minimal()+ # xlim(-10,110) + ylim(-10,110)+
      coord_cartesian(xlim=c(-9,109), ylim=c(-9,109))+
      scale_x_continuous(breaks=seq(0,100,25)) + scale_y_continuous(breaks=seq(0,100,25)) +
      # theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))+
      geom_point(size=0.2)+
      geom_point(aes(size=ifelse(cup_epidem_mean==0, 0.3, cup_epidem_mean),   # size=2.2*cup_epidem_mean
                     color=(cup_epidem_mean / cup_not_epidem_mean)), 
                 position=position_jitter(h=3.5, w=3.5, seed = 2),
                 alpha=0.7) +
      scale_size_area(max_size=7)+
      geom_text_repel(aes(label=Parent), position=position_jitter(h=3.5, w=3.5, seed = 2), 
                      max.overlaps =30, size=2.7, segment.size=0.4,
                      max.time = 1, max.iter = 20000, force=1.1) +
      scale_color_gradient2(midpoint = 1, low = "blue", mid = "lightgrey",
                            high = "red") + # , space = "Lab"
      geom_abline(intercept=0, slope = 1, colour='grey')+
      ggtitle(paste(discats_fnames[k1],"-",discats_fnames[k2],sep=" "), subtitle = "Common overexpressed pathways") +
      xlab(paste0("% interactions not in epidemiology (",ninter_not_epidem,")")) + 
      ylab(paste0("% interactions in epidemiology (",ninter_epidem,")"))+
      theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right", legend.box = "vertical", legend.title=element_text(size=52)) + 
      guides(size = guide_legend(order = 1),
             colour = guide_colourbar(order = 2))+
      labs(size = "   ", color="  ")
    print(p)
    
    # Underexpressed common pathways
    p <- ggplot(cdfdown, aes(cdown_not_epidem_perc, cdown_epidem_perc)) +   # aes(not_epidem_down_perc, epidem_down_perc)
      theme_minimal() + #xlim(-10,110) + ylim(-10,110)+
      coord_cartesian(xlim=c(-9,109), ylim=c(-9,109))+
      scale_x_continuous(breaks=seq(0,100,25)) + scale_y_continuous(breaks=seq(0,100,25)) +
      # theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))+
      geom_point(size=0.2)+
      geom_point(aes(size=ifelse(cdown_epidem_mean==0, 0.3, cdown_epidem_mean),
                     color=(cdown_epidem_mean / cdown_not_epidem_mean)), 
                 position=position_jitter(h=3.5, w=3.5, seed = 2),
                 alpha=0.7) +
      scale_size_area(max_size=7)+
      geom_text_repel(aes(label=Parent), position=position_jitter(h=3.5, w=3.5, seed = 2),
                      max.overlaps =30, size=2.7, segment.size=0.4,
                      max.time = 1, max.iter = 20000, force=1.1) +
      scale_color_gradient2(midpoint = 1, low = "blue", mid = "lightgrey",
                            high = "red") + # , space = "Lab"
      geom_abline(intercept=0, slope = 1, colour='grey')+
      ggtitle(paste(discats_fnames[k1],"-",discats_fnames[k2],sep=" "), subtitle = "Common underexpressed pathways") +
      xlab(paste0("% interactions not in epidemiology (",ninter_not_epidem,")")) + 
      ylab(paste0("% interactions in epidemiology (",ninter_epidem,")"))+
      theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right", legend.box = "vertical", legend.title=element_text(size=52)) + 
      # guide_legend(keyheight=6)+
      guides(size = guide_legend(order = 1),
             colour = guide_colourbar(order = 2))+
      labs(size = "   ", color="   ")
    print(p)
    
  }
}

dev.off()

discat_pairsdf = as.data.frame(discat_pairsdf)
colnames(discat_pairsdf) = c("dis1_category", "dis2_category", "dsn_int_epidem", "dsn_not_int_epidem")
discat_pairsdf$dsn_int_epidem = as.numeric(as.character(discat_pairsdf$dsn_int_epidem))
discat_pairsdf$dsn_not_int_epidem = as.numeric(as.character(discat_pairsdf$dsn_not_int_epidem))
sum(discat_pairsdf$dsn_int_epidem)        # 146
sum(discat_pairsdf$dsn_not_int_epidem)    # 184
discat_pairsdf$in_dsn = discat_pairsdf$dsn_int_epidem + discat_pairsdf$dsn_not_int_epidem
discat_pairsdf$precision = discat_pairsdf$dsn_int_epidem / discat_pairsdf$in_dsn

discat_pairsdf <- discat_pairsdf[order(-discat_pairsdf$precision), ]
# g <- ggplot(discat_pairsdf, aes(x=dis1_category, y=dis2_category, fill=precision)) + geom_tile()
# g

library(igraph)
# 'infectious and parasitic diseases' has cero links in our DSN. 
filt_discat_pairsdf = discat_pairsdf[((discat_pairsdf$dis1_category != 'infectious and parasitic diseases') & (discat_pairsdf$dis2_category != 'infectious and parasitic diseases')), ]
graph = graph_from_data_frame(filt_discat_pairsdf, directed=FALSE)
matrix = as.matrix(get.adjacency(graph,attr = "precision"))
matrix[1:5, 1:2]

library(gplots)

g <- heatmap.2(matrix, scale="none",dendrogram = "row",
               cellnote = matrix, notecol="black", notecex = 0.7,
               trace="none",margins=c(18,18),
               cexRow = 0.8,cexCol = 0.8,
               main = "Precision of the DSN")


#### OBTAIN OVERLAP BY DISEASE CATEGORY ####
dev.off()

# Read the comparable Barabasi & DSN (only positive and significant links)
# DISEASE LEVEL NETWORKS
barabasi_filename = "Network_building/Overlapping_results/Shuffling_labels/pairwise_union_spearman_distance_sDEGs_pos_consistent_B_TRUE_final_barabasi.txt"
dsn_filename = "Network_building/Overlapping_results/Shuffling_labels/pairwise_union_spearman_distance_sDEGs_pos_consistent_B_TRUE_final_network.txt"

# ICD9 LEVEL NETWORKS
barabasi_filename = "../Desktop_20200508/Network_building/Overlapping_results/Shuffling_labels/icd9_pairwise_union_spearman_distance_sDEGs_pos_B_TRUE_final_barabasi.txt"
dsn_filename = "../Desktop_20200508/Network_building/Overlapping_results/Shuffling_labels/icd9_pairwise_union_spearman_distance_sDEGs_pos_B_TRUE_final_network.txt"

barabasi =  read.csv2(barabasi_filename,stringsAsFactors = F,sep="\t",header=T)
dsn =  read.csv2(dsn_filename,stringsAsFactors = F,sep="\t",header=T)

if(grepl("icd9", dsn_filename, fixed=TRUE)){
  outputdir = "PLOTS/Pathway_counts/icd9_level/icd9_"
}else{
  outputdir = "PLOTS/Pathway_counts/dis_level/dis_"
}

barabasi$RR_99perc_left_bound = as.numeric(barabasi$RR_99perc_left_bound)
dsn$Distance = as.numeric(dsn$Distance); 

dim(barabasi)     # 331 --> 329
dim(dsn)      # 304     --> 347 # Already contains only positive interactions. 

# Add the disease category
icd9_discat = metadata[, c("icd9", "disease_cat")]; dim(icd9_discat)           # 45 icd9
icd9_discat = icd9_discat[!duplicated(icd9_discat$icd9), ]; dim(icd9_discat)   # 41 icd9

dsn = merge(dsn, icd9_discat, by.x = "Dis1", by.y = "icd9", all.x = FALSE, all.y = FALSE); dim(dsn)
dsn = merge(dsn, icd9_discat, by.x = "Dis2", by.y = "icd9", all.x = FALSE, all.y = FALSE); dim(dsn)
colnames(dsn)[4:5] <- c("Dis1_cat", "Dis2_cat")
dsn = dsn[, c('Dis1', 'Dis2', 'Distance', 'Dis1_cat', 'Dis2_cat')]

barabasi = merge(barabasi, icd9_discat, by.x = "Dis1", by.y = "icd9", all.x = FALSE, all.y = FALSE); dim(barabasi)
barabasi = merge(barabasi, icd9_discat, by.x = "Dis2", by.y = "icd9", all.x = FALSE, all.y = FALSE); dim(barabasi)
colnames(barabasi)[4:5] <- c("Dis1_cat", "Dis2_cat")
barabasi = barabasi[, c('Dis1', 'Dis2', 'RR_99perc_left_bound', 'Dis1_cat', 'Dis2_cat')]

# Remove uler
# dim(dsn) # 347
# dsn = dsn[(dsn$Dis1 != 569) & (dsn$Dis2 != 569), ]; dim(dsn) # 321
# dim(barabasi)
# barabasi = barabasi[(barabasi$Dis1 != 569) & (barabasi$Dis2 != 569), ]; dim(barabasi)


# Compute entire overlap
barabasi_g = graph_from_data_frame(barabasi, directed = FALSE)
dsn_g = graph_from_data_frame(dsn, directed = FALSE)
common_g = intersection(barabasi_g, dsn_g)
intersect = gsize(common_g); intersect
gperc_de_barabasi = (intersect/gsize(barabasi_g))*100; gperc_de_barabasi  # 46.20061
gperc_de_dsn = (intersect/gsize(dsn_g))*100; gperc_de_dsn                 # 43.80403

# PARENTHESIS: COMPUTE CORRELATION between Barabasi's RR and the Distance in our DSN
library("ggpubr")

get_dsn_barabasi_correlation <- function(df,specifics="all data"){
  print(ggscatter(df, x = "RR_99perc_left_bound", y = "Distance", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "Barabasi's RR", ylab = "DSN's correlation",
                  title=paste0("Pearson - ",specifics)))
  print(ggscatter(df, x = "RR_99perc_left_bound", y = "Distance", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "Barabasi's RR", ylab = "DSN's correlation",
                  title=paste0("Spearman - ",specifics)))
}

commondf = as_data_frame(common_g)
cor(commondf$RR_99perc_left_bound, commondf$Distance, method='pearson')         # DL: 0.5369963 | ICD9 L: -0.00954203
cor(commondf$RR_99perc_left_bound, commondf$Distance, method='spearman')        # DL: 0.3040407 | ICD9 L: 0.1997043

metadata_icd9 = read.csv2('../Desktop_20200508/ICD9_RESULTS_GREIN/dl_general_info.txt',stringsAsFactors = F,sep="\t",header=T)
metadata_icd9[metadata_icd9$disease_name == 707, ]$disease_name <- 569
dis_wo_sdegs = metadata_icd9[metadata_icd9$Number_sDEGs == 0, ]$disease_name; length(dis_wo_sdegs); 41-length(dis_wo_sdegs)
commondf_wsdegs = commondf[(!commondf$from %in% dis_wo_sdegs) & (!commondf$to %in% dis_wo_sdegs),  ]

pdf(file=paste0(outputdir,"correlation_barabasi_dsn.pdf"), width = 8, height = 8)
cor.test(commondf$RR_99perc_left_bound, commondf$Distance, method='pearson', na.action = "na.exclude")    # DL: 0.5369963, pv: 4.509e-11 | ICD9 L: -0.00954203, pv:0.9065
cor.test(commondf$RR_99perc_left_bound, commondf$Distance, method='spearman', na.action = "na.exclude")   # DL: 0.3040407, pv: 0.0004628 | ICD9 L: 0.1997043, pv: 0.01314 
get_dsn_barabasi_correlation(commondf)
get_dsn_barabasi_correlation(commondf_wsdegs, specifics = "with sDEGs")

# Without neoplasms
commondf_woneoplasms = commondf[(commondf$Dis1_cat_1 != 'neoplasms') & (commondf$Dis2_cat_1 != 'neoplasms'), ]
get_dsn_barabasi_correlation(commondf_woneoplasms, specifics = "wo neoplasms")

commondf_wsdegs_woneoplasms = commondf_wsdegs[(commondf_wsdegs$Dis1_cat_1 != 'neoplasms') & (commondf_wsdegs$Dis2_cat_1 != 'neoplasms'), ]
get_dsn_barabasi_correlation(commondf_wsdegs_woneoplasms, specifics = "with sDEGs, wo neoplasms")

corr_wo_categories = data.frame()

# Remove outliers 
if(outputdir == 'PLOTS/Pathway_counts/icd9_level/icd9_'){
  # Remove outlier 300 # QUIZÁ este es el bueno!
  commondf = commondf[commondf$RR_99perc_left_bound < 300, ]
  commondf_woneoplasms = commondf_woneoplasms[commondf_woneoplasms$RR_99perc_left_bound < 300, ]
  
  commondf_wsdegs = commondf_wsdegs[commondf_wsdegs$RR_99perc_left_bound < 300, ]
  commondf_wsdegs_woneoplasms = commondf_wsdegs_woneoplasms[commondf_wsdegs_woneoplasms$RR_99perc_left_bound < 300, ]
  
  cor.test(commondf$RR_99perc_left_bound, commondf$Distance, method='pearson', na.action = "na.exclude")   # DL: 0.5369963, pv: 4.509e-11 | ICD9 L: 0.488465, pv:1.5e-10 
  cor.test(commondf$RR_99perc_left_bound, commondf$Distance, method='spearman', na.action = "na.exclude")   # DL: 0.3040407, pv: 0.0004628 | ICD9 L: 0.2112232, pv: 0.008876
  get_dsn_barabasi_correlation(commondf, specifics = "without RR=300 point" )
  get_dsn_barabasi_correlation(commondf_woneoplasms, specifics = "without RR=300 point, wo neoplasms" )
  # w sDEGs
  get_dsn_barabasi_correlation(commondf_wsdegs, specifics = "with sDEGs, without RR=300 point" )
  get_dsn_barabasi_correlation(commondf_wsdegs_woneoplasms, specifics = "with sDEGs, without RR=300 point, wo neoplasms" )
  for (cat in discats){
    ccommondf = commondf[(commondf$Dis1_cat_1 != cat) & (commondf$Dis2_cat_1 != cat), ]
    ccorr1 = cor.test(ccommondf$RR_99perc_left_bound, ccommondf$Distance, method='pearson', na.action = "na.exclude")
    ccommondf_wo_27 = ccommondf[ccommondf$RR_99perc_left_bound < 25, ]
    ccorr2 = cor.test(ccommondf_wo_27$RR_99perc_left_bound, ccommondf_wo_27$Distance, method='pearson', na.action = "na.exclude")
    corr_wo_categories = rbind(corr_wo_categories, data.frame(cat,as.numeric(ccorr1$estimate), ccorr1$p.value, as.numeric(ccorr2$estimate), ccorr2$p.value))
  }
  
  # Remove outlier 27 # SOLO para comprobar qué pasa. 
  commondf = commondf[commondf$RR_99perc_left_bound < 25, ]
  commondf_woneoplasms = commondf_woneoplasms[commondf_woneoplasms$RR_99perc_left_bound < 25, ]
  
  commondf_wsdegs = commondf_wsdegs[commondf_wsdegs$RR_99perc_left_bound < 25, ]
  commondf_wsdegs_woneoplasms = commondf_wsdegs_woneoplasms[commondf_wsdegs_woneoplasms$RR_99perc_left_bound < 25, ]
  
  cor.test(commondf$RR_99perc_left_bound, commondf$Distance, method='pearson', na.action = "na.exclude")    # DL: 0.5369963, pv: 4.509e-11 | ICD9 L: 0.3302223, pv:3.256e-05
  cor.test(commondf$RR_99perc_left_bound, commondf$Distance, method='spearman', na.action = "na.exclude")   # DL: 0.3040407, pv: 0.0004628 | ICD9 L: 0.1955522, pv: 0.01588
  get_dsn_barabasi_correlation(commondf, specifics = "without RR=27 & 300 point" )
  get_dsn_barabasi_correlation(commondf_woneoplasms, specifics = "without RR=27 & 300 point, wo neoplasms" )
  
  get_dsn_barabasi_correlation(commondf_wsdegs, specifics = "with sDEGs, without RR=27 & 300 point" )
  get_dsn_barabasi_correlation(commondf_wsdegs_woneoplasms, specifics = "with sDEGs, without RR=27 & 300 point, wo neoplasms" )
  
}else{
  commondf = commondf[commondf$RR_99perc_left_bound < 20, ]
  commondf_woneoplasms = commondf_woneoplasms[commondf_woneoplasms$RR_99perc_left_bound < 20, ]
  cor.test(commondf$RR_99perc_left_bound, commondf$Distance, method='pearson', na.action = "na.exclude")   # DL: 0.5369963, pv: 4.509e-11 | ICD9 L: 0.488465, pv:1.5e-10 
  cor.test(commondf$RR_99perc_left_bound, commondf$Distance, method='spearman', na.action = "na.exclude")   # DL: 0.3040407, pv: 0.0004628 | ICD9 L: 0.2112232, pv: 0.008876
  get_dsn_barabasi_correlation(commondf, specifics = "without RR=20 point" )
  get_dsn_barabasi_correlation(commondf_woneoplasms, specifics = "without RR=20 point, wo neoplasms" )
  
  commondf = commondf[commondf$RR_99perc_left_bound < 10, ]
  commondf_woneoplasms = commondf_woneoplasms[commondf_woneoplasms$RR_99perc_left_bound < 20, ]
  cor.test(commondf$RR_99perc_left_bound, commondf$Distance, method='pearson', na.action = "na.exclude")   # DL: 0.5369963, pv: 4.509e-11 | ICD9 L: 0.488465, pv:1.5e-10 
  cor.test(commondf$RR_99perc_left_bound, commondf$Distance, method='spearman', na.action = "na.exclude")   # DL: 0.3040407, pv: 0.0004628 | ICD9 L: 0.2112232, pv: 0.008876
  get_dsn_barabasi_correlation(commondf, specifics = "without RR=20 & 10 point" )
  get_dsn_barabasi_correlation(commondf_woneoplasms, specifics = "without RR=20 & 10 point, wo neoplasms" )
}

head(corr_wo_categories)
colnames(corr_wo_categories) = c("Discats", "Corr1", "pvalue1", "Corr2", "pvalue2")
corr_wo_categories$abbrev = c("Congenital","Circulatory","Digestive","Musculoskeletal","Nervous",
                              "Respiratory","Infectious and parasitic",
                              "Mental disorders","Neoplasms")
orange = "#C9A624"
green = "#28A894"

g = ggplot(corr_wo_categories, aes(reorder(abbrev, Corr1), Corr1, fill=Corr1)) + 
  theme_minimal() +
  geom_bar(stat="identity", color="lightgrey", size=0.05)+
  scale_fill_gradient2(low = orange, high = green, mid = "white", 
                       midpoint = 0.49, space = "Lab", #limit=c(0.35, 0.6),
                       name="Correlation")+
  geom_abline(intercept=0.49, slope = 0, colour='grey')+
  # scale_y_continuous(breaks = sort(c(seq(min(corr_wo_categories$Corr1), max(corr_wo_categories$Corr1), length.out=5), 0.49)))+
  scale_y_continuous(breaks = sort(c(seq(0, 0.7, length.out=8), 0.49)))+
  # annotate("text", x=0, y=0.49, label="0.49", size=3)+
  geom_text(aes(reorder(abbrev, Corr1), Corr1+0.01, label=signif(pvalue1,digits=2)), size=3)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("Correlation between the epidemiological networks and the DSN weights") +
  xlab("Disease category not included") + 
  ylab("Pearson correlation")+
  theme(plot.title = element_text(hjust = 0.5))
g

dev.off()

# COMPUTE THE OVERLAP without neoplasms
# cbarabasi = barabasi[(barabasi$Dis1_cat == 'neoplasms') | (barabasi$Dis2_cat == 'neoplasms'), ]
# cdsn = dsn[(dsn$Dis1_cat == 'neoplasms') | (dsn$Dis2_cat == 'neoplasms'), ]
# woneo_barabasi = graph_from_data_frame(cbarabasi, directed = FALSE)
# woneo_dsn = graph_from_data_frame(cdsn, directed = FALSE)
# woneo_intersect = gsize(intersection(woneo_barabasi, woneo_dsn)); woneo_intersect
# woneo_recall = (woneo_intersect/gsize(woneo_barabasi))*100; woneo_recall
# woneo_precision = (woneo_intersect/gsize(woneo_dsn))*100; woneo_precision

# COMPUTE THE OVERLAP BY DISEASE CATEGORY
overlap_by_cat = data.frame()
for (cat in discats){
  cbarabasi = barabasi[(barabasi$Dis1_cat == cat) | (barabasi$Dis2_cat == cat), ]
  cdsn = dsn[(dsn$Dis1_cat == cat) | (dsn$Dis2_cat == cat), ]
  barabasi_g = graph_from_data_frame(cbarabasi, directed = FALSE)
  dsn_g = graph_from_data_frame(cdsn, directed = FALSE)
  intersect = gsize(intersection(barabasi_g, dsn_g))
  overlap_by_cat = rbind(overlap_by_cat, data.frame(cat, intersect, gsize(dsn_g), gsize(barabasi_g)))
}

colnames(overlap_by_cat)[3:4] = c('size_dsn', 'size_barabasi')
overlap_by_cat$perc_de_barabasi = (overlap_by_cat$intersect / overlap_by_cat$size_barabasi)*100
overlap_by_cat$perc_de_dsn = (overlap_by_cat$intersect / overlap_by_cat$size_dsn)*100

# CONTINUE TO COMPUTE THE OVERLAP BY DISEASE CATEGORY PAIR
# head(icd9_discat)
overlap_df = data.frame()
for (k1 in 1:ncats){
  cat1 = discats[k1]; # cat1 = discats[5]; cat2 = discats[10]; cat1 = discats[2]; cat2 = discats[4]
  for(k2 in k1:ncats){
    cat2 = discats[k2]
    cbarabasi = barabasi[((barabasi$Dis1_cat == cat1 & barabasi$Dis2_cat == cat2) | (barabasi$Dis1_cat == cat2 & barabasi$Dis2_cat == cat1)),  ]
    cdsn = dsn[((dsn$Dis1_cat == cat1 & dsn$Dis2_cat == cat2) | (dsn$Dis1_cat == cat2 & dsn$Dis2_cat == cat1)),  ]
    barabasi_g = graph_from_data_frame(cbarabasi, directed = FALSE)
    dsn_g = graph_from_data_frame(cdsn, directed = FALSE)
    intersect = gsize(intersection(barabasi_g, dsn_g))
    overlap_df = rbind(overlap_df, data.frame(cat1, cat2, intersect, gsize(dsn_g), gsize(barabasi_g)))
  }
}

colnames(overlap_df)[4:5] = c('size_dsn', 'size_barabasi')

overlap_df$perc_de_barabasi = (overlap_df$intersect / overlap_df$size_barabasi)*100
overlap_df$perc_de_dsn = (overlap_df$intersect / overlap_df$size_dsn)*100
# head(overlap_df)

# 'infectious and parasitic diseases' has cero links in our DSN. 
if(outputdir == "PLOTS/Pathway_counts/dis_level/dis_"){  # At the disease level, remove infectious and parasitic diseases
  filt_overlap_df = overlap_df[((overlap_df$cat1 != 'infectious and parasitic diseases') & (overlap_df$cat2 != 'infectious and parasitic diseases')), ]
  filt_discat = discats[discats != 'infectious and parasitic diseases']
  abbreviatures = data.frame(dis_cat = filt_discat, abbrev= c("Congenital","Circulatory","Digestive",
                                                          "Musculoskeletal","Nervous","Respiratory",
                                                          "Mental disorders","Neoplasms"))
}else{
  filt_overlap_df = overlap_df
  abbreviatures = data.frame(dis_cat = discats, abbrev= c("Congenital","Circulatory","Digestive",
                                                          "Musculoskeletal","Nervous","Respiratory",
                                                          "Infectious and parasitic",
                                                          "Mental disorders","Neoplasms"))
}

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
write.table(meandf,file="PLOTS/Pathway_counts/mean_precision_recall_dsn.txt",sep="\t",row.names=F, quote = FALSE)

median1 = apply(matrix1, 1, median, na.rm=TRUE)
median2 = apply(matrix2, 1, median, na.rm=TRUE)
mediandf = data.frame(median1, median2)
mediandf$discat = names(median1); rownames(mediandf) = NULL
write.table(meandf,file="PLOTS/Pathway_counts/median_precision_recall_dsn.txt",sep="\t",row.names=F, quote = FALSE)


library(gplots)
pdf(file=paste0(outputdir,"overlap_by_category_pair.pdf"), width = 8, height = 8)

palette.breaks=seq(0,100,length=1000)
color.palette<-colorRampPalette(c("#C81E17","#FFFFFF","#405191"))(length(palette.breaks)-1)

g1 <- heatmap.2(matrix1, scale="none",dendrogram = "none",
                cellnote = m1, notecol="black", notecex = 0.7,
                trace="none",margins=c(22,22),#revC = TRUE,
                cexRow = 0.8,cexCol = 0.8,
                col=color.palette, 
                main = "% Overlap of the Epidemiological Networks")
print(g1)

g2 <- heatmap.2(matrix2, scale="none",dendrogram = "none",
                cellnote = m2, notecol="black", notecex = 0.7,
                trace="none",margins=c(22,22),#revC = TRUE,
                cexRow = 0.8,cexCol = 0.8,
                col=color.palette,
                main = "% Overlap of the DSN")
print(g2)
dev.off()


# Getting triangular matrices
library(reshape2)

# Get triangular matrix
# matrix1[lower.tri(matrix1)] = NA; m1[lower.tri(m1)] = NA
# matrix2[lower.tri(matrix2)] = NA; m2[lower.tri(m2)] = NA

df1 = melt(matrix1); rounded1 = melt(m1)
df2 = melt(matrix2); rounded2 = melt(m2)


pdf(file=paste0(outputdir,"overlap_by_category_pair_triangular.pdf"), width =6.5, height = 6)

# (orange, green): ("#D6341C","#1CD67E")
# orange = "#BD2C13"
# green = "#13BD73"
orange = "#C9A624"
green = "#28A894"


# Respecto a Barabasi: RECALL (df1)
h1 <- ggplot(data=df1, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color='white', na.rm = TRUE)+
        theme_minimal()+ 
        scale_fill_gradient2(low = orange, high = green, mid = "white", 
                       midpoint = gperc_de_barabasi, space = "Lab",limit=c(0,100),
                       name="Recall")+
        geom_text(data=rounded1, aes(x=Var1, y=Var2, label=value), size=2)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
        ggtitle(paste0('DSN Recall by disease category pairs (',round(gperc_de_barabasi,2), " %)"), subtitle = "(Percentage of the epimiological network explained by the DSN)") +
        xlab('ICD9 disease categories')+ylab("")
print(h1)

# Respecto a la DSN: PRECISION (df2)
h2 <- ggplot(data=df2, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color='white', na.rm = TRUE)+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white", 
                       midpoint = gperc_de_dsn, space = "Lab",limit=c(0,100),
                       name="Precision")+
  geom_text(data=rounded2, aes(x=Var1, y=Var2, label=value), size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  ggtitle(paste0('DSN Precision by disease category pairs (',round(gperc_de_dsn, 2)," %)"), subtitle = "(Percentage of the DSN described in the epimiological network)") +
  xlab('ICD9 disease categories')+ylab("")
print(h2)


# Respecto a Barabasi + puntos de DSN
h1 <- ggplot(data=df1, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color='white', na.rm = TRUE)+ # df1 RECALL:COLOR
  theme_minimal()+ 
  geom_point(data=df2, aes(Var1, Var2), size=df2$value/9, color='black', alpha=0.1)+   # df2 size
  scale_fill_gradient2(low = orange, high = green, mid = "white", 
                       midpoint = gperc_de_barabasi, space = "Lab",limit=c(0,100),
                       name="Recall")+
  # geom_text(data=rounded1, aes(x=Var1, y=Var2, label=value), size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  ggtitle('DSN Precision and Recall by disease category pairs') +
  xlab('ICD9 disease categories')+ylab("")
print(h1)

# Respecto a Barabasi + puntos de DSN AL REVÉS
h1 <- ggplot(data=df2, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color='white', na.rm = TRUE)+ # df2 precision COLOR
  theme_minimal()+ 
  geom_point(data=df1, aes(Var1, Var2), size=df1$value/8, color='black', alpha=0.1)+    # df1 recall size
  scale_fill_gradient2(low = orange, high = green, mid = "white", 
                       midpoint = gperc_de_dsn, space = "Lab",limit=c(0,100),
                       name="Precision")+
  # geom_text(data=rounded1, aes(x=Var1, y=Var2, label=value), size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  ggtitle('DSN Precision and Recall by disease category pairs') +
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
  scale_size_area(max_size=11.5)+
  # scale_size_continuous(range = c(0.5, 10))+
  geom_point(aes(Var1, Var2), size=0.05)+
  scale_colour_gradient2(low = orange, high = green, mid = "white", 
                         midpoint = gperc_de_dsn, space = "Lab",limit=c(0,100),
                         name="Precision")+
  # geom_text(data=rounded1, aes(x=Var1, y=Var2, label=value), size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  ggtitle('DSN Precision and Recall by disease category pairs', subtitle = "(% overlap from Barabasi)") +
  xlab('ICD9 disease categories')+ylab("")
print(h1)

with_na_recall = subset(all_df, subset=((all_df$Var1 == 'Congenital' & all_df$Var2 == 'Infectious and parasitic') | (all_df$Var2 == 'Congenital' & all_df$Var1 == 'Infectious and parasitic') | (all_df$Var1 == 'Circulatory' & all_df$Var2 == 'Infectious and parasitic') | (all_df$Var2 == 'Circulatory' & all_df$Var1 == 'Infectious and parasitic') | (all_df$Var1 == 'Infectious and parasitic' & all_df$Var2 == 'Infectious and parasitic')))

h1 <- ggplot(data=all_df, aes(x=Var1, y=Var2)) + geom_tile(color='lightgrey',fill="white", na.rm = TRUE)+   # df2 precision COLOR
  theme_minimal()+ 
  geom_point(aes(Var1, Var2, size=recall), colour="grey", stroke=0.6)+  # 
  # scale_size_area()+
  geom_point(aes(Var1, Var2, colour=precision, size=recall))+   # , colour=df1$value
  # geom_point(data=with_na_recall, shape=95, size=4)+
  scale_size_area(max_size=14.4)+
  # scale_size_continuous(range = c(0.5, 10))+
  geom_point(aes(Var1, Var2), size=0.05)+
  scale_colour_gradient2(low = orange, high = green, mid = "white", 
                         midpoint = gperc_de_dsn, space = "Lab",limit=c(0,100),
                         name="Precision")+
  # geom_text(data=rounded1, aes(x=Var1, y=Var2, label=value), size=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle('DSN Precision and Recall by disease category pairs') +
  xlab('ICD9 disease categories')+ylab("")+ labs(size= "Recall")
print(h1)
ggsave(paste0(outputdir,"overlap_by_category_pair_triangular_combined.pdf"), width =7.03, height = 6.08)

meang1 <- ggplot(meandf, aes(reorder(discat, mean1), mean1, fill=mean1)) + 
  geom_bar(stat = "identity", color="lightgrey", size=0.05)+ylim(0,80)+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white",
                       midpoint = gperc_de_barabasi, space = "Lab",limit=c(0,100),
                       name="Mean recall")+
  geom_text(aes(x=reorder(discat, mean1), y=mean1+2, label=round(mean1,2)), size=2.5)+
  geom_abline(intercept=gperc_de_barabasi, slope = 0, colour='grey', size=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))+
  ggtitle('Mean recall') +
  xlab('ICD9 disease categories')+ylab("Mean recall")
meang1

meang2 <- ggplot(meandf, aes(reorder(discat, mean2), mean2, fill=mean2)) + 
  geom_bar(stat = "identity", color="lightgrey", size=0.05)+ylim(0,90)+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white",
                       midpoint = gperc_de_dsn, space = "Lab",limit=c(0,100),
                       name="Mean precision")+
  geom_text(aes(x=reorder(discat, mean2), y=mean2+2, label=round(mean2,2)), size=2.5)+
  geom_abline(intercept=gperc_de_dsn, slope = 0, colour='grey', size=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))+
  ggtitle('Mean precision') +
  xlab('ICD9 disease categories')+ylab("Mean precision")
meang2

mediang1 <- ggplot(mediandf, aes(reorder(discat, median1), median1, fill=median1)) + 
  geom_bar(stat = "identity", color="lightgrey", size=0.05)+ylim(0,80)+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white",
                       midpoint = gperc_de_barabasi, space = "Lab",limit=c(0,100),
                       name="Median recall")+
  geom_text(aes(x=reorder(discat, median1), y=median1+2, label=round(median1,2)), size=2.5)+
  geom_abline(intercept=gperc_de_barabasi, slope = 0, colour='grey', size=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))+
  ggtitle('Median recall') +
  xlab('ICD9 disease categories')+ylab("Median recall")
mediang1

mediang2 <- ggplot(meandf, aes(reorder(discat, median2), median2, fill=median2)) + 
  geom_bar(stat = "identity", color="lightgrey", size=0.05)+ylim(0,90)+
  theme_minimal()+ 
  scale_fill_gradient2(low = orange, high = green, mid = "white",
                       midpoint = gperc_de_dsn, space = "Lab",limit=c(0,100),
                       name="Median precision")+
  geom_text(aes(x=reorder(discat, median2), y=median2+2, label=round(median2,2)), size=2.5)+
  geom_abline(intercept=gperc_de_dsn, slope = 0, colour='grey', size=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))+
  ggtitle('Median precision') +
  xlab('ICD9 disease categories')+ylab("Median precision")
mediang2

dev.off()



col_to_interpret <- c("Parent", "epidem_up", "epidem_down", "not_epidem_up", "not_epidem_down",
          "epidem_up_mean", "epidem_down_mean", "not_epidem_up_mean", "not_epidem_down_mean",
          "ratio_up", "ratio_down")
df_to_interpret = parent_level[, col_to_interpret]


##### SANKEY PLOT: INSPECTING THE RESULTS OF THE PATHWAY COUNTS AT DISEASE LEVEL
pdf(file=paste0(outputdir,"sankey_diagram.pdf"), width =10, height = 6)

test = parent_level[, c("Parent", "ratio_up","ratio_down", "epidem_up","epidem_down","not_epidem_up",
                        "not_epidem_down","epidem_up_mean","epidem_down_mean","not_epidem_up_mean",
                        "not_epidem_down_mean")]

ggscatter(test, x = "ratio_up", y = "ratio_down", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "Ratio overexpressed", ylab = "Ratio underexpressed",
                title="Correlation between the ratios of over and underexpressed common pathways")


### Sankey PLOT OF THE RATE: # PATHWAYS IN EPIDEM / # PATHWAYS NOT IN EPIDEM
tosankey = test[,1:2]; colnames(tosankey)[2] = "Ratio"; tosankey$Misregul = rep("overxpressed",length(tosankey$Parent))
head(tosankey)
to_add_to_sankey = data.frame("Parent" = test$Parent, "Ratio" = test$ratio_down, 
                              "Misregul"=rep("underexpressed",length(tosankey$Parent)))
head(to_add_to_sankey)
tosankey = rbind(tosankey, to_add_to_sankey)


library(ggalluvial)

parent_colors = readRDS(file="Reactome/parent_colors.rds")
pie(rep(1,length(parent_colors$Parent)), labels=parent_colors$Parent, col=as.character(parent_colors$Color))
dim(tosankey); dim(parent_colors)
tosankey = merge(tosankey, parent_colors, by= 'Parent', all.x=TRUE, all.y = FALSE)
parent_labels = unique(tosankey$Parent)
tosankey$Rounded_ratio = round(tosankey$Ratio, 4)


ggplot(tosankey, aes(x=Misregul, stratum=reorder(Ratio, -Ratio), alluvium=Parent, fill=Color, label=Parent))+   #, alpha=0.5
  theme_classic()+ 
  scale_fill_identity()+
  geom_flow(stat = "alluvium", lode.guidance = "rightleft", color = "darkgray", na.rm=TRUE) +
  geom_stratum(size=0.2, na.rm=TRUE) +
  # scale_linetype_manual(values = "blank") +
  # scale_linetype_manual(values = c("blank", "solid")) +scale_color_manual(values=c("white", "white"))+
  ggrepel::geom_text_repel(
    aes(label = ifelse(Misregul == "overxpressed", as.character(Parent), NA)),
    stat = "alluvium", size = 2.6, direction = "y", nudge_x = -1, segment.alpha=0
  ) +

  # geom_text(stat="alluvium", aes(label=after_stat(alluvium)), size=2.5)+
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=2.6)+
  # geom_text(aes(y=c(c(1:27),c(1:27))-0.5,label=Parent), size=2.5)+
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+   # do not plot y axis
  xlab('Misregulated pathways')+ylab("Reactome pathathway categories")+
  ggtitle("Pathway enrichment in epidemiological interactions")



ggplot(tosankey, aes(x=Misregul, stratum=reorder(Rounded_ratio, -Rounded_ratio), alluvium=Parent, fill=Color, label=Parent))+   #, alpha=0.5
  theme_classic()+ 
  scale_fill_identity()+
  geom_flow(stat = "alluvium", lode.guidance = "rightleft", color = "darkgray", na.rm=TRUE) +
  geom_stratum(size=0.2, na.rm=TRUE) +
  # scale_linetype_manual(values = "blank") +
  # scale_linetype_manual(values = c("blank", "solid")) +scale_color_manual(values=c("white", "white"))+
  ggrepel::geom_text_repel(
    aes(label = ifelse(Misregul == "overxpressed", as.character(Parent), NA)),
    stat = "alluvium", size = 2.6, direction = "y", nudge_x = -1, segment.alpha=0
  ) +
  
  # geom_text(stat="alluvium", aes(label=after_stat(alluvium)), size=2.5)+
  geom_text(stat="stratum", aes(label=after_stat(stratum)), size=2.6)+
  # geom_text(aes(y=c(c(1:27),c(1:27))-0.5,label=Parent), size=2.5)+
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+   # do not plot y axis
  xlab('Misregulated pathways')+ylab("Reactome pathathway categories")+
  ggtitle("Pathway enrichment in epidemiological interactions")

dev.off()

