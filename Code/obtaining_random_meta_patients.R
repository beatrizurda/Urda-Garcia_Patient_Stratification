#####################################################################################
# OBTAINING RANDOM META-PATIENTS
#
#   We use PAM algorithm to find clusters for each disease based on their
#   normalized (and if necessary) batch effect removal counts
# 
# Beatriz Urda García 2020
######################################################################################

library("randomizr")
setwd("~/Desktop/ANALYSIS")

##### Create runs for the building the network
# runs_iters <- c(501:1000)
# runs <- paste("python","build_metapatient_dis_network_randomization.py",runs_iters, sep="\t")
# write.table(runs,file="Randomization/Network_building/jobs_mare/runs_build_network_random_metapatients_501_1000.txt",
#             row.names=F, quote=FALSE, col.names=FALSE)


# Disease metadata
metadata =  read.csv2('new_disease_metadata_final_names.txt',stringsAsFactors = F,sep="\t",header=T)

# Number of clusters for each disease
nclustersdf =  read.csv2("Metapatients/with_entire_count_matrix/Summary_best_silhouette_metapatients.txt",stringsAsFactors = F,sep="\t",header=T)

# Sample metadata (e.g., patients, controls)
sample_meta = readRDS("new_meta_jon_grein_good_quality.rds")

# Int list with the number of the cluster and names with the sample name. 
example = readRDS("Metapatients/with_entire_count_matrix/final_clusters/AdenomatousPolyps_normalized_counts.rds")

combat_files <- list.files(path = "DL_RESULTS_GREIN/after_combat_counts/")
combat_diseases <- gsub("_combat_counts.rds","",combat_files)

dislist = nclustersdf$Disease

clusterpath <- "Metapatients/with_entire_count_matrix/Randomization/final_clusters/"

# Create a folder by disease  for the final clusters
for (dis in dislist){
  dir.create(paste0(clusterpath,dis))
}

# Create a folder by disease  for the DEA results 
for (dis in dislist){
  dir.create(paste0("Metapatients/with_entire_count_matrix/Randomization/DEA_results/",dis))
}


nramdomizations = 1000

# Create runs for the DEA
runs <- paste("Rscript","DEA_random_metapatients.R","1000",dislist, sep="\t")
write.table(runs,file="Metapatients/with_entire_count_matrix/Randomization/jobs_mare/runs_DEA_random_metapatients.txt",
            row.names=F, quote=FALSE, col.names=FALSE)


for (dis in dislist){
  # Generate Random meta-patient
  # dis = dislist[1]; dis   # To comment
  nclusters = nclustersdf[nclustersdf$Disease == dis, ]$Best_silhouette; nclusters
  
  if(dis %in% combat_diseases){
    cclusters = readRDS(paste0("Metapatients/with_entire_count_matrix/final_clusters/",dis,"_combat_counts",".rds"))
  }else{
    cclusters = readRDS(paste0("Metapatients/with_entire_count_matrix/final_clusters/",dis,"_normalized_counts",".rds"))
  }
  set.seed(5)
  sample_names = names(cclusters); sample_names
  nsamples = length(sample_names)
  ctable = as.data.frame(table(cclusters)); #print(ctable)
  # ctable=data.frame(cclusters=c(1,2), Freq=c(5,5)); nclusters=2  # To test the randomization
  
  for(niter in 1:nramdomizations){
    randomized_order = sample(c(1:nsamples),nsamples, replace=FALSE); randomized_order  # SET THE SEED!
    # print(randomized_order)
    
    randomized_samples = c()
    index=1
    for(k in 1:nclusters){
      # print(" ")
      # print(k)
      # print(randomized_order[c(index:(index+ctable[k,2]-1))])
      cindexes = randomized_order[c(index:(index+ctable[k,2]-1))]
      cvector = rep(k, length(cindexes)); names(cvector) <- sample_names[cindexes]
      randomized_samples <- append(randomized_samples, cvector)
      # radomized_samples = append(radomized_samples, randomized_order[c(index:ctable[k,2])])
      index=index+ctable[k,2]
    }
    # print(randomized_samples)
    saveRDS(randomized_samples, file=paste0(clusterpath,dis,"/","clusters","_",niter))
    
  }
}


