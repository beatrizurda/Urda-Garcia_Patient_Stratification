
setwd("~/Desktop/ANALYSIS/Metapatients/with_entire_count_matrix/final_clusters/")

hey <- readRDS("BreastCancer_combat_counts.rds")

cl1 <- names(hey[hey == 1])
cl2 <- names(hey[hey == 2])
cl3 <- names(hey[hey == 3])

length(cl1)
length(cl2)
length(cl3)

cl1
cl2
cl3

