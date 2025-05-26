#########################################################################################
# Exploring the molecular mechanisms under the negative link between 
# Huntington disease and cancer
#
# Beatriz Urda García 2025
#########################################################################################

library(ggplot2)
library(tidyr) # Reshape to long format
library(dplyr)
library(ggplot2)
library(scales)


# DATA LOADING
meta <- read.csv("shiny_network_app_pnas/new_disease_metadata_final_names.txt",header=T, sep="\t",stringsAsFactors = F); dim(meta)
discats <- unique(meta$disease_cat)

dis_colors <- meta[, c("final_disease_name","new_dis_cat_colors","disease_cat")]; colnames(dis_colors) = c("Dis", "Color","Dis_category")

discat_abbv <- data.frame(discats = sort(unique(meta$disease_cat)), 
                          abbv = c("congenital","circulatory","digestive","musculoskeletal",
                                   "nervous","respiratory","infectious","mentaldisorders","neoplasms"))

reactome_cats <- read.csv("shiny_network_app_pnas/additional_data/reactome_pathway_categories.txt",header=F, sep="\t",stringsAsFactors = F)
reactome_cats <- reactome_cats$V1; dim(reactome_cats)
reactome_cats <- reactome_cats[!reactome_cats == "Digestion and absorption"]   # It is not present. 

diseases_dic <- readRDS("shiny_network_app_pnas/disease_dic.rds")
extended_node_names <- read.csv("shiny_network_app_pnas/additional_data/extended_node_names.txt",header=T, sep="\t",stringsAsFactors = F); dim(extended_node_names)

genesdf <- read.csv("shiny_network_app_pnas/additional_data/gene_list_df.txt",header=T, sep="\t",stringsAsFactors = F); dim(genesdf) # Contains ensemble - gene symbol correspondance. 
pathwaysdf <- read.csv("shiny_network_app_pnas/additional_data/pcats.txt",header=T, sep="\t",stringsAsFactors = F); dim(genesdf)
all_pathways = sort(readRDS("shiny_network_app_pnas/additional_data/all_pathways.rds"))

dsn_dic = readRDS("shiny_network_app_pnas/additional_data/DSN_genes_pathways_dics.rds")
ssn_dic = readRDS("shiny_network_app_pnas/additional_data/SSN_genes_pathways_dics.rds")

# HD
HD_up = diseases_dic[["Huntington's disease"]]@genes_up; length(HD_up)
HD_down = diseases_dic[["Huntington's disease"]]@genes_down; length(HD_down)

HD_pup = diseases_dic[["Huntington's disease"]]@pathways_up; length(HD_pup)
HD_pdown = diseases_dic[["Huntington's disease"]]@pathways_down; length(HD_pdown)

# Cancer
cancer_list = c("Lung cancer", "Chronic lymphocytic leukemia", "Breast cancer", "Liver cancer", "Hyperplastic polyps", "Sessile serrated polyps", "Papillary thyroid cancer")
cancer_list = c("Liver cancer","Lung cancer", "Breast cancer", "Chronic lymphocytic leukemia")
ccancer = "Lung cancer"
C_up = diseases_dic[[ccancer]]@genes_up; length(C_up)
C_down = diseases_dic[[ccancer]]@genes_down; length(C_down)

C_pup = diseases_dic[[ccancer]]@pathways_up; length(C_down)
C_pdown = diseases_dic[[ccancer]]@pathways_down; length(C_pdown)

# Genes
up_down=intersect(HD_up,C_down); length(up_down) # 377
down_up=intersect(HD_down, C_up); length(down_up) # 231

# Pathways
pup_pdown=intersect(HD_pup,C_pdown); length(pup_pdown) # 79
pdown_pup=intersect(HD_pdown, C_pup); length(pdown_pup) # 17


# Initialize results data frame
results <- data.frame(
  Disease = character(),
  
  # Raw counts
  Gene_up_down = integer(),
  Gene_down_up = integer(),
  Total_gene_overlap = integer(),
  Path_up_down = integer(),
  Path_down_up = integer(),
  Total_pathway_overlap = integer(),
  
  # Percentages relative to HD
  Perc_gene_up_down_HD = numeric(),
  Perc_gene_down_up_HD = numeric(),
  Perc_path_up_down_HD = numeric(),
  Perc_path_down_up_HD = numeric(),
  
  # Percentages relative to cancer
  Perc_gene_up_down_Cancer = numeric(),
  Perc_gene_down_up_Cancer = numeric(),
  Perc_path_up_down_Cancer = numeric(),
  Perc_path_down_up_Cancer = numeric(),
  
  # Totals as % of HD and cancer
  Perc_total_genes_HD = numeric(),
  Perc_total_genes_Cancer = numeric(),
  Perc_total_pathways_HD = numeric(),
  Perc_total_pathways_Cancer = numeric(),
  
  stringsAsFactors = FALSE
)

# Loop through each cancer
for (ccancer in cancer_list) {
  # Get cancer gene/pathway lists
  C_up <- diseases_dic[[ccancer]]@genes_up
  C_down <- diseases_dic[[ccancer]]@genes_down
  C_pup <- diseases_dic[[ccancer]]@pathways_up
  C_pdown <- diseases_dic[[ccancer]]@pathways_down
  
  # Intersections
  up_down <- intersect(HD_up, C_down)
  down_up <- intersect(HD_down, C_up)
  pup_pdown <- intersect(HD_pup, C_pdown)
  pdown_pup <- intersect(HD_pdown, C_pup)
  
  # Counts
  n_up_down <- length(up_down)
  n_down_up <- length(down_up)
  n_pup_pdown <- length(pup_pdown)
  n_pdown_pup <- length(pdown_pup)
  
  total_gene_overlap <- n_up_down + n_down_up
  total_pathway_overlap <- n_pup_pdown + n_pdown_pup
  
  # Total counts
  total_HD_genes <- length(HD_up) + length(HD_down)
  total_HD_pathways <- length(HD_pup) + length(HD_pdown)
  total_Cancer_genes <- length(C_up) + length(C_down)
  total_Cancer_pathways <- length(C_pup) + length(C_pdown)
  
  # Percentages relative to HD
  perc_up_down_HD <- ifelse(length(HD_up) > 0, n_up_down / length(HD_up) * 100, NA)
  perc_down_up_HD <- ifelse(length(HD_down) > 0, n_down_up / length(HD_down) * 100, NA)
  perc_pup_pdown_HD <- ifelse(length(HD_pup) > 0, n_pup_pdown / length(HD_pup) * 100, NA)
  perc_pdown_pup_HD <- ifelse(length(HD_pdown) > 0, n_pdown_pup / length(HD_pdown) * 100, NA)
  
  # Percentages relative to cancer
  perc_up_down_Cancer <- ifelse(length(C_down) > 0, n_up_down / length(C_down) * 100, NA)
  perc_down_up_Cancer <- ifelse(length(C_up) > 0, n_down_up / length(C_up) * 100, NA)
  perc_pup_pdown_Cancer <- ifelse(length(C_pdown) > 0, n_pup_pdown / length(C_pdown) * 100, NA)
  perc_pdown_pup_Cancer <- ifelse(length(C_pup) > 0, n_pdown_pup / length(C_pup) * 100, NA)
  
  # Total overlaps as % of total HD and cancer
  perc_total_genes_HD <- ifelse(total_HD_genes > 0, total_gene_overlap / total_HD_genes * 100, NA)
  perc_total_genes_Cancer <- ifelse(total_Cancer_genes > 0, total_gene_overlap / total_Cancer_genes * 100, NA)
  perc_total_pathways_HD <- ifelse(total_HD_pathways > 0, total_pathway_overlap / total_HD_pathways * 100, NA)
  perc_total_pathways_Cancer <- ifelse(total_Cancer_pathways > 0, total_pathway_overlap / total_Cancer_pathways * 100, NA)
  
  # Store all values
  results <- rbind(results, data.frame(
    Disease = ccancer,
    
    Gene_up_down = n_up_down,
    Gene_down_up = n_down_up,
    Total_gene_overlap = total_gene_overlap,
    Path_up_down = n_pup_pdown,
    Path_down_up = n_pdown_pup,
    Total_pathway_overlap = total_pathway_overlap,
    
    Perc_gene_up_down_HD = perc_up_down_HD,
    Perc_gene_down_up_HD = perc_down_up_HD,
    Perc_path_up_down_HD = perc_pup_pdown_HD,
    Perc_path_down_up_HD = perc_pdown_pup_HD,
    
    Perc_gene_up_down_Cancer = perc_up_down_Cancer,
    Perc_gene_down_up_Cancer = perc_down_up_Cancer,
    Perc_path_up_down_Cancer = perc_pup_pdown_Cancer,
    Perc_path_down_up_Cancer = perc_pdown_pup_Cancer,
    
    Perc_total_genes_HD = perc_total_genes_HD,
    Perc_total_genes_Cancer = perc_total_genes_Cancer,
    Perc_total_pathways_HD = perc_total_pathways_HD,
    Perc_total_pathways_Cancer = perc_total_pathways_Cancer
  ))
}

# Round percentages for clarity
results[ , grepl("Perc_", names(results))] <- round(results[ , grepl("Perc_", names(results))], 2)

# View final results
print(results)

## Plots based on percentages
# Genes
ggplot(results, aes(x = Disease, y = Perc_total_genes_HD)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Total Gene Overlap with HD (% of HD genes)",
       y = "% overlap", x = "Disease") +
  theme_minimal()

# Pathways
ggplot(results, aes(x = Disease, y = Perc_total_pathways_HD)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(title = "Total Pathway Overlap with HD (% of HD pathways)",
       y = "% overlap", x = "Disease") +
  theme_minimal()

# Plots based on numbers
# Genes
ggplot(results, aes(x = Disease, y = Total_gene_overlap)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Total Gene Overlap with HD (% of HD genes)",
       y = "% overlap", x = "Disease") +
  theme_minimal()



ggplot(results, aes(x = Disease, y = Total_pathway_overlap)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(title = "Total Pathway Overlap with HD (% of HD pathways)",
       y = "% overlap", x = "Disease") +
  theme_minimal()



# Final version
# Genes
p1 = ggplot(results, aes(x = reorder(Disease, -Total_gene_overlap), y = Total_gene_overlap)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Total_gene_overlap), 
            vjust = -0.4, size = 3.5) +
  labs(y = "# overlapping genes", x = "Disease") +
  theme_minimal()
p1

# Pathways
p2 = ggplot(results, aes(x = reorder(Disease, -Total_pathway_overlap), y = Total_pathway_overlap)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  geom_text(aes(label = Total_pathway_overlap), 
            vjust = -0.4, size = 3.5) +
  labs(y = "# overlapping pathways", x = "Disease") +
  theme_minimal()
p2 

# ggsave()

# Stacked for over and under
# Long format for gene overlaps
genes_stack <- results %>%
  select(Disease, Gene_up_down, Gene_down_up) %>%
  pivot_longer(cols = c(Gene_up_down, Gene_down_up),
               names_to = "Direction",
               values_to = "Count") %>%
  mutate(Direction = recode(Direction,
                            Gene_up_down = "HD up ∩ Cancer down",
                            Gene_down_up = "HD down ∩ Cancer up"))

# Add total gene overlaps
gene_totals <- genes_stack %>%
  group_by(Disease) %>%
  summarise(Total = sum(Count))

ggplot(genes_stack, aes(x = reorder(Disease, -Count), y = Count, fill = Direction)) +
  geom_bar(stat = "identity") +
  geom_text(data = gene_totals,
            aes(x = Disease, y = Total, label = Total),
            inherit.aes = FALSE,  # prevents ggplot from expecting Direction here
            vjust = -0.4, size = 3.5) +
  labs(title = "Gene Overlap with HD (Stacked)",
       y = "# overlapping genes", x = "Disease", fill = "Direction") +
  scale_fill_manual(values = c("HD up ∩ Cancer down" = "steelblue", 
                               "HD down ∩ Cancer up" = "firebrick")) +
  theme_minimal()


# Prepare stacked pathway data
pathways_stack <- results %>%
  select(Disease, Path_up_down, Path_down_up) %>%
  pivot_longer(cols = c(Path_up_down, Path_down_up),
               names_to = "Direction",
               values_to = "Count") %>%
  mutate(Direction = recode(Direction,
                            Path_up_down = "HD up ∩ Cancer down",
                            Path_down_up = "HD down ∩ Cancer up"))

# Add total pathway overlaps
pathway_totals <- pathways_stack %>%
  group_by(Disease) %>%
  summarise(Total = sum(Count))

# Plot
ggplot(pathways_stack, aes(x = reorder(Disease, -Count), y = Count, fill = Direction)) +
  geom_bar(stat = "identity") +
  geom_text(data = pathway_totals,
            aes(x = Disease, y = Total, label = Total),
            inherit.aes = FALSE,
            vjust = -0.4, size = 3.5) +
  labs(title = "Pathway Overlap with HD (Stacked)",
       y = "# overlapping pathways", x = "Disease", fill = "Direction") +
  scale_fill_manual(values = c("HD up ∩ Cancer down" = "orange", 
                               "HD down ∩ Cancer up" = "darkred")) +
  theme_minimal()



library(UpSetR)

# Initialize lists
gene_up_down_list <- list()
gene_down_up_list <- list()
path_up_down_list <- list()
path_down_up_list <- list()

# Build the overlap lists
for (ccancer in cancer_list) {
  # Gene sets
  C_up <- diseases_dic[[ccancer]]@genes_up
  C_down <- diseases_dic[[ccancer]]@genes_down
  
  # Pathway sets
  C_pup <- diseases_dic[[ccancer]]@pathways_up
  C_pdown <- diseases_dic[[ccancer]]@pathways_down
  
  # Genes
  gene_up_down_list[[ccancer]] <- intersect(HD_up, C_down)     # HD up ∩ Cancer down
  gene_down_up_list[[ccancer]] <- intersect(HD_down, C_up)     # HD down ∩ Cancer up
  
  # Pathways
  path_up_down_list[[ccancer]] <- intersect(HD_pup, C_pdown)   # HD up ∩ Cancer down
  path_down_up_list[[ccancer]] <- intersect(HD_pdown, C_pup)   # HD down ∩ Cancer up
}

# Convert to binary matrices
gene_up_down_matrix <- fromList(gene_up_down_list)
gene_down_up_matrix <- fromList(gene_down_up_list)
path_up_down_matrix <- fromList(path_up_down_list)
path_down_up_matrix <- fromList(path_down_up_list)

# Plot UpSet for genes (HD up ∩ Cancer down)
upset1 = upset(gene_up_down_matrix,
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "steelblue",
      mainbar.y.label = "HD up ∩ Cancer down (genes)",
      sets.x.label = "Genes per cancer")
print(upset1)

# Plot UpSet for genes (HD down ∩ Cancer up)
upset2 = upset(gene_down_up_matrix,
      order.by = "freq",
      main.bar.color = "firebrick",
      sets.bar.color = "firebrick",
      mainbar.y.label = "HD down ∩ Cancer up (genes)",
      sets.x.label = "Genes per cancer")
print(upset2)

# Plot UpSet for pathways (HD up ∩ Cancer down)
upset3 = upset(path_up_down_matrix,
      order.by = "freq",
      main.bar.color = "orange",
      sets.bar.color = "orange",
      mainbar.y.label = "HD up ∩ Cancer down (pathways)",
      sets.x.label = "Pathways per cancer")
print(upset3)

# Plot UpSet for pathways (HD down ∩ Cancer up)
upset4 = upset(path_down_up_matrix,
      order.by = "freq",
      main.bar.color = "darkred",
      sets.bar.color = "darkred",
      mainbar.y.label = "HD down ∩ Cancer up (pathways)",
      sets.x.label = "Pathways per cancer")
print(upset4)

# Helper function: count frequency across lists
get_top_elements <- function(list_of_sets, top_n = 20) {
  unlisted <- unlist(list_of_sets)
  as.data.frame(table(unlisted)) %>%
    arrange(desc(Freq)) %>%
    slice(1:top_n) %>%
    rename(Element = unlisted, NumDiseases = Freq)
}

# Get top 10 for each
top_gene_up_down <- get_top_elements(gene_up_down_list) # Up in HD and down in cancer
top_gene_down_up <- get_top_elements(gene_down_up_list) # Down in HD and up in cancer
top_path_up_down <- get_top_elements(path_up_down_list)
top_path_down_up <- get_top_elements(path_down_up_list)


# Show them
cat(" Top 10 genes: HD up ∩ Cancer down\n")
print(top_gene_up_down)

cat("\n Top 10 genes: HD down ∩ Cancer up\n")
print(top_gene_down_up)

cat("\n Top 10 pathways: HD up ∩ Cancer down\n")
print(top_path_up_down)

cat("\n Top 10 pathways: HD down ∩ Cancer up\n")
print(top_path_down_up)

#### Proportion of genes and pathways dysregulated in opposite direction
# Initialize output table
opposite_proportions <- data.frame(
  Disease = character(),
  Total_shared_genes = integer(),
  Opposite_genes = integer(),
  Prop_opposite_genes = numeric(),
  Total_shared_pathways = integer(),
  Opposite_pathways = integer(),
  Prop_opposite_pathways = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each cancer
for (ccancer in cancer_list) {
  # Get cancer signatures
  C_up <- diseases_dic[[ccancer]]@genes_up
  C_down <- diseases_dic[[ccancer]]@genes_down
  C_pup <- diseases_dic[[ccancer]]@pathways_up
  C_pdown <- diseases_dic[[ccancer]]@pathways_down
  
  # Shared genes
  shared_genes <- intersect(union(HD_up, HD_down), union(C_up, C_down))
  up_down_genes <- intersect(HD_up, C_down)
  down_up_genes <- intersect(HD_down, C_up)
  opposite_genes <- union(up_down_genes, down_up_genes)
  
  # Shared pathways
  shared_pathways <- intersect(union(HD_pup, HD_pdown), union(C_pup, C_pdown))
  up_down_pathways <- intersect(HD_pup, C_pdown)
  down_up_pathways <- intersect(HD_pdown, C_pup)
  opposite_pathways <- union(up_down_pathways, down_up_pathways)
  
  # Proportions
  prop_genes <- ifelse(length(shared_genes) > 0,
                       length(opposite_genes) / length(shared_genes), NA)
  prop_pathways <- ifelse(length(shared_pathways) > 0,
                          length(opposite_pathways) / length(shared_pathways), NA)
  
  # Add to output
  opposite_proportions <- rbind(opposite_proportions, data.frame(
    Disease = ccancer,
    Total_shared_genes = length(shared_genes),
    Opposite_genes = length(opposite_genes),
    Prop_opposite_genes = round(prop_genes, 3),
    Total_shared_pathways = length(shared_pathways),
    Opposite_pathways = length(opposite_pathways),
    Prop_opposite_pathways = round(prop_pathways, 3)
  ))
}

# View result
print(opposite_proportions)

#                        Disease Total_shared_genes Opposite_genes Prop_opposite_genes Total_shared_pathways Opposite_pathways Prop_opposite_pathways
# 1                 Liver cancer               1367            791               0.579                   252               179                  0.710
# 2                  Lung cancer               1451            843               0.581                   412               316                  0.767
# 3                Breast cancer               1714           1028               0.600                   671               617                  0.920
# 4 Chronic lymphocytic leukemia                207            144               0.696                    20                20                  1.000

mean(opposite_proportions$Prop_opposite_genes) # 0.614
mean(opposite_proportions$Prop_opposite_pathways) # 0.84925

#### Only for genes< ####
# Initialize named lists
gene_up_down_list <- list()
gene_down_up_list <- list()
path_up_down_list <- list()
path_down_up_list <- list()

for (ccancer in cancer_list) {
  C_up <- diseases_dic[[ccancer]]@genes_up
  C_down <- diseases_dic[[ccancer]]@genes_down
  
  # HD up ∩ Cancer down
  up_down_genes <- intersect(HD_up, C_down)
  gene_up_down_list[[ccancer]] <- up_down_genes
  
  # HD down ∩ Cancer up
  down_up_genes <- intersect(HD_down, C_up)
  gene_down_up_list[[ccancer]] <- down_up_genes
}

#  Convert the gene lists into a binary presence matrix
# For up_down
up_down_matrix <- fromList(gene_up_down_list)

# For down_up
down_up_matrix <- fromList(gene_down_up_list)

# UpSet for HD up ∩ Cancer down
upset(up_down_matrix,
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "steelblue",
      mainbar.y.label = "HD up ∩ Cancer down gene overlaps",
      sets.x.label = "Genes per cancer")

# UpSet for HD down ∩ Cancer up
upset(down_up_matrix,
      order.by = "freq",
      main.bar.color = "firebrick",
      sets.bar.color = "firebrick",
      mainbar.y.label = "HD down ∩ Cancer up gene overlaps",
      sets.x.label = "Genes per cancer")
#### >Only for genes #####

####### OLD PLOTS
ggplot(genes_stack, aes(x = reorder(Disease, -Count), y = Count, fill = Direction)) +
  geom_bar(stat = "identity") +
  labs(title = "Gene Overlap with HD (Stacked)",
       y = "# overlapping genes", x = "Disease", fill = "Direction") +
  scale_fill_manual(values = c("HD up ∩ Cancer down" = "steelblue", 
                               "HD down ∩ Cancer up" = "firebrick")) +
  theme_minimal()


# Transform for gene overlap
genes_long <- results %>%
  select(Disease, Perc_total_genes_HD, Perc_total_genes_Cancer) %>%
  pivot_longer(
    cols = c(Perc_total_genes_HD, Perc_total_genes_Cancer),
    names_to = "Perspective",
    values_to = "Percentage"
  ) %>%
  mutate(Perspective = gsub("Perc_total_genes_", "", Perspective))

# Plot
ggplot(genes_long, aes(x = Disease, y = Percentage, fill = Perspective)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Total Gene Overlap: HD vs Cancer Perspective",
       y = "% overlap", x = "Disease") +
  scale_fill_manual(values = c("HD" = "steelblue", "Cancer" = "darkred")) +
  theme_minimal()

# Transform for pathway overlap
pathways_long <- results %>%
  select(Disease, Perc_total_pathways_HD, Perc_total_pathways_Cancer) %>%
  pivot_longer(
    cols = c(Perc_total_pathways_HD, Perc_total_pathways_Cancer),
    names_to = "Perspective",
    values_to = "Percentage"
  ) %>%
  mutate(Perspective = gsub("Perc_total_pathways_", "", Perspective))

# Plot
ggplot(pathways_long, aes(x = Disease, y = Percentage, fill = Perspective)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Total Pathway Overlap: HD vs Cancer Perspective",
       y = "% overlap", x = "Disease") +
  scale_fill_manual(values = c("HD" = "darkorange", "Cancer" = "forestgreen")) +
  theme_minimal()


##### COLORED BY NUMBERS 
# Create normalized intensity values (0–1 scale) for coloring
genes_long <- results %>%
  select(Disease, 
         Perc_total_genes_HD, Perc_total_genes_Cancer, 
         Total_gene_overlap) %>%
  pivot_longer(
    cols = c(Perc_total_genes_HD, Perc_total_genes_Cancer),
    names_to = "Perspective",
    values_to = "Percentage"
  ) %>%
  mutate(
    Perspective = gsub("Perc_total_genes_", "", Perspective),
    OverlapCount = Total_gene_overlap
  ) %>%
  group_by(Perspective) %>%
  mutate(
    NormOverlap = rescale(OverlapCount, to = c(0.3, 1))  # avoid pure white
  ) %>%
  ungroup()

# Create the color column manually
genes_long <- results %>%
  select(Disease, 
         Perc_total_genes_HD, Perc_total_genes_Cancer, 
         Total_gene_overlap) %>%
  pivot_longer(
    cols = c(Perc_total_genes_HD, Perc_total_genes_Cancer),
    names_to = "Perspective",
    values_to = "Percentage"
  ) %>%
  mutate(
    Perspective = gsub("Perc_total_genes_", "", Perspective),
    OverlapCount = Total_gene_overlap,
    NormOverlap = rescale(OverlapCount, to = c(0.3, 1))  # GLOBAL rescale
  )

# Now define color shades manually using both group and intensity
genes_long$FillColor <- mapply(function(p, alpha) {
  if (p == "HD") {
    rgb(0, 0, 1, alpha = alpha)  # blue
  } else {
    rgb(1, 0, 0, alpha = alpha)  # red
  }
}, genes_long$Perspective, genes_long$NormOverlap)


# Plot
ggplot(genes_long, aes(x = Disease, y = Percentage, fill = FillColor)) +
  geom_col(position = position_dodge()) +
  scale_fill_identity() +
  labs(title = "Total Gene Overlap: HD vs Cancer",
       subtitle = "Bar height = % overlap, Color intensity = # overlapping genes",
       y = "% overlap", x = "Disease") +
  theme_minimal() +
  theme(legend.position = "none")



#### COLORED BY NUMBERS - no funciona. 
# Prepare long data for genes
genes_long <- results %>%
  select(Disease, 
         Perc_total_genes_HD, Perc_total_genes_Cancer, 
         Total_gene_overlap) %>%
  pivot_longer(
    cols = c(Perc_total_genes_HD, Perc_total_genes_Cancer),
    names_to = "Perspective",
    values_to = "Percentage"
  ) %>%
  mutate(
    Perspective = gsub("Perc_total_genes_", "", Perspective),
    OverlapCount = Total_gene_overlap  # use for color fill
  )

# Plot
ggplot(genes_long, aes(x = Disease, y = Percentage, fill = OverlapCount)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Gene Overlap: HD vs Cancer (% overlap, color = overlap count)",
       y = "% overlap", x = "Disease", fill = "# overlapping\ngenes") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal()

# Prepare long data for pathways
pathways_long <- results %>%
  select(Disease, 
         Perc_total_pathways_HD, Perc_total_pathways_Cancer, 
         Total_pathway_overlap) %>%
  pivot_longer(
    cols = c(Perc_total_pathways_HD, Perc_total_pathways_Cancer),
    names_to = "Perspective",
    values_to = "Percentage"
  ) %>%
  mutate(
    Perspective = gsub("Perc_total_pathways_", "", Perspective),
    OverlapCount = Total_pathway_overlap
  )

# Plot
ggplot(pathways_long, aes(x = Disease, y = Percentage, fill = OverlapCount)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Pathway Overlap: HD vs Cancer (% overlap, color = overlap count)",
       y = "% overlap", x = "Disease", fill = "# overlapping\npathways") +
  scale_fill_gradient(low = "moccasin", high = "darkorange4") +
  theme_minimal()



