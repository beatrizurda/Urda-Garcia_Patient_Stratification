#########################################################################################
# PLOT BARPLOT OF THE OVERLAP OF DISEASE-DISEASE MOLECULAR NETWORKS
# 
# 
#
# Beatriz Urda García 2021
#########################################################################################



overlaps <- c(46.2,36.09,36.03,16,8.71)
names(overlaps) <- c("RNA-seq", "Microbiome", "miRNAs", "microarrays", "PPIs")
barplot(overlaps,ylim=c(0,50),las=2,col=parent_colors[,2],
        legend=T)

colors <- c("#621D7D","#568AB8","#D6C85E")
colors <- c("#36396B","#5C62B8","#BEC2EE")
colors <- c("#4874B5","#E69B4C","#E3C542")
colors <- c("#4849B5","#AB9430","#E3C542")
colors <- c("#364F70","#4E73A3","#BBD2F0")

overlaps <- c(46.2,16,8.71)
names(overlaps) <- c("RNA-seq", "Microarrays", "PPIs")
par(mar=c(6,5,5,5))
barplot(overlaps,ylim=c(0,60),las=2,col=colors,
        legend=T)

df <- as.data.frame(overlaps)
df$label <- c("RNA-seq", "Microarrays", "PPIs")
g <- ggplot(df, aes(x=reorder(label, -overlaps), y=overlaps, fill=reorder(label, -overlaps))) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("#364F70","#4E73A3","#BBD2F0")) + ylim(0,50) +
  xlab("Disease-disease networks based on molecular information") +
  ylab("Percentage of overlap with the epidemiology") +
  labs(fill="")+
  theme_classic()+
  # labs(fill="Disease-disease networks based on molecular information")+
  theme(axis.text.x =element_text(size=11), axis.title = element_text(size=12)) +
  theme(legend.position="bottom")
  # theme(axis.text.x =element_text(size=11, angle=45, hjust=1))
g

# Save with size 5.5x8.0

### WITH THE SSN
colors =  c("#23344A","#364F70","#4E73A3","#BBD2F0")
overlaps <- c(64.1,46.2,16,8.71)
names(overlaps) <- c("RNA-seq", "Microarrays", "PPIs")
par(mar=c(6,5,5,5))
barplot(overlaps,ylim=c(0,70),las=2,col=colors,
        legend=T)

df <- as.data.frame(overlaps)
df$label <- c("RNA-seq\n(meta-patients)", "RNA-seq\n(diseases)", "Microarrays", "PPIs")
g <- ggplot(df, aes(x=reorder(label, -overlaps), y=overlaps, fill=reorder(label, -overlaps))) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("#23344A","#364F70","#4E73A3","#BBD2F0")) + ylim(0,70) +
  xlab("Disease-disease networks based on molecular information") +
  ylab("Percentage of overlap with the epidemiology") +
  labs(fill="")+
  theme_classic()+
  # labs(fill="Disease-disease networks based on molecular information")+
  theme(axis.text.x =element_text(size=11), axis.title = element_text(size=12)) +
  theme(legend.position="bottom")
# theme(axis.text.x =element_text(size=11, angle=45, hjust=1))
g


c("#621D7D","#568AB8","#D6C85E")
c("#621D7D","#568AB8","#D6C85E")