##Creating heatmap for selected dataset
##load libraries
library("phyloseq")
library("microbiome")
library("ggplot2")
library("RColorBrewer")

###Filtering table removing taxa in low abundance, using filter_taxa()
##subset dominant taxa
##Transform count data to square root
#Jangadeiros
psO_perifiton16S_jangadeiros_filt_annotation_sqr <- microbiome::transform(psO_perifiton16S_jangadeiros_filt_annotation, "hellinger")

df_psO_perifiton16S_jangadeiros_filt_annotation_sqr <- data.frame(tax_table(psO_perifiton16S_jangadeiros_filt_annotation_sqr), otu_table(psO_perifiton16S_jangadeiros_filt_annotation_sqr))
write.csv(df_psO_perifiton16S_jangadeiros_filt_annotation_sqr, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_jangadeiros_filt_annotation_sqr.csv")

#Veleiros
psO_perifiton16S_veleiros_filt_annotation_sqr <- microbiome::transform(psO_perifiton16S_veleiros_filt_annotation, "hellinger")

df_psO_perifiton16S_veleiros_filt_annotation_sqr <- data.frame(tax_table(psO_perifiton16S_veleiros_filt_annotation_sqr), otu_table(psO_perifiton16S_veleiros_filt_annotation_sqr))
write.csv(df_psO_perifiton16S_veleiros_filt_annotation_sqr, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_veleiros_filt_annotation_sqr.csv")

#subset the 50 more abundant otus at annotation level
psO_perifiton16S_jangadeiros_filt_annotation_sqr_50 <- prune_taxa(taxa_names(psO_perifiton16S_jangadeiros_filt_annotation_sqr)[1:50], psO_perifiton16S_jangadeiros_filt_annotation_sqr)
psO_perifiton16S_jangadeiros_filt_annotation_sqr_50

psO_perifiton16S_veleiros_filt_annotation_sqr_50 <- prune_taxa(taxa_names(psO_perifiton16S_veleiros_filt_annotation_sqr)[1:50], psO_perifiton16S_veleiros_filt_annotation_sqr)
psO_perifiton16S_veleiros_filt_annotation_sqr_50

###Heatmap plot
##Using plot_Heatmap()
#Heatmap for bacterial otus observed on psO_perifiton16S_jangadeiros_rel facet on time
heatmap_jangadeiros_filt_annotation_sqr_50_time <- plot_heatmap(psO_perifiton16S_jangadeiros_filt_annotation_sqr_50,sample.order = "Sample_id", taxa.label = "Annotation",taxa.order="Phylum",low="white",high="#03018C",na.value="white", trans = NULL) +
  facet_grid(~Time, scales = "free") +
  geom_tile(colour= "gray",size= 0.25) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size= 8, angle= 90, vjust= 0.5, hjust= 1, colour= "black"), axis.text.y = element_text(size = 8, colour= "black", face= "italic")) +
  theme(legend.text = element_text(size= 6), legend.title= element_text(size= 6)) +
  scale_fill_continuous(low="white",high="#03018C", labels = scales::percent) +
  theme(strip.text.x = element_text(size= 8, face= "bold", colour= "black"))
heatmap_jangadeiros_filt_annotation_sqr_50_time

#ggsave("heatmap_jangadeiros_filt_annotation_sqr_50_time.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Heatmap_16S/", width = 23, height = 38, units = "cm",dpi = 600)
ggsave("heatmap_jangadeiros_filt_annotation_sqr_50_time.png", path = "~/Desktop/", width = 10, height = 14, units = "cm",dpi = 600, device = "png")
ggsave("heatmap_jangadeiros_filt_annotation_sqr_50_time.eps", path = "~/Desktop/", width = 10, height = 14, units = "cm",dpi = 600, device = "eps")

#Heatmap for bacterial otus observed on psO_perifiton16S_veleiros_rel facet on time
heatmap_veleiros_filt_annotation_sqr_50_time <- plot_heatmap(psO_perifiton16S_veleiros_filt_annotation_sqr_50,sample.order = "Sample_id", taxa.label = "Annotation",taxa.order="Phylum",low="white",high="#03018C",na.value="white", trans = NULL) +
  facet_grid(~Time, scales = "free") +
  geom_tile(colour= "gray",size= 0.25) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size= 8, angle= 90, vjust= 0.5, hjust= 1, colour = "black" ), axis.text.y = element_text(size= 8, colour= "black", face= "italic")) +
  #theme(legend.text = element_text(size = 6), legend.title = element_text(size = 6)) +
  theme(legend.text = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size= 8, face= "bold", colour= "black"))
heatmap_veleiros_filt_annotation_sqr_50_time

#ggsave("heatmap_veleiros_filt_annotation_sqr_50_time.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Heatmap_16S/", width = 29, height = 38.5, units = "cm",dpi = 600)
ggsave("heatmap_veleiros_filt_annotation_sqr_50_time.png", path = "~/Desktop/", width = 11, height = 14, units = "cm",dpi = 600, device = "png")
ggsave("heatmap_veleiros_filt_annotation_sqr_50_time.eps", path = "~/Desktop/", width = 11, height = 14, units = "cm",dpi = 600, device = "eps")

#############################

