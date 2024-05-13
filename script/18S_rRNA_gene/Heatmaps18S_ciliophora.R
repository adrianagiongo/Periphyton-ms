##Creating heatmap for selected dataset

##load libraries
library("phyloseq")
library("microbiome")
library("ggplot2")
library("RColorBrewer")

## Jangadeiros & Veleiros
#Aglomerate taxa using tax glom function and removing NAs
psO_perifiton18S_ciliophora_annotation = tax_glom(psO_perifiton18S_ciliophora, "Annotation_curated", NArm= TRUE)

#Transfom data to square root abundance
psO_perifiton18S_ciliophora_sqr <- microbiome::transform(psO_perifiton18S_ciliophora_annotation, "hellinger")

psO_perifiton18S_ciliophora_rel <- microbiome::transform(psO_perifiton18S_ciliophora_annotation, "compositional")
df_psO_perifiton18S_ciliophora_rel <- data.frame(tax_table(psO_perifiton18S_ciliophora_rel), otu_table(psO_perifiton18S_ciliophora_rel))
write.csv(df_psO_perifiton18S_ciliophora_rel, "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Tables_18S_Ciliophora/df_psO_perifiton18S_ciliophora_rel.csv")

###Heatmap plot
##Using plot_Heatmap()
#Heatmap for ciliophora otus observed on psO_perifiton18S_jangadeiros_sqr facet on time
psO_perifiton18S_ciliophora_sqr_time = plot_heatmap(psO_perifiton18S_ciliophora_sqr,sample.order = "Time", taxa.label = "Annotation_curated", taxa.order="Annotation_curated", low="white", high="#03018C", na.value="white", trans = NULL) +
  facet_grid(~Place, scales = "free") +
  geom_tile(colour = "gray", size = 0.25) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 7, angle = 90,vjust = 0.5, hjust = 1, colour = "black" ), axis.text.y = element_text(size= 7,colour = "black", face = "italic")) +
  #theme(legend.text = element_text(size = 6), legend.title = element_text(size = 6), legend.key.size = unit(0.3, 'cm')) +
  theme(legend.position = "none") +
  scale_fill_continuous(low = "white",high = "#03018C", labels = scales::percent) +
  theme(strip.text.x = element_text(size = 8, face ="bold", colour = "black"))
psO_perifiton18S_ciliophora_sqr_time

#ggsave("psO_perifiton18S_ciliophora_sqr_time.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Heatmap_18S_Ciliophora/", width = 28, height = 34, units = "cm",dpi = 600)
ggsave("psO_perifiton18S_ciliophora_sqr_time.png", path = "~/Desktop/", width = 10, height = 14, units = "cm", dpi = 600, device = "png")
ggsave("psO_perifiton18S_ciliophora_sqr_time.eps", path = "~/Desktop/", width = 10, height = 14, units = "cm", dpi = 600, device = "eps")

