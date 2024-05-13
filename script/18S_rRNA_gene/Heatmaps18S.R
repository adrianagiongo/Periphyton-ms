##Creating heatmap for selected extended dataset

##load libraries
library("phyloseq")
library("microbiome")
library("ggplot2")
library("RColorBrewer")

#Agglomerate taxa using tax glom function and removing NAs
psO_perifiton18S_extended_clade <- tax_glom(psO_perifiton18S_extended, "Clade", NArm= TRUE)
psO_perifiton18S_extended_clade

psO_perifiton18S_extended_annotation <- tax_glom(psO_perifiton18S_extended, "Annotation", NArm= TRUE)
psO_perifiton18S_extended_annotation

#Transform data to square root abundance
psO_perifiton18S_extended_clade_rel <- microbiome::transform(psO_perifiton18S_extended_clade, "compositional")

df_psO_perifiton18S_extended_clade_rel <- data.frame(tax_table(psO_perifiton18S_extended_clade_rel), otu_table(psO_perifiton18S_extended_clade_rel))
write.csv(df_psO_perifiton18S_extended_clade_rel, "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/Tables_Perifiton_18S/df_psO_perifiton18S_extended_clade_rel.csv")

psO_perifiton18S_extended_clade_sqr <- microbiome::transform(psO_perifiton18S_extended_clade, "hellinger")

df_psO_perifiton18S_extended_clade_sqr <- data.frame(tax_table(psO_perifiton18S_extended_clade_sqr), otu_table(psO_perifiton18S_extended_clade_sqr))
write.csv(df_psO_perifiton18S_extended_clade_sqr, "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/Tables_Perifiton_18S/df_psO_perifiton18S_extended_clade_sqr.csv")

psO_perifiton18S_extended_annotation_sqr <- microbiome::transform(psO_perifiton18S_extended_annotation, "hellinger")

df_psO_perifiton18S_extended_annotation_sqr <- data.frame(tax_table(psO_perifiton18S_extended_annotation_sqr), otu_table(psO_perifiton18S_extended_annotation_sqr))
write.csv(df_psO_perifiton18S_extended_annotation_sqr, "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/Tables_Perifiton_18S/df_psO_perifiton18S_extended_annotation_sqr.csv")

###Heatmap plot
##Using plot_Heatmap()
#Heatmap for taxa observed on psO_perifiton18S_extended_annotation_sqr facet on place
heatmap_psO_perifiton18S_extended_clade_sqr_place <- plot_heatmap(psO_perifiton18S_extended_clade_sqr,sample.order = "Sample_id", taxa.label = "Clade",taxa.order="Clade",low="white",high="#03018C",na.value="white", trans = NULL) +
  facet_grid(~Place, scales = "free") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size=7, angle = 90,vjust =0.5, hjust = 1, colour = "black" ), axis.text.y = element_text(size=8,colour = "black", face = "italic")) +
  theme(legend.text = element_text(size = 7), legend.title = element_text(size = 7), legend.key.size = unit (0.3, "cm")) +
 #theme(legend.position = "none") +
  geom_tile(colour="gray",size=0.25) +
  theme(strip.text.x = element_text(size = 10, face="bold")) +
  scale_fill_continuous(low="white",high="#03018C", labels = scales::percent)
heatmap_psO_perifiton18S_extended_clade_sqr_place

#ggsave("heatmap_psO_perifiton18S_extended_clade_sqr_place.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/Heatmap_18S/", width = 29, height = 12, units = "cm",dpi = 600)
ggsave("heatmap_psO_perifiton18S_extended_clade_sqr_place.eps", path = "~/Desktop/", width = 11, height = 5, units = "cm",dpi = 600, device = "eps")
ggsave("heatmap_psO_perifiton18S_extended_clade_sqr_place.png", path = "~/Desktop/", width = 11, height = 5, units = "cm",dpi = 600, device = "png")

