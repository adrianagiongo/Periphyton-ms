##After the clean dataset script, this is an unsupervised filtering

#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")
library("gridExtra")
library("ggpubr")

#Calculate number of genera present after filtration
length(get_taxa_unique(psO_perifiton18S, taxonomic.rank = "Phylum"))
length(get_taxa_unique(psO_perifiton18S, taxonomic.rank = "Annotation"))

#Agglomerate taxa using tax glom function and removing NAs
psO_perifiton18S_phylum <- tax_glom(psO_perifiton18S, "Phylum", NArm= TRUE)
psO_perifiton18S_phylum
psO_perifiton18S_annotation <- tax_glom(psO_perifiton18S, "Annotation", NArm= TRUE)
psO_perifiton18S_annotation

#Transform data to relative abundance
psO_perifiton18S_phylum_rel <- transform_sample_counts(psO_perifiton18S_phylum, function(x) {x / sum(x)})
psO_perifiton18S_annotation_rel <- transform_sample_counts(psO_perifiton18S_annotation, function(x) {x / sum(x)})

##Create a function to plot abundances for a subset of samples (Subset KINGDOM and facet PHYLUM)
#for Place
plot_abundance_place_phy = function(phyloseq, title = "",
                                   Facet ="Phylum", Color = "Place") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Eukaryota"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Place", y="Abundance", color = Color, fill = Color)) +
    geom_point(size = 1, alpha = 0.9,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top" ) +
    theme(legend.position = "bottom")
}


##plot abundance relative of OTUs (level phylum) using the same plot_abundance function
#for Place
p_abund_rel_place_phylum <- plot_abundance_place_phy(psO_perifiton18S_phylum_rel,"Eukaryota") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=8, colour = "black")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=8, colour = "black")) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size=7, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.position = "bottom") +
  theme(strip.text.x = element_text(size = 7, colour = "black")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), label = c("0", "25", "50", "75", "100")) +
  scale_color_manual(values = c("#26a4b5", "#014560"))
p_abund_rel_place_phylum

#ggsave("p_abund_rel_place_phylum.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/Taxa_Composition/", width = 20, height = 8, units = "cm",dpi = 600)
ggsave("p_abund_rel_place_phylum.eps", path = "~/Desktop/", width = 8, height = 6, units = "cm",dpi = 600)
ggsave("p_abund_rel_place_phylum.png", path = "~/Desktop/", width = 8, height = 6, units = "cm",dpi = 600)
