##Creating MDS plot for selected dataset
#load package
library("phyloseq")
library("microbiome")
library("vegan")
library("ggplot2")
library("dplyr")

#Palette ("#23207C", "#00B5ED", "#F08C01", "#D6000A") or ("#969a2a", "#566017", "e1d167", "4c3e17") or ("#1B9E77", "#D95F02", "#7570B3", "#E7298A") or ("#9E0142", "#F46D43", "#FEE08B", "#E6F598")

#set seed
set.seed(2021)

#Transform abundances to square root
psO_perifiton16S_sqr <- microbiome::transform(psO_perifiton16S, "hellinger")

#Multivariate analysis based on Bray-Curtis distance and MDS ordination method
MDS_bray_psO_perifiton16S_sqr<-ordinate(psO_perifiton16S_sqr, "MDS","bray")

#Print stress data, dimensions and number of tries
head(MDS_bray_psO_perifiton16S_sqr)

####################################### 2D plot
##Create a MDS plot (all)
#Place-Time
plot_MDS_bray_psO_perifiton16S_sqr_all<-plot_ordination(psO_perifiton16S_sqr,MDS_bray_psO_perifiton16S_sqr, type="sample",color="Time",shape="Place") +
  geom_point(size=3.0) +
  theme_bw() +
  theme(legend.title = element_text(size = 8, color = "black")) +
  theme(legend.text = element_text(size=8, color = "black")) +
  theme(axis.title = element_text(size = 10, color = "black")) +
  theme(axis.text = element_text(size = 10, colour = "black")) +
  scale_color_manual(values = c("#23207C", "#00B5ED", "#F08C01", "#D6000A"))
plot_MDS_bray_psO_perifiton16S_sqr_all

#ggsave("plot_MDS_bray_psO_perifiton16S_sqr_all.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Ordination/", width = 18, height = 15, units = "cm",dpi = 600)
ggsave("plot_MDS_bray_psO_perifiton16S_sqr_all.png", path = "~/Desktop/", width = 18, height = 15, units = "cm",dpi = 600, device = "png")
ggsave("plot_MDS_bray_psO_perifiton16S_sqr_all.eps", path = "~/Desktop/", width = 10, height = 7, units = "cm",dpi = 600, device = "eps")

####################################### STATISTICS
##Multivariate Anova for ordinate files
set.seed(2021)

#Convert phyloseq object to dataframe using abundances function and meta function from microbiome package
psO_perifiton16S_sqr_abundances <- abundances(psO_perifiton16S_sqr)
psO_perifiton16S_sqr_meta <- meta(psO_perifiton16S_sqr)

##Permanova using adonis function from vegan and print p value
#all data
Permanova_psO_perifiton16S_place <- adonis(t(psO_perifiton16S_sqr_abundances) ~ Place, data = psO_perifiton16S_sqr_meta, permutations = 10000, method = "bray")
Permanova_psO_perifiton16S_place

Permanova_psO_perifiton16S_time <- adonis(t(psO_perifiton16S_sqr_abundances) ~ Time, data = psO_perifiton16S_sqr_meta, permutations = 10000, method = "bray")
Permanova_psO_perifiton16S_time

Permanova_psO_perifiton16S_place_time <- adonis(t(psO_perifiton16S_sqr_abundances) ~ Place*Time, data = psO_perifiton16S_sqr_meta, permutations = 10000, method = "bray")
Permanova_psO_perifiton16S_place_time

print(as.data.frame(Permanova_psO_perifiton16S_place_time$aov.tab)["Place", "Pr(>F)"])
print(as.data.frame(Permanova_psO_perifiton16S_place_time$aov.tab)["Time", "Pr(>F)"])
