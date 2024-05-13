#Create rarefied data from phyloseq object using vegan package or microbiome package
#Loading package
library("phyloseq")
library("vegan")
library("microbiome")
library("ggpubr")

#Palette color code ("#26a4b5", "#3c809b", "#05b6fc","#8bd8f4", "#014560")

##Creating rarefaction curve of non-rarefied samples for prokaryotes or only bacteria
#Using rarecurve()
rarecurve(t(otu_table(psO_perifiton16S)), step=20, cex=0.5, col = "blue")
rarecurve(t(otu_table(psO_perifiton16S_jangadeiros_filt)), step=20, cex=0.5, col = "blue")
rarecurve(t(otu_table(psO_perifiton16S_veleiros_filt)), step=20, cex=0.5, col = "blue")

##Rarefying samples to the minimum number of reads among samples (only bacteria)
#Using rarefy_even_depth() rarefy sample reads to a specific number of total sequences
psO_perifiton16S_rarefied<-rarefy_even_depth(psO_perifiton16S, rngseed=2021,sample.size = 4200, trimOTUs=TRUE)
psO_perifiton16S_rarefied
sample_sums(psO_perifiton16S_rarefied)

##Creating rarefaction curve of rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_perifiton16S_rarefied)), step=20, cex=0.5, col = "black")


##Calculates the alfa diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_perifiton16S_rarefied_richness <- richness(psO_perifiton16S_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_perifiton16S_rarefied_eveness <- evenness(psO_perifiton16S_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_perifiton16S_rarefied_diversity <- microbiome::diversity(psO_perifiton16S_rarefied, index = 'shannon', zeroes = TRUE)

##################################################
##Boxplot alpha diversity using ggviolin function from ggpubr
#Prepare file using meta function, and combn(seq_along)
psO_perifiton16S_rarefied.meta <- meta(psO_perifiton16S_rarefied)

#select collumn from output file to be plotted and tested
#for all
psO_perifiton16S_rarefied.meta$chao1 <- alfa_div_psO_perifiton16S_rarefied_richness$chao1
psO_perifiton16S_rarefied.meta$pielou <- alfa_div_psO_perifiton16S_rarefied_eveness$pielou
psO_perifiton16S_rarefied.meta$shannon <- alfa_div_psO_perifiton16S_rarefied_diversity$shannon
psO_perifiton16S_rarefied.meta$observed <- alfa_div_psO_perifiton16S_rarefied_richness$observed

head(psO_perifiton16S_rarefied.meta)
write.csv(psO_perifiton16S_rarefied.meta, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_meta_psO_perifiton16S_rarefied.csv")

##check for the distribution of the diversity using "hist" function to plot, "shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05) then the null hypotesis (the population is normally distributed) is rejected
#for data with only two variable

hist(psO_perifiton16S_rarefied.meta$chao1)
shapiro.test(psO_perifiton16S_rarefied.meta$chao1)
qqnorm(psO_perifiton16S_rarefied.meta$chao1)

hist(psO_perifiton16S_rarefied.meta$pielou)
shapiro.test(psO_perifiton16S_rarefied.meta$pielou)
qqnorm(psO_perifiton16S_rarefied.meta$pielou)

hist(psO_perifiton16S_rarefied.meta$shannon)
shapiro.test(psO_perifiton16S_rarefied.meta$shannon)
qqnorm(psO_perifiton16S_rarefied.meta$shannon)

#Combine files data and defining comparisons
alfa_div_comparison_place <- list(c("Jangadeiros","Veleiros"))
alfa_div_comparison_time<- list(c("T1","T2"),c("T1","T3"),c("T1","T4"),c("T2","T3"),c("T2","T4"),c("T3","T4"))


##plot using ggviolin
#test kruskall_wallis to compare places
violin_plot_psO_perifiton16S_rarefied_shannon_place <- ggviolin(psO_perifiton16S_rarefied.meta, x="Place", y= "shannon", add = "jitter",fill="Place", palette = c("#26a4b5", "#014560"), ylab = "Shannon index", ylim = c(0, 6)) +
  geom_jitter(width=0.1, size = 1) +
  stat_compare_means(method = "wilcox.test", label ="p", label.y = 5.5, label.x = 2, size = 3) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 10, color = "black")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 10, colour = "black"))
violin_plot_psO_perifiton16S_rarefied_shannon_place

#ggsave("violin_plot_psO_perifiton16S_rarefied_shannon_place.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Alpha_diversity/", width = 15, height = 15, units = "cm",dpi = 600)
ggsave("violin_plot_psO_perifiton16S_rarefied_shannon_place.png", path = "~/Desktop/", width = 15, height = 15, units = "cm",dpi = 600, device = "png")
ggsave("violin_plot_psO_perifiton16S_rarefied_shannon_place.eps", path = "~/Desktop/", width = 5, height = 5, units = "cm",dpi = 600, device = "eps")

violin_plot_psO_perifiton16S_rarefied_richness_place <- ggviolin(psO_perifiton16S_rarefied.meta, x = "Place", y = "chao1", add = "jitter", fill = "Place", palette = c("#26a4b5", "#014560"), ylab = "Chao1 index", ylim = c(0, 400)) +
  geom_jitter(width = 0.1, size = 1) +
  stat_compare_means(method = "wilcox.test", label = "p", label.y = 340, label.x = 2, size = 3) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 10, color = "black")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 10, colour = "black"))
violin_plot_psO_perifiton16S_rarefied_richness_place

#ggsave("violin_plot_psO_perifiton16S_rarefied_richness_place.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Alpha_diversity/", width = 15, height = 15, units = "cm",dpi = 600)
ggsave("violin_plot_psO_perifiton16S_rarefied_richness_place.png", path = "~/Desktop/", width = 15, height = 15, units = "cm",dpi = 600, device = "png")
ggsave("violin_plot_psO_perifiton16S_rarefied_richness_place.eps", path = "~/Desktop/", width = 6.3, height = 5, units = "cm",dpi = 600, device = "eps")

violin_plot_psO_perifiton16S_rarefied_pielou_place <- ggviolin(psO_perifiton16S_rarefied.meta, x="Place", y= "pielou", add = "jitter",fill="Place", palette = c("#26a4b5", "#014560"),ylab = "Pielou", ylim = c(0.4, 1.2)) +
  geom_jitter(width=0.1, size = 1) +
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 1.15, label.x = 2, size=3) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_blank()) +
  theme(axis.title = element_text(size = 10, color = "black")) +
  theme(axis.text = element_text(size = 10, colour = "black"))
violin_plot_psO_perifiton16S_rarefied_pielou_place

#ggsave("violin_plot_psO_perifiton16S_rarefied_pielou_place.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Alpha_diversity/", width = 15, height = 15, units = "cm",dpi = 600)
ggsave("violin_plot_psO_perifiton16S_rarefied_pielou_place.png", path = "~/Desktop/", width = 15, height = 15, units = "cm",dpi = 600, device = "png")
ggsave("violin_plot_psO_perifiton16S_rarefied_pielou_place.eps", path = "~/Desktop/", width = 6.3, height = 5, units = "cm",dpi = 600, device = "eps")

#test kruskall_wallis to compare time points
#colored by place
violin_plot_psO_perifiton16S_rarefied_shannon_time <- ggviolin(psO_perifiton16S_rarefied.meta, x="Time", y= "shannon",ylab = "Shannon index", ylim = c(0, 6)) +
  geom_jitter(aes(color=Place), width=0.1, size = 1) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 10, colour = "black")) +
  scale_color_manual(values =c("#26a4b5", "#014560")) +
  stat_compare_means(method = "kruskal.test", label = "p", label.y = 6.0, label.x = 3, size = 3) 
violin_plot_psO_perifiton16S_rarefied_shannon_time

#ggsave("violin_plot_psO_perifiton16S_rarefied_shannon_time.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Alpha_diversity/", width = 15, height = 15, units = "cm",dpi = 600)
ggsave("violin_plot_psO_perifiton16S_rarefied_shannon_time.png", path = "~/Desktop/", width = 15, height = 15, units = "cm",dpi = 600, device = "png")
ggsave("violin_plot_psO_perifiton16S_rarefied_shannon_time.eps", path = "~/Desktop/", width = 6, height = 5, units = "cm",dpi = 600, device = "eps")

violin_plot_psO_perifiton16S_rarefied_richness_time <- ggviolin(psO_perifiton16S_rarefied.meta, x="Time", y= "chao1",ylab = "Chao1 index", ylim = c(0, 400)) +
  geom_jitter(aes(color=Place), width=0.1, size = 1) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 10, colour = "black")) +
  scale_color_manual(values =c("#26a4b5", "#014560")) +
  stat_compare_means(method = "kruskal.test", label= "p", label.y = 340, label.x = 3, size=3)
violin_plot_psO_perifiton16S_rarefied_richness_time

#ggsave("violin_plot_psO_perifiton16S_rarefied_richness_time.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Alpha_diversity/", width = 15, height = 15, units = "cm",dpi = 600)
ggsave("violin_plot_psO_perifiton16S_rarefied_richness_time.png", path = "~/Desktop/", width = 15, height = 15, units = "cm",dpi = 600, device = "png")
ggsave("violin_plot_psO_perifiton16S_rarefied_richness_time.eps", path = "~/Desktop/", width = 6.2, height = 5, units = "cm",dpi = 600, device = "eps")

violin_plot_psO_perifiton16S_rarefied_pielou_time <- ggviolin(psO_perifiton16S_rarefied.meta, x ="Time", y = "pielou", ylab = "Pielou", ylim = c(0.4, 1.2)) +
  geom_jitter(aes(color=Place), width=0.1, size = 1) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 10, color = "black")) +
  theme(axis.text = element_text(size = 10, colour = "black")) +
  scale_color_manual(values =c("#26a4b5", "#014560")) +
  stat_compare_means(comparisons = alfa_div_comparison_time, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(1.02, 1.1)) +
  stat_compare_means(method = "kruskal.test", label= "p", label.y = 1.15, label.x = 3, size=3)
violin_plot_psO_perifiton16S_rarefied_pielou_time

#ggsave("violin_plot_psO_perifiton16S_rarefied_pielou_time.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Alpha_diversity/", width = 15, height = 15, units = "cm",dpi = 600)
ggsave("violin_plot_psO_perifiton16S_rarefied_pielou_time.png", path = "~/Desktop/", width = 15, height = 15, units = "cm",dpi = 600, device = "png")
ggsave("violin_plot_psO_perifiton16S_rarefied_pielou_time.eps", path = "~/Desktop/", width = 6.2, height = 5, units = "cm",dpi = 600, device = "eps")

###########################################
#arrange multiple ggplots in one single figure 
#ggarrange(violin_plot_chao1_C, violin_plot_pielou_C, violin_plot_shannon_C, violin_plot_dominance_C, labels = c ("A", "B","C","D"), ncol =2, nrow =2)
#ggarrange(violin_plot_chao1_L, violin_plot_shannon_L, violin_plot_dominance_L, labels = c ("A", "B","C"), ncol =3, nrow =1)
alpha_plots <- ggarrange(violin_plot_psO_perifiton16S_rarefied_shannon_place, violin_plot_psO_perifiton16S_rarefied_richness_place, violin_plot_psO_perifiton16S_rarefied_pielou_place, violin_plot_psO_perifiton16S_rarefied_shannon_time, violin_plot_psO_perifiton16S_rarefied_richness_time, violin_plot_psO_perifiton16S_rarefied_pielou_time, ncol =3, nrow =2)
ggsave("alpha_plots.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Alpha_diversity/", width = 40, height = 35, units = "cm",dpi = 600)

