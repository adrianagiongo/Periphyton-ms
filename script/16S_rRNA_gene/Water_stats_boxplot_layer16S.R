# Stats qPCR
#install.packages("multcompView", repos="http://R-Forge.R-project.org")
#Load packages
library("readxl")
library("ggplot2")
library("emmeans")
library("multcompView")
library("multcomp")
library("tidyverse")
library("ggpubr")

###### Load table (absolute abundance)
perifiton_water_stat <- read_excel("~/Desktop/perifiton_water.xlsx")
perifiton_water_stat
str(perifiton_water_stat)

################
#Stats
# Encode vector as factor. Explanatory variable (Sample), RA_16S is dependent
perifiton_water_stat$Site=as.factor(perifiton_water_stat$Site)
levels(perifiton_water_stat$Site)

# Fitting Linear Model to produce ANOVA table 
perifiton_water_stat_temp <- lm(Temp ~ Site, data = perifiton_water_stat)
anova(perifiton_water_stat_temp)
summary(perifiton_water_stat_temp)

perifiton_water_stat_pH <- lm(pH ~ Site, data = perifiton_water_stat)
anova(perifiton_water_stat_pH)
summary(perifiton_water_stat_pH)

perifiton_water_stat_solids <- lm(Solids ~ Site, data = perifiton_water_stat)
anova(perifiton_water_stat_solids)
summary(perifiton_water_stat_solids)

perifiton_water_stat_conduct <- lm(Conduct ~ Site, data = perifiton_water_stat)
anova(perifiton_water_stat_conduct)
summary(perifiton_water_stat_conduct)

perifiton_water_stat_dissolvO2 <- lm(DissolvO2 ~ Site, data = perifiton_water_stat)
anova(perifiton_water_stat_dissolvO2)
summary(perifiton_water_stat_dissolvO2)

# Post-hoc multicomparison Test, using emmeans package(methods: pairwise (pairwise-test), dunnet (dunn-test), tukey (Tukey-Test) - adjust: BH, bonferroni, fdr, BY, hommel, hockberg, holm)
#contrast(emmeans(root_weight_stat,~Layer*Rotation), method="pairwise", adjust = "none")
contrast(emmeans(root_weight_JKI_stat,~Layer*Rotation), method="pairwise", adjust = "bonferroni")
contrast(emmeans(root_weight_JKI_stat,~Rotation), method="pairwise", adjust = "bonferroni")
summary(emmeans(root_weight_JKI_stat,~Layer*Rotation))
plot(emmeans(root_weight_JKI_stat,~Layer*Rotation))

#contrast(emmeans(root_length_IFZ_stat,~Layer*Rotation), method="pairwise", adjust = "none")
contrast(emmeans(root_length_IFZ_stat,~Layer*Rotation), method="pairwise", adjust = "bonferroni")
contrast(emmeans(root_length_IFZ_stat,~Rotation), method="pairwise", adjust = "bonferroni")
summary(emmeans(root_length_IFZ_stat,~Layer*Rotation))
plot(emmeans(root_length_IFZ_stat,~Layer*Rotation))

#contrast(emmeans(root_surface_IFZ_stat,~Layer*Rotation), method="pairwise", adjust = "none")
contrast(emmeans(root_surface_IFZ_stat,~Layer*Rotation), method="pairwise", adjust = "bonferroni")
contrast(emmeans(root_surface_IFZ_stat,~Rotation), method="pairwise", adjust = "bonferroni")
summary(emmeans(root_surface_IFZ_stat,~Layer*Rotation))
plot(emmeans(root_surface_IFZ_stat,~Layer*Rotation))

#contrast(emmeans(root_weight_IFZ_stat,~Layer*Rotation), method="pairwise", adjust = "none")
contrast(emmeans(root_weight_IFZ_stat,~Layer*Rotation), method="pairwise", adjust = "bonferroni")
contrast(emmeans(root_weight_IFZ_stat,~Rotation), method="pairwise", adjust = "bonferroni")
summary(emmeans(root_weight_IFZ_stat,~Layer*Rotation))
plot(emmeans(root_weight_IFZ_stat,~Layer*Rotation))

#contrast(emmeans(root_weight_m2_IFZ_stat,~Layer*Rotation), method="pairwise", adjust = "none")
contrast(emmeans(root_weight_m2_IFZ_stat,~Layer*Rotation), method="pairwise", adjust = "bonferroni")
contrast(emmeans(root_weight_m2_IFZ_stat,~Rotation), method="pairwise", adjust = "bonferroni")
summary(emmeans(root_weight_m2_IFZ_stat,~Layer*Rotation))
plot(emmeans(root_weight_m2_IFZ_stat,~Layer*Rotation))

#contrast(emmeans(root_weight_m3_IFZ_stat,~Layer*Rotation), method="pairwise", adjust = "none")
contrast(emmeans(root_weight_m3_IFZ_stat,~Layer*Rotation), method="pairwise", adjust = "bonferroni")
contrast(emmeans(root_weight_m3_IFZ_stat,~Rotation), method="pairwise", adjust = "bonferroni")
summary(emmeans(root_weight_m3_IFZ_stat,~Layer*Rotation))
plot(emmeans(root_weight_m3_IFZ_stat,~Layer*Rotation))



######## Plots

## Weight JKI
finalplot_root_weight_JKI_stat <- ggplot(data=root_weight_JKI_stat, aes(x=Layer, y=Roots, fill=Layer)) + 
  theme_bw() +
  facet_wrap(~Rotation, scale = "free_x") +
  geom_boxplot(alpha=0.8, outlier.shape=NA, width = 0.6, aes_string(x ="Layer", color = "Layer")) +
  #geom_jitter(width=0.1, size = 1) +
  scale_color_manual(values= c("#dcc8ba", "#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) + 
  scale_fill_manual(values= c("#dcc8ba", "#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) +
  xlab("Soil layer") + ylab("Root weight (g fw)") +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by = 0.5), label = c("0", "0.5", "1.0", "1.5", "2.0", "2.5")) +
  theme(legend.position = 'right') +
  theme(axis.text.x = element_text(angle = 90, size=18, color="black", vjust = 0.5)) +
  theme(axis.text.y = element_text(size=18, color="black")) +
  theme(axis.title = element_text(size=20, color="black")) +
  theme(strip.text.x = element_text(size = 14, face="bold"))
finalplot_root_weight_JKI_stat

ggsave("finalplot_root_weight_JKI_stat.png", path = "~/Documents/R_analysis_2/jki_seq1/Stats/output/", width = 32, height = 12, units = "cm",dpi = 300)

## Length IFZ
finalplot_root_length_IFZ_stat <- ggplot(data=root_length_IFZ_stat, aes(x=Layer, y=Length, fill=Layer)) + 
  theme_bw() +
  facet_wrap(~Rotation, scale = "free_x") +
  geom_boxplot(alpha=0.8, outlier.shape=NA, width = 0.6, aes_string(x ="Layer", color = "Layer")) +
  #geom_jitter(width=0.1, size = 1) +
  scale_color_manual(values= c("#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) + 
  scale_fill_manual(values= c("#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) +
  xlab("Soil layer") + ylab("Root length (cm/cm³)") +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 0.5), label = c("0", "0.5", "1.0", "1.5", "2.0")) +
  theme(legend.position = 'right') +
  theme(axis.text.x = element_text(angle = 90, size=18, color="black", vjust = 0.5)) +
  theme(axis.text.y = element_text(size=18, color="black")) +
  theme(axis.title = element_text(size=20, color="black")) +
  theme(strip.text.x = element_text(size = 14, face="bold"))
finalplot_root_length_IFZ_stat

ggsave("finalplot_root_length_IFZ_stat.png", path = "~/Documents/R_analysis_2/jki_seq1/Stats/output/", width = 32, height = 12, units = "cm",dpi = 300)

## Surface IFZ 
finalplot_root_surface_IFZ_stat <- ggplot(data=root_surface_IFZ_stat, aes(x=Layer, y=Surface, fill=Layer)) + 
  theme_bw() +
  facet_wrap(~Rotation, scale = "free_x") +
  geom_boxplot(alpha=0.8, outlier.shape=NA, width = 0.6, aes_string(x ="Layer", color = "Layer")) +
  #geom_jitter(width=0.1, size = 1) +
  #geom_text(data=CIs_root_length_IFZ_stat, aes(y = 2.1, label = trimws(.group)), size = 5) +
  scale_color_manual(values= c("#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) + 
  scale_fill_manual(values= c("#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) +
  xlab("Soil layer") + ylab("Root surface (cm²/cm³)") +
  scale_y_continuous(limits = c(0, 0.06), breaks = seq(0, 0.06, by = 0.02), label = c("0", "0.02", "0.04", "0.06")) +
  theme(legend.position = '') +
  theme(axis.text.x = element_text(angle = 90, size=18, color="black", vjust = 0.5)) +
  theme(axis.text.y = element_text(size=18, color="black")) +
  theme(axis.title = element_text(size=20, color="black")) +
  theme(strip.text.x = element_text(size = 14, face="bold"))
finalplot_root_surface_IFZ_stat

ggsave("finalplot_root_surface_IFZ_stat.png", path = "~/Documents/R_analysis_2/jki_seq1/Stats/output/", width = 32, height = 12, units = "cm",dpi = 300)

## Weight IFZ 
finalplot_root_weight_IFZ_stat <- ggplot(data=root_weight_IFZ_stat, aes(x=Layer, y=Weight, fill=Layer)) + 
  theme_bw() +
  facet_wrap(~Rotation, scale = "free_x") +
  geom_boxplot(alpha=0.8, outlier.shape=NA, width = 0.6, aes_string(x ="Layer", color = "Layer")) +
  #geom_jitter(width=0.1, size = 1) +
  #geom_text(data=CIs_root_length_IFZ_stat, aes(y = 2.1, label = trimws(.group)), size = 5) +
  scale_color_manual(values= c("#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) + 
  scale_fill_manual(values= c("#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) +
  xlab("Soil layer") + ylab("Root weight (g dw)") +
  scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, by = 0.02), label = c("0", "0.02", "0.04", "0.06", "0.08")) +
  theme(legend.position = '') +
  theme(axis.text.x = element_text(angle = 90, size=18, color="black", vjust = 0.5)) +
  theme(axis.text.y = element_text(size=18, color="black")) +
  theme(axis.title = element_text(size=20, color="black")) +
  theme(strip.text.x = element_text(size = 14, face="bold"))
finalplot_root_weight_IFZ_stat

ggsave("finalplot_root_weight_IFZ_stat.png", path = "~/Documents/R_analysis_2/jki_seq1/Stats/output/", width = 32, height = 12, units = "cm",dpi = 300)

## Weight_m2 IFZ 
finalplot_root_weight_m2_IFZ_stat <- ggplot(data=root_weight_m2_IFZ_stat, aes(x=Layer, y=Weight_m2, fill=Layer)) + 
  theme_bw() +
  facet_wrap(~Rotation, scale = "free_x") +
  geom_boxplot(alpha=0.8, outlier.shape=NA, width = 0.6, aes_string(x ="Layer", color = "Layer")) +
  #geom_jitter(width=0.1, size = 1) +
  #geom_text(data=CIs_root_length_IFZ_stat, aes(y = 2.1, label = trimws(.group)), size = 5) +
  scale_color_manual(values= c("#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) + 
  scale_fill_manual(values= c("#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) +
  xlab("Soil layer") + ylab("Root weight (g dw/m²)") +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5), label = c("0", "5", "10", "15", "20", "25")) +
  theme(legend.position = '') +
  theme(axis.text.x = element_text(angle = 90, size=18, color="black", vjust = 0.5)) +
  theme(axis.text.y = element_text(size=18, color="black")) +
  theme(axis.title = element_text(size=20, color="black")) +
  theme(strip.text.x = element_text(size = 14, face="bold"))
finalplot_root_weight_m2_IFZ_stat

ggsave("finalplot_root_weight_m2_IFZ_stat.png", path = "~/Documents/R_analysis_2/jki_seq1/Stats/output/", width = 32, height = 12, units = "cm",dpi = 300)

## Weight_m3 IFZ 
finalplot_root_weight_m3_IFZ_stat <- ggplot(data=root_weight_m3_IFZ_stat, aes(x=Layer, y=Weight_m3, fill=Layer)) + 
  theme_bw() +
  facet_wrap(~Rotation, scale = "free_x") +
  geom_boxplot(alpha=0.8, outlier.shape=NA, width = 0.6, aes_string(x ="Layer", color = "Layer")) +
  #geom_jitter(width=0.1, size = 1) +
  #geom_text(data=CIs_root_length_IFZ_stat, aes(y = 2.1, label = trimws(.group)), size = 5) +
  scale_color_manual(values= c("#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) + 
  scale_fill_manual(values= c("#DEB3AD", "#DE847B", "#B95C50", "#3B0404")) +
  xlab("Soil layer") + ylab("Root weight (g dw/m³)") +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10), label = c("0", "10", "20", "35", "40", "50", "60")) +
  theme(legend.position = '') +
  theme(axis.text.x = element_text(angle = 90, size=18, color="black", vjust = 0.5)) +
  theme(axis.text.y = element_text(size=18, color="black")) +
  theme(axis.title = element_text(size=20, color="black")) +
  theme(strip.text.x = element_text(size = 14, face="bold"))
finalplot_root_weight_m3_IFZ_stat

ggsave("finalplot_root_weight_m3_IFZ_stat.png", path = "~/Documents/R_analysis_2/jki_seq1/Stats/output/", width = 32, height = 12, units = "cm",dpi = 300)


#arrange multiple ggplots in one single figure ALL
###########################################
roots_data <- ggarrange(finalplot_root_weight_JKI_stat,
                        finalplot_root_weight_IFZ_stat,
                        finalplot_root_length_IFZ_stat,
                        finalplot_root_surface_IFZ_stat,
                        ncol =4, nrow =1, widths = c(1.25, 1, 1, 1), common.legend = TRUE, legend="none", label.y = 1, align = c("v"))
ggsave("roots_data.png", path = "~/Documents/R_analysis_2/jki_seq1/Stats/output/", width = 46, height = 12, units = "cm",dpi = 300)
###########################################
