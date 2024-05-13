##Create a pie chart using ggplot2
#Load package
library("ggplot2")
library("ggtext")
library("RColorBrewer")
library("dplyr")

#Load data
counts_jangadeiros <- read.csv("~/Documents/R_analysis_2/Perifiton/Perifiton_18S/data/data_Ciliophora/data_Ciliophora_pie/Perifiton_counts_ciliophora_jangadeiros.csv",header = T,row.names = 1)
counts_veleiros <- read.csv("~/Documents/R_analysis_2/Perifiton/Perifiton_18S/data/data_Ciliophora/data_Ciliophora_pie/Perifiton_counts_ciliophora_veleiros.csv",header = T,row.names = 1)

#Turn relative abundance to numeric
counts_jangadeiros$Relative_abundance <- as.numeric(as.character(counts_jangadeiros$Relative_abundance))
counts_jangadeiros

counts_veleiros$Relative_abundance <- as.numeric(as.character(counts_veleiros$Relative_abundance))
counts_veleiros

##Plot a pie chart
#Concentric piechart jangadeiros
ggplot(counts_jangadeiros, aes(x= Time, y= Relative_abundance, fill=Species)) +
  geom_col() +
  coord_polar(theta = "y") +
  scale_fill_manual(name=NULL,
                    breaks = c("*Carchesium polypinum*", "Epistylis plicatilis", "Epistylis portoalegrensis", "Epistylis smalli", "Epistylis spp.", "*Opercularia* spp.", "Vaginicola tincta", "Vorticella campanula", "Vorticella convalaria", "Vorticella spp.", "Zoothamnium arbuscula", "Myoschistom duplicatum"),
                    values = c("#B2DF8A", "#85191E", "#B22128", "#DF2A33", "#E6565D", "#999999", "#6A3D9A", "#AE9E06", "#E0CB08", "#F9E209", "#A6CEE3", "#000000")
) +
  scale_x_discrete(breaks=c("T1", "T2", "T3", "T4")) +
  scale_y_continuous(breaks = NULL) +
  labs(x=NULL, y=NULL) +
  theme_void() +
  theme(legend.text = element_markdown(), legend.key.size = unit(10, "pt"), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_markdown(size = 8, hjust = 0, face = "bold"))

ggsave("piechart_concentric_ciliophora_counts_jangadeiros.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Piechart_18S_Ciliophora/", width = 12, height = 8, units = "cm", dpi = 300)


#Concentric piechart veleiros
ggplot(counts_veleiros, aes(x= Time, y= Relative_abundance, fill=Species)) +
  geom_col() +
  coord_polar(theta = "y") +
  scale_fill_manual(name=NULL,
                    breaks = c("*Carchesium polypinum*", "Epistylis plicatilis", "Epistylis portoalegrensis", "Epistylis smalli", "*Epistylis* spp.", "Opercularia spp.", "Vorticella campanula", "Vorticella convalaria", "Vorticella spp.", "Zoothamnium spp."),
                    values = c("#B2DF8A", "#85191E", "#B22128", "#DF2A33", "#E6565D", "#999999", "#AE9E06", "#E0CB08", "#F9E209", "#1F78B4")
  ) +
  scale_x_discrete(breaks=c("T1", "T2", "T3", "T4")) +
  scale_y_continuous(breaks = NULL) +
  labs(x=NULL, y=NULL) +
  theme_void() +
  theme(legend.text = element_markdown(), legend.key.size = unit(10, "pt"), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_markdown(size = 8, hjust = 0, face = "bold"))

ggsave("piechart_concentric_ciliophora_counts_veleiros.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Piechart_18S_Ciliophora/", width = 12, height = 8, units = "cm", dpi = 300)

#vertical piechart jangadeiros
ggplot(counts_jangadeiros, aes(x= 1, y= Relative_abundance, fill=Species)) +
  geom_col() +
  coord_polar(theta = "y") +
  facet_wrap(~Time, nrow = 4, strip.position = "left") +
  scale_fill_manual(name=NULL,
                    breaks = c("Carchesium polypinum", "Epistylis plicatilis", "Epistylis portoalegrensis", "Epistylis smalli", "Epistylis spp.", "Opercularia spp.", "Vaginicola tincta", "Vorticella campanula", "Vorticella convalaria", "Vorticella spp.", "Zoothamnium arbuscula", "Myoschistom duplicatum"),
                    values = c("#B2DF8A", "#85191E", "#B22128", "#DF2A33", "#E6565D", "#999999", "#6A3D9A", "#AE9E06", "#E0CB08", "#F9E209", "#A6CEE3", "#000000")
  ) +
  scale_x_discrete(breaks=c("T1", "T2", "T3", "T4")) +
  scale_y_continuous(breaks = NULL) +
  labs(x=NULL, y=NULL) +
  theme_void() +
  theme(legend.text = element_markdown(size = 16), legend.key.size = unit(20, "pt"), axis.line = element_blank(), axis.ticks = element_blank(), strip.text.y.left = element_markdown(angle = 0, size = 16, hjust = 0, face = "bold"))

ggsave("piechart_ciliophora_counts_jangadeirosa.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Piechart_18S_Ciliophora/", width = 12, height = 15, units = "cm", dpi = 600)

#vertical piechart veleiros
ggplot(counts_veleiros, aes(x= 1, y= Relative_abundance, fill=Species)) +
  geom_col() +
  coord_polar(theta = "y") +
  facet_wrap(~Time, nrow = 4, strip.position = "left") +
  scale_fill_manual(name=NULL,
                    breaks = c("Carchesium polypinum", "Epistylis plicatilis", "Epistylis portoalegrensis", "Epistylis smalli", "Epistylis spp.", "Opercularia spp.", "Vorticella campanula", "Vorticella convalaria", "Vorticella spp.", "Zoothamnium spp."),
                    values = c("#B2DF8A", "#85191E", "#B22128", "#DF2A33", "#E6565D", "#999999", "#AE9E06", "#E0CB08", "#F9E209", "#1F78B4")
  ) +
  scale_x_discrete(breaks=c("T1", "T2", "T3", "T4")) +
  scale_y_continuous(breaks = NULL) +
  labs(x=NULL, y=NULL) +
  theme_void() +
  theme(legend.text = element_markdown(size = 16), legend.key.size = unit(20, "pt"), axis.line = element_blank(), axis.ticks = element_blank(), strip.text.y.left = element_markdown(angle = 0, size = 16, hjust = 0, face = "bold"))

ggsave("piechart_ciliophora_counts_veleiros.png", path = "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Piechart_18S_Ciliophora/", width = 12, height = 15, units = "cm", dpi = 600)
