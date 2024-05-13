#load libraries

library("ggplot2")
library("tidyr")
library("readxl")
library("dplyr")
library("RColorBrewer")
#palette Paired - ("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

### JANGADEIROS
#load data 
ciliophora_jangadeiros_correl<-read_excel("~/Desktop/R_analysis/Perifiton_project/Perifiton_18S/data/Spearman_16S_18S_correlation_jangadeiros_table_abund1_heatmap.xlsx")
ciliophora_jangadeiros_correl

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
ciliophora_jangadeiros_correl_processed<-gather(data=ciliophora_jangadeiros_correl,key = Ciliophora, value = Result,Amphileptus:Uronema, -ID,-Isolate_id)
head(ciliophora_jangadeiros_correl_processed)

###Constructing the ggplot heatmap
#Convert variables on the x-axis to factors with desired order
#ciliophora_jangadeiros_correl_processed$Ciliophora <- factor(ciliophora_jangadeiros_correl_processed$Ciliophora, levels = c("Spores2Y", "Spores15Y", "Roots2Y","Roots15Y"))

#Convert variables on the x-axis to factors with desired order
ciliophora_jangadeiros_correl_processed$Isolate_id <- factor(ciliophora_jangadeiros_correl_processed$Isolate_id, levels = rev(sort(unique(ciliophora_jangadeiros_correl_processed$Isolate_id))))

heatmap_ciliophora_jangadeiros_correl_processed<-ggplot(data=ciliophora_jangadeiros_correl_processed,mapping=aes(x=Ciliophora, y=Isolate_id, fill=Result))
heatmap_ciliophora_jangadeiros_correl_processed + 
  geom_tile(colour="grey",size=0.25) +
  xlab(label="") +
  scale_fill_gradient2(name="Correlation",low="#E31A1C", mid="white", high="#1F78B4") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size=16, colour="black", angle = 90, hjust=1, vjust=0.5)) +
  theme(axis.text.y = element_text(size=16, colour="black")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size = 16, colour = "black")) +
  theme(axis.title.y = element_blank()) +
  ggtitle(label="")

ggsave("heatmap_ciliophora_jangadeiros_correl_processed.png", path = "~/Desktop/R_analysis/Perifiton_project/Perifiton_18S/output/Heatmap/", width = 25, height = 31, units = "cm",dpi = 600)


### VELEIROS
#load data 
ciliophora_veleiros_correl<-read_excel("~/Desktop/R_analysis/Perifiton_project/Perifiton_18S/data/Spearman_16S_18S_correlation_veleiros_table_abund1_heatmap.xlsx")
ciliophora_veleiros_correl

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
ciliophora_veleiros_correl_processed<-gather(data=ciliophora_veleiros_correl,key = Ciliophora, value = Result,Dexiotricha:Oxitricha, -ID,-Isolate_id)
head(ciliophora_veleiros_correl_processed)

###Constructing the ggplot heatmap
#Convert variables on the x-axis to factors with desired order
#ciliophora_veleiros_correl_processed$Ciliophora <- factor(ciliophora_veleiros_correl_processed$Ciliophora, levels = c("Spores2Y", "Spores15Y", "Roots2Y","Roots15Y"))

#Convert variables on the x-axis to factors with desired order
ciliophora_veleiros_correl_processed$Isolate_id <- factor(ciliophora_veleiros_correl_processed$Isolate_id, levels = rev(sort(unique(ciliophora_veleiros_correl_processed$Isolate_id))))

heatmap_ciliophora_veleiros_correl_processed<-ggplot(data=ciliophora_veleiros_correl_processed,mapping=aes(x=Ciliophora, y=Isolate_id, fill=Result))
heatmap_ciliophora_veleiros_correl_processed + 
  geom_tile(colour="grey",size=0.25) +
  xlab(label="") +
  scale_fill_gradient2(name="Correlation",low="#E31A1C", mid="white", high="#1F78B4") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size=14, colour="black", angle = 90, hjust=1, vjust=0.5)) +
  theme(axis.text.y = element_text(size=14, colour="black")) +
  theme(legend.title = element_text(size=14), legend.text = element_text(size = 14, colour = "black")) +
  theme(axis.title.y = element_blank()) +
  ggtitle(label="")

ggsave("heatmap_ciliophora_veleiros_correl_processed.png", path = "~/Desktop/R_analysis/Perifiton_project/Perifiton_18S/output/Heatmap/", width = 13, height = 12, units = "cm",dpi = 600)
