#Create a venndiagram 

#loading packages
library("phyloseq")
library("limma")
library("VennDiagram")
library("ggplot2")

###############
#Venndiagram JANGADEIROS

#verify name of variables
sample_variables(psO_perifiton16S_jangadeiros_filt_annotation)

#merge samples by time
venndiagram_psO_perifiton16S_jangadeiros_filt_annotation_time <- merge_samples(psO_perifiton16S_jangadeiros_filt_annotation, "Time")
sample_data(venndiagram_psO_perifiton16S_jangadeiros_filt_annotation_time)
venndiagram_psO_perifiton16S_jangadeiros_filt_annotation_time

#create the object to calculate the variable intersections
table_venndiagram_psO_perifiton16S_jangadeiros_filt_annotation_time <- t(otu_table(venndiagram_psO_perifiton16S_jangadeiros_filt_annotation_time))
table_venndiagram_psO_perifiton16S_jangadeiros_filt_annotation_time

venn_counts_jangadeiros_16S_time <- vennCounts(table_venndiagram_psO_perifiton16S_jangadeiros_filt_annotation_time)
venn_counts_jangadeiros_16S_time

#plot Venndiagram
venn_plot_jangadeiros_16S_time <- draw.quad.venn(
  area1 = 176,
  area2 = 147,
  area3 = 215,
  area4 = 109,
  n12 = 72,
  n13 = 89,
  n14 = 53,
  n23 = 91,
  n24 = 60,
  n34 = 73,
  n123 = 60,
  n124 = 44,
  n134 = 48,
  n234 = 53,
  n1234 = 43,
  category = c("T1", "T2", "T3", "T4"),
  fill = c("#23207C", "#00B5ED", "#F08C01", "#D6000A"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#23207C", "#00B5ED", "#F08C01", "#D6000A")
);

tiff(filename = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/VennDiagram/venn_plot_jangadeiros_16S_time.tiff", width = 12, height = 12, units = "cm", res = 300, compression = "lzw");
grid.draw(venn_plot_jangadeiros_16S_time);
dev.off()

###############
#Venndiagram VELEIROS

#verify name of variables
sample_variables(psO_perifiton16S_veleiros_filt_annotation)

#merge samples by time
venndiagram_psO_perifiton16S_veleiros_filt_annotation_time = merge_samples(psO_perifiton16S_veleiros_filt_annotation, "Time")
sample_data(venndiagram_psO_perifiton16S_veleiros_filt_annotation_time)
venndiagram_psO_perifiton16S_veleiros_filt_annotation_time

#create the object to calculate the variable intersections
table_venndiagram_psO_perifiton16S_veleiros_filt_annotation_time <- t(otu_table(venndiagram_psO_perifiton16S_veleiros_filt_annotation_time))
table_venndiagram_psO_perifiton16S_veleiros_filt_annotation_time

venn_counts_veleiros_16S_time <- vennCounts(table_venndiagram_psO_perifiton16S_veleiros_filt_annotation_time)
venn_counts_veleiros_16S_time

#plot Venndiagram
venn_plot_veleiros_16S_time <- draw.quad.venn(
  area1 = 118,
  area2 = 110,
  area3 = 73,
  area4 = 117,
  n12 = 54,
  n13 = 48,
  n14 = 65,
  n23 = 36,
  n24 = 50,
  n34 = 56,
  n123 = 30,
  n124 = 41,
  n134 = 43,
  n234 = 31,
  n1234 = 28,
  category = c("T1", "T2", "T3", "T4"),
  fill = c("#23207C", "#00B5ED", "#F08C01", "#D6000A"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#23207C", "#00B5ED", "#F08C01", "#D6000A")
);

tiff(filename = "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/VennDiagram/venn_plot_veleiros_16S_time.tiff", width = 12, height = 12, units = "cm", res = 300, compression = "lzw");
grid.draw(venn_plot_veleiros_16S_time);
dev.off()

####################
psO_perifiton16S_jangadeiros_filt_annotation_rel <- microbiome::transform(psO_perifiton16S_jangadeiros_filt_annotation, "compositional")
df_psO_perifiton16S_jangadeiros_filt_annotation_rel <- data.frame(tax_table(psO_perifiton16S_jangadeiros_filt_annotation_rel), otu_table(psO_perifiton16S_jangadeiros_filt_annotation_rel))
write.csv(df_psO_perifiton16S_jangadeiros_filt_annotation_rel, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_jangadeiros_filt_annotation_rel.csv")

psO_perifiton16S_veleiros_filt_annotation_rel <- microbiome::transform(psO_perifiton16S_veleiros_filt_annotation, "compositional")
df_psO_perifiton16S_veleiros_filt_annotation_rel <- data.frame(tax_table(psO_perifiton16S_veleiros_filt_annotation_rel), otu_table(psO_perifiton16S_veleiros_filt_annotation_rel))
write.csv(df_psO_perifiton16S_veleiros_filt_annotation_rel, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_veleiros_filt_annotation_rel.csv")

