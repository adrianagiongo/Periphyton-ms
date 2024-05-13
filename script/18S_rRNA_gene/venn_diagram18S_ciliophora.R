#Create a venndiagram 

#loading packages
library("phyloseq")
library("limma")
library("VennDiagram")
library("ggplot2")

###############
#Venndiagram JANGADEIROS

#verify name of variables
sample_variables(psO_perifiton18S_ciliophora_jangadeiros_filt_annotation)

#merge samples by time
venndiagram_psO_perifiton18S_ciliophora_jangadeiros_filt_annotation_time <- merge_samples(psO_perifiton18S_ciliophora_jangadeiros_filt_annotation, "Time")
sample_data(venndiagram_psO_perifiton18S_ciliophora_jangadeiros_filt_annotation_time)
venndiagram_psO_perifiton18S_ciliophora_jangadeiros_filt_annotation_time

#create the object to calculate the variable intersections
table_venndiagram_psO_perifiton18S_ciliophora_jangadeiros_filt_annotation_time <- t(otu_table(venndiagram_psO_perifiton18S_ciliophora_jangadeiros_filt_annotation_time))
table_venndiagram_psO_perifiton18S_ciliophora_jangadeiros_filt_annotation_time

venn_counts_ciliophora_jangadeiros_18S_time <- vennCounts(table_venndiagram_psO_perifiton18S_ciliophora_jangadeiros_filt_annotation_time)
venn_counts_ciliophora_jangadeiros_18S_time

#plot Venndiagram
#vennDiagram(venn_counts_ciliophora_jangadeiros_18S_time, cex=c(1,1.2,0.8,), names = c("T1", "T2", "T3", "T4"),circle.col = c("#26a4b5", "#014560"))

venn_plot_ciliophora_jangadeiros_18S_time <- draw.quad.venn(
  area1 = 24,
  area2 = 23,
  area3 = 36,
  area4 = 24,
  n12 = 11,
  n13 = 20,
  n14 = 15,
  n23 = 18,
  n24 = 9,
  n34 = 19,
  n123 = 11,
  n124 = 8,
  n134 = 14,
  n234 = 9,
  n1234 = 8,
  category = c("T1", "T2", "T3", "T4"),
  fill = c("#23207C", "#00B5ED", "#F08C01", "#D6000A"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#23207C", "#00B5ED", "#F08C01", "#D6000A")
);

tiff(filename = "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Venn_diagram_18S_Ciliophora/venn_plot_ciliophora_jangadeiros_18S_time.tiff", width = 12, height = 12, units = "cm", res = 300, compression = "lzw");
grid.draw(venn_plot_ciliophora_jangadeiros_18S_time);
dev.off()

###############
#Venndiagram VELEIROS

#verify name of variables
sample_variables(psO_perifiton18S_ciliophora_veleiros_filt_annotation)

#merge samples by time
venndiagram_psO_perifiton18S_ciliophora_veleiros_filt_annotation_time = merge_samples(psO_perifiton18S_ciliophora_veleiros_filt_annotation, "Time")
sample_data(venndiagram_psO_perifiton18S_ciliophora_veleiros_filt_annotation_time)
venndiagram_psO_perifiton18S_ciliophora_veleiros_filt_annotation_time

#create the object to calculate the variable intersections
table_venndiagram_psO_perifiton18S_ciliophora_veleiros_filt_annotation_time <- t(otu_table(venndiagram_psO_perifiton18S_ciliophora_veleiros_filt_annotation_time))
table_venndiagram_psO_perifiton18S_ciliophora_veleiros_filt_annotation_time

venn_counts_ciliophora_veleiros_18S_time <- vennCounts(table_venndiagram_psO_perifiton18S_ciliophora_veleiros_filt_annotation_time)
venn_counts_ciliophora_veleiros_18S_time

#plot Venndiagram
#vennDiagram(venn_counts_ciliophora_veleiros_18S_time, cex=c(1,1.2), names = c("JANGADEIROS", "VELEIROS"), circle.col = c("#26a4b5", "#014560"))

venn_plot_ciliophora_veleiros_18S_time <- draw.quad.venn(
  area1 = 14,
  area2 = 22,
  area3 = 13,
  area4 = 18,
  n12 = 12,
  n13 = 10,
  n14 = 11,
  n23 = 12,
  n24 = 15,
  n34 = 11,
  n123 = 10,
  n124 = 11,
  n134 = 9,
  n234 = 11,
  n1234 = 9,
  category = c("T1", "T2", "T3", "T4"),
  fill = c("#23207C", "#00B5ED", "#F08C01", "#D6000A"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#23207C", "#00B5ED", "#F08C01", "#D6000A")
);

tiff(filename = "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Venn_diagram_18S_Ciliophora/venn_plot_ciliophora_veleiros_18S_time.tiff" ,width = 12, height = 12, units = "cm", res = 300, compression = "lzw");
grid.draw(venn_plot_ciliophora_veleiros_18S_time);
dev.off()
