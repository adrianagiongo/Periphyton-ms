###Select groups and subset data from the cleaned dataset
##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")

##Transform count data to relative abundance
#Jangadeiros
psO_perifiton16S_jangadeiros_filt_annotation_rel <- microbiome::transform(psO_perifiton16S_jangadeiros_filt_annotation, "composition")

df_psO_perifiton16S_jangadeiros_filt_annotation_rel <- data.frame(tax_table(psO_perifiton16S_jangadeiros_filt_annotation_rel), otu_table(psO_perifiton16S_jangadeiros_filt_annotation_rel))
write.csv(df_psO_perifiton16S_jangadeiros_filt_annotation_rel, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_jangadeiros_filt_annotation_rel.csv")

#Veleiros
psO_perifiton16S_veleiros_filt_annotation_rel <- microbiome::transform(psO_perifiton16S_veleiros_filt_annotation, "composition")

df_psO_perifiton16S_veleiros_filt_annotation_rel <- data.frame(tax_table(psO_perifiton16S_veleiros_filt_annotation_rel), otu_table(psO_perifiton16S_veleiros_filt_annotation_rel))
write.csv(df_psO_perifiton16S_veleiros_filt_annotation_rel, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_veleiros_filt_annotation_rel.csv")

