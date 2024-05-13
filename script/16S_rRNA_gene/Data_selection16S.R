###Select groups and subset data from the cleaned dataset
##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")

###Subseting samples from cleaned dataset to keep only samples that represent every case on a specific variable, using subset_samples()
#Place
sample_variables(psO_perifiton16S)
psO_perifiton16S_jangadeiros <- subset_samples(psO_perifiton16S, Place=="Jangadeiros")
psO_perifiton16S_jangadeiros

psO_perifiton16S_veleiros <- subset_samples(psO_perifiton16S, Place=="Veleiros")
psO_perifiton16S_veleiros

#Time
psO_perifiton16S_time1 <- subset_samples(psO_perifiton16S, Time=="T1")
psO_perifiton16S_time1

psO_perifiton16S_time2 <- subset_samples(psO_perifiton16S, Time=="T2")
psO_perifiton16S_time2

psO_perifiton16S_time3 <- subset_samples(psO_perifiton16S, Time=="T3")
psO_perifiton16S_time3

psO_perifiton16S_time4 <- subset_samples(psO_perifiton16S, Time=="T4")
psO_perifiton16S_time4

##Filtering table removing ASV that do not have counts in any sample on the group subset, using prune_taxa()
#for psO_perifiton16S_jangadeiros (coded group jangadeiros)
psO_perifiton16S_jangadeiros_filt<-prune_taxa(taxa_sums(psO_perifiton16S_jangadeiros) > 0, psO_perifiton16S_jangadeiros)
psO_perifiton16S_jangadeiros_filt

#for psO_perifiton16S_veleiros (coded group veleiros)
psO_perifiton16S_veleiros_filt<-prune_taxa(taxa_sums(psO_perifiton16S_veleiros) > 0, psO_perifiton16S_veleiros)
psO_perifiton16S_veleiros_filt

#for psO_perifiton16S_time1 (coded group T1)
psO_perifiton16S_time1_filt<-prune_taxa(taxa_sums(psO_perifiton16S_time1) > 0, psO_perifiton16S_time1)
psO_perifiton16S_time1_filt

#for psO_perifiton16S_time2 (coded group T2)
psO_perifiton16S_time2_filt<-prune_taxa(taxa_sums(psO_perifiton16S_time2) > 0, psO_perifiton16S_time2)
psO_perifiton16S_time2_filt

#for psO_perifiton16S_time3 (coded group T3)
psO_perifiton16S_time3_filt<-prune_taxa(taxa_sums(psO_perifiton16S_time3) > 0, psO_perifiton16S_time3)
psO_perifiton16S_time3_filt

#for psO_perifiton16S_time4 (coded group T4)
psO_perifiton16S_time4_filt<-prune_taxa(taxa_sums(psO_perifiton16S_time4) > 0, psO_perifiton16S_time4)
psO_perifiton16S_time4_filt

##Agglomerating taxa per sample at the appropriated taxonomic rank, using tax_glom()
#For Jangadeiros
colnames(tax_table(psO_perifiton16S_jangadeiros_filt))
psO_perifiton16S_jangadeiros_filt_annotation <- tax_glom(psO_perifiton16S_jangadeiros_filt, taxrank = "Annotation")
ntaxa(psO_perifiton16S_jangadeiros_filt); ntaxa(psO_perifiton16S_jangadeiros_filt_annotation)

df_psO_perifiton16S_jangadeiros_filt_annotation <- data.frame(tax_table(psO_perifiton16S_jangadeiros_filt_annotation), otu_table(psO_perifiton16S_jangadeiros_filt_annotation))
write.csv(df_psO_perifiton16S_jangadeiros_filt_annotation, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_jangadeiros_filt_annotation.csv")

#For Veleiros
colnames(tax_table(psO_perifiton16S_veleiros_filt))
psO_perifiton16S_veleiros_filt_annotation <- tax_glom(psO_perifiton16S_veleiros_filt, taxrank = "Annotation")
ntaxa(psO_perifiton16S_veleiros_filt); ntaxa(psO_perifiton16S_veleiros_filt_annotation)

df_psO_perifiton16S_veleiros_filt_annotation <- data.frame(tax_table(psO_perifiton16S_veleiros_filt_annotation), otu_table(psO_perifiton16S_veleiros_filt_annotation))
write.csv(df_psO_perifiton16S_veleiros_filt_annotation, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_veleiros_filt_annotation.csv")

#For T1
colnames(tax_table(psO_perifiton16S_time1_filt))
psO_perifiton16S_time1_filt_annotation <- tax_glom(psO_perifiton16S_time1_filt, taxrank = "Annotation")
ntaxa(psO_perifiton16S_time1_filt); ntaxa(psO_perifiton16S_time1_filt_annotation)

df_psO_perifiton16S_time1_filt_annotation <- data.frame(tax_table(psO_perifiton16S_time1_filt_annotation), otu_table(psO_perifiton16S_time1_filt_annotation))
write.csv(df_psO_perifiton16S_time1_filt_annotation, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_time1_filt_annotation.csv")

#For T2
colnames(tax_table(psO_perifiton16S_time2_filt))
psO_perifiton16S_time2_filt_annotation <- tax_glom(psO_perifiton16S_time2_filt, taxrank = "Annotation")
ntaxa(psO_perifiton16S_time2_filt); ntaxa(psO_perifiton16S_time2_filt_annotation)

df_psO_perifiton16S_time2_filt_annotation <- data.frame(tax_table(psO_perifiton16S_time2_filt_annotation), otu_table(psO_perifiton16S_time2_filt_annotation))
write.csv(df_psO_perifiton16S_time2_filt_annotation, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_time2_filt_annotation.csv")

#For T3
colnames(tax_table(psO_perifiton16S_time3_filt))
psO_perifiton16S_time3_filt_annotation <- tax_glom(psO_perifiton16S_time3_filt, taxrank = "Annotation")
ntaxa(psO_perifiton16S_time3_filt); ntaxa(psO_perifiton16S_time3_filt_annotation)

df_psO_perifiton16S_time3_filt_annotation <- data.frame(tax_table(psO_perifiton16S_time3_filt_annotation), otu_table(psO_perifiton16S_time3_filt_annotation))
write.csv(df_psO_perifiton16S_time3_filt_annotation, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_time3_filt_annotation.csv")

#For T4
colnames(tax_table(psO_perifiton16S_time4_filt))
psO_perifiton16S_time4_filt_annotation <- tax_glom(psO_perifiton16S_time4_filt, taxrank = "Annotation")
ntaxa(psO_perifiton16S_time4_filt); ntaxa(psO_perifiton16S_time4_filt_annotation)

df_psO_perifiton16S_time4_filt_annotation <- data.frame(tax_table(psO_perifiton16S_time4_filt_annotation), otu_table(psO_perifiton16S_time4_filt_annotation))
write.csv(df_psO_perifiton16S_time4_filt_annotation, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_time4_filt_annotation.csv")


