###Select groups and subset data from the cleaned dataset
##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")

###Subseting samples from cleaned dataset to keep only samples that represent every case on a specific variable, using subset_samples()
#Place
sample_variables(psO_perifiton18S_ciliophora)
psO_perifiton18S_ciliophora_jangadeiros <- subset_samples(psO_perifiton18S_ciliophora, Place=="Jangadeiros")
psO_perifiton18S_ciliophora_jangadeiros

psO_perifiton18S_ciliophora_veleiros <- subset_samples(psO_perifiton18S_ciliophora, Place=="Veleiros")
psO_perifiton18S_ciliophora_veleiros

#for psO_perifiton18S_jangadeiros (coded group jangadeiros)
psO_perifiton18S_ciliophora_jangadeiros_filt<-prune_taxa(taxa_sums(psO_perifiton18S_ciliophora_jangadeiros) > 0, psO_perifiton18S_ciliophora_jangadeiros)
psO_perifiton18S_ciliophora_jangadeiros_filt

#for psO_perifiton18S_jangadeiros (coded group jangadeiros)
psO_perifiton18S_ciliophora_veleiros_filt<-prune_taxa(taxa_sums(psO_perifiton18S_ciliophora_veleiros) > 0, psO_perifiton18S_ciliophora_veleiros)
psO_perifiton18S_ciliophora_veleiros_filt

####JANGADEIROS
#Aglomerate taxa using tax glom function and removing NAs
psO_perifiton18S_ciliophora_jangadeiros_filt_annotation = tax_glom(psO_perifiton18S_ciliophora_jangadeiros_filt, "Annotation_curated", NArm= TRUE)

#Transfom data to square root abundance
psO_perifiton18S_ciliophora_jangadeiros_sqr <- microbiome::transform(psO_perifiton18S_ciliophora_jangadeiros_filt_annotation, "hellinger")

psO_perifiton18S_ciliophora_jangadeiros_rel <- microbiome::transform(psO_perifiton18S_ciliophora_jangadeiros_filt_annotation, "compositional")
df_psO_perifiton18S_ciliophora_jangadeiros_rel <- data.frame(tax_table(psO_perifiton18S_ciliophora_jangadeiros_rel), otu_table(psO_perifiton18S_ciliophora_jangadeiros_rel))
write.csv(df_psO_perifiton18S_ciliophora_jangadeiros_rel, "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Tables_18S_Ciliophora/df_psO_perifiton18S_ciliophora_jangadeiros_rel.csv")

####VELEIROS
#Aglomerate taxa using tax glom function and removing NAs
psO_perifiton18S_ciliophora_veleiros_filt_annotation = tax_glom(psO_perifiton18S_ciliophora_veleiros_filt, "Annotation_curated", NArm= TRUE)

#Transfom data to square root abundance
psO_perifiton18S_ciliophora_veleiros_sqr <- microbiome::transform(psO_perifiton18S_ciliophora_veleiros_filt_annotation, "hellinger")

psO_perifiton18S_ciliophora_veleiros_rel <- microbiome::transform(psO_perifiton18S_ciliophora_veleiros_filt_annotation, "compositional")
df_psO_perifiton18S_ciliophora_veleiros_rel <- data.frame(tax_table(psO_perifiton18S_ciliophora_veleiros_rel), otu_table(psO_perifiton18S_ciliophora_veleiros_rel))
write.csv(df_psO_perifiton18S_ciliophora_veleiros_rel, "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Tables_18S_Ciliophora/df_psO_perifiton18S_ciliophora_veleiros_rel.csv")



