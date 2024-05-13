###Select groups and subset data from the cleaned dataset
##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")

###Subseting samples from cleaned dataset to keep only samples that represent every case on a specific variable, using subset_samples()
#Place
sample_variables(psO_perifiton18S)
psO_perifiton18S_alveolata_jangadeiros <- subset_samples(psO_perifiton18S_alveolata, Place=="Jangadeiros")
psO_perifiton18S_alveolata_jangadeiros

psO_perifiton18S_alveolata_veleiros <- subset_samples(psO_perifiton18S_alveolata, Place=="Veleiros")
psO_perifiton18S_alveolata_veleiros

#Time
psO_perifiton18S_alveolata_time1 <- subset_samples(psO_perifiton18S_alveolata, Time=="T1")
psO_perifiton18S_alveolata_time1

psO_perifiton18S_alveolata_time2 <- subset_samples(psO_perifiton18S_alveolata, Time=="T2")
psO_perifiton18S_alveolata_time2

psO_perifiton18S_alveolata_time3 <- subset_samples(psO_perifiton18S_alveolata, Time=="T3")
psO_perifiton18S_alveolata_time3

psO_perifiton18S_alveolata_time4 <- subset_samples(psO_perifiton18S_alveolata, Time=="T4")
psO_perifiton18S_alveolata_time4

##Filtering table removing ASV that do not have counts in any sample on the group subset, using prune_taxa()
#for psO_perifiton18S_jangadeiros (coded group jangadeiros)
psO_perifiton18S_alveolata_jangadeiros_filt<-prune_taxa(taxa_sums(psO_perifiton18S_alveolata_jangadeiros) > 0, psO_perifiton18S_alveolata_jangadeiros)
psO_perifiton18S_alveolata_jangadeiros_filt

#for psO_perifiton18S_veleiros (coded group veleiros)
psO_perifiton18S_alveolata_veleiros_filt<-prune_taxa(taxa_sums(psO_perifiton18S_alveolata_veleiros) > 0, psO_perifiton18S_alveolata_veleiros)
psO_perifiton18S_alveolata_veleiros_filt

#for psO_perifiton18S_time1 (coded group T1)
psO_perifiton18S_alveolata_time1_filt<-prune_taxa(taxa_sums(psO_perifiton18S_alveolata_time1) > 0, psO_perifiton18S_alveolata_time1)
psO_perifiton18S_alveolata_time1_filt

#for psO_perifiton18S_time2 (coded group T2)
psO_perifiton18S_alveolata_time2_filt<-prune_taxa(taxa_sums(psO_perifiton18S_alveolata_time2) > 0, psO_perifiton18S_alveolata_time2)
psO_perifiton18S_alveolata_time2_filt

#for psO_perifiton18S_time3 (coded group T3)
psO_perifiton18S_alveolata_time3_filt<-prune_taxa(taxa_sums(psO_perifiton18S_alveolata_time3) > 0, psO_perifiton18S_alveolata_time3)
psO_perifiton18S_alveolata_time3_filt

#for psO_perifiton18S_time4 (coded group T4)
psO_perifiton18S_alveolata_time4_filt<-prune_taxa(taxa_sums(psO_perifiton18S_alveolata_time4) > 0, psO_perifiton18S_alveolata_time4)
psO_perifiton18S_alveolata_time4_filt
