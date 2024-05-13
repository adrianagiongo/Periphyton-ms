##Supervised cleaning process to clean the raw dataset

#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")

#Show available ranks in the dataset
rank_names(Perifiton18S_data)

##Create table with number of features for each phylum
##Remove respective artefacts on a specific taxonomic level
#Kingdom
table(tax_table(Perifiton18S_data)[,"Kingdom"], exclude = NULL)
psO_perifiton18S <- subset_taxa(Perifiton18S_data, !is.na(Kingdom) & !Kingdom %in% c("Bacteria", "Archaea"))
table(tax_table(psO_perifiton18S)[,"Kingdom"], exclude = NULL)

#Clade
table(tax_table(psO_perifiton18S)[,"Clade"], exclude = NULL)
psO_perifiton18S_extended <- subset_taxa(psO_perifiton18S, !Clade %in% c("NA"))
table(tax_table(psO_perifiton18S_extended)[,"Clade"], exclude = NULL)
table(tax_table(psO_perifiton18S_extended)[,"Kingdom"], exclude = NULL)
table(tax_table(psO_perifiton18S_extended)[,"Phylum"], exclude = NULL)

#Phylum
psO_perifiton18S_alveolata <- subset_taxa(psO_perifiton18S, Phylum %in% c("Acavomonidia", "Apicomplexa", "Ciliophora", "Dinoflagellata","Protalveolata"))
table(tax_table(psO_perifiton18S_alveolata)[,"Phylum"], exclude = NULL)
table(tax_table(psO_perifiton18S_alveolata)[,"Kingdom"], exclude = NULL)

#Order
table(tax_table(psO_perifiton18S_alveolata)[,"Order"], exclude = NULL)

#Family
table(tax_table(psO_perifiton18S_alveolata)[,"Family"], exclude = NULL)

#Genus
table(tax_table(psO_perifiton18S_alveolata)[,"Genus"], exclude = NULL)

#Annotation
table(tax_table(psO_perifiton18S_alveolata)[,"Annotation"], exclude = NULL)

##Observe psO after clean spurious taxa
Perifiton18S_data
psO_perifiton18S_alveolata
psO_perifiton18S_extended

#Compute prevalence of each feature and store as data.frame
prevdf_extended = apply(X = otu_table(psO_perifiton18S_extended),
               MARGIN = ifelse(taxa_are_rows(psO_perifiton18S_extended), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

prevdf_alveolata = apply(X = otu_table(psO_perifiton18S_alveolata),
               MARGIN = ifelse(taxa_are_rows(psO_perifiton18S_alveolata), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to the data frame
prevdf_extended = data.frame(Prevalence=prevdf_extended, TotalAbundance = taxa_sums(psO_perifiton18S_extended), tax_table(psO_perifiton18S_extended))

prevdf_alveolata = data.frame(Prevalence=prevdf_alveolata, TotalAbundance = taxa_sums(psO_perifiton18S_alveolata), tax_table(psO_perifiton18S_alveolata))

#Write data frame in csv format
write.csv(prevdf_extended, "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/Tables_Perifiton_18S/psO_perifiton18S_prevdf_extended.csv")

write.csv(prevdf_alveolata, "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/Tables_Perifiton_18S/psO_perifiton18S_prevdf_alveolata.csv")
