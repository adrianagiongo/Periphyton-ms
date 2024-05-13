##Supervised cleaning process to clean the raw dataset

#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")

#Show available ranks in the dataset
rank_names(Perifiton16S_data)

##Create table with number of features for each phylum
##Remove respective artefacts on a specific taxonomic level
#Kingdom
table(tax_table(Perifiton16S_data)[,"Kingdom"], exclude = NULL)
psO_perifiton16S <- subset_taxa(Perifiton16S_data, !is.na(Kingdom) & !Kingdom %in% c("Eukaryota", "Archaea"))
table(tax_table(psO_perifiton16S)[,"Kingdom"], exclude = NULL)

#Phylum
table(tax_table(psO_perifiton16S)[,"Phylum"], exclude = NULL)

#Order
table(tax_table(psO_perifiton16S)[,"Order"], exclude = NULL)
psO_perifiton16S <- subset_taxa(psO_perifiton16S, !is.na(Order) & !Order %in% c("Chloroplast"))
table(tax_table(psO_perifiton16S)[,"Order"], exclude = NULL)

#Family
table(tax_table(psO_perifiton16S)[,"Family"], exclude = NULL)
psO_perifiton16S <- subset_taxa(psO_perifiton16S, !is.na(Family) & !Family %in% c("Mitochondria"))
table(tax_table(psO_perifiton16S)[,"Family"], exclude = NULL)

#Genus
table(tax_table(psO_perifiton16S)[,"Genus"], exclude = NULL)

#Annotation
table(tax_table(psO_perifiton16S)[,"Annotation"], exclude = NULL)

##Observe psO after clean expurious taxa
Perifiton16S_data
psO_perifiton16S
table(tax_table(psO_perifiton16S)[,"Kingdom"], exclude = NULL)

#Compute prevalence of each feature and store as data.frame
prevdf = apply(X = otu_table(psO_perifiton16S),
               MARGIN = ifelse(taxa_are_rows(psO_perifiton16S), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to the data frame
prevdf_16S = data.frame(Prevalence=prevdf, TotalAbundance = taxa_sums(psO_perifiton16S), tax_table(psO_perifiton16S), otu_table(psO_perifiton16S))

#Write data frame in csv format
write.csv(prevdf_16S, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/psO_perifiton16S_prevdf.csv")

####################
#Calculate number of genera present after filtration
length(get_taxa_unique(psO_perifiton16S, taxonomic.rank = "Phylum"))
length(get_taxa_unique(psO_perifiton16S, taxonomic.rank = "Annotation"))

##Check number of taxa after aglomerate different taxa-levels using tax glom function and removing NAs
#for phylum
psO_perifiton16S_phylum = tax_glom(psO_perifiton16S, "Phylum", NArm= TRUE)
psO_perifiton16S_phylum
ntaxa(psO_perifiton16S); ntaxa(psO_perifiton16S_phylum)

df_psO_perifiton16S_phylum <- data.frame(tax_table(psO_perifiton16S_phylum), otu_table(psO_perifiton16S_phylum))
write.csv(df_psO_perifiton16S_phylum, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_phylum.csv")

#for annotation
psO_perifiton16S_annotation = tax_glom(psO_perifiton16S, "Annotation", NArm= TRUE)
psO_perifiton16S_annotation
ntaxa(psO_perifiton16S); ntaxa(psO_perifiton16S_annotation)

df_psO_perifiton16S_annotation <- data.frame(tax_table(psO_perifiton16S_annotation), otu_table(psO_perifiton16S_annotation))
write.csv(df_psO_perifiton16S_annotation, "~/Documents/R_analysis_2/Perifiton/Perifiton_16S/output/Tables_Perifiton_16S/df_psO_perifiton16S_annotation.csv")


