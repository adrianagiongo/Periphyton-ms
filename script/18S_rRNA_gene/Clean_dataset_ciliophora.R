##Supervised cleaning process to clean the raw dataset

#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")

#Show available ranks in the dataset
rank_names(Perifiton18S_ciliophora_data)

##Create table with number of features for each taxonomic level
##Remove respective artefacts on a specific taxonomic level

#Kingdom
table(tax_table(Perifiton18S_ciliophora_data)[,"Kingdom"], exclude = NULL)

#Phylum
table(tax_table(Perifiton18S_ciliophora_data)[,"Phylum"], exclude = NULL)

#Class_Subclass_curated
table(tax_table(Perifiton18S_ciliophora_data)[,"Class_subclass_curated"], exclude = NULL)
psO_perifiton18S_ciliophora <- subset_taxa(Perifiton18S_ciliophora_data, !is.na(Class_subclass_curated) & !Class_subclass_curated %in% c("Amoebozoa", "Bacillariophycea", "Bacteria", "Chlorophyta", "Cryptocaryon", "Diatomacea", "Fungi", "Haptophyta", "Heterolobosea", "Microsporidium", "Mollusca", "Ochrophyta", "NA", "Uncultured eukaryote"))
table(tax_table(psO_perifiton18S_ciliophora)[,"Class_subclass_curated"], exclude = NULL)

#Annotation_curated
table(tax_table(psO_perifiton18S_ciliophora)[,"Annotation_curated"], exclude = NULL)

##Observe psO after clean expurious taxa
Perifiton18S_ciliophora_data
psO_perifiton18S_ciliophora

#Compute prevalence of each feature and store as data.frame
prevdf_ciliophora = apply(X = otu_table(psO_perifiton18S_ciliophora),
               MARGIN = ifelse(taxa_are_rows(psO_perifiton18S_ciliophora), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to the data frame
prevdf_ciliophora <- data.frame(Prevalence=prevdf_ciliophora, TotalAbundance = taxa_sums(psO_perifiton18S_ciliophora), tax_table(psO_perifiton18S_ciliophora))

#Write data frame in csv format
write.csv(prevdf_ciliophora, "~/Documents/R_analysis_2/Perifiton/Perifiton_18S/output/output_Ciliophora/Tables_18S_Ciliophora/psO_perifiton18S_ciliophora_prevdf.csv")

