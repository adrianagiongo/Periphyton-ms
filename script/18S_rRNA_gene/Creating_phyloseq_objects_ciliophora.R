#Transform dataset to phyloseq format

#load libraries
library("phyloseq")
library("ggplot2")
library("readxl")
library("dplyr")

#move excel files to R objects and metadata on .csv format with names on row 1. 
otu_ciliophora_mat<-read_excel("~/Documents/R_analysis_2/Perifiton/Perifiton_18S/data/data_Ciliophora/Perifiton_Ciliophora_curated_otu.xlsx")
taxo_ciliophora_mat<-read_excel("~/Documents/R_analysis_2/Perifiton/Perifiton_18S/data/data_Ciliophora/Perifiton_Ciliophora_curated_taxa.xlsx")
metadata_ciliophora<-read.csv("~/Documents/R_analysis_2/Perifiton/Perifiton_18S/data/data_Ciliophora/Perifiton_metadata.csv", row.names = 1)

#define the row names from the otu column under header "OTU"
#remove the column otu from otu, taxo and metadata tables since it is now used as a row name
row.names(otu_ciliophora_mat)<-otu_ciliophora_mat$OTU
otu_ciliophora_mat_unamed<-otu_ciliophora_mat %>% select (-OTU)

row.names(taxo_ciliophora_mat)<-taxo_ciliophora_mat$OTU
taxo_ciliophora_mat_unamed<-taxo_ciliophora_mat %>% select (-OTU)

#transform otu and taxo tables into matrixes
otu_ciliophora_matrix<-as.matrix(otu_ciliophora_mat_unamed)
taxo_ciliophora_matrix<-as.matrix(taxo_ciliophora_mat_unamed)

#transform to phyloseq objects using functions otu_table, tax_table and sample_data
OTU_ciliophora=otu_table(otu_ciliophora_matrix, taxa_are_rows = TRUE)
TAX_ciliophora=tax_table(taxo_ciliophora_matrix)
SAMPLE_ciliophora=sample_data(metadata_ciliophora, errorIfNULL = T)

Perifiton18S_ciliophora_data<-phyloseq(OTU_ciliophora, TAX_ciliophora, SAMPLE_ciliophora)

#visualise phyloseq-class
Perifiton18S_ciliophora_data

#visualise data
sample_names(Perifiton18S_ciliophora_data)
rank_names(Perifiton18S_ciliophora_data)
sample_variables(Perifiton18S_ciliophora_data)

#go through the next scripts
