### Github repository for 

## Dataset periphyton (16S rRNA gene amplicon sequencing) pipeline analyses
Performed in R v.4.1.3 [R Core Team](https://www.r-project.org)

#### Color code
- Sites
  - ![#26a4b5](https://placehold.co/15x15/26a4b5/26a4b5.png) `#26a4b5` (Jangadeiros)
  - ![#014560](https://placehold.co/15x15/014560/014560.png) `#014560` (Veleiros)
- Time points
  - ![#23207C](https://placehold.co/15x15/23207C/23207C.png) `#23207C` (T1)
  - ![#00B5ED](https://placehold.co/15x15/00B5ED/00B5ED.png) `#00B5ED` (T2)
  - ![#F08C01](https://placehold.co/15x15/F08C01/F08C01.png) `#F08C01` (T3)
  - ![#D6000A](https://placehold.co/15x15/D6000A/D6000A.png) `#D6000A` (T4)
   


### Packages
library("dada2")\
library("ShortRead")\
library("Biostrings")\
library("phyloseq")\
library("microbiome")\
library("vegan")\
library("DESeq2") \
library("dplyr")\
library("stringr")\
library("ggpubr")\
library("tidyr")\
library("ggplot2")\
library("readxl")\
library("RColorBrewer")
  
### 1. Dada2
This script uses the raw data obtained from the BioProject SRA.\
A R server is required. \
Database used: [SILVA 138 SSU](https://www.arb-silva.de/documentation/release-138/) 

### 2. Creating phyloseq object
This script creates a phyloseq object based on these files:

- jki_seq1_metadata.csv
- jki_seq1_otu.xlsx
- jki_seq1_seqs.xlsx
- jki_seq1_taxa.xlsx

### 3. Rename NA
This script replaces NA for the latest taxonomy found for an ASV.

### 4. Clean dataset
This script removes unwanted taxonomic groups from the dataset.
- Root NAs (Domain, Phylum)
- Eukaryotes (Domain)
- Chloroplasts (Order)
- Mitochondria (Family)

### 5. Data selection
This script selects a group of samples to be analyzed separately.

### 6. Rarefaction
This script performs rarefaction based on the minimum sequences.

### 7. Alpha diversity
This script calculates alpha diversity based on the rarefied data.

### 8. Ordination 
This script creates MDS plots and calculates PERMANOVA and ANOSIM based on the rarefied data.

### 9. DESeq2
Based on the rarefied data, this script performs differential abundance (DA) between two groups.
