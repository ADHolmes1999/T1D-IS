### Data Tidying Code
## Antonio Holmes
# 06/14/2024

# Install & Load Packages -------------------------------------------------

install.packages(c("dplyr", "ggplot2","readr","readxl","tidyverse"))

# Load libraries (adding as I edit the code with the data)
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(tidyverse)

# Loading & Tidying Data --------------------------------------------------

# Second
## THIS WILL HAVE TO BE DONE ON NON-DIABs TOO FOR COMPARING (SAME CODE JUST DIFF DATASETS)

# Load gene expression data 
gene_expected_count <- read_csv("project/Data/gene_expected_count.annot.csv") # all in one dataset

#Filtering only the genes of interest
 gene_expected<-gene_expected_count %>% 
  filter(external_gene_name%in%c("Foxp3","Ptpn2","Tigit","Tgfb1", "Il1b", "Il17a", "Ifng","Serpine1", "Ptgs2","Nos1", "Nos2", "Nos3","Sod1","Sod2","Sod3","Cat", "Gpx1","Gpx2","Gclc", "Gclm","Gsr","Gstm1", "Gstt1","Gstp1", "Gzma","Prf1","Klrg1")) ## add the rest

#Free up some space and RAM
rm(gene_expected_count)

# NOD Samples -------------------------------------------------------------

# NOD sample 1
nod_sample_1 <-gene_expected %>% select(entrezgene_id, external_gene_name,`6872-JK-1`,`6872-JK-6`,`6872-JK-11`,`6872-JK-16`,`6872-JK-21`)

# NOD Sample 2
nod_sample_2 <-gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-2`,`6872-JK-7`,`6872-JK-12`,`6872-JK-17`,`6872-JK-22`)

# NOD Sample 3
nod_sample_3 <-gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-3`,`6872-JK-8`,`6872-JK-13`,`6872-JK-18`,`6872-JK-23`)

# NOD Sample 4
nod_sample_4 <-gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-4`,`6872-JK-9`,`6872-JK-14`,`6872-JK-19`,`6872-JK-24`)

# NOD Sample 5
nod_sample_5 <-gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-77`,`6872-JK-10`,`6872-JK-15`,`6872-JK-20`,`6872-JK-25`)

# Combining Data sets
NOD_dataset<- nod_sample_1 %>% 
  full_join(nod_sample_2) %>% full_join(nod_sample_3) %>%
  full_join(nod_sample_4) %>% full_join(nod_sample_5)

#removing the NOD sample datasets for now
rm(nod_sample_1,nod_sample_2,nod_sample_3,nod_sample_4,nod_sample_5)

# NOD (AHM) Samples [IGNORE] -------------------------------------------------------

# AHM NOD Sample 1
gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-71`,`6872-JK-73`,`6872-JK-75`)

# AHM NOD Sample 2
gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-72`,`6872-JK-74`,`6872-JK-76`)


# NOR Samples -------------------------------------------------------------

# NOR Sample 1
gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-26`,`6872-JK-31`,`6872-JK-36`,`6872-JK-41`,`6872-JK-46`)

# NOR Sample 2 == Seq 42 is missing from Gene_Expected
gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-27`,`6872-JK-32`,`6872-JK-37`,#`6872-JK-42`,
                               `6872-JK-47`)

# NOR Sample 3 == Seq 28 missing from Gene_Expected
gene_expected %>% select(entrezgene_id,external_gene_name,#`6872-JK-28`,
                               `6872-JK-33`,`6872-JK-38`,`6872-JK-43`,`6872-JK-48`)

# NOR Sample 4
gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-29`,`6872-JK-34`,`6872-JK-39`,`6872-JK-44`,`6872-JK-49`)

# NOR Sample 5
gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-30`,`6872-JK-35`,`6872-JK-40`,`6872-JK-45`,`6872-JK-50`)


# NS Samples --------------------------------------------------------------

# NS Sample 1
ns_sample_1 <-gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-51`,`6872-JK-55`,`6872-JK-59`,`6872-JK-63`,`6872-JK-67`)

# NS Sample 2 == Seq 52 missing or the real column is 52_REPREP
ns_sample_2 <-gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-52-REPREP`,`6872-JK-56`,`6872-JK-60`,`6872-JK-64`,`6872-JK-68`)

# NS Sample 3
ns_sample_3 <-gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-53`,`6872-JK-57`,`6872-JK-61`,`6872-JK-65`,`6872-JK-69`)

# NS Sample 4
ns_sample_4 <-gene_expected %>% select(entrezgene_id,external_gene_name,`6872-JK-54`,`6872-JK-58`,`6872-JK-62`,`6872-JK-66`,`6872-JK-70`)

# NS Sample 5 = none

NS_dataset<-ns_sample_1 %>% full_join(ns_sample_2) %>%
  full_join(ns_sample_3) %>%full_join(ns_sample_4)

#
rm(ns_sample_1,ns_sample_2,ns_sample_3,ns_sample_4)

# Free up some more space and RAM
rm(gene_expected)

# [WORK ON THIS, PURE TRASH] Results -----------------------------------------------------------------
 
# # Visualize inflammation scores (making a new dataset so I can see better)
# scores_df <- data.frame(Sample = names(inflammation_score), Score = inflammation_score)
# 
# ggplot(scores_df, aes(x = Sample, y = Score)) +
#   geom_bar(stat = "identity") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "Inflammation Scores", x = "Sample", y = "Inflammation Score")
