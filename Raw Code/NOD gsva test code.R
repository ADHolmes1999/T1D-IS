### NOD sample 1 GSVA Test Code
## Antonio Holmes
# 06/18/2024

# Install & Load Packages -------------------------------------------------


if (!("GSVA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("GSVA", update = FALSE)
}

# DPLYR Package
install.packages("dplyr")

#TIDYVERSE Package
install.packages("tidyverse")

# Load necessary libraries
library(GSVA)
library(Biobase)
library(GSEABase)
library(tidyverse,dplyr)
library(readr)
library(readxl)

#[DOWNLOADING AND LOADING GSVA, NOT WORKING]-----------------------------------------


install.packages("BiocManager")
BiocManager::install("GSVA", version="3.19", force=TRUE)

BiocManager::install("beachmat", version="3.19", force=TRUE)


install.packages("remotes")
library(remotes)
install_github("rcastelo/GSVA")


# Load & Process Data -----------------------------------------------------

## USING THE TIDY DATASET (GENE_EXPECTED) FROM DATA_TIDYING CODE 

# GSVA requires a matrix so make sure there are no duplicates
sum(duplicated(nod_sample_1$entrezgene_id)) # = 34235

#Mutate duplicate column for Entrezgene IDs
nod_sample_1<-nod_sample_1 %>% mutate("Duplicate"=duplicated(entrezgene_id))

#Filtering out the duplcates
nod_sample_1<-nod_sample_1 %>% filter(Duplicate==FALSE)

# Checking once again
sum(duplicated(nod_sample_1$entrezgene_id)) # should = 0

# Matrix for GSVA ---------------------------------------------------------

# IF YOU RUN PREVIOUS CHUNK, change gene_expected to filtered_gene_expected
gene_expected_matrix <- nod_sample_1 %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Duplicate) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("entrezgene_id") %>%
  # Now we can convert our object into a matrix
  as.matrix()

# Results ---------------------------------------------------------

# Make hallmarks list
hallmarks_list <- split(
  gene_expected_count$entrezgene_id, # The genes we want split into pathways
  gene_expected_count$external_gene_name # The pathways made as the higher levels of the list
)

# The gsva() function documentation says use kcdf = "Gaussian" if expression values are continuous such as log-CPMs, log-RPKMs or log-TPMs. Use kcdf = "Poisson" on integer counts. 
# What is our vst() transformed data? I am guessing log2-like scale, so Gaussian works for us.
#gsva_results <-
  gsva(
  gene_expected_matrix,
  hallmarks_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size, may need to change
  min.sz = 15,
  # Maximum gene set size, may need to change
  max.sz = 500,
  # Compute Gaussian-distributed scores, play around with this
  mx.diff = TRUE,
  # Don't print out the progress bar, play around with this
  verbose = FALSE
)

# Print 10 rows (checking the data out)
head(gsva_results[, 1:10])

# [WORK ON THIS, PURE TRASH] Plots (Visuals) ---------------------------------------------------------
# 
# ## DATA MAY NEED TO BE SPREAD() OR GATHER()
# # 
# annot_df <- gene_expected %>%
#   # We need the sample IDs and the main column that contains the metadata info
#   dplyr::select(
#     refinebio_accession_code, #
#     refinebio_title #
#   ) %>%
#   # Create our `time_point` variable based on `refinebio_title`
#   dplyr::mutate(
#     time_point = dplyr::case_when(
#       # Create our new variable based whether the refinebio_title column
#       # contains _AV_ or _CV_
#       stringr::str_detect(refinebio_title, "_AV_") ~ "acute illness",
#       stringr::str_detect(refinebio_title, "_CV_") ~ "recovering"
#     )
#   ) %>%
#   # We don't need the older version of the variable anymore
#   dplyr::select(-refinebio_title)
# 
# 
# #
# annot_df <- annot_df %>%
#   # pheatmap will want our sample names that match our data to
#   tibble::column_to_rownames("refinebio_accession_code")
# 
# pathway_heatmap <- pheatmap::pheatmap(gsva_results,
# annotation_col = annot_df, # Add metadata labels!
# show_colnames = FALSE, # Don't show sample labels
# fontsize_row = 6 # Shrink the pathway labels a tad
# )
# 
# # Print out heatmap here
# pathway_heatmap
