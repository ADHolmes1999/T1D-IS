### GSVA Test Code
## Antonio Holmes
# 05/29/2024

# Install & Load Packages -------------------------------------------------

# GSVA Package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSVA")

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


# Load & Process Data -----------------------------------------------------

## USING THE TIDY DATASET (GENE_EXPECTED) FROM DATA_TIDYING CODE 

# GSVA requires a matrix so make sure there are no duplicates
sum(duplicated(gene_expected_count$entrezgene_id)) # = 34235

#Mutate duplicate column for Entrezgene IDs
gene_expected_count<-gene_expected_count %>% mutate("Duplicate1"=duplicated(entrezgene_id),
                                                    "Duplicate2"=duplicated(external_gene_name))

#Filtering out the duplcates
gene_expected_count<-gene_expected_count %>% filter(Duplicate1==FALSE,
                                                    Duplicate2==FALSE)

# Checking once again
sum(duplicated(gene_expected_count$entrezgene_id)) # should = 0

# Matrix for GSVA ---------------------------------------------------------

# IF YOU RUN PREVIOUS CHUNK, change gene_expected to filtered_gene_expected
gene_expected_matrix <- gene_expected_count %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Duplicate1,-Duplicate2,-entrezgene_id, -gene_id,-description) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("external_gene_name") %>%
  # Now we can convert our object into a matrix
  as.matrix()

# Results ---------------------------------------------------------

# Makes hallmarks list
hallmarks_list <- list(gene_expected_count$external_gene_name)

# Param used in gsva()
gsva_param<-GSVA::gsvaParam(gene_expected_matrix,hallmarks_list)

#GSVA Function Results
gsva_results<-GSVA::gsva(gsva_param)


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
