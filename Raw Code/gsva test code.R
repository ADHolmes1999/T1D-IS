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

 # First let's determine the gene means
 gene_means <- rowMeans(gene_expected %>% dplyr::select(-Ensembl, -entrez_id))

 # Let's add this as a column in our `mapped_df`.
 gene_expected <- gene_expected %>%
   # Add gene_means as a column called gene_means
   dplyr::mutate(gene_means) %>%
   # Reorder the columns so `gene_means` column is upfront
   dplyr::select(Ensembl, entrez_id, gene_means, dplyr::everything())

 # This is the new data set without duplication
 filtered_gene_expected <- gene_expected %>%
   # Sort so that the highest mean expression values are at the top
   dplyr::arrange(dplyr::desc(gene_means)) %>%
   # Filter out the duplicated rows using `dplyr::distinct()`
   dplyr::distinct(entrez_id, .keep_all = TRUE)

# Checking once again
sum(duplicated(filtered_mapped_df$entrez_id)) # should = 0

# Matrix for GSVA ---------------------------------------------------------

# IF YOU RUN PREVIOUS CHUNK, change gene_expected to filtered_gene_expected
gene_expected_matrix <- gene_expected %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Ensembl, -gene_means) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("entrez_id") %>%
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
gsva_results <- gsva(
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
