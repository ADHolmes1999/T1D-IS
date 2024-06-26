### NOD sample 1 GSVA Test Code
## Antonio Holmes
# 06/18/2024

#Load Data & Req Packages  -------------------------------------------------

#Running Data Tidying Script
source("/Users/antonioholmes/project/Code/Raw Code/data_tidying.R")

#Installing GSVA for this script
install.packages("BiocManager")
BiocManager::install("GSVA", version="3.19", force=TRUE)

#Loading the required packages
library(GSVA)

# Load & Process Data -----------------------------------------------------

## USING THE TIDY DATASET (GENE_EXPECTED) FROM DATA_TIDYING CODE 

# GSVA requires a matrix so make sure there are no duplicates
sum(duplicated(NOD_dataset$entrezgene_id)) # = 34235

#Mutate duplicate column for Entrezgene IDs
NOD_dataset<-NOD_dataset %>% mutate("Duplicate"=duplicated(entrezgene_id))

#Filtering out the duplcates
NOD_dataset<-NOD_dataset %>% filter(Duplicate==FALSE)

# Checking once again
sum(duplicated(NOD_dataset$entrezgene_id)) # should = 0

# Matrix for GSVA ---------------------------------------------------------

# IF YOU RUN PREVIOUS CHUNK, change gene_expected to filtered_gene_expected
gene_expected_matrix <- NOD_dataset %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Duplicate,-entrezgene_id) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("external_gene_name") %>%
  # Now we can convert our object into a matrix
  as.matrix()

# Results ---------------------------------------------------------

# Make hallmarks list
hallmarks_list <- list(NOD_dataset$external_gene_name)

# Param used in gsva()
nod_gsva_paramt<-GSVA::gsvaParam(gene_expected_matrix,hallmarks_list)

#GSVA Function Results
nod_gsva_results<-GSVA::gsva(gsva_test)

# NOT WORKING (POSSIBLE DELETION)
#gsva_results <-
  # GSVA::gsva(
  #   gsvaParam(gene_expected_matrix, hallmarks_list),
  # # Appropriate for our vst transformed data
  # kcdf = "Gaussian",
  # # Minimum gene set size, may need to change
  # min.sz = 15,
  # # Maximum gene set size, may need to change
  # max.sz = 500,
  # # Compute Gaussian-distributed scores, play around with this
  # mx.diff = TRUE,
  # # Don't print out the progress bar, play around with this
  # verbose = FALSE)
      

# Print 10 rows (checking the data out)
head(gsva_results[, 1:10])

# [WORK ON THIS, PURE TRASH] Plots (Visuals) ---------------------------------------------------------

# Shows gene expression as a heat map
pheatmap()

# # Multiple line graph
# x <- #timepoints
# y1 <- #NOD GSVA Results
# y2 <- #NOR GSVA Results
# y3 <- #NS GSVA Results
# 
# # Plot multiple lines using matplot
# matplot(x, cbind(y1, y2, y3), type = "l", lty = 1,
#         col = c("red", "blue", "green"), xlab = "X",
#         ylab = "Y", main = "Multiple Lines Plot")
# legend("topright", legend = c("NOD", "NOR", "NS"),
#        col = c("red", "blue", "green"),
#        lty = 1)