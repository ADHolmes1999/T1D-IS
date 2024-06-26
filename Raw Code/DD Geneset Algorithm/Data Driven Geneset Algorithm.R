### Data Driven GS Algorithm (PathfindR)
## Antonio Holmes
# 06/07/2024


## WHAT NEEDS TO HAPPEN? ---------------------------------------------------

# I NEED A DATASET WITH JUST GENE NAME, LOG FC(optional), AND P-VALUE
# then
# INSERT THIS TIABLE INTO run_pathfindR(example_pathfindR_input, p_val_threshold = 0.01)


#PATHFINR WEBSITE:  https://cran.r-project.org/web/packages/pathfindR/vignettes/intro_vignette.html

# Install & Load Packages -------------------------------------------------

install.packages("pathfindR")
BiocManager::install("org.Hs.eg.db")

library(pathfindR)

# Loading & Tidying Data --------------------------------------------------

# DATASET W/ ALL SAMPLES WITH HLA MUTATION
mutation_dataset<- NOD_dataset %>% full_join(NS_dataset)

# GSVA requires a matrix so make sure there are no duplicates
sum(duplicated(mutation_dataset$entrezgene_id)) # = 0

# Matrix for Log FC & P-Val ---------------------------------------------------------

# IF YOU RUN PREVIOUS CHUNK, change gene_expected to filtered_gene_expected
mutation_matrix <- mutation_dataset %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-entrezgene_id) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("external_gene_name") %>%
  # Now we can convert our object into a matrix
  as.matrix()


# Creating Mutation Calculation Table ------------------------------------------------------------------

n_samples <- 45 ## number of samples
nod_Grp1 <- 25 ## number of samples in group 1
ns_Grp2 <- n_samples - nod_Grp1 ## number of samples in group 2

## not sure if I need this yet
#gene_expected_matrix[hallmarks_list$set1, (nGrp1+1):n] <- gene_expected_matrix[hallmarks_list$set1, (nGrp1+1):n] + 2

## build design matrix
gene_design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nod_Grp1), rep(1, ns_Grp2)))

## fit linear model
fit_mut <- lmFit(mutation_matrix, gene_design)

## estimate moderated t-statistics
fit_mut <- eBayes(fit_mut)

## Making the results into a dataset, Mutation Calc
mut_calc<-topTable(fit_mut, coef="sampleGroup2vs1",number = 27) #10 vs 27 looks weird

# Selecting the columns I need
mut_calc<- mut_calc %>% select(logFC,P.Value)

# Plugging Results into PathfindR ---------------------------------------------------------------

pathfindR::run_pathfindR(mut_calc,p_val_threshold = 0.5) #threshold should be 0.05
        # ## Testing input
        # The input looks OK
        # ## Processing input. Converting gene symbols,
        # if necessary (and if human gene symbols provided)
        # Number of genes provided in input: 10
        # Number of genes in input after p-value filtering: 10
        ################
        # Error in input_processing(input, p_val_threshold, pin_path, convert2alias) : 
        #   None of the genes were in the PIN
        # Please check your gene symbols
              #THIS CODE NEEDS GENE SYMBOLS



