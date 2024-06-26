


# Why it does not matter if there are more than one row of genes ----------

#https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html


# GSVA Code that works (EX CODE) ----------------------------------------------------

BiocManager::install("limma")

library(limma)

p <- 10 ## number of genes
n <- 30 ## number of samples
nGrp1 <- 15 ## number of samples in group 1
nGrp2 <- n - nGrp1 ## number of samples in group 2

## consider three disjoint gene sets
geneSets <- list(set1=paste("g", 1:3, sep=""),
                 set2=paste("g", 4:6, sep=""),
                 set3=paste("g", 7:10, sep=""))

## sample data from a normal distribution with mean 0 and st.dev. 1
y <- matrix(rnorm(n*p), nrow=p, ncol=n,
            dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2

## build GSVA parameter object
gsvapar <- gsvaParam(y, geneSets, maxDiff=TRUE)

## estimate GSVA enrichment scores for the three sets
gsva_es <- gsva(gsvapar)

## fit the same linear model now to the GSVA enrichment scores
fit <- lmFit(gsva_es, design)

## estimate moderated t-statistics
fit <- eBayes(fit)

## set1 is differentially expressed
topTable(fit, coef="sampleGroup2vs1")

# Making Hallmarks List like GeneSets (FUNCTIONING)-------------------------------------

hal_test <- list(nod_sample_1$external_gene_name)

gsvapa_test <- gsvaParam(gene_expected_matrix, hal_test, maxDiff=FALSE)

gsva_test <- gsva(gsvapa_test)




# Log FC & P-Val (EX CODE) -----------------------------------------

## sample data from a normal distribution with mean 0 and st.dev. 1
y <- matrix(rnorm(n*p), nrow=p, ncol=n,
            dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2

## build design matrix
design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nGrp1), rep(1, nGrp2)))

## fit linear model
fit <- lmFit(y, design)

## estimate moderated t-statistics
fit <- eBayes(fit)

## genes in set1 are differentially expressed
topTable(fit, coef="sampleGroup2vs1")


# Making Log FC & P-Val (FUNCTIONING W/IN CORRECT SCRIPT) ---------------------------------

## I TRANSFERED THIS OVER TO DD DENESET ALGORITHM, I NEED TO MAKE A NEW MATRIX

# DATASET W/ ALL SAMPLES WITH HLA MUTATION
mutation_dataset<- NOD_dataset %>% full_join(NS_dataset)

n_samples <- 45 ## number of samples
nod_Grp1 <- 25 ## number of samples in group 1
ns_Grp2 <- n_samples - nod_Grp1 ## number of samples in group 2


## the matrix needed for lmFit()
gene_expected_matrix

## not sure if I need this yet
#gene_expected_matrix[hallmarks_list$set1, (nGrp1+1):n] <- gene_expected_matrix[hallmarks_list$set1, (nGrp1+1):n] + 2

## build design matrix
gene_design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nod_Grp1), rep(1, ns_Grp2)))

## fit linear model
#fit <- 
  lmFit(gene_expected_matrix, gene_design)














# Interaction Graph (FUNCTIONING) -------------------------------------

install.packages(c("igraph", "ggplot2"))
BiocManager::install(c("clusterProfiler", "ReactomePA"))
library(igraph)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)

# Example gene list
genes <- c("11171", "8243", "112464", "2194",
          "9318", "79026", "1654", "65003",
          "6240", "3476", "6238", "3836",
          "4176", "1017", "249")

# Pathway enrichment analysis
enriched_pathways <- enrichPathway(gene = genes, organism = "human")

# Extract pathway and gene information
pathways <- as.data.frame(enriched_pathways)

# CAN I REUSE THE CODE BELOW PIPED INTO THE PATHFINDR FUNX

`# Create an edge list for the graph
edge_list <- data.frame(
  from = rep(pathways$Description, each = length(genes)),
  to = rep(genes, times = nrow(pathways)))

# Create the graph
g <- graph_from_data_frame(edge_list, directed = FALSE)

plot(g, vertex.label.cex = 0.8, vertex.size = 5)



  