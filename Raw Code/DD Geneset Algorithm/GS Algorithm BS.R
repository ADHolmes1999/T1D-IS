### Data Driven GS Algorithm (Bootstrapping)
## Antonio Holmes
# 06/27/2024


# DECODING THE FUNX ------------------------------------------------------------------

## What does the function require to produce the results I am looking for 

# Nat Log, Log FC. P-Val (adj P-val)
# Pre-defined gene-set (w/ gene dictionary and annot)
# modeling gene expression
# pseudolilikehood funx

# Install & Load Packages ------------------------------------------------------------------

install.packages(c("igraph", "ggplot2"))
BiocManager::install("clusterProfiler")
BiocManager::install("ReactomePA",force = TRUE)

library(igraph)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)


# Gene Set & Pathways ------------------------------------------------------------------

#Example gene list
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

# Graph & Results ------------------------------------------------------------------

# Create the graph
g <- graph_from_data_frame(edge_list, directed = FALSE)

plot(g, vertex.label.cex = 0.8, vertex.size = 5)