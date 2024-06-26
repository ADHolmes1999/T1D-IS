# All samples GSVA
library(msigdbr)
msigdbr_collections() # Take a look at all pathway groups in msigdbr database
sets_hallmark <- msigdbr(species="Mus musculus", category="H")
pwl_hallmark <- split(sets_hallmark$gene_symbol,sets_hallmark$gs_name)
sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME")
pwl_reactome <- split(sets_reactome$gene_symbol, sets_reactome$gs_name)
kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG")
pwl_kegg <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)
biocarta_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:BIOCARTA")
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, biocarta_gene_sets$gs_name)
pwl_msigdbr <- c(pwl_hallmark, pwl_reactome, pwl_kegg, pwl_biocarta) # Compile them all
length(pwl_msigdbr)
filtered_mapped_matrix <- r_var # GSVA requires data in matrix w/ genes as row names
filtered_mapped_matrix <- as.matrix(sapply(filtered_mapped_matrix, as.numeric))
rownames(filtered_mapped_matrix) <- rownames(r_var)
library(GSVA)
gsva_msig <- gsva(filtered_mapped_matrix, c(pwl_msigdbr), method = "gsva", # GSVA conversion
                  kcdf = "Gaussian", # "Gaussian" for continuous counts; "Poisson" on integers
                  min.sz = 15, max.sz = 500, # Minimum and max gene set size
                  mx.diff = TRUE, # Compute Gaussian-distributed scores
                  verbose = TRUE) # Progress bar

# Applying the GSVA package and pathway list to my data matrix r_var0

filtered_mapped_matrix <- r_var0 # GSVA requires data in matrix w/ genes as row names
filtered_mapped_matrix <- as.matrix(sapply(filtered_mapped_matrix, as.numeric))
rownames(filtered_mapped_matrix) <- rownames(r_var0)
gsva_msig0 <- gsva(filtered_mapped_matrix, c(pwl_msigdbr), method = "gsva", 
                   kcdf = "Gaussian", # "Gaussian" for continuous counts; "Poisson" on integers
                   min.sz = 15, max.sz = 500, # Minimum & max gene set size
                   mx.diff = TRUE, # Compute Gaussian-distributed scores
                   verbose = TRUE) # Progress bar, Output is matrix of gene set scores by samples
# dev.off()

