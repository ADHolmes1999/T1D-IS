### Scoring System Test Code
## Antonio Holmes
# 06/07/2024

# Install & Load Packages -------------------------------------------------

# Cutpointr
install.packages("cutpointr")

# load library
library(cutpointr)

# Loading & Tidying Data --------------------------------------------------

# Load gene expression data 
gene_expected_count <- read_csv("Data/gene_expected_count.annot.csv") # all in one data set

# Load seq key (for reference)
#sequencing_key <- read_excel("Data/sequencing key.xlsx", col_names = FALSE)


# RENAME THIS SECTION -----------------------------------------------------

# Try whenever

# make yes or no t1d column

#cutpointr(gene_data, gene_count, t1d_yes_no_col, 
        #  method = maximize_metric, metric = sum_sens_spec) #change to minimum

# summary(data_from_code_above)

