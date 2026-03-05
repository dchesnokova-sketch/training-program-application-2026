# Load libraries -------------------
# You may use base R or tidyverse for this exercise

library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)

# Load data here ----------------------
# Load each file with a meaningful variable name.
gene_lvl <- read_csv("data/GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv")
metadata <- read_csv("data/GSE60450_filtered_metadata.csv")

# Inspect the data -------------------------

# What are the dimensions of each data set? (How many rows/columns in each?)
# Keep the code here for each file.

## Expression data
dim(gene_lvl)
# 23735 rows 14 columns

## Metadata
dim(metadata)
# 12 rows 4 columns

# Prepare/combine the data for plotting ------------------------
# How can you combine this data into one data.frame?
colnames(metadata)[1] <- "sample_id"
colnames(gene_lvl)[1] <- "gene_id"

metadata <- metadata %>%
  mutate(
    # Creating cell_type column from immunophenotype
    cell_type = case_when(
      str_detect(tolower(immunophenotype), "luminal") ~ "Luminal",
      str_detect(tolower(immunophenotype), "basal") ~ "Basal"))

# Identify the columns in gene_lvl that contain expression data (those starting with "GSM")
gsm_cols <- names(gene_lvl)[grepl("^GSM", names(gene_lvl))]

# Long-format expression: one row per (gene, sample)
expr_long <- gene_lvl %>%
  pivot_longer(
    cols = all_of(gsm_cols),
    names_to = "sample_id",
    values_to = "expression"
  ) %>%
  mutate(
    sample_id = as.character(sample_id),
    expression = suppressWarnings(as.numeric(expression))
  )

# Joining with metadata
data_long <- expr_long %>%
  left_join(metadata %>% mutate(sample_id = as.character(sample_id)),
            by = "sample_id")

# Plot the data --------------------------
## Plot the expression by cell type
## Can use boxplot() or geom_boxplot() in ggplot2
# To overcome the issue with expression skewing to 0 and to make the plot more interpretable, I log-transformed the expression values.

p <- ggplot(data_long, aes(x = cell_type, y = log2(expression + 1), color = cell_type, fill = cell_type)) +
  geom_boxplot(width = 0.4, alpha = 0.6, 
               outlier.shape = 21,
               outlier.size = 1.5) +
  theme_cowplot() +
  labs(
    title = "Gene Expression by Cell Type",
    x = "Cell type",
    y = "log2(Expression + 1)",
    color = "Cell type",
    fill = "Cell type"
  )

print(p)

## Save the plot
### Show code for saving the plot with ggsave() or a similar function

ggsave("results/gene_expression_by_cell_type.png", plot = p, width = 8, height = 6)
