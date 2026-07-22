# ==========================================================
# Part 1: Identify top TFs by positive-contribution gene associations
# and generate all pairwise TF combinations
# ==========================================================


library(dplyr)
library(tidyr)
library(purrr)

# Load original dataset
# Choose: TF_gene_centric_with_zscores.csv
data <- read.csv(file.choose())

# Convert contribution and expression to long format

tf_binding_timepoint <- bind_rows(
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_0",
      contrib = mean_contrib_day_0,
      z_expression = z_expr_0d
    ),
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_1",
      contrib = mean_contrib_day_1,
      z_expression = z_expr_1d
    ),
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_7",
      contrib = mean_contrib_day_7,
      z_expression = z_expr_7d
    )
)

View(tf_binding_timepoint)


# Keep only positive contribution cases
positive_contrib_data <- tf_binding_timepoint %>%
  filter(contrib > 0)

View(positive_contrib_data)

# =============================================
# Create unique TF-gene associations
# =============================================
# A TF is associated with a gene if it has at least one positive contribution case
positive_TF_gene_associations <- positive_contrib_data %>%
  distinct(TFName, gene_name)

View(positive_TF_gene_associations)

# =============================================
# Count unique associated genes per TF
# =============================================

TF_gene_counts <- positive_TF_gene_associations %>%
  group_by(TFName) %>%
  summarise(
    n_associated_genes = n_distinct(gene_name),
    .groups = "drop"
  ) %>%
  arrange(desc(n_associated_genes))

View(TF_gene_counts)

write.csv(
  TF_gene_counts,
  "~/Desktop/TF_gene_counts_positive_contribution.csv",
  row.names = FALSE
)

# =============================================
# Select top 30 TFs with the most associated genes
# =============================================

top30_TFs <- TF_gene_counts %>%
  slice_max(
    order_by = n_associated_genes,
    n = 30,
    with_ties = FALSE
  )

View(top30_TFs)

top30_TF_names <- top30_TFs$TFName

write.csv(
  top30_TFs,
  "~/Desktop/top30_TFs_by_positive_associated_genes.csv",
  row.names = FALSE
)

# =====================================================
#  Keep TF-gene associations only for top 30 TFs
# =====================================================

top30_TF_gene_associations <- positive_TF_gene_associations %>%
  filter(TFName %in% top30_TF_names)

View(top30_TF_gene_associations)

# =============================================
# Create all unique TF pairs from top 30 TFs
# =============================================

TF_pairs <- combn(top30_TF_names, 2, simplify = FALSE) %>%
  map_dfr(~ data.frame(
    TF1 = .x[1],
    TF2 = .x[2]
  ))

View(TF_pairs)

write.csv(
  TF_pairs,
  "~/Desktop/top30_TF_pairs.csv",
  row.names = FALSE
)

nrow(TF_pairs)
