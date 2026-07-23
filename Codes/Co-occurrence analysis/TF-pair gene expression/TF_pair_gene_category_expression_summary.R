# ==========================================================
# Positive contribution TF pair analysis
# Part 2: Gene category classification and expression summary (continue from part1)
# ==========================================================

# =============================================
# Function to classify genes for one TF pair
# =============================================

classify_genes_for_pair <- function(tf1, tf2, association_data) {
  
  # Genes associated with TF1
  genes_tf1 <- association_data %>%
    filter(TFName == tf1) %>%
    pull(gene_name) %>%
    unique()
  
  # Genes associated with TF2
  genes_tf2 <- association_data %>%
    filter(TFName == tf2) %>%
    pull(gene_name) %>%
    unique()
  
  # Genes only associated with TF1
  genes_tf1_only <- setdiff(genes_tf1, genes_tf2)
  
  # Genes only associated with TF2
  genes_tf2_only <- setdiff(genes_tf2, genes_tf1)
  
  # Genes associated with both TF1 and TF2
  genes_both <- intersect(genes_tf1, genes_tf2)
  
  # extra step so empty categories do not create errors
  make_category_df <- function(genes, category_name) {
    
    if (length(genes) == 0) {
      return(data.frame(
        TF1 = character(0),
        TF2 = character(0),
        gene_name = character(0),
        gene_category = character(0)
      ))
    }
    
    data.frame(
      TF1 = rep(tf1, length(genes)),
      TF2 = rep(tf2, length(genes)),
      gene_name = genes,
      gene_category = rep(category_name, length(genes))
    )
  }
  
  bind_rows(
    make_category_df(genes_tf1_only, "TF1_only"),
    make_category_df(genes_tf2_only, "TF2_only"),
    make_category_df(genes_both, "Both_TFs")
  )
}

# =============================================
# Classify genes for all TF pairs
# =============================================

TF_pair_gene_categories <- TF_pairs %>%
  mutate(
    pair_data = purrr::map2(
      TF1,
      TF2,
      ~ classify_genes_for_pair(
        tf1 = .x,
        tf2 = .y,
        association_data = top30_TF_gene_associations
      )
    )
  ) %>%
  select(pair_data) %>%
  tidyr::unnest(pair_data)

View(TF_pair_gene_categories)

write.csv(
  TF_pair_gene_categories,
  "~/Desktop/TF_pair_gene_categories_positive_contribution.csv",
  row.names = FALSE
)

# =============================================
#  Count genes per TF pair and category
# =============================================

TF_pair_category_counts <- TF_pair_gene_categories %>%
  group_by(TF1, TF2, gene_category) %>%
  summarise(
    n_genes = n_distinct(gene_name),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = gene_category,
    values_from = n_genes,
    values_fill = 0
  ) %>%
  mutate(
    total_genes_in_pair = TF1_only + TF2_only + Both_TFs
  ) %>%
  arrange(desc(Both_TFs))

View(TF_pair_category_counts)

write.csv(
  TF_pair_category_counts,
  "~/Desktop/TF_pair_category_counts.csv",
  row.names = FALSE
)

# =============================================
#  Create gene expression table
# =============================================

gene_expression <- data %>%
  group_by(gene_name) %>%
  summarise(
    z_expr_0d = first(z_expr_0d),
    z_expr_1d = first(z_expr_1d),
    z_expr_7d = first(z_expr_7d),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c(z_expr_0d, z_expr_1d, z_expr_7d),
    names_to = "timepoint",
    values_to = "z_expression"
  ) %>%
  mutate(
    timepoint = recode(
      timepoint,
      z_expr_0d = "day_0",
      z_expr_1d = "day_1",
      z_expr_7d = "day_7"
    )
  )

View(gene_expression)

# =============================================
# Add expression to TF pair gene categories
# =============================================

TF_pair_gene_expression <- TF_pair_gene_categories %>%
  left_join(
    gene_expression,
    by = "gene_name"
  )

View(TF_pair_gene_expression)

write.csv(
  TF_pair_gene_expression,
  "~/Desktop/TF_pair_gene_expression_positive_contribution.csv",
  row.names = FALSE
)

# =============================================
# Summarise expression per TF pair, category and timepoint
# =============================================

TF_pair_expression_summary <- TF_pair_gene_expression %>%
  group_by(TF1, TF2, timepoint, gene_category) %>%
  summarise(
    n_genes = n_distinct(gene_name),
    mean_z_expression = mean(z_expression, na.rm = TRUE),
    median_z_expression = median(z_expression, na.rm = TRUE),
    sd_z_expression = sd(z_expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(TF1, TF2, timepoint, gene_category)

View(TF_pair_expression_summary)

write.csv(
  TF_pair_expression_summary,
  "~/Desktop/TF_pair_expression_summary.csv",
  row.names = FALSE
)