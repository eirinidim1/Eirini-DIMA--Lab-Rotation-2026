#=================================
# Pairwise Wilcoxon tests for each TF pair and each timepoint
#==================================
library(dplyr)
library(tidyr)
library(purrr)

# Load final TF pair gene expression dataset
# Choose: TF_pair_gene_expression_positive_contribution.csv
TF_pair_gene_expression <- read.csv(file.choose())

View(TF_pair_gene_expression)

# =============================================
# Count genes per TF pair, timepoint and category
# =============================================

TF_pair_group_counts <- TF_pair_gene_expression %>%
  group_by(TF1, TF2, timepoint, gene_category) %>%
  summarise(
    n_genes = n_distinct(gene_name),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = gene_category,
    values_from = n_genes,
    values_fill = 0
  )

View(TF_pair_group_counts)

# =============================================
# Keep valid TF pair-timepoint combinations with at least 30 genes in all three categories
# =============================================

valid_pair_timepoints <- TF_pair_group_counts %>%
  filter(
    TF1_only >= 30,
    TF2_only >= 30,
    Both_TFs >= 30
  ) %>%
  select(TF1, TF2, timepoint)

View(valid_pair_timepoints)

TF_pair_gene_expression_valid <- TF_pair_gene_expression %>%
  semi_join(
    valid_pair_timepoints,
    by = c("TF1", "TF2", "timepoint")
  )

View(TF_pair_gene_expression_valid)

# =============================================
# Pairwise Wilcoxon tests per TF pair and timepoint
# =============================================

TF_pair_pairwise_wilcoxon <- TF_pair_gene_expression_valid %>%
  group_by(TF1, TF2, timepoint) %>%
  summarise(
    n_TF1_only = n_distinct(gene_name[gene_category == "TF1_only"]),
    n_TF2_only = n_distinct(gene_name[gene_category == "TF2_only"]),
    n_Both_TFs = n_distinct(gene_name[gene_category == "Both_TFs"]),
    
    median_TF1_only = median(
      z_expression[gene_category == "TF1_only"],
      na.rm = TRUE
    ),
    
    median_TF2_only = median(
      z_expression[gene_category == "TF2_only"],
      na.rm = TRUE
    ),
    
    median_Both_TFs = median(
      z_expression[gene_category == "Both_TFs"],
      na.rm = TRUE
    ),
    
    mean_TF1_only = mean(
      z_expression[gene_category == "TF1_only"],
      na.rm = TRUE
    ),
    
    mean_TF2_only = mean(
      z_expression[gene_category == "TF2_only"],
      na.rm = TRUE
    ),
    
    mean_Both_TFs = mean(
      z_expression[gene_category == "Both_TFs"],
      na.rm = TRUE
    ),
    
    # Both_TFs vs TF1_only
    p_Both_vs_TF1 = wilcox.test(
      z_expression[gene_category == "Both_TFs"],
      z_expression[gene_category == "TF1_only"],
      exact = FALSE
    )$p.value,
    
    # Both_TFs vs TF2_only
    p_Both_vs_TF2 = wilcox.test(
      z_expression[gene_category == "Both_TFs"],
      z_expression[gene_category == "TF2_only"],
      exact = FALSE
    )$p.value,
    
    # TF1_only vs TF2_only
    p_TF1_vs_TF2 = wilcox.test(
      z_expression[gene_category == "TF1_only"],
      z_expression[gene_category == "TF2_only"],
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  )

View(TF_pair_pairwise_wilcoxon)

# =============================================
# BH correction per timepoint and comparison type
# =============================================

TF_pair_pairwise_wilcoxon_corrected <- TF_pair_pairwise_wilcoxon %>%
  group_by(timepoint) %>%
  mutate(
    p_adj_Both_vs_TF1 = p.adjust(p_Both_vs_TF1, method = "BH"),
    p_adj_Both_vs_TF2 = p.adjust(p_Both_vs_TF2, method = "BH"),
    p_adj_TF1_vs_TF2 = p.adjust(p_TF1_vs_TF2, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    diff_Both_minus_TF1 = median_Both_TFs - median_TF1_only,
    diff_Both_minus_TF2 = median_Both_TFs - median_TF2_only,
    diff_TF1_minus_TF2 = median_TF1_only - median_TF2_only
  )

View(TF_pair_pairwise_wilcoxon_corrected)

write.csv(
  TF_pair_pairwise_wilcoxon_corrected,
  "~/Desktop/TF_pair_pairwise_wilcoxon_corrected.csv",
  row.names = FALSE
)

# Find significant pairwise results

#Significant Both vs TF1
significant_Both_vs_TF1 <- TF_pair_pairwise_wilcoxon_corrected %>%
  filter(p_adj_Both_vs_TF1 < 0.05) %>%
  arrange(p_adj_Both_vs_TF1)

View(significant_Both_vs_TF1)

write.csv(
  significant_Both_vs_TF1,
  "~/Desktop/significant_Both_vs_TF1.csv",
  row.names = FALSE
)

#Significant Both vs TF2
significant_Both_vs_TF2 <- TF_pair_pairwise_wilcoxon_corrected %>%
  filter(p_adj_Both_vs_TF2 < 0.05) %>%
  arrange(p_adj_Both_vs_TF2)

View(significant_Both_vs_TF2)

write.csv(
  significant_Both_vs_TF2,
  "~/Desktop/significant_Both_vs_TF2.csv",
  row.names = FALSE
)

#Significant TF1 vs TF2
significant_TF1_vs_TF2 <- TF_pair_pairwise_wilcoxon_corrected %>%
  filter(p_adj_TF1_vs_TF2 < 0.05) %>%
  arrange(p_adj_TF1_vs_TF2)

View(significant_TF1_vs_TF2)

write.csv(
  significant_TF1_vs_TF2,
  "~/Desktop/significant_TF1_vs_TF2.csv",
  row.names = FALSE
)

# =============================================
# Summary of significant pairwise results
# =============================================

significant_summary <- data.frame(
  comparison = c("Both_vs_TF1", "Both_vs_TF2", "TF1_vs_TF2"),
  n_significant = c(
    nrow(significant_Both_vs_TF1),
    nrow(significant_Both_vs_TF2),
    nrow(significant_TF1_vs_TF2)
  )
)

View(significant_summary)

write.csv(
  significant_summary,
  "~/Desktop/significant_pairwise_summary.csv",
  row.names = FALSE
)


# =============================================
# Final significant summary table:
# mean expression, raw p-values and adjusted p-values
# Keep only rows with at least one significant adjusted p-value
# =============================================

TF_pair_final_significant_summary_table <- TF_pair_pairwise_wilcoxon_corrected %>%
  filter(
    p_adj_Both_vs_TF1 < 0.05 |
      p_adj_Both_vs_TF2 < 0.05 |
      p_adj_TF1_vs_TF2 < 0.05
  ) %>%
  select(
    TF1,
    TF2,
    timepoint,
    
    n_TF1_only,
    n_TF2_only,
    n_Both_TFs,
    
    mean_TF1_only,
    mean_TF2_only,
    mean_Both_TFs,
    
    p_Both_vs_TF1,
    p_Both_vs_TF2,
    p_TF1_vs_TF2,
    
    p_adj_Both_vs_TF1,
    p_adj_Both_vs_TF2,
    p_adj_TF1_vs_TF2,
    
    diff_Both_minus_TF1,
    diff_Both_minus_TF2,
    diff_TF1_minus_TF2
  ) %>%
  arrange(timepoint, p_adj_Both_vs_TF1, p_adj_Both_vs_TF2, p_adj_TF1_vs_TF2)

View(TF_pair_final_significant_summary_table)

write.csv(
  TF_pair_final_significant_summary_table,
  "~/Desktop/TF_pair_final_significant_summary_table.csv",
  row.names = FALSE
)

# =============================================
# Direction of significant pairwise effects
# =============================================

TF_pair_final_interpreted <- TF_pair_final_significant_summary_table %>%
  mutate(
    direction_Both_vs_TF1 = case_when(
      p_adj_Both_vs_TF1 < 0.05 & mean_Both_TFs > mean_TF1_only ~ "Both_higher_than_TF1_only",
      p_adj_Both_vs_TF1 < 0.05 & mean_Both_TFs < mean_TF1_only ~ "Both_lower_than_TF1_only",
      TRUE ~ "not_significant"
    ),
    
    direction_Both_vs_TF2 = case_when(
      p_adj_Both_vs_TF2 < 0.05 & mean_Both_TFs > mean_TF2_only ~ "Both_higher_than_TF2_only",
      p_adj_Both_vs_TF2 < 0.05 & mean_Both_TFs < mean_TF2_only ~ "Both_lower_than_TF2_only",
      TRUE ~ "not_significant"
    ),
    
    direction_TF1_vs_TF2 = case_when(
      p_adj_TF1_vs_TF2 < 0.05 & mean_TF1_only > mean_TF2_only ~ "TF1_only_higher_than_TF2_only",
      p_adj_TF1_vs_TF2 < 0.05 & mean_TF1_only < mean_TF2_only ~ "TF1_only_lower_than_TF2_only",
      TRUE ~ "not_significant"
    )
  )

View(TF_pair_final_interpreted)

write.csv(
  TF_pair_final_interpreted,
  "~/Desktop/TF_pair_final_interpreted.csv",
  row.names = FALSE
)

#summary
direction_summary_Both_vs_TF1 <- TF_pair_final_interpreted %>%
  count(direction_Both_vs_TF1)

direction_summary_Both_vs_TF2 <- TF_pair_final_interpreted %>%
  count(direction_Both_vs_TF2)

direction_summary_TF1_vs_TF2 <- TF_pair_final_interpreted %>%
  count(direction_TF1_vs_TF2)

View(direction_summary_Both_vs_TF1)
View(direction_summary_Both_vs_TF2)
View(direction_summary_TF1_vs_TF2)

# =============================================
# Rank strongest significant TF pair hits
# =============================================

top_candidate_hits <- TF_pair_final_interpreted %>%
  mutate(
    min_adjusted_p = pmin(
      p_adj_Both_vs_TF1,
      p_adj_Both_vs_TF2,
      p_adj_TF1_vs_TF2,
      na.rm = TRUE
    ),
    
    max_abs_mean_difference = pmax(
      abs(diff_Both_minus_TF1),
      abs(diff_Both_minus_TF2),
      abs(diff_TF1_minus_TF2),
      na.rm = TRUE
    )
  ) %>%
  arrange(min_adjusted_p, desc(max_abs_mean_difference))

View(top_candidate_hits)

write.csv(
  top_candidate_hits,
  "~/Desktop/top_candidate_TF_pair_hits.csv",
  row.names = FALSE
)

#=====================================================
#Bar plot for significant TF pair cases by timepoint
#=====================================================

p_TF_pair_timepoint <- ggplot(
  TF_pair_significant_by_timepoint,
  aes(x = timepoint, y = n_significant_cases)
) +
  geom_col(
    width = 0.65,
    fill = "grey50"
  ) +
  geom_text(
    aes(label = n_significant_cases),
    vjust = -0.4,
    size = 5
  ) +
  labs(
    title = "Significant TF pair cases by timepoint",
    x = "Timepoint",
    y = "Number of significant TF pair-timepoint cases"
  ) +
  ylim(0, 40) +
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_TF_pair_timepoint)

# =============================================
# Simple Table with expression pattern for report
# =============================================

Table_top10_TF_pairs_simple <- top_candidate_hits %>%
  slice_head(n = 10) %>%
  mutate(
    min_adjusted_p = signif(min_adjusted_p, 3),
    expression_pattern = case_when(
      mean_Both_TFs < mean_TF1_only & mean_Both_TFs > mean_TF2_only ~
        "Both_TFs intermediate",
      mean_Both_TFs < mean_TF1_only ~
        "Both_TFs lower than TF1_only",
      mean_Both_TFs > mean_TF2_only ~
        "Both_TFs higher than TF2_only",
      TRUE ~ "Distinct expression pattern"
    )
  ) %>%
  select(
    TF1,
    TF2,
    timepoint,
    expression_pattern,
    min_adjusted_p
  )

View(Table_top10_TF_pairs_simple)

write.csv(
  Table_top10_TF_pairs_simple,
  "~/Desktop/top10_TF_pairs_simpl.csv",
  row.names = FALSE
)