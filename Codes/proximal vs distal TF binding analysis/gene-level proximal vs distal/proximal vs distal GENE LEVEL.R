# ========================================================================
# Gene-level proximal vs distal analysis-Positive contribution scores only
# ========================================================================

library(dplyr)
library(tidyr)

# Load original dataset
# file name: TF_gene_centric_with_zscores.csv
data <- read.csv(file.choose())

# ==============================================================================
# Convert contribution scores to long format and keep only positive contribution events
# ==============================================================================

positive_contribution_long <- data %>%
  mutate(
    TSS_region = ifelse(abs(dist_to_TSS) <= 2000, "proximal", "distal")
  ) %>%
  select(
    TFName,
    gene_name,
    dist_to_TSS,
    TSS_region,
    mean_contrib_day_0,
    mean_contrib_day_1,
    mean_contrib_day_7,
    z_expr_0d,
    z_expr_1d,
    z_expr_7d
  ) %>%
  pivot_longer(
    cols = c(mean_contrib_day_0, mean_contrib_day_1, mean_contrib_day_7),
    names_to = "timepoint",
    values_to = "contribution_score"
  ) %>%
  mutate(
    timepoint = recode(
      timepoint,
      mean_contrib_day_0 = "day_0",
      mean_contrib_day_1 = "day_1",
      mean_contrib_day_7 = "day_7"
    )
  ) %>%
  filter(
    !is.na(contribution_score),
    contribution_score > 0
  )

# ========================================================
# Collapse positive contribution events into TF-gene pairs
# ========================================================
#Positive contribution events were collapsed at the TF–gene pair level
#to summarize whether each pair contained proximal, distal, or both types
#of positive contribution events

tf_gene_TSS_category_positive <- positive_contribution_long %>%
  group_by(TFName, gene_name) %>%
  summarise(
    n_positive_events = n(),
    n_proximal_events = sum(TSS_region == "proximal"),
    n_distal_events = sum(TSS_region == "distal"),
    
    z_expr_0d = first(z_expr_0d),
    z_expr_1d = first(z_expr_1d),
    z_expr_7d = first(z_expr_7d),
    
    .groups = "drop"
  ) %>%
  mutate(
    TF_gene_TSS_category = case_when(
      n_proximal_events > 0 & n_distal_events == 0 ~ "proximal_only",
      n_proximal_events == 0 & n_distal_events > 0 ~ "distal_only",
      n_proximal_events > 0 & n_distal_events > 0 ~ "both",
      TRUE ~ NA_character_
    )
  )

View(tf_gene_TSS_category_positive)

write.csv(
  tf_gene_TSS_category_positive,
  "~/Desktop/gene_level_TF_gene_TSS_category_positive_contribution.csv",
  row.names = FALSE
)

# =============================================
# Count TF-gene pair categories
# =============================================

TF_gene_TSS_category_counts_positive <- tf_gene_TSS_category_positive %>%
  count(TF_gene_TSS_category)

View(TF_gene_TSS_category_counts_positive)

write.csv(
  TF_gene_TSS_category_counts_positive,
  "~/Desktop/gene_level_TF_gene_TSS_category_counts_positive_contribution.csv",
  row.names = FALSE
)

# =============================================
# Convert expression to long format so that each TF-gene pair is represented separately at each timepoint
# =============================================

tf_gene_expression_long_positive <- tf_gene_TSS_category_positive %>%
  select(
    TFName,
    gene_name,
    TF_gene_TSS_category,
    n_positive_events,
    n_proximal_events,
    n_distal_events,
    z_expr_0d,
    z_expr_1d,
    z_expr_7d
  ) %>%
  pivot_longer(
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

# =============================================
# Wilcoxon test: proximal_only vs distal_only per TF/timepoint
# =============================================

TF_proximal_vs_distal_tests_positive <- tf_gene_expression_long_positive %>%
  filter(
    TF_gene_TSS_category %in% c("proximal_only", "distal_only")
  ) %>%
  group_by(TFName, timepoint) %>%
  summarise(
    n_proximal_only = sum(TF_gene_TSS_category == "proximal_only"),
    n_distal_only = sum(TF_gene_TSS_category == "distal_only"),
    
    median_proximal_only = median(
      z_expression[TF_gene_TSS_category == "proximal_only"],
      na.rm = TRUE
    ),
    
    median_distal_only = median(
      z_expression[TF_gene_TSS_category == "distal_only"],
      na.rm = TRUE
    ),
    
    difference_median_proximal_minus_distal =
      median_proximal_only - median_distal_only,
    
    wilcox_p_value = ifelse(
      n_proximal_only >= 3 & n_distal_only >= 3,
      wilcox.test(
        z_expression ~ TF_gene_TSS_category,
        exact = FALSE
      )$p.value,
      NA
    ),
    
    .groups = "drop"
  ) %>%
  group_by(timepoint) %>%
  mutate(
    wilcox_p_adj = p.adjust(wilcox_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(wilcox_p_adj, wilcox_p_value)

View(TF_proximal_vs_distal_tests_positive)

write.csv(
  TF_proximal_vs_distal_tests_positive,
  "~/Desktop/gene_level_TF_proximal_vs_distal_expression_tests_positive_contribution.csv",
  row.names = FALSE
)

# =============================================
# Clean top 10 table for report
# =============================================

TF_proximal_vs_distal_top10_for_report <- TF_proximal_vs_distal_tests_positive %>%
  filter(
    !is.na(wilcox_p_value),
    n_proximal_only >= 3,
    n_distal_only >= 3
  ) %>%
  select(
    TFName,
    timepoint,
    n_proximal_only,
    n_distal_only,
    median_proximal_only,
    median_distal_only,
    difference_median_proximal_minus_distal,
    wilcox_p_value,
    wilcox_p_adj
  ) %>%
  mutate(
    median_proximal_only = round(median_proximal_only, 3),
    median_distal_only = round(median_distal_only, 3),
    difference_median_proximal_minus_distal =
      round(difference_median_proximal_minus_distal, 3),
    wilcox_p_value = signif(wilcox_p_value, 3),
    wilcox_p_adj = signif(wilcox_p_adj, 3)
  ) %>%
  arrange(wilcox_p_adj, wilcox_p_value) %>%
  slice_head(n = 10)

View(TF_proximal_vs_distal_top10_for_report)

write.csv(
  TF_proximal_vs_distal_top10_for_report,
  "~/Desktop/gene_level_TF_proximal_vs_distal_top10_for_report.csv",
  row.names = FALSE
)