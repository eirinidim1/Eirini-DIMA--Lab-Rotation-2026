# ==========================================================================
# Positive contribution score vs expression correlation analysis PER TIMEPOINT
# ==========================================================================

library(dplyr)
library(tidyr)

# Load original dataset
# file: TF_gene_centric_with_zscores.csv
data <- read.csv(file.choose())

View(data)

# =============================================
# Convert dataset to long format
# =============================================

corr_long_data <- bind_rows(
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_0",
      contribution_score = mean_contrib_day_0,
      expression_score = z_expr_0d
    ),
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_1",
      contribution_score = mean_contrib_day_1,
      expression_score = z_expr_1d
    ),
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_7",
      contribution_score = mean_contrib_day_7,
      expression_score = z_expr_7d
    )
)

View(corr_long_data)

write.csv(
  corr_long_data,
  "~/Desktop/corr_long_data.csv",
  row.names = FALSE
)

# =============================================
# Keep only positive contribution scores
# and classify expression direction
# =============================================

corr_positive_data <- corr_long_data %>%
  filter(
    contribution_score > 0,
    !is.na(contribution_score),
    !is.na(expression_score),
    expression_score != 0
  ) %>%
  mutate(
    expression_direction = case_when(
      expression_score > 0 ~ "positive_expression",
      expression_score < 0 ~ "negative_expression"
    ),
    timepoint = factor(timepoint, levels = c("day_0", "day_1", "day_7"))
  )

View(corr_positive_data)

write.csv(
  corr_positive_data,
  "~/Desktop/corr_positive_data.csv",
  row.names = FALSE
)

# =============================================
# Count observations by timepoint and expression direction
# =============================================

corr_counts_timepoint <- corr_positive_data %>%
  count(timepoint, expression_direction, name = "n_events") %>%
  arrange(timepoint, expression_direction)

View(corr_counts_timepoint)

write.csv(
  corr_counts_timepoint,
  "~/Desktop/corr_counts_timepoint.csv",
  row.names = FALSE
)

# ===========================================
# Global correlation by timepoint and expression direction
# ===========================================

corr_global_timepoint <- corr_positive_data %>%
  group_by(timepoint, expression_direction) %>%
  summarise(
    n_observations = n(),
    
    pearson_r = cor(
      contribution_score,
      expression_score,
      method = "pearson",
      use = "complete.obs"
    ),
    
    pearson_p_value = cor.test(
      contribution_score,
      expression_score,
      method = "pearson"
    )$p.value,
    
    spearman_rho = cor(
      contribution_score,
      expression_score,
      method = "spearman",
      use = "complete.obs"
    ),
    
    spearman_p_value = cor.test(
      contribution_score,
      expression_score,
      method = "spearman",
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  arrange(timepoint, expression_direction)

View(corr_global_timepoint)

write.csv(
  corr_global_timepoint,
  "~/Desktop/corr_global_timepoint.csv",
  row.names = FALSE
)

# ===========================================
# Per-TF correlation by timepoint and expression direction
# ===========================================

corr_TF_timepoint <- corr_positive_data %>%
  group_by(TFName, timepoint, expression_direction) %>%
  filter(
    n() >= 3,
    sd(contribution_score, na.rm = TRUE) > 0,
    sd(expression_score, na.rm = TRUE) > 0
  ) %>%
  summarise(
    n_observations = n(),
    
    pearson_r = cor(
      contribution_score,
      expression_score,
      method = "pearson",
      use = "complete.obs"
    ),
    
    pearson_p_value = cor.test(
      contribution_score,
      expression_score,
      method = "pearson"
    )$p.value,
    
    spearman_rho = cor(
      contribution_score,
      expression_score,
      method = "spearman",
      use = "complete.obs"
    ),
    
    spearman_p_value = cor.test(
      contribution_score,
      expression_score,
      method = "spearman",
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  arrange(timepoint, expression_direction, TFName)

View(corr_TF_timepoint)

write.csv(
  corr_TF_timepoint,
  "~/Desktop/corr_TF_timepoint.csv",
  row.names = FALSE
)

# ==============
# BH correction
# ==============

corr_TF_timepoint_corrected <- corr_TF_timepoint %>%
  group_by(timepoint, expression_direction) %>%
  mutate(
    pearson_p_adj = p.adjust(pearson_p_value, method = "BH"),
    spearman_p_adj = p.adjust(spearman_p_value, method = "BH"),
    abs_pearson_r = abs(pearson_r),
    abs_spearman_rho = abs(spearman_rho)
  ) %>%
  ungroup() %>%
  arrange(timepoint, expression_direction, spearman_p_adj)

View(corr_TF_timepoint_corrected)

write.csv(
  corr_TF_timepoint_corrected,
  "~/Desktop/corr_TF_timepoint_corrected.csv",
  row.names = FALSE
)

# ====================
# Significant TF-timepoint correlation hits
# =====================

# Keep TF-timepoint-expression_direction groups with:
# n_observations >= 50
# spearman_p_adj < 0.05
# abs_spearman_rho >= 0.2

corr_significant_hits <- corr_TF_timepoint_corrected %>%
  filter(
    n_observations >= 50,
    spearman_p_adj < 0.05,
    abs_spearman_rho >= 0.2
  ) %>%
  arrange(timepoint, expression_direction, spearman_p_adj)

View(corr_significant_hits)

write.csv(
  corr_significant_hits,
  "~/Desktop/corr_significant_hits.csv",
  row.names = FALSE
)

# ====================
# Top 10 significant hits for report
# Ranked by adjusted Spearman p-value
# =====================

corr_top10_report <- corr_significant_hits %>%
  arrange(spearman_p_adj) %>%
  select(
    TFName,
    timepoint,
    expression_direction,
    n_observations,
    pearson_r,
    spearman_rho,
    spearman_p_adj
  ) %>%
  slice_head(n = 10)

View(corr_top10_report)

write.csv(
  corr_top10_report,
  "~/Desktop/corr_top10_report.csv",
  row.names = FALSE
)

#===================
#Bar plot for report
#===================

library(ggplot2)
library(dplyr)

# Make expression direction labels prettier
corr_global_plot_data <- corr_global_timepoint %>%
  mutate(
    expression_direction = recode(
      expression_direction,
      "positive_expression" = "Positive expression",
      "negative_expression" = "Negative expression"
    ),
    timepoint = factor(timepoint, levels = c("day_0", "day_1", "day_7"),
                       labels = c("Day 0", "Day 1", "Day 7"))
  )

# Plot global Spearman correlations per timepoint
global_corr_plot <- ggplot(
  corr_global_plot_data,
  aes(
    x = timepoint,
    y = spearman_rho,
    fill = expression_direction
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = c(
      "Negative expression" = "#E6D2C3",
      "Positive expression" = "#CFE3C2"
    )
  ) +
  labs(
    title = "Global contribution-expression correlations across timepoints",
    x = "Timepoint",
    y = "Spearman correlation coefficient",
    fill = "Expression direction"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

global_corr_plot

ggsave(
  "~/Desktop/Figure_global_correlation_timepoint.png",
  global_corr_plot,
  width = 7,
  height = 5,
  dpi = 300
)




