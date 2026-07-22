# ==========================================================
# TF-GENE TIMECOURSE CONTRIBUTION-EXPRESSION CORRELATION
# Input: corr_positive_data.csv
# ==========================================================

library(dplyr)
library(ggplot2)


#==========================================================
#Load positive-contribution event-level data
#==========================================================

# filename: corr_positive_data.csv
corr_positive_data <- read.csv(file.choose())

colnames(corr_positive_data)


# ==========================================================
# Summarise positive contribution events per TF-gene pair and timepoint
# ==========================================================
#
# More than one positive contribution event may exist for
# the same TF-gene pair at the same timepoint.
#
# Therefore, we calculate one mean positive contribution
# score per TF-gene-timepoint combination.

TF_gene_timepoint_summary <- corr_positive_data %>%
  mutate(
    timepoint = as.character(timepoint)
  ) %>%
  group_by(
    TFName,
    gene_name,
    timepoint
  ) %>%
  summarise(
    n_positive_events = n(),
    
    mean_positive_contribution = mean(
      contribution_score,
      na.rm = TRUE
    ),
    
    expression_score = first(expression_score),
    
    .groups = "drop"
  )

View(TF_gene_timepoint_summary)


write.csv(
  TF_gene_timepoint_summary,
  "~/Desktop/TF_gene_timepoint_summary.csv",
  row.names = FALSE
)


# ==========================================================
# Keep TF-gene pairs represented at all three timepoints
# ==========================================================

complete_TF_gene_timecourse <- TF_gene_timepoint_summary %>%
  group_by(
    TFName,
    gene_name
  ) %>%
  filter(
    n_distinct(timepoint) == 3,
    all(
      c("day_0", "day_1", "day_7") %in% timepoint
    )
  ) %>%
  ungroup()

View(complete_TF_gene_timecourse)


n_complete_TF_gene_pairs <- complete_TF_gene_timecourse %>%
  distinct(
    TFName,
    gene_name
  ) %>%
  nrow()

print(
  paste(
    "Complete TF-gene pairs:",
    n_complete_TF_gene_pairs
  )
)


write.csv(
  complete_TF_gene_timecourse,
  "~/Desktop/complete_TF_gene_timecourse.csv",
  row.names = FALSE
)


# ==========================================================
# Calculate contribution-expression correlation across Day 0,1 and 7 for each TF-gene pair
# ==========================================================

TF_gene_timecourse_correlations <- complete_TF_gene_timecourse %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = c(
        "day_0",
        "day_1",
        "day_7"
      )
    )
  ) %>%
  arrange(
    TFName,
    gene_name,
    timepoint
  ) %>%
  group_by(
    TFName,
    gene_name
  ) %>%
  summarise(
    n_timepoints = n_distinct(timepoint),
    
    contribution_day_0 =
      mean_positive_contribution[
        timepoint == "day_0"
      ],
    
    contribution_day_1 =
      mean_positive_contribution[
        timepoint == "day_1"
      ],
    
    contribution_day_7 =
      mean_positive_contribution[
        timepoint == "day_7"
      ],
    
    expression_day_0 =
      expression_score[
        timepoint == "day_0"
      ],
    
    expression_day_1 =
      expression_score[
        timepoint == "day_1"
      ],
    
    expression_day_7 =
      expression_score[
        timepoint == "day_7"
      ],
    
    contribution_sd = sd(
      mean_positive_contribution,
      na.rm = TRUE
    ),
    
    expression_sd = sd(
      expression_score,
      na.rm = TRUE
    ),
    
    pearson_r = if (
      sd(mean_positive_contribution, na.rm = TRUE) > 0 &&
      sd(expression_score, na.rm = TRUE) > 0
    ) {
      cor(
        mean_positive_contribution,
        expression_score,
        method = "pearson",
        use = "complete.obs"
      )
    } else {
      NA_real_
    },
    
    spearman_rho = if (
      sd(mean_positive_contribution, na.rm = TRUE) > 0 &&
      sd(expression_score, na.rm = TRUE) > 0
    ) {
      cor(
        mean_positive_contribution,
        expression_score,
        method = "spearman",
        use = "complete.obs"
      )
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    absolute_spearman_rho = abs(spearman_rho),
    
    correlation_direction = case_when(
      is.na(spearman_rho) ~ "not_testable",
      spearman_rho > 0 ~ "positive",
      spearman_rho < 0 ~ "negative",
      TRUE ~ "zero"
    )
  ) %>%
  arrange(
    desc(absolute_spearman_rho)
  )

View(TF_gene_timecourse_correlations)


write.csv(
  TF_gene_timecourse_correlations,
  "~/Desktop/TF_gene_timecourse_correlations.csv",
  row.names = FALSE
)


# ===================================
# Summarise correlation results
# ===================================

TF_gene_timecourse_correlation_summary <-
  TF_gene_timecourse_correlations %>%
  summarise(
    n_complete_TF_gene_pairs = n(),
    
    n_testable_TF_gene_pairs = sum(
      !is.na(spearman_rho)
    ),
    
    n_not_testable_TF_gene_pairs = sum(
      is.na(spearman_rho)
    ),
    
    n_positive_correlations = sum(
      spearman_rho > 0,
      na.rm = TRUE
    ),
    
    n_negative_correlations = sum(
      spearman_rho < 0,
      na.rm = TRUE
    ),
    
    n_zero_correlations = sum(
      spearman_rho == 0,
      na.rm = TRUE
    ),
    
    mean_spearman_rho = mean(
      spearman_rho,
      na.rm = TRUE
    ),
    
    median_spearman_rho = median(
      spearman_rho,
      na.rm = TRUE
    ),
    
    mean_pearson_r = mean(
      pearson_r,
      na.rm = TRUE
    ),
    
    median_pearson_r = median(
      pearson_r,
      na.rm = TRUE
    )
  )

View(TF_gene_timecourse_correlation_summary)


write.csv(
  TF_gene_timecourse_correlation_summary,
  "~/Desktop/TF_gene_timecourse_correlation_summary.csv",
  row.names = FALSE
)


# ==========================================================
#  Count positive and negative correlation patterns
# ==========================================================

TF_gene_correlation_direction_counts <-
  TF_gene_timecourse_correlations %>%
  filter(
    correlation_direction != "not_testable"
  ) %>%
  count(
    correlation_direction,
    name = "n_TF_gene_pairs"
  )

View(TF_gene_correlation_direction_counts)


write.csv(
  TF_gene_correlation_direction_counts,
  "~/Desktop/TF_gene_correlation_direction_counts.csv",
  row.names = FALSE
)
# ==========================================================
# Distribution of exact Spearman correlation values across all testable TF-gene pairs
# ==========================================================

Spearman_distribution_data <- TF_gene_timecourse_correlations %>%
  filter(
    !is.na(spearman_rho)
  ) %>%
  mutate(
    spearman_rho_group = factor(
      spearman_rho,
      levels = sort(unique(spearman_rho))
    )
  ) %>%
  count(
    spearman_rho_group,
    name = "n_TF_gene_pairs"
  )

View(Spearman_distribution_data)

write.csv(
  Spearman_distribution_data,
  "~/Desktop/TF_gene_spearman_distribution_counts.csv",
  row.names = FALSE
)


Spearman_distribution_plot <- ggplot(
  Spearman_distribution_data,
  aes(
    x = spearman_rho_group,
    y = n_TF_gene_pairs
  )
) +
  geom_col(
    width = 0.7,
    fill = "grey50"
  ) +
  geom_text(
    aes(label = n_TF_gene_pairs),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    title = "Distribution of TF-gene timecourse correlations",
    x = "Spearman correlation coefficient",
    y = "Number of TF-gene pairs"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

Spearman_distribution_plot


ggsave(
  "~/Desktop/Figure_TF_gene_spearman_distribution.png",
  plot = Spearman_distribution_plot,
  width = 7,
  height = 5,
  dpi = 300
)