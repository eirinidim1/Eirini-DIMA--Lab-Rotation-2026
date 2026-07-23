# ==========================================================
# BINOMIAL TEST PER TF PER TIMEPOINT AGAINST TIMEPOINT BACKGROUND
# ==========================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# Load per TF per timepoint summary
positive_contrib_per_TF_timepoint <- read.csv(file.choose())
# file name: positive_contrib_per_TF_timepoint.csv


# ==========================================================
# Calculate global background per timepoint
# ==========================================================

global_background_per_timepoint <- positive_contrib_per_TF_timepoint %>%
  group_by(timepoint) %>%
  summarise(
    total_positive_expression = sum(n_positive_expression),
    total_negative_expression = sum(n_negative_expression),
    total_positive_contribution_cases = sum(n_total_positive_contrib),
    global_positive_fraction = total_positive_expression / total_positive_contribution_cases,
    global_negative_fraction = total_negative_expression / total_positive_contribution_cases,
    .groups = "drop"
  )

View(global_background_per_timepoint)

write.csv(
  global_background_per_timepoint,
  "~/Desktop/global_background_per_timepoint.csv",
  row.names = FALSE
)


# ==========================================================
# Binomial test per TF per timepoint
# ==========================================================

positive_contrib_binomial_results_per_TF_timepoint <- positive_contrib_per_TF_timepoint %>%
  left_join(
    global_background_per_timepoint %>%
      select(timepoint, global_positive_fraction),
    by = "timepoint"
  ) %>%
  rowwise() %>%
  mutate(
    binomial_p_value = binom.test(
      x = n_positive_expression,
      n = n_total_positive_contrib,
      p = global_positive_fraction,
      alternative = "two.sided"
    )$p.value
  ) %>%
  ungroup() %>%
  group_by(timepoint) %>%
  mutate(
    binomial_p_adj = p.adjust(binomial_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    difference_from_timepoint_global =
      positive_expression_fraction - global_positive_fraction,
    
    expression_bias = case_when(
      positive_expression_fraction > global_positive_fraction ~ "positive_expression_enriched",
      positive_expression_fraction < global_positive_fraction ~ "negative_expression_enriched",
      TRUE ~ "balanced"
    )
  ) %>%
  arrange(timepoint, binomial_p_adj, binomial_p_value)

View(positive_contrib_binomial_results_per_TF_timepoint)

write.csv(
  positive_contrib_binomial_results_per_TF_timepoint,
  "~/Desktop/positive_contrib_binomial_results_per_TF_timepoint.csv",
  row.names = FALSE
)


# ==========================================================
# Candidate TF-timepoint cases
# ==========================================================

positive_contrib_candidate_TFs_per_timepoint <- positive_contrib_binomial_results_per_TF_timepoint %>%
  filter(
    n_total_positive_contrib >= 50,
    binomial_p_value < 0.05
  ) %>%
  arrange(timepoint, binomial_p_value)

View(positive_contrib_candidate_TFs_per_timepoint)

write.csv(
  positive_contrib_candidate_TFs_per_timepoint,
  "~/Desktop/positive_contrib_candidate_TFs_per_timepoint.csv",
  row.names = FALSE
)

#Significant
positive_contrib_significant_TFs_per_timepoint <- positive_contrib_binomial_results_per_TF_timepoint %>%
  filter(
    n_total_positive_contrib >= 50,
    binomial_p_adj < 0.05
  ) %>%
  arrange(timepoint, binomial_p_adj)

View(positive_contrib_significant_TFs_per_timepoint)

write.csv(
  positive_contrib_significant_TFs_per_timepoint,
  "~/Desktop/positive_contrib_significant_TFs_per_timepoint.csv",
  row.names = FALSE
)


# ==========================================================
#  Figure: Expression direction per timepoint
# ==========================================================

expression_direction_per_timepoint <- global_background_per_timepoint %>%
  select(
    timepoint,
    total_positive_expression,
    total_negative_expression
  ) %>%
  pivot_longer(
    cols = c(total_positive_expression, total_negative_expression),
    names_to = "expression_direction",
    values_to = "n_events"
  ) %>%
  mutate(
    expression_direction = case_when(
      expression_direction == "total_positive_expression" ~ "Positive expression",
      expression_direction == "total_negative_expression" ~ "Negative expression"
    ),
    timepoint = factor(timepoint, levels = c("day_0", "day_1", "day_7")),
    expression_direction = factor(
      expression_direction,
      levels = c("Negative expression", "Positive expression")
    )
  )

View(expression_direction_per_timepoint)

write.csv(
  expression_direction_per_timepoint,
  "~/Desktop/Figure2_expression_direction_per_timepoint_data.csv",
  row.names = FALSE
)


p_expression_direction_per_timepoint <- ggplot(
  expression_direction_per_timepoint,
  aes(
    x = timepoint,
    y = n_events,
    fill = expression_direction
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  scale_fill_manual(
    values = c(
      "Positive expression" = "#CFE3C2",
      "Negative expression" = "#E6D2C3"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "day_0" = "Day 0",
      "day_1" = "Day 1",
      "day_7" = "Day 7"
    )
  ) +
  labs(
    title = "Expression direction among positive contribution events across timepoints",
    x = "Timepoint",
    y = "Number of positive contribution events",
    fill = "Expression direction"
  ) +
  theme_bw() +
  theme(
    legend.position = "top"
  )

print(p_expression_direction_per_timepoint)

ggsave(
  "~/Desktop/Figure2_expression_direction_per_timepoint.png",
  plot = p_expression_direction_per_timepoint,
  width = 7,
  height = 4.5,
  dpi = 300
)

