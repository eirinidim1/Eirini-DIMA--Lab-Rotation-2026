# ==========================================================================
# Compare positive contribution scores with absolute negative contribution scores
# ==========================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# Load original dataset
# Select: TF_gene_centric_with_zscores.csv
data <- read.csv(file.choose())

# Convert dataset to long format
tf_binding_timepoint <- bind_rows(
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_0",
      contribution_score = mean_contrib_day_0
    ),
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_1",
      contribution_score = mean_contrib_day_1
    ),
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_7",
      contribution_score = mean_contrib_day_7
    )
)

# Classify positive and negative contribution scores
contribution_magnitude_data <- tf_binding_timepoint %>%
  filter(
    !is.na(contribution_score),
    contribution_score != 0
  ) %>%
  mutate(
    contribution_direction = case_when(
      contribution_score > 0 ~ "positive_contribution",
      contribution_score < 0 ~ "negative_contribution"
    ),
    contribution_magnitude = abs(contribution_score)
  )

View(contribution_magnitude_data)

write.csv(
  contribution_magnitude_data,
  "~/Desktop/contribution_magnitude_data.csv",
  row.names = FALSE
)

# ==========================================================================
# Summary statistics
# ==========================================================================

contribution_magnitude_summary <- contribution_magnitude_data %>%
  group_by(contribution_direction) %>%
  summarise(
    n_observations = n(),
    mean_magnitude = mean(contribution_magnitude, na.rm = TRUE),
    median_magnitude = median(contribution_magnitude, na.rm = TRUE),
    q25_magnitude = quantile(contribution_magnitude, 0.25, na.rm = TRUE),
    q75_magnitude = quantile(contribution_magnitude, 0.75, na.rm = TRUE),
    q90_magnitude = quantile(contribution_magnitude, 0.90, na.rm = TRUE),
    q95_magnitude = quantile(contribution_magnitude, 0.95, na.rm = TRUE),
    max_magnitude = max(contribution_magnitude, na.rm = TRUE),
    .groups = "drop"
  )

View(contribution_magnitude_summary)

write.csv(
  contribution_magnitude_summary,
  "~/Desktop/contribution_magnitude_summary_positive_vs_abs_negative.csv",
  row.names = FALSE
)

# ==========================================================================
# Boxplot
# ==========================================================================

p_boxplot <- ggplot(
  contribution_magnitude_data,
  aes(
    x = contribution_direction,
    y = contribution_magnitude
  )
) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(
    ylim = c(
      0,
      quantile(contribution_magnitude_data$contribution_magnitude, 0.99, na.rm = TRUE)
    )
  ) +
  labs(
    title = "Magnitude comparison of positive and negative contribution scores",
    x = "Contribution direction",
    y = "Absolute contribution score"
  ) +
  theme_bw()

print(p_boxplot)

ggsave(
  "~/Desktop/contribution_magnitude_boxplot_positive_vs_abs_negative.png",
  plot = p_boxplot,
  width = 7,
  height = 5,
  dpi = 300
)

