# ==========================================================
# BOXPLOT:
# Positive contribution scores according to expression direction
# ==========================================================

library(dplyr)
library(ggplot2)

# ==============================
# Select the original dataset
# ==============================

data <- read.csv(file.choose())
# TF_gene_centric_with_zscores.csv


# ===============================================================
# Convert contribution scores and expression z-scores  to long format
# ===============================================================

expression_boxplot_data <- bind_rows(
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_0",
      contribution_score = mean_contrib_day_0,
      z_expression = z_expr_0d
    ),
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_1",
      contribution_score = mean_contrib_day_1,
      z_expression = z_expr_1d
    ),
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_7",
      contribution_score = mean_contrib_day_7,
      z_expression = z_expr_7d
    )
)


# =================================================
# Keep only positive contribution events and classify expression direction
# =================================================

expression_boxplot_data <- expression_boxplot_data %>%
  
  filter(
    !is.na(contribution_score),
    !is.na(z_expression),
    contribution_score > 0
  ) %>%
  
  mutate(
    
    expression_direction = case_when(
      z_expression > 0 ~ "Positive expression",
      z_expression < 0 ~ "Negative expression",
      TRUE ~ "Zero expression"
    )
    
  ) %>%
  
  # events with z-expression exactly equal to zero are excluded
  filter(expression_direction != "Zero expression") %>%
  
  mutate(
    
    timepoint = factor(
      timepoint,
      levels = c("day_0", "day_1", "day_7"),
      labels = c("Day 0", "Day 1", "Day 7")
    ),
    
    expression_direction = factor(
      expression_direction,
      levels = c(
        "Negative expression",
        "Positive expression"
      )
    )
  )


# ================================
# Calculate n for each boxplot
# ================================

boxplot_labels <- expression_boxplot_data %>%
  
  group_by(timepoint, expression_direction) %>%
  
  summarise(
    
    n = n(),
    
    # Upper whisker of each boxplot
    upper_whisker = boxplot.stats(
      contribution_score
    )$stats[5],
    
    .groups = "drop"
  )


#====================================
# Define the visible y-axis range
#====================================

# Contribution scores may contain a small number of very large values.
#The 99th percentile is used only to zoom the figure.
# No observations are removed from the boxplot calculations.

y_limit <- quantile(
  expression_boxplot_data$contribution_score,
  probs = 0.99,
  na.rm = TRUE
)

# Place each n label above its corresponding boxplot
boxplot_labels <- boxplot_labels %>%
  mutate(
    label_y = pmin(
      upper_whisker + 0.035 * y_limit,
      1.02 * y_limit
    )
  )


# =====================
# Create boxplot
# =====================

p_expression_direction_boxplot <- ggplot(
  expression_boxplot_data,
  aes(
    x = timepoint,
    y = contribution_score,
    fill = expression_direction
  )
) +
  
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.shape = NA
  ) +
  
  geom_text(
    data = boxplot_labels,
    aes(
      x = timepoint,
      y = label_y,
      label = paste0("n = ", n),
      group = expression_direction
    ),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = 0
  ) +
  
  
  scale_fill_manual(
    values = c(
      "Negative expression" = "#E6D2C3",
      "Positive expression" = "#CFE3C2"
    )
  ) +
  
  coord_cartesian(
    ylim = c(0, 1.10 * y_limit)
  ) +
  
  labs(
    title = paste(
      "Positive contribution scores according to",
      "expression direction across timepoints"
    ),
    x = "Timepoint",
    y = "Positive contribution score",
    fill = "Expression direction"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "top",
    plot.title = element_text(
      hjust = 0.5
    )
  )


# Display figure
print(p_expression_direction_boxplot)


# ==========================
# Save figure and the data used in the figure
# ==========================

#figure
ggsave(
  "~/Desktop/expression_direction_boxplot.png",
  plot = p_expression_direction_boxplot,
  width = 7,
  height = 4.5,
  dpi = 300
)

#data
write.csv(
  expression_boxplot_data,
  "~/Desktop/expression_direction_boxplot_data.csv",
  row.names = FALSE
)

write.csv(
  boxplot_labels,
  "~/Desktop/expression_direction_boxplot_sample_sizes.csv",
  row.names = FALSE
)