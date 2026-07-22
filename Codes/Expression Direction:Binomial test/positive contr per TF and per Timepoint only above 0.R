# =============================================
# Positive contribution/expression analysis
# Overall per TF + per TF per timepoint
# =============================================

library(dplyr)

# Load original dataset
data <- read.csv(file.choose())

# Convert to long format
tf_binding_timepoint <- bind_rows(
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_0",
      contrib = mean_contrib_day_0,
      z_expr = z_expr_0d
    ),
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_1",
      contrib = mean_contrib_day_1,
      z_expr = z_expr_1d
    ),
  
  data %>%
    transmute(
      gene_name,
      TFName,
      timepoint = "day_7",
      contrib = mean_contrib_day_7,
      z_expr = z_expr_7d
    )
)

# Keep only positive contributions and classify expression direction
positive_contrib_data <- tf_binding_timepoint %>%
  filter(contrib > 0) %>%
  mutate(
    expression_sign = case_when(
      z_expr > 0 ~ "positive_expression",
      z_expr < 0 ~ "negative_expression",
      TRUE ~ "zero_expression"
    )
  ) %>%
  filter(expression_sign != "zero_expression")

# ==================
# Overall per TF
# ==================

positive_contrib_overall_per_TF <- positive_contrib_data %>%
  group_by(TFName) %>%
  summarise(
    n_positive_expression = sum(expression_sign == "positive_expression"),
    n_negative_expression = sum(expression_sign == "negative_expression"),
    n_total_positive_contrib = n(),
    positive_expression_fraction = n_positive_expression / n_total_positive_contrib,
    negative_expression_fraction = n_negative_expression / n_total_positive_contrib,
    mean_contribution = mean(contrib, na.rm = TRUE),
    median_contribution = median(contrib, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(positive_expression_fraction))

View(positive_contrib_overall_per_TF)

write.csv(
  positive_contrib_overall_per_TF,
  "~/Desktop/positive_contrib_overall_per_TF.csv",
  row.names = FALSE
)

# =============================================
# Per TF per timepoint
# =============================================

positive_contrib_per_TF_timepoint <- positive_contrib_data %>%
  group_by(TFName, timepoint) %>%
  summarise(
    n_positive_expression = sum(expression_sign == "positive_expression"),
    n_negative_expression = sum(expression_sign == "negative_expression"),
    n_total_positive_contrib = n(),
    positive_expression_fraction = n_positive_expression / n_total_positive_contrib,
    negative_expression_fraction = n_negative_expression / n_total_positive_contrib,
    mean_contribution = mean(contrib, na.rm = TRUE),
    median_contribution = median(contrib, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(timepoint, desc(positive_expression_fraction))

View(positive_contrib_per_TF_timepoint)

write.csv(
  positive_contrib_per_TF_timepoint,
  "~/Desktop/positive_contrib_per_TF_timepoint.csv",
  row.names = FALSE
)