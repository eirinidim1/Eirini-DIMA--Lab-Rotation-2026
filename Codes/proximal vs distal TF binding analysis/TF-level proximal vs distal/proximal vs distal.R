# ==========================================================================
# TF-level mean TSS distance and expression comparison
# Positive contribution scores only
# ==========================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Load original dataset
# Choose: TF_gene_centric_with_zscores.csv
data <- read.csv(file.choose())

# Check column names
colnames(data)

# =============================================
# Convert contribution scores to long format
# =============================================

contribution_long <- data %>%
  mutate(
    abs_dist_to_TSS = abs(dist_to_TSS)
  ) %>%
  select(
    TFName,
    gene_name,
    dist_to_TSS,
    abs_dist_to_TSS,
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
  )

View(contribution_long)

# =============================================
# Keep only positive contribution scores
# =============================================

contribution_positive <- contribution_long %>%
  filter(
    !is.na(contribution_score),
    contribution_score > 0
  )

View(contribution_positive)

# =============================================
# Mean TSS distance per TF-gene pair
# =============================================

TF_gene_mean_distance <- contribution_positive %>%
  group_by(TFName, gene_name) %>%
  summarise(
    n_positive_events = n(),
    mean_abs_TSS_distance_per_gene = mean(abs_dist_to_TSS, na.rm = TRUE),
    median_abs_TSS_distance_per_gene = median(abs_dist_to_TSS, na.rm = TRUE),
    .groups = "drop"
  )

View(TF_gene_mean_distance)

# =============================================
# Mean of gene-level mean distances per TF
# =============================================

TF_level_mean_distance <- TF_gene_mean_distance %>%
  group_by(TFName) %>%
  summarise(
    n_genes_targeted = n(),
    total_positive_events = sum(n_positive_events),
    mean_of_gene_mean_abs_TSS_distance = mean(
      mean_abs_TSS_distance_per_gene,
      na.rm = TRUE
    ),
    median_of_gene_mean_abs_TSS_distance = median(
      mean_abs_TSS_distance_per_gene,
      na.rm = TRUE
    ),
    sd_gene_mean_abs_TSS_distance = sd(
      mean_abs_TSS_distance_per_gene,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(mean_of_gene_mean_abs_TSS_distance)

View(TF_level_mean_distance)

write.csv(
  TF_level_mean_distance,
  "~/Desktop/TF_level_mean_TSS_distance.csv",
  row.names = FALSE
)

# =============================================
# Full distribution of TF-level mean TSS distances
# =============================================

p_full <- ggplot(
  TF_level_mean_distance,
  aes(x = mean_of_gene_mean_abs_TSS_distance)
) +
  geom_histogram(bins = 50) +
  scale_x_continuous(
    breaks = c(0, 100000, 200000, 400000, 600000),
    labels = comma
  ) +
  labs(
    title = "TF-level mean TSS distances",
    x = "Mean absolute distance from TSS per TF (bp)",
    y = "Number of TFs"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank()
  )

print(p_full)

ggsave(
  filename = "~/Desktop/TF_level_mean_TSS_distance_distribution_full.png",
  plot = p_full,
  width = 10,
  height = 6,
  dpi = 300
)


# =============================================
# Zoomed distribution of TF-level mean TSS distances
# =============================================

p_zoom <- ggplot(
  TF_level_mean_distance,
  aes(x = mean_of_gene_mean_abs_TSS_distance)
) +
  geom_histogram(bins = 50) +
  geom_vline(
    xintercept = 10000,
    color = "red",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 10000,
    y = Inf,
    label = "10 kb cutoff",
    color = "red",
    angle = 90,
    vjust = -0.5,
    hjust = 1.1,
    size = 4
  ) +
  scale_x_continuous(
    breaks = c(0, 50000, 100000, 150000, 200000),
    labels = comma
  ) +
  coord_cartesian(xlim = c(0, 200000)) +
  labs(
    title = "TF-level mean TSS distances",
    x = "Mean absolute distance from TSS per TF (bp)",
    y = "Number of TFs"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank()
  )

print(p_zoom)

ggsave(
  filename = "~/Desktop/TF_level_mean_TSS_distance_distribution_zoomed.png",
  plot = p_zoom,
  width = 10,
  height = 6,
  dpi = 300
)

# =============================================
# Classify TFs as proximal or distal
# Cutoff = 10,000 bp
# =============================================

TF_level_class <- TF_level_mean_distance %>%
  mutate(
    TF_distance_class = ifelse(
      mean_of_gene_mean_abs_TSS_distance <= 10000,
      "proximal_TF",
      "distal_TF"
    )
  ) %>%
  select(
    TFName,
    mean_of_gene_mean_abs_TSS_distance,
    TF_distance_class
  )

View(TF_level_class)

# Check how many TFs are proximal/distal
table(TF_level_class$TF_distance_class)

write.csv(
  TF_level_class,
  "~/Desktop/TF_level_class.csv",
  row.names = FALSE
)

# =============================================
# Collapse expression to TF-gene pair level
# =============================================

TF_gene_expression <- contribution_positive %>%
  group_by(TFName, gene_name) %>%
  summarise(
    z_expr_0d = first(z_expr_0d),
    z_expr_1d = first(z_expr_1d),
    z_expr_7d = first(z_expr_7d),
    .groups = "drop"
  )

View(TF_gene_expression)

# =============================================
#  Add TF proximal/distal class
# =============================================

TF_gene_expression_class <- TF_gene_expression %>%
  left_join(
    TF_level_class,
    by = "TFName"
  )

View(TF_gene_expression_class)

# Check if any TFs did not get a class
sum(is.na(TF_gene_expression_class$TF_distance_class))

# =============================================
#  Convert expression to long format
# =============================================

TF_gene_expression_long <- TF_gene_expression_class %>%
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

View(TF_gene_expression_long)

# =============================================
# Summarise expression by TF class and timepoint
# =============================================

TF_class_expression_summary <- TF_gene_expression_long %>%
  group_by(timepoint, TF_distance_class) %>%
  summarise(
    n_TF_gene_pairs = n(),
    mean_z_expression = mean(z_expression, na.rm = TRUE),
    median_z_expression = median(z_expression, na.rm = TRUE),
    sd_z_expression = sd(z_expression, na.rm = TRUE),
    .groups = "drop"
  )

View(TF_class_expression_summary)

write.csv(
  TF_class_expression_summary,
  "~/Desktop/TF_class_expression_summary.csv",
  row.names = FALSE
)

# =============================================
# Test expression difference between proximal_TF and distal_TF
# =============================================

TF_class_expression_tests <- TF_gene_expression_long %>%
  group_by(timepoint) %>%
  summarise(
    n_proximal_TF_pairs = sum(TF_distance_class == "proximal_TF"),
    n_distal_TF_pairs = sum(TF_distance_class == "distal_TF"),
    
    median_proximal_TF_expression = median(
      z_expression[TF_distance_class == "proximal_TF"],
      na.rm = TRUE
    ),
    
    median_distal_TF_expression = median(
      z_expression[TF_distance_class == "distal_TF"],
      na.rm = TRUE
    ),
    
    difference_median_proximal_minus_distal =
      median_proximal_TF_expression - median_distal_TF_expression,
    
    wilcox_p_value = ifelse(
      n_proximal_TF_pairs >= 3 & n_distal_TF_pairs >= 3,
      wilcox.test(
        z_expression ~ TF_distance_class,
        exact = FALSE
      )$p.value,
      NA
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    wilcox_p_adj = p.adjust(wilcox_p_value, method = "BH")
  )

View(TF_class_expression_tests)

write.csv(
  TF_class_expression_tests,
  "~/Desktop/TF_class_expression_tests.csv",
  row.names = FALSE
)



















