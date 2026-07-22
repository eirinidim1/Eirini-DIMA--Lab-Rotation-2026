# =========================================
# Upstream vs downstream analysis
# Using only positive contribution scores
# =========================================

library(dplyr)
library(tidyr)

# Choose the original dataset: TF_gene_centric_with_zscores.csv
data <- read.csv(file.choose())

# ================================================
# Define upstream/downstream and convert contribution to long format
# Keep only positive contribution scores
# ===============================================

up_down_long <- data %>%
  mutate(
    TSS_side = case_when(
      dist_to_TSS < 0 ~ "upstream",
      dist_to_TSS > 0 ~ "downstream",
      dist_to_TSS == 0 ~ "at_TSS"
    )
  ) %>%
  filter(TSS_side %in% c("upstream", "downstream")) %>%
  select(
    TFName,
    gene_name,
    dist_to_TSS,
    TSS_side,
    mean_contrib_day_0,
    mean_contrib_day_1,
    mean_contrib_day_7
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

View(up_down_long)

write.csv(
  up_down_long,
  "~/Desktop/up_down_long.csv",
  row.names = FALSE
)


# Summary per TF, timepoint and TSS side
up_down_summary <- up_down_long %>%
  group_by(TFName, timepoint, TSS_side) %>%
  summarise(
    n_binding_sites = n(),
    mean_contribution = mean(contribution_score, na.rm = TRUE),
    median_contribution = median(contribution_score, na.rm = TRUE),
    sd_contribution = sd(contribution_score, na.rm = TRUE),
    .groups = "drop"
  )

View(up_down_summary)

write.csv(
  up_down_summary,
  "~/Desktop/up_down_summary.csv",
  row.names = FALSE
)

# =========================================================
# Wilcoxon test upstream vs downstream per TF and timepoint
# =========================================================

up_down_tests <- up_down_long %>%
  group_by(TFName, timepoint) %>%
  summarise(
    n_upstream = sum(TSS_side == "upstream"),
    n_downstream = sum(TSS_side == "downstream"),
    
    median_upstream = median(
      contribution_score[TSS_side == "upstream"],
      na.rm = TRUE
    ),
    
    median_downstream = median(
      contribution_score[TSS_side == "downstream"],
      na.rm = TRUE
    ),
    
    difference_median_upstream_minus_downstream =
      median_upstream - median_downstream,
    
    wilcox_p_value = ifelse(
      n_upstream >= 3 & n_downstream >= 3,
      wilcox.test(
        contribution_score ~ TSS_side,
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

View(up_down_tests)

write.csv(
  up_down_tests,
  "~/Desktop/up_down_tests.csv",
  row.names = FALSE
)

# =============================================
# Significant upstream/downstream hits
# =============================================

up_down_hits <- up_down_tests %>%
  filter(
    wilcox_p_adj < 0.05,
    n_upstream >= 3,
    n_downstream >= 3
  ) %>%
  mutate(
    direction = case_when(
      difference_median_upstream_minus_downstream > 0 ~ "upstream_higher",
      difference_median_upstream_minus_downstream < 0 ~ "downstream_higher",
      TRUE ~ "no_difference"
    )
  ) %>%
  arrange(wilcox_p_adj)

View(up_down_hits)

write.csv(
  up_down_hits,
  "~/Desktop/up_down_hits.csv",
  row.names = FALSE
)

# =============================================
# Count direction of significant hits
# =============================================

up_down_direction_summary <- up_down_hits %>%
  count(direction)

View(up_down_direction_summary)

write.csv(
  up_down_direction_summary,
  "~/Desktop/up_down_direction_summary.csv",
  row.names = FALSE
)

# =================================
# Count direction per timepoint
# =================================

up_down_direction_by_timepoint <- up_down_hits %>%
  count(timepoint, direction)

View(up_down_direction_by_timepoint)

write.csv(
  up_down_direction_by_timepoint,
  "~/Desktop/up_down_direction_by_timepoint.csv",
  row.names = FALSE
)

# =============================================
# Figure: Significant upstream/downstream contribution score differences
# =============================================

library(ggplot2)

up_down_direction_by_timepoint_plot <- up_down_direction_by_timepoint %>%
  tidyr::complete(
    timepoint = c("day_0", "day_1", "day_7"),
    direction = c("upstream_higher", "downstream_higher"),
    fill = list(n = 0)
  ) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("day_0", "day_1", "day_7")),
    direction = recode(
      direction,
      upstream_higher = "Upstream higher",
      downstream_higher = "Downstream higher"
    )
  )

p_up_down_direction <- ggplot(
  up_down_direction_by_timepoint_plot,
  aes(
    x = timepoint,
    y = n,
    fill = direction
  )
) +
  geom_col(position = "dodge", width = 0.7) +
  labs(
    title = "Direction of significant upstream/downstream contribution score differences",
    x = "Timepoint",
    y = "Number of significant TF-timepoint cases",
    fill = "Higher contribution score"
  ) +
  theme_bw()

print(p_up_down_direction)

ggsave(
  "~/Desktop/Figure4_up_down_direction_by_timepoint.png",
  plot = p_up_down_direction,
  width = 7,
  height = 5,
  dpi = 300
)

# =============================================
# Table for report: Top upstream/downstream contribution score differences 
# =============================================

up_down_hits_top10_for_report <- up_down_hits %>%
  select(
    TFName,
    timepoint,
    n_upstream,
    n_downstream,
    difference_median_upstream_minus_downstream,
    wilcox_p_adj,
    direction
  ) %>%
  mutate(
    difference_median_upstream_minus_downstream =
      round(difference_median_upstream_minus_downstream, 4),
    wilcox_p_adj = signif(wilcox_p_adj, 3),
    direction = recode(
      direction,
      upstream_higher = "Upstream higher",
      downstream_higher = "Downstream higher"
    )
  ) %>%
  rename(
    TF = TFName,
    Timepoint = timepoint,
    `n upstream` = n_upstream,
    `n downstream` = n_downstream,
    `Median difference (upstream - downstream)` =
      difference_median_upstream_minus_downstream,
    `BH-adjusted p-value` = wilcox_p_adj,
    `Higher positive contribution score` = direction
  ) %>%
  arrange(`BH-adjusted p-value`) %>%
  slice_head(n = 10)

View(up_down_hits_top10_for_report)

write.csv(
  up_down_hits_top10_for_report,
  "~/Desktop/Table_up_down_top10_contribution_score_differences.csv",
  row.names = FALSE
)

