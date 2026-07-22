# ==========================================================
# HEATMAP OF SIGNIFICANT TF PAIR-TIMEPOINT CASES
# ==========================================================

library(dplyr)
library(pheatmap)

# ==========================================================
# Load summary table of significant TF-pair cases
# ==========================================================

# file name selected: TF_pair_final_significant_summary_table.csv
tf_pair_sig <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE
)

colnames(tf_pair_sig)
View(tf_pair_sig)


# ==========================================================
#  Create row labels for each TF pair-timepoint case
# ==========================================================

tf_pair_sig <- tf_pair_sig %>%
  mutate(
    case_id = paste(
      TF1,
      TF2,
      timepoint,
      sep = "_"
    )
  )

View(tf_pair_sig)


# ==========================================================
# Keep the three mean expression columns for heatmap
# ==========================================================

heatmap_data <- tf_pair_sig %>%
  select(
    case_id,
    mean_TF1_only,
    mean_TF2_only,
    mean_Both_TFs
  )

# Convert to data frame and set row names
heatmap_matrix <- as.data.frame(heatmap_data)

rownames(heatmap_matrix) <- heatmap_matrix$case_id

heatmap_matrix$case_id <- NULL

# Convert to numeric matrix
heatmap_matrix <- as.matrix(heatmap_matrix)

View(heatmap_matrix)


# ========================================
# Create row annotation for timepoint
# ========================================

annotation_row <- tf_pair_sig %>%
  select(
    case_id,
    timepoint
  ) %>%
  as.data.frame()

rownames(annotation_row) <- annotation_row$case_id

annotation_row$case_id <- NULL

# Make timepoint a factor
annotation_row$timepoint <- factor(
  annotation_row$timepoint,
  levels = c(
    "day_1",
    "day_7"
  )
)

View(annotation_row)



# Define annotation colours
annotation_colors <- list(
  timepoint = c(
    day_1 = "#29C6D1",
    day_7 = "#57A773"
  )
)


# ============================
# Create clustered heatmap
#=============================

pheatmap(
  heatmap_matrix,
  
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  
  color = colorRampPalette(
    c(
      "navy",
      "white",
      "firebrick3"
    )
  )(100),
  
  border_color = "grey80",
  
  fontsize_row = 7,
  fontsize_col = 11,
  
  labels_col = c(
    "TF1 only",
    "TF2 only",
    "Both TFs"
  ),
  
  angle_col = 45,
  
  main =
    "Expression profiles of significant TF pair-timepoint cases",
  
  cellwidth = 45,
  cellheight = 10,
  
  filename =
    "~/Desktop/Figure_TF_pair_heatmap.png",
  
  width = 9,
  height = 12
)