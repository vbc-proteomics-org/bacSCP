# =============================================================================
# Supplemental Figure 1, data of Figure 2
# Protein Quantity vs C score, Q value, # Precursors, Data Completeness
# =============================================================================

# install.packages(c("tidyverse", "patchwork"))  # uncomment if needed
library(tidyverse)
library(patchwork)

# --- Set working directory (adjust for local RStudio) ---
setwd("C:\\Users\\manuel.matzinger\\OneDrive - VBC\\Share_Julia_Manuel\\GitHub\\Supplemental_Figure_1")

# --- Read data ---
protein_data <- read.delim("proteins_MM_Astral1_TC003_colispike_MBRpergroup_Report.tsv",
                           sep = "\t", header = TRUE, check.names = FALSE)
file_mapping <- read.delim("file_mapping_fig2.txt",
                            sep = "\t", header = TRUE, check.names = FALSE)

# --- Identify column types by suffix pattern ---
all_cols <- colnames(protein_data)

# Extract raw file name from column headers like:
# [1] 20250908_...0cell_P1.raw.PG.Quantity
extract_raw_file <- function(col_name) {
  m <- regmatches(col_name, regexpr("\\d+_.*\\.raw", col_name))
  if (length(m) == 0) return(NA_character_)
  return(m)
}

# Classify per-run columns
quan_cols   <- all_cols[grepl("\\.raw\\.PG\\.Quantity$", all_cols)]
cscore_cols <- all_cols[grepl("\\.raw\\.PG\\.Cscore", all_cols)]
nprec_cols  <- all_cols[grepl("\\.raw\\.PG\\.NrOfPrecursorsUsedForQuantification$", all_cols)]
qval_cols   <- all_cols[grepl("\\.raw\\.PG\\.QValue", all_cols)]

cat("Found", length(quan_cols), "Quantity columns,",
    length(cscore_cols), "Cscore columns,",
    length(nprec_cols), "NrPrecursors columns,",
    length(qval_cols), "QValue columns\n")

# --- Build mapping: raw file -> condition ---
file_to_condition <- setNames(file_mapping$condition, file_mapping$FileName)

# Map each per-run column to its condition
map_col_to_condition <- function(col_name) {
  raw_file <- extract_raw_file(col_name)
  if (is.na(raw_file)) return(NA_character_)
  cond <- file_to_condition[raw_file]
  return(unname(cond))
}

quan_conditions   <- sapply(quan_cols, map_col_to_condition)
cscore_conditions <- sapply(cscore_cols, map_col_to_condition)
nprec_conditions  <- sapply(nprec_cols, map_col_to_condition)
qval_conditions   <- sapply(qval_cols, map_col_to_condition)

# --- Filter for 3 conditions of interest ---
conditions_of_interest <- c("0cell", "1x_intact_st", "200pg")

# --- Contaminant flag ---
protein_data$is_contaminant <- grepl("^cont_", protein_data$PG.ProteinGroups)

# --- Helper: compute per-protein average across runs of a condition ---
condition_avg <- function(data, cols, cond_map, condition) {
  sel_cols <- cols[cond_map == condition]
  if (length(sel_cols) == 0) return(rep(NA_real_, nrow(data)))
  sub <- data[, sel_cols, drop = FALSE]
  sub <- as.data.frame(lapply(sub, function(x) as.numeric(gsub(",", ".", as.character(x)))))
  rowMeans(sub, na.rm = TRUE)
}

# --- Helper: compute per-protein average, excluding 0 and NaN ---
condition_avg_nonzero <- function(data, cols, cond_map, condition) {
  sel_cols <- cols[cond_map == condition]
  if (length(sel_cols) == 0) return(rep(NA_real_, nrow(data)))
  sub <- data[, sel_cols, drop = FALSE]
  sub <- as.data.frame(lapply(sub, function(x) as.numeric(gsub(",", ".", as.character(x)))))
  sub[sub == 0] <- NA
  rowMeans(sub, na.rm = TRUE)
}

# --- Helper: compute per-protein MINIMUM across runs of a condition ---
condition_min <- function(data, cols, cond_map, condition) {
  sel_cols <- cols[cond_map == condition]
  if (length(sel_cols) == 0) return(rep(NA_real_, nrow(data)))
  sub <- data[, sel_cols, drop = FALSE]
  sub <- as.data.frame(lapply(sub, function(x) as.numeric(gsub(",", ".", as.character(x)))))
  apply(sub, 1, function(row) {
    vals <- row[!is.na(row)]
    if (length(vals) == 0) return(NaN)
    min(vals)
  })
}

# --- Helper: compute data completeness per condition ---
condition_completeness <- function(data, cols, cond_map, condition) {
  sel_cols <- cols[cond_map == condition]
  if (length(sel_cols) == 0) return(rep(NA_real_, nrow(data)))
  sub <- data[, sel_cols, drop = FALSE]
  sub <- as.data.frame(lapply(sub, function(x) as.numeric(gsub(",", ".", as.character(x)))))
  n_total <- ncol(sub)
  n_present <- rowSums(!is.na(sub))
  (n_present / n_total) * 100
}

# --- Build long-format data for plotting ---
plot_data_list <- list()

for (cond in conditions_of_interest) {
  avg_quan   <- condition_avg(protein_data, quan_cols, quan_conditions, cond)
  avg_cscore <- condition_avg(protein_data, cscore_cols, cscore_conditions, cond)
  avg_nprec  <- condition_avg_nonzero(protein_data, nprec_cols, nprec_conditions, cond)
  avg_qval   <- condition_min(protein_data, qval_cols, qval_conditions, cond)
  completeness <- condition_completeness(protein_data, quan_cols, quan_conditions, cond)

  df_cond <- data.frame(
    ProteinGroup   = protein_data$PG.ProteinGroups,
    Condition      = cond,
    is_contaminant = protein_data$is_contaminant,
    Avg_Quantity   = avg_quan,
    Avg_Cscore     = avg_cscore,
    Avg_NrPrecursors = avg_nprec,
    Min_Qvalue     = avg_qval,
    Completeness   = completeness,
    stringsAsFactors = FALSE
  )
  plot_data_list[[cond]] <- df_cond
}

plot_df <- bind_rows(plot_data_list)

# Remove rows where Avg_Quantity is NaN (protein not detected at all in condition)
plot_df <- plot_df %>% filter(!is.nan(Avg_Quantity) & !is.na(Avg_Quantity))

# --- Color mapping ---
plot_df$ProteinType <- ifelse(plot_df$is_contaminant, "Contaminant", "Protein")
color_map <- c("Contaminant" = "grey60", "Protein" = "#2ca02c")

# --- Condition labels for facets ---
plot_df$Condition <- factor(plot_df$Condition,
                            levels = c("0cell", "1x_intact_st", "200pg"),
                            labels = c("0 cell", "1x intact st", "200 pg"))

# --- Shared y-axis limits across all 4 plots ---
y_vals <- log2(plot_df$Avg_Quantity)
y_range <- range(y_vals, na.rm = TRUE)
y_lim <- c(floor(y_range[1]), ceiling(y_range[2]))

# --- Common theme ---
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5)
  )

# --- Plot 1: Protein Quantity vs C score ---
p1 <- ggplot(plot_df, aes(x = Avg_Cscore, y = log2(Avg_Quantity), color = ProteinType)) +
  geom_point(alpha = 0.5, size = 1) +
  facet_wrap(~ Condition, scales = "free_x") +
  coord_cartesian(ylim = y_lim) +
  scale_color_manual(values = color_map) +
  labs(x = "Average C Score", y = expression(log[2]("Protein Quantity")),
       title = "Protein Quantity vs C Score") +
  theme_pub

# --- Plot 2: Protein Quantity vs Q value (lowest/best per condition) ---
p2 <- ggplot(plot_df, aes(x = Min_Qvalue, y = log2(Avg_Quantity), color = ProteinType)) +
  geom_point(alpha = 0.5, size = 1) +
  facet_wrap(~ Condition, scales = "free_x") +
  coord_cartesian(ylim = y_lim) +
  scale_color_manual(values = color_map) +
  labs(x = "Best (Min) Q Value (Run-Wise)", y = expression(log[2]("Protein Quantity")),
       title = "Protein Quantity vs Q Value") +
  theme_pub

# --- Plot 3: Protein Quantity vs # Precursors (split x-axis) ---
# Custom piecewise-linear transformation:
#   [0, 5]    maps to [0, 0.5]   (50% of axis width)
#   (5, max]  maps to (0.5, 1.0] (50% of axis width)
split_breakpoint <- 5
nprec_max <- ceiling(max(plot_df$Avg_NrPrecursors, na.rm = TRUE))

split_transform <- function(x) {
  ifelse(x <= split_breakpoint,
         x / split_breakpoint * 0.5,
         0.5 + (x - split_breakpoint) / (nprec_max - split_breakpoint) * 0.5)
}
split_inverse <- function(x) {
  ifelse(x <= 0.5,
         x / 0.5 * split_breakpoint,
         split_breakpoint + (x - 0.5) / 0.5 * (nprec_max - split_breakpoint))
}

split_trans <- scales::trans_new(
  name      = "split_axis",
  transform = split_transform,
  inverse   = split_inverse
)

# Tick marks: 1-5 in the left half, sensible ticks in the right half
right_ticks <- seq(10, nprec_max, by = 10)
right_ticks <- right_ticks[right_ticks > split_breakpoint]
split_breaks <- sort(unique(c(1:split_breakpoint, right_ticks)))

p3 <- ggplot(plot_df, aes(x = Avg_NrPrecursors, y = log2(Avg_Quantity), color = ProteinType)) +
  geom_point(alpha = 0.5, size = 1) +
  # Vertical line at the split boundary
  geom_vline(xintercept = split_breakpoint, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(~ Condition, scales = "free_x") +
  coord_cartesian(ylim = y_lim) +
  scale_x_continuous(trans = split_trans, breaks = split_breaks) +
  scale_color_manual(values = color_map) +
  labs(x = "Avg # Precursors Used for Quantification",
       y = expression(log[2]("Protein Quantity")),
       title = "Protein Quantity vs # Precursors") +
  theme_pub

# --- Plot 4: Protein Quantity vs Data Completeness ---
p4 <- ggplot(plot_df, aes(x = Completeness, y = log2(Avg_Quantity), color = ProteinType)) +
  geom_point(alpha = 0.5, size = 1) +
  facet_wrap(~ Condition, scales = "free_x") +
  coord_cartesian(ylim = y_lim) +
  scale_color_manual(values = color_map) +
  labs(x = "Data Completeness (%)",
       y = expression(log[2]("Protein Quantity")),
       title = "Protein Quantity vs Data Completeness") +
  theme_pub

# --- Combine and save ---
combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave("proteomics_4panel_figure.svg", combined, width = 16, height = 8)

cat("\nPlots saved successfully!\n")
cat("Proteins in plot data:", nrow(plot_df), "\n")
cat("Conditions:", levels(plot_df$Condition), "\n")
