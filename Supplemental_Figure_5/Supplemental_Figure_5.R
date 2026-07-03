library(tidyverse)
library(patchwork)

# Define the subdirectory where plots should be saved
plots_dir <- "plots" 

# 1) define files -------------------------------------------------------------------------
dpp_file          <- "data/DPP_sup_fig_2.txt"
precursorNr_file  <- "data/protein_table_sup_fig_2.tsv"
file_mapping_file <- "data/file_mapping_sup_fig_2.tsv"

# 2) FIXING THE DPP FILE (BROKEN HEADER) ----------------------------------------------------------------------
header_line <- readLines(dpp_file, n = 1)
filenames   <- str_extract_all(header_line, "[^ ]+?\\.raw")[[1]]

DPP_matrix <- read_delim(
  dpp_file,
  delim     = "\t",
  skip      = 1,
  col_names = FALSE
)

if (ncol(DPP_matrix) != length(filenames)) {
  warning("Column count mismatch — applying truncation to smallest length.")
  min_cols   <- min(ncol(DPP_matrix), length(filenames))
  DPP_matrix <- DPP_matrix[, 1:min_cols]
  filenames  <- filenames[1:min_cols]
}

colnames(DPP_matrix) <- filenames
DPP <- as.data.frame(DPP_matrix)

# 3) LOAD MAPPING + PRECURSOR FILES-------------------------------------------------------------------------
file_mapping <- read_tsv(file_mapping_file)
precursorNr  <- read_tsv(precursorNr_file)

# 4) CLEAN MAPPING & DEFINE ORDER-------------------------------------------------------------------------
mapping <- file_mapping %>%
  mutate(FileName = as.character(FileName))

all_conditions_present <- unique(mapping$condition)

desired_order <- c(
  "0cell",
  sort(all_conditions_present[str_starts(all_conditions_present, "1x")]),
  sort(all_conditions_present[str_starts(all_conditions_present, "10x")]),
  "1pg", "5pg", "10pg", "20pg", "50pg", "100pg", "200pg"
)
desired_order_present  <- intersect(desired_order, all_conditions_present)
final_condition_levels <- c(desired_order_present, sort(setdiff(all_conditions_present, desired_order_present)))

# Global annotation reference
table_row_names <- c("#cells", "intact", "protoplast", "stained")
annotation_cols <- c("condition", "#cells", "intact", "protoplast", "stained")

annotation_base <- mapping %>%
  select(any_of(annotation_cols)) %>%
  distinct(condition, .keep_all = TRUE)

# 5) SHARED THEME & DYNAMIC ANNOTATION FUNCTIONS
shared_theme <- function() {
  list(
    theme_classic(base_size = 12, base_family = "Arial"),
    theme(
      axis.title.y  = element_text(size = 13, face = "bold", family = "Arial"),
      axis.title.x  = element_blank(),
      axis.text     = element_text(size = 12, color = "black", family = "Arial"),
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      legend.title  = element_text(size = 13, face = "bold", family = "Arial"),
      legend.text   = element_text(size = 12, family = "Arial"),
      axis.line     = element_line(linewidth = 0.8),
      axis.ticks    = element_line(linewidth = 0.8),
      plot.margin   = margin(t = 5, r = 5, b = 0, l = 5)
    )
  )
}

# Function to build the annotation table plot based on conditions actually present in the data
make_dynamic_annotation_plot <- function(current_data_conditions) {
  
  # Identify which conditions from the global order are actually in this specific plot
  present_conds <- intersect(final_condition_levels, unique(as.character(current_data_conditions)))
  
  ann_df <- annotation_base %>%
    filter(condition %in% present_conds) %>%
    mutate(condition = factor(condition, levels = present_conds)) %>%
    arrange(condition)
  
  table_long <- ann_df %>%
    mutate(across(-condition, as.character)) %>%
    pivot_longer(cols = -condition, names_to = "row_label", values_to = "cell_value") %>%
    mutate(
      cell_value = ifelse(is.na(cell_value) | cell_value == "NA", "", cell_value),
      row_label  = factor(row_label, levels = rev(table_row_names)),
      x_pos      = as.integer(condition)
    )
  
  n_conds <- length(present_conds)
  
  ggplot(table_long, aes(x = x_pos, y = row_label)) +
    geom_hline(yintercept = seq(0.5, length(table_row_names) + 0.5, by = 1),
               colour = "grey70", linewidth = 0.4) +
    geom_text(aes(label = cell_value), size = 5, family = "Arial", hjust = 0.5, vjust = 0.5) +
    scale_x_continuous(limits = c(0.5, n_conds + 0.5),
                       breaks = seq_len(n_conds), expand = c(0, 0)) +
    theme_minimal(base_family = "Arial") +
    theme(
      axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
      axis.title   = element_blank(),
      axis.text.y  = element_text(size = 15, family = "Arial", hjust = 1),
      panel.grid   = element_blank(),
      plot.margin  = margin(t = 0, r = 5, b = 5, l = 5)
    )
}


# 6) Supplemental Figure 2 B: DPP per condition-------------------------------------------------------------------------
dpp_long <- DPP %>%
  pivot_longer(cols = everything(), names_to = "FileName", values_to = "DPP_Value") %>%
  left_join(mapping, by = "FileName") %>%
  filter(!is.na(condition)) %>%
  mutate(condition = factor(condition, levels = final_condition_levels))

dpp_medians <- dpp_long %>%
  group_by(condition) %>%
  summarise(median_val = median(DPP_Value, na.rm = TRUE), .groups = "drop")

p_dpp <- ggplot(dpp_long, aes(x = condition, y = DPP_Value)) +
  geom_boxplot(fill = "#56B4E9", colour = "black", alpha = 0.85, outlier.size = 1) +
  geom_text(data = dpp_medians, aes(x = condition, y = median_val, label = round(median_val, 1)),
            size = 3.5, family = "Arial", fontface = "bold", vjust = -0.5) +
  labs(y = "Datapoints Per Peak") +
  shared_theme()

combined1 <- p_dpp / make_dynamic_annotation_plot(dpp_long$condition) + plot_layout(heights = c(4, 1))
combined1

# save the DPP plot
ggsave(file.path(plots_dir, "sup_fig2_B.svg"),
       plot = combined1, device = "svg", width = 9.48, height = 4.68, units = "in")


# 7) Supplemental Figure 2 A: Precursor distribution per condition -------------------------------------------------------------------------

precursor_long <- precursorNr %>%
  select(contains("PG.NrOfPrecursorsUsedForQuantification")) %>%
  pivot_longer(everything(), names_to = "RawColumn", values_to = "PrecursorNr") %>%
  mutate(FileName = str_extract(RawColumn, "[^ ]+?\\.raw")) %>%
  left_join(mapping, by = "FileName") %>%
  filter(!is.na(PrecursorNr), !is.na(condition)) %>%
  mutate(condition = factor(condition, levels = final_condition_levels))

precursor_medians <- precursor_long %>%
  group_by(condition) %>%
  summarise(median_val = median(PrecursorNr, na.rm = TRUE), .groups = "drop")

p_precursors <- ggplot(precursor_long, aes(x = condition, y = PrecursorNr)) +
  geom_boxplot(fill = "#56B4E9", colour = "black", alpha = 0.85, outlier.size = 1) +
  geom_text(data = precursor_medians, aes(x = condition, y = median_val, label = round(median_val, 1)),
            size = 3.5, family = "Arial", fontface = "bold", vjust = -0.5) +
  labs(y = "Nr of Precursors Used per Protein") +
  shared_theme()

combined2 <- p_precursors / make_dynamic_annotation_plot(precursor_long$condition) + plot_layout(heights = c(4, 1))
combined2

# save the precursor per protein plot
ggsave(file.path(plots_dir, "sup_fig2_A.svg"),
       plot = combined2, device = "svg", width = 9.48, height = 4.68, units = "in")


# 8) CV PLOTS: Preparation & Logic -------------------------------------------------------------------------

# Build master long table for CV calculations
pg_long <- precursorNr %>%
  select(PG.ProteinGroups, contains("PG.Quantity")) %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "RawColumn", values_to = "Quantity") %>%
  mutate(FileName = str_extract(RawColumn, "[^ ]+?\\.raw")) %>%
  left_join(
    precursorNr %>%
      select(PG.ProteinGroups, contains("PG.NrOfPrecursorsUsedForQuantification")) %>%
      pivot_longer(cols = -PG.ProteinGroups, names_to = "RawColumn", values_to = "PrecursorNr") %>%
      mutate(FileName = str_extract(RawColumn, "[^ ]+?\\.raw")) %>%
      select(PG.ProteinGroups, FileName, PrecursorNr),
    by = c("PG.ProteinGroups", "FileName")
  ) %>%
  left_join(mapping, by = "FileName") %>%
  filter(!is.na(condition))

calc_cv <- function(x) {
  x_use <- x[!is.na(x)]
  if (length(x_use) < 3) return(NA_real_)
  m <- mean(x_use)
  if (m == 0) return(NA_real_)
  (sd(x_use) / m) * 100
}

# Helper to plot CV with dynamic annotation
plot_cv_with_ann <- function(data, ylabel) {
  # Ensure factor levels match only what is present
  present_levels <- intersect(final_condition_levels, unique(as.character(data$condition)))
  data <- data %>% mutate(condition = factor(condition, levels = present_levels))
  
  medians <- data %>%
    group_by(condition) %>%
    summarise(median_val = median(CV, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(data, aes(x = condition, y = CV)) +
    geom_boxplot(fill = "#56B4E9", colour = "black", alpha = 0.85, outlier.size = 1) +
    geom_text(data = medians, aes(x = condition, y = median_val, label = paste0(round(median_val, 1), "%")),
              size = 3.5, family = "Arial", fontface = "bold", vjust = -0.5) +
    labs(y = ylabel) +
    shared_theme()
  
  print(p / make_dynamic_annotation_plot(data$condition) + plot_layout(heights = c(4, 1)))
}


# 9) Supplemental Figure 2 C: All Proteins (n >= 3 replicates) -------------------------------------------------------------------------
cv_all <- pg_long %>%
  filter(!is.na(Quantity)) %>%
  group_by(PG.ProteinGroups, condition) %>%
  summarise(CV = calc_cv(Quantity), .groups = "drop") %>%
  filter(!is.na(CV))

combined3 <- plot_cv_with_ann(cv_all, "CV% (all proteins, ≥3 replicates)")
combined3

# save the CV% all proteins plot
ggsave(file.path(plots_dir, "sup_fig2_C.svg"),
       plot = combined3, device = "svg", width = 9.48, height = 4.68, units = "in")


# 10) Supplemental Figure 2 D: Non-Contaminants + >= 3 Precursors -------------------------------------------------------------------------

cv_nocont_3prec <- pg_long %>%
  filter(!str_starts(PG.ProteinGroups, "cont_"),
         !is.na(Quantity), !is.na(PrecursorNr), PrecursorNr >= 3) %>%
  group_by(PG.ProteinGroups, condition) %>%
  summarise(CV = calc_cv(Quantity), .groups = "drop") %>%
  filter(!is.na(CV))

combined4 <- plot_cv_with_ann(cv_nocont_3prec, "CV% (non-contaminant, ≥3 precursors & replicates)")
combined4

# save the CV% non-contaminant plot
ggsave(file.path(plots_dir, "sup_fig2_D.svg"),
       plot = combined4, device = "svg", width = 9.48, height = 4.68, units = "in")
