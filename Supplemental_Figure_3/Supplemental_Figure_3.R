library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

# Define the subdirectory where plots should be saved "plots"
plots_dir <- "plots"  # Subdirectory with plots

# read data
BacSCP_data <- read_delim(
  "data/SC_protein_table_sup_fig_3.tsv",
  delim = "\t"
)

bulk200ng_data <- read_delim(
  "data/bulk_protein_table_sup_fig_3.tsv",
  delim = "\t"
)

# process BacSCP_data  (keep ALL proteins incl. contaminants)
protein_to_gene <- BacSCP_data %>%
  select(PG.ProteinGroups, PG.Genes) %>%
  distinct()

# Keep the full single-cell table (with precursor columns) BEFORE we reduce
# BacSCP_data to quantity columns only. We need the precursor information for
# splitting the matched single-cell overlay by precursor count.
BacSCP_raw <- BacSCP_data

BacSCP_data <- BacSCP_data %>%
  select(PG.ProteinGroups, contains(".raw.PG.Quantity"))

protein_data_t <- as.data.frame(t(BacSCP_data))
colnames(protein_data_t) <- BacSCP_data$PG.ProteinGroups
colnames(protein_data_t) <- make.names(colnames(protein_data_t), unique = TRUE)
protein_data_t <- rownames_to_column(protein_data_t, var = "Filename")
protein_data_t <- protein_data_t %>% slice(-1)

protein_data_t <- protein_data_t %>%
  mutate(condition = case_when(
    grepl("THIDTC006_0cell_", Filename) ~ "0cell",
    grepl("BS_dMcsB_chil_10x_", Filename) ~ "10x_BS_chil",
    grepl("BS_dMcsB_HS_10x_", Filename)   ~ "10x_BS_HS",
    grepl("BS_dMcsB_chil_", Filename) & !grepl("10x_", Filename) ~ "1x_BS_chil",
    grepl("BS_dMcsB_HS_",   Filename) & !grepl("10x_", Filename) ~ "1x_BS_HS",
    TRUE ~ NA_character_
  ))

# Identify contaminant vs. non-contaminant protein columns
cont_cols <- grep("^cont_", colnames(protein_data_t), value = TRUE)
prot_cols <- setdiff(colnames(protein_data_t), c("Filename", "condition", cont_cols))
all_prot_cols <- c(prot_cols, cont_cols)  # all proteins including contaminants


# process bulk200ng_data  (keep ALL proteins incl. contaminants)

bulk200ng_data <- bulk200ng_data %>%
  select(PG.ProteinGroups, contains("PG.Quantity"))

# Long format WITH contaminants
bulk200ng_long_all <- bulk200ng_data %>%
  pivot_longer(-PG.ProteinGroups, names_to = "sample_col", values_to = "quantity") %>%
  mutate(condition = case_when(
    grepl("_dMscB_BS_hs", sample_col) ~ "200ng dMscB 50°C",
    grepl("_dMscB_BS_st", sample_col) ~ "200ng dMscB 37°C",
    TRUE ~ "other"
  ))

# Long format WITHOUT contaminants
bulk200ng_long_nocont <- bulk200ng_long_all %>%
  filter(!grepl("^cont_", PG.ProteinGroups))

# HELPER: BacSCP long format (two versions)

make_bacscp_long <- function(include_contaminants) {
  cols <- if (include_contaminants) all_prot_cols else prot_cols
  protein_data_t %>%
    select(Filename, condition, all_of(cols)) %>%
    pivot_longer(-c(Filename, condition), names_to = "protein_col", values_to = "quantity") %>%
    mutate(quantity = suppressWarnings(as.numeric(quantity))) %>%
    filter(!is.na(quantity) & !is.nan(quantity) & quantity > 0)
}

bacscp_long_all    <- make_bacscp_long(include_contaminants = TRUE)
bacscp_long_nocont <- make_bacscp_long(include_contaminants = FALSE)

# HELPER: single-cell precursor category per protein --------------------------
# For a given single-cell condition (identified by sample_pattern), compute the
# mean number of precursors used for quantification across the replicate samples
# in which the protein is quantified (precursors > 0). Each protein is then
# binned into "1", "2" or "3+" precursors. Keyed by make.names(PG.ProteinGroups)
# so it can be joined onto the matched bulk proteins.
get_sc_precursor_cat <- function(raw_df, sample_pattern) {
  prec_cols <- grep("PG.NrOfPrecursorsUsedForQuantification",
                    colnames(raw_df), value = TRUE)
  prec_cols <- prec_cols[grepl(sample_pattern, prec_cols) & !grepl("10x", prec_cols)]

  if (length(prec_cols) == 0)
    stop(paste("No precursor columns found for sample pattern:", sample_pattern))

  m <- as.matrix(raw_df[, prec_cols])
  # coerce to numeric (handle possible comma decimals defensively)
  m <- matrix(suppressWarnings(as.numeric(gsub(",", ".", m))),
              nrow = nrow(m), ncol = ncol(m))
  m[is.na(m) | m <= 0] <- NA

  mean_prec <- rowMeans(m, na.rm = TRUE)

  tibble(
    protein_id     = make.names(raw_df$PG.ProteinGroups),
    mean_precursors = mean_prec
  ) %>%
    filter(!is.na(mean_precursors) & !is.nan(mean_precursors)) %>%
    # if make.names collapses duplicates, keep the max mean (defensive)
    group_by(protein_id) %>%
    summarise(mean_precursors = mean(mean_precursors, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(prec_cat = case_when(
      mean_precursors < 1.5 ~ "1",
      mean_precursors < 2.5 ~ "2",
      TRUE                  ~ "3+"
    ))
}

# Precursor categories per single-cell condition
sc_prec_cat_37 <- get_sc_precursor_cat(BacSCP_raw, "BS_dMcsB_chil")  # 37 °C
sc_prec_cat_50 <- get_sc_precursor_cat(BacSCP_raw, "BS_dMcsB_HS")    # 50 °C

# HELPER FUNCTION: histogram plot (matched overlay split by precursor count) ---
make_quantity_hist <- function(bulk_long,
                    bacscp_long,
                    bulk_condition,
                    bacscp_condition,
                    sc_prec_cat,
                    highlight_colors,   # named vector c("1"=, "2"=, "3+"=)
                    bulk_label,
                    sc_label_prefix     # e.g. "Single cell dMcsB 37 °C"
                    ) {

  quant_bulk <- bulk_long %>%
    filter(condition == bulk_condition,
           !is.na(quantity) & !is.nan(quantity) & quantity > 0) %>%
    select(PG.ProteinGroups, quantity) %>%
    mutate(source = "bulk_all")

  proteins_bacscp <- bacscp_long %>%
    filter(condition == bacscp_condition) %>%
    pull(protein_col) %>%
    unique()

  matched_bulk <- quant_bulk %>%
    filter(make.names(PG.ProteinGroups) %in% proteins_bacscp) %>%
    mutate(source = "matched",
           protein_id = make.names(PG.ProteinGroups)) %>%
    # attach the single-cell precursor category
    left_join(sc_prec_cat %>% select(protein_id, prec_cat), by = "protein_id")

  if (nrow(matched_bulk) == 0)
    warning(paste("No matched proteins found for BacSCP condition:", bacscp_condition))

  n_missing_cat <- sum(is.na(matched_bulk$prec_cat))
  if (n_missing_cat > 0)
    warning(paste(n_missing_cat,
                  "matched protein(s) had no single-cell precursor info and were dropped for",
                  bacscp_condition))

  matched_bulk <- matched_bulk %>% filter(!is.na(prec_cat))

  # human readable category labels
  cat_labels <- c(
    "1"  = paste0(sc_label_prefix, ", 1 precursor"),
    "2"  = paste0(sc_label_prefix, ", 2 precursors"),
    "3+" = paste0(sc_label_prefix, ", \u22653 precursors")
  )

  # build a single shared fill scale: grey base + 3 highlight tones
  fill_levels <- c(bulk_label, cat_labels[["1"]], cat_labels[["2"]], cat_labels[["3+"]])
  fill_values <- setNames(
    c("#808080", highlight_colors[["1"]], highlight_colors[["2"]], highlight_colors[["3+"]]),
    fill_levels
  )

  base_df <- quant_bulk %>%
    transmute(quantity, fill_grp = bulk_label)

  overlay_df <- matched_bulk %>%
    transmute(quantity, fill_grp = cat_labels[prec_cat])

  plot_df <- bind_rows(base_df, overlay_df) %>%
    mutate(fill_grp = factor(fill_grp, levels = fill_levels))

  # common binwidth so base and overlay bars align exactly
  bin_width <- diff(range(log10(quant_bulk$quantity))) / 50

  ggplot(plot_df, aes(x = log10(quantity), fill = fill_grp)) +
    geom_histogram(
      data     = subset(plot_df, fill_grp == bulk_label),
      binwidth = bin_width, boundary = 0, closed = "left",
      alpha = 0.6, color = "black"
    ) +
    geom_histogram(
      data     = subset(plot_df, fill_grp != bulk_label),
      binwidth = bin_width, boundary = 0, closed = "left",
      alpha = 0.95, color = "black", position = "stack"
    ) +
    scale_fill_manual(values = fill_values, breaks = fill_levels, drop = FALSE) +
    labs(
      x     = expression(bold(log[10]~"(Protein quantity)")),
      y     = "Count",
      fill  = ""
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_classic(base_size = 14, base_family = "Arial") +
    theme(
      axis.title      = element_text(size = 14, face = "bold", family = "Arial"),
      axis.text       = element_text(size = 12, color = "black", family = "Arial"),
      legend.text     = element_text(size = 11, family = "Arial"),
      legend.position = "top"
    )
}

# Supplemental Figure 3 B: 50°C — WITHOUT contaminants ----
# Matched single-cell overlay split into 3 violet tones by precursor count.

p1 <- make_quantity_hist(
  bulk_long        = bulk200ng_long_nocont,
  bacscp_long      = bacscp_long_nocont,
  bulk_condition   = "200ng dMscB 50°C",
  bacscp_condition = "1x_BS_HS",
  sc_prec_cat      = sc_prec_cat_50,
  highlight_colors = c("1" = "#BCBDDC", "2" = "#9E9AC8", "3+" = "#756BB1"),
  bulk_label       = "200ng dMcsB 50 °C (contaminants removed)",
  sc_label_prefix  = "Matched single cell dMcsB 50 °C (n=20)"
)
p1

# save the histogram
ggsave(file.path(plots_dir, "sup_fig3_B_hist_50°C.svg"),
       plot = p1, device = "svg", width = 10.5, height = 5, units = "in")


# Supplemental Figure 3 A: 37 °C — WITHOUT contaminants ----
# Matched single-cell overlay split into 3 orange tones by precursor count.

p2 <- make_quantity_hist(
  bulk_long        = bulk200ng_long_nocont,
  bacscp_long      = bacscp_long_nocont,
  bulk_condition   = "200ng dMscB 37°C",
  bacscp_condition = "1x_BS_chil",
  sc_prec_cat      = sc_prec_cat_37,
  highlight_colors = c("1" = "#FED976", "2" = "#FD8D3C", "3+" = "#E31A1C"),
  bulk_label       = "200ng dMcsB 37 °C (contaminants removed)",
  sc_label_prefix  = "Matched single cell dMcsB 37 °C (n=20)"
)
p2

# save the histogram
ggsave(file.path(plots_dir, "sup_fig3_A_hist_37°C.svg"),
       plot = p2, device = "svg", width = 11, height = 5.5, units = "in")
