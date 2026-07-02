library(tidyverse)
library(extrafont)
library(patchwork)
library(gridExtra)
library(grid)
library(ggrepel)

loadfonts(device = "win", quiet = TRUE)

# ---- User-definable parameters ----
# Minimum number of precursors a protein must have (per file) to be kept.
# Values are taken from the "PG.NrOfPrecursorsUsedForQuantification" columns.
# (Used for Figure 2 B histogram.)
min_precursors <- 1
# Maximum Q value a protein may have (per file) to be kept.
# Values are taken from the "PG.QValue" columns.
# (Used for Figure 2 B histogram.)
max_qvalue <- 0.05

# ---- Proteins to annotate on the rank plots ----
# Add protein accession numbers here (as they appear in the PG.ProteinGroups column,
# e.g. "P0A6Y8" or "cont_P00761"). Each listed protein will be labelled on BOTH rank
# plots using its gene name from the PG.Genes column (via ggrepel). Leave empty (c())
# to draw the rank plots with no labels.
proteins_to_annotate <- c("P0A6T1", "P0A6Y8","P0A9B2","P0ABU2","P0AC53","P0ADY7","P21362","P33360","P64634", "P76165")

# define directories, read in data and prepare data ------------------------------------
# Define the subdirectory where plots should be saved.
# A relative "plots" folder is used so the script is portable; change this back to your
# own absolute path if you prefer (e.g. the original OneDrive/GitHub location).
plots_dir <- "plots"
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Read input files
coli_data <- read_delim("data/protein_table_fig2_revision.tsv", delim = "\t")
file_mapping <- read_delim("data/file_mapping_fig2.txt", delim = "\t")

# Columns containing PG.Quantity
protein_cols <- grep("PG.Quantity", colnames(coli_data), value = TRUE)

# Helper to extract base filename
extract_base_name <- function(colname) {
  colname_clean <- str_remove(colname, "^\\[\\d+\\]\\s*")
  pattern <- "([^.]+\\.raw)"
  match_result <- regmatches(colname_clean, regexec(pattern, colname_clean))
  if(length(match_result[[1]]) > 0) return(match_result[[1]][2])
  NA
}
base_names <- sapply(protein_cols, extract_base_name, USE.NAMES = FALSE)

# Link protein columns to condition from mapping
col_condition_df <- tibble(
  protein_col = protein_cols,
  base_name = base_names
) %>%
  left_join(file_mapping %>% select(FileName, condition),
            by = c("base_name" = "FileName"))

if(any(is.na(col_condition_df$condition))) {
  warning("Some filenames could not be matched to conditions. See printed rows:")
  print(col_condition_df[is.na(col_condition_df$condition), ])
}

# Mark contaminant vs E. coli
coli_data <- coli_data %>%
  mutate(protein_type = if_else(str_starts(PG.ProteinGroups, "cont_"), "Contaminant", "Ecoli"))

# Count quantified proteins per replicate column
plot_data <- map_dfr(protein_cols, function(col) {
  condition <- col_condition_df$condition[col_condition_df$protein_col == col]
  base_name <- col_condition_df$base_name[col_condition_df$protein_col == col]
  data_col <- coli_data[[col]]
  contaminant_count <- sum(!is.na(data_col[coli_data$protein_type == "Contaminant"]) &
                    !is.nan(data_col[coli_data$protein_type == "Contaminant"]))
  ecoli_count <- sum(!is.na(data_col[coli_data$protein_type == "Ecoli"]) &
                    !is.nan(data_col[coli_data$protein_type == "Ecoli"]))
  total_count <- contaminant_count + ecoli_count
  tibble(
    protein_col = col, base_name = base_name, condition = condition,
    contaminant_count = contaminant_count, ecoli_count = ecoli_count, total_count = total_count
  )
})

plot_data <- plot_data %>% filter(!is.na(condition))
if(nrow(plot_data) == 0) stop("No matched data left after linking filenames. Check mapping.")

# Filter conditions
exclude_conditions <- c("10x_intact_st", "10x_intact_us")
plot_data <- plot_data %>% filter(!condition %in% exclude_conditions)

# Bar summary
bar_data <- plot_data %>%
  group_by(condition) %>%
  summarise(contaminant_count = mean(contaminant_count), ecoli_count = mean(ecoli_count), .groups = "drop") %>%
  pivot_longer(cols = c(contaminant_count, ecoli_count), names_to = "protein_type", values_to = "count") %>%
  mutate(protein_type = recode(protein_type, contaminant_count = "Contaminant", ecoli_count = "E. coli"))

# Factor ordering
all_conditions_present <- unique(plot_data$condition)
desired_order <- c(
  "0cell", "1x_intact_us", "1x_intact_st", "1x_protop_ceph_st", "10x_protop_ceph_st", "10x_protop_st",
  "1pg", "5pg", "10pg", "20pg", "50pg", "100pg", "200pg"
)
desired_order_present <- intersect(desired_order, all_conditions_present)
leftover <- setdiff(all_conditions_present, desired_order_present)
final_condition_levels <- c(desired_order_present, sort(leftover))

plot_data <- plot_data %>% mutate(condition = factor(condition, levels = final_condition_levels))
bar_data  <- bar_data  %>% mutate(condition = factor(condition, levels = final_condition_levels))

# ---- Build annotation table (shared across plots) ----
table_row_names <- c("#cells", "bulk [pg]", "intact", "spheroplast", "stained")
annotation_cols <- c("condition", "#cells", "bulk [pg]", "intact", "spheroplast", "stained")

annotation_df <- file_mapping %>%
  select(any_of(annotation_cols)) %>%
  distinct(condition, .keep_all = TRUE) %>%
  filter(condition %in% final_condition_levels) %>%
  mutate(condition = factor(condition, levels = final_condition_levels)) %>%
  arrange(condition)

table_long <- annotation_df %>%
  mutate(across(-condition, as.character)) %>%
  pivot_longer(cols = -condition, names_to = "row_label", values_to = "cell_value") %>%
  mutate(
    row_label = factor(row_label, levels = rev(table_row_names)),
    x_pos = as.integer(condition)
  )

n_conditions <- length(final_condition_levels)

# ---- Shared function to build annotation table plot ----
make_table_plot <- function() {
  ggplot(table_long, aes(x = x_pos, y = row_label)) +
    geom_hline(
      yintercept = seq(0.5, length(table_row_names) + 0.5, by = 1),
      colour = "grey70", linewidth = 0.4
    ) +
    geom_text(aes(label = cell_value), size = 5, family = "Arial", hjust = 0.5, vjust = 0.5) +
    scale_x_continuous(
      limits = c(0.5, n_conditions + 0.5),
      breaks = seq_len(n_conditions),
      expand = c(0, 0)
    ) +
    scale_y_discrete() +
    theme_minimal(base_family = "Arial") +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 15, family = "Arial", hjust = 1),
      panel.grid   = element_blank(),
      plot.margin  = margin(t = 0, r = 5, b = 5, l = 5)
    )
}

# ---- Shared theme (reduced font sizes) ----
shared_theme <- function() {
  list(
    theme_classic(base_size = 12, base_family = "Arial"),
    theme(
      axis.title.y = element_text(size = 13, face = "bold", family = "Arial"),
      axis.title.x = element_blank(),
      axis.text    = element_text(size = 12, color = "black", family = "Arial"),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_text(size = 13, face = "bold", family = "Arial"),
      legend.text  = element_text(size = 12, family = "Arial"),
      axis.line    = element_line(linewidth = 0.8),
      axis.ticks   = element_line(linewidth = 0.8),
      plot.margin  = margin(t = 5, r = 5, b = 0, l = 5)
    )
  )
}



# Figure 2 A: Stacked bar plot with replicate dots ----

bar_annotation <- bar_data %>%
  group_by(condition) %>%
  summarise(total = sum(count), ecoli = sum(count[protein_type == "E. coli"]), .groups = "drop")

p1 <- ggplot() +
  geom_bar(
    data = bar_data,
    aes(x = condition, y = count, fill = protein_type),
    stat = "identity", position = "stack", width = 0.75, alpha = 0.85
  ) +
  geom_point(
    data = plot_data,
    aes(x = condition, y = total_count),
    position = position_jitter(width = 0.15, seed = 42),
    shape = 21, size = 3.2, stroke = 1, colour = "black", fill = "white"
  ) +
  geom_text(
    data = bar_annotation,
    aes(x = condition, y = total + max(total) * 0.05, label = round(ecoli)),
    family = "Arial", size = 4, vjust = 0
  ) +
  scale_fill_manual(
    values = c("Contaminant" = "#8F8989", "E. coli" = "#148509"),
    name = "Protein type",
    labels = c("Contaminant", expression(italic("E. coli")))
  ) +
  labs(x = NULL, y = "#quantified proteins") +
  shared_theme()

combined1 <- p1 / make_table_plot() + plot_layout(heights = c(4, 1))
combined1

# save the Fig. 2 A barplot
ggsave(file.path(plots_dir, "fig2_A_barplot.svg"),
       plot = combined1, device = "svg", width = 9.48, height = 4.68, units = "in")



# Figure 2 C: Stacked bar chart of total protein quantity (sum) ----

sum_data <- map_dfr(protein_cols, function(col) {
  condition <- col_condition_df$condition[col_condition_df$protein_col == col]
  base_name <- col_condition_df$base_name[col_condition_df$protein_col == col]
  if(is.na(condition) || condition %in% exclude_conditions) return(NULL)
  data_col <- coli_data[[col]]
  tibble(
    condition = condition, base_name = base_name,
    contaminant_quantity = sum(data_col[coli_data$protein_type == "Contaminant"], na.rm = TRUE),
    ecoli_quantity       = sum(data_col[coli_data$protein_type == "Ecoli"],        na.rm = TRUE),
    total_quantity       = sum(data_col[!is.na(data_col) & !is.nan(data_col)],     na.rm = TRUE)
  )
})

sum_data <- sum_data %>% mutate(condition = factor(condition, levels = final_condition_levels))

sum_summary <- sum_data %>%
  group_by(condition) %>%
  summarise(
    mean_contaminant = mean(contaminant_quantity), sd_contaminant = sd(contaminant_quantity),
    mean_ecoli       = mean(ecoli_quantity),       sd_ecoli       = sd(ecoli_quantity),
    mean_total       = mean(total_quantity),       sd_total       = sd(total_quantity),
    n = n(), .groups = "drop"
  ) %>%
  mutate(across(starts_with("sd_"), ~ replace_na(.x, 0)))

sum_summary_long <- sum_summary %>%
  pivot_longer(cols = c(mean_contaminant, mean_ecoli), names_to = "protein_type", values_to = "mean_quantity") %>%
  mutate(
    sd_quantity  = if_else(protein_type == "mean_contaminant", sd_contaminant, sd_ecoli),
    protein_type = recode(protein_type, mean_contaminant = "Contaminant", mean_ecoli = "E. coli")
  )

p3 <- ggplot() +
  geom_bar(
    data = sum_summary_long,
    aes(x = condition, y = mean_quantity, fill = protein_type),
    stat = "identity", position = "stack", width = 0.75, alpha = 0.85
  ) +
  geom_errorbar(
    data = sum_summary,
    aes(x = condition, ymin = mean_total - sd_total, ymax = mean_total + sd_total),
    width = 0.3, linewidth = 0.6
  ) +
  # Individual replicate dots (total quantity per replicate)
  geom_point(
    data = sum_data,
    aes(x = condition, y = total_quantity),
    position = position_jitter(width = 0.15, seed = 42),
    shape = 21, size = 3.2, stroke = 1, colour = "black", fill = "white"
  ) +
  scale_fill_manual(
    values = c("Contaminant" = "#8F8989", "E. coli" = "#148509"),
    name = "Protein type",
    labels = c("Contaminant", expression(italic("E. coli")))
  ) +
  labs(x = NULL, y = "Total protein quantity (sum)") +
  shared_theme()

combined3 <- p3 / make_table_plot() + plot_layout(heights = c(4, 1))
combined3

# save the Fig. 2 C barplot
ggsave(file.path(plots_dir, "fig2_C_barplot.svg"),
       plot = combined3, device = "svg", width = 9.48, height = 4.68, units = "in")


# Figure 2 zoom in from C: E. coli protein quantity only with reference line ----

reference_value <- sum_summary %>% filter(condition == "1x_intact_st") %>% pull(mean_ecoli)

p4 <- ggplot() +
  geom_bar(
    data = sum_summary,
    aes(x = condition, y = mean_ecoli),
    stat = "identity", width = 0.75, alpha = 0.85, fill = "#148509"
  ) +
  geom_errorbar(
    data = sum_summary,
    aes(x = condition, ymin = mean_ecoli - sd_ecoli, ymax = mean_ecoli + sd_ecoli),
    width = 0.3, linewidth = 0.6
  ) +
  # Individual replicate dots (E. coli quantity per replicate)
  geom_point(
    data = sum_data,
    aes(x = condition, y = ecoli_quantity),
    position = position_jitter(width = 0.15, seed = 42),
    shape = 21, size = 3.2, stroke = 1, colour = "black", fill = "white"
  ) +
  geom_hline(yintercept = reference_value, linetype = "dashed", linewidth = 0.5, color = "black") +
  labs(x = NULL, y = expression(bold(italic("E. coli")~"protein quantity (sum)"))) +
  coord_cartesian(ylim = c(0, 20000)) +
  shared_theme() +
  theme(legend.position = "none")

combined4 <- p4 / make_table_plot() + plot_layout(heights = c(4, 1))
combined4

# save the Fig. 2 zoom in barplot
ggsave(file.path(plots_dir, "fig2_zoom-in_barplot.svg"),
       plot = combined4, device = "svg", width = 9.48, height = 4.68, units = "in")



# Figure 2 B: Quantity histogram coloured by precursor count (NO annotation table) ----

# Contaminant proteins (PG.ProteinGroups starting with "cont_") are removed here:
# protein_type == "Ecoli" keeps only non-contaminant (E. coli) proteins.
coli_ecoli <- coli_data %>% filter(protein_type == "Ecoli")

cols_200pg <- col_condition_df %>% filter(condition == "200pg") %>% pull(protein_col) %>% unique()

if(length(cols_200pg) == 0) stop("No columns found for 200pg condition.")

# ---- Gather quantity + precursor count per protein per replicate (Figure 2 B) ----
# For each PG.Quantity column we look up the matching
# PG.NrOfPrecursorsUsedForQuantification and PG.QValue columns of the SAME raw file.
# We then keep only observations that pass the precursor (>= min_precursors) and
# Q value (<= max_qvalue) thresholds, and tag each surviving observation with its
# precursor count so the histogram can be coloured by precursor support level.
gather_quant_with_precursors <- function(df, quant_cols, min_precursors, max_qvalue) {
  all_cols <- colnames(df)
  map_dfr(quant_cols, function(qcol) {
    base        <- str_remove(qcol, "\\.PG\\.Quantity$")
    base_prefix <- paste0(base, ".PG.")

    prec_col <- all_cols[startsWith(all_cols, base_prefix) &
                    grepl("PG.NrOfPrecursorsUsedForQuantification", all_cols, fixed = TRUE)]
    qval_col <- all_cols[startsWith(all_cols, base_prefix) &
                    grepl("PG.QValue", all_cols, fixed = TRUE)]

    quantity   <- suppressWarnings(as.numeric(gsub(",", ".", as.character(df[[qcol]]))))
    precursors <- if (length(prec_col) == 1) suppressWarnings(as.numeric(df[[prec_col]])) else rep(NA_real_, nrow(df))
    qvalue     <- if (length(qval_col) == 1) suppressWarnings(as.numeric(df[[qval_col]])) else rep(NA_real_, nrow(df))

    if (length(prec_col) != 1)
      warning(sprintf("Could not uniquely match a precursor column for '%s' (found %d).", qcol, length(prec_col)))
    if (length(qval_col) != 1)
      warning(sprintf("Could not uniquely match a Q value column for '%s' (found %d).", qcol, length(qval_col)))

    tibble(
      PG.ProteinGroups = df$PG.ProteinGroups,
      protein_col      = qcol,
      quantity         = quantity,
      precursors       = precursors,
      qvalue           = qvalue
    )
  }) %>%
    filter(!is.na(quantity) & !is.nan(quantity) & quantity > 0) %>%
    filter(!is.na(precursors) & precursors >= min_precursors) %>%
    filter(!is.na(qvalue) & qvalue <= max_qvalue)
}

# ---- Base layer: all 200pg non-contaminant proteins (grey) ----
# Per (protein x replicate) observation of the 200pg sample, on the x = log10(200pg
# quantity) axis. This is the full reference distribution drawn in grey.
base_200pg <- gather_quant_with_precursors(
  coli_ecoli,
  quant_cols     = cols_200pg,
  min_precursors = min_precursors,
  max_qvalue     = max_qvalue
)

if(nrow(base_200pg) == 0) stop("No 200pg proteins passed the precursor / Q value filters for Figure 2 B.")

# ---- Overlay layer: matched single cell proteins (from 1x_intact_st) ----
# 1) Find proteins detected in the 1x_intact_st single cell condition (after filters)
#    and bin each one by its single cell precursor support (mean across detected
#    replicates) into 1 / 2 / 3+ precursors.
cols_1x_intact_st_B <- col_condition_df %>% filter(condition == "1x_intact_st") %>% pull(protein_col) %>% unique()
if(length(cols_1x_intact_st_B) == 0) stop("No columns found for 1x_intact_st condition (Figure 2 B overlay).")

sc_1x_obs <- gather_quant_with_precursors(
  coli_ecoli,
  quant_cols     = cols_1x_intact_st_B,
  min_precursors = min_precursors,
  max_qvalue     = max_qvalue
)

sc_1x_protein <- sc_1x_obs %>%
  group_by(PG.ProteinGroups) %>%
  summarise(mean_precursors = mean(precursors, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    prec_cat = case_when(
      mean_precursors < 1.5 ~ "1",
      mean_precursors < 2.5 ~ "2",
      TRUE                  ~ "3+"
    )
  )

# 2) The overlay shows those matched proteins AT THEIR 200pg quantity (same x axis as
#    the grey base), coloured by their single cell precursor category.
matched_overlay <- base_200pg %>%
  inner_join(sc_1x_protein %>% select(PG.ProteinGroups, prec_cat), by = "PG.ProteinGroups")

# ---- Shared fill scale (one combined legend for grey base + 3 green categories) ----
hist_levels <- c("200pg (all non contaminant proteins)", "n = 1", "n = 2", "n \u2265 3")
hist_colors <- c(
  "200pg (all non contaminant proteins)" = "#B5B5B5",
  "n = 1"                                = "#A1D99B",
  "n = 2"                                = "#41AB5D",
  "n \u2265 3"                           = "#006D2C"
)

base_plot_df <- base_200pg %>%
  mutate(fill_cat = factor("200pg (all non contaminant proteins)", levels = hist_levels))

overlay_plot_df <- matched_overlay %>%
  mutate(
    fill_cat = recode(prec_cat, "1" = "n = 1", "2" = "n = 2", "3+" = "n \u2265 3"),
    fill_cat = factor(fill_cat, levels = hist_levels)
  )

# Common bin width / origin so the green overlay bars line up exactly with the grey bars.
log_range <- range(log10(base_200pg$quantity))
bin_width  <- diff(log_range) / 50

p_hist <- ggplot() +
  # grey base: all 200pg non-contaminant proteins
  geom_histogram(
    data = base_plot_df,
    aes(x = log10(quantity), fill = fill_cat),
    binwidth = bin_width, boundary = 0, closed = "left",
    color = "black", linewidth = 0.2
  ) +
  # green overlay: matched single cell proteins, stacked by precursor category
  geom_histogram(
    data = overlay_plot_df,
    aes(x = log10(quantity), fill = fill_cat),
    binwidth = bin_width, boundary = 0, closed = "left",
    color = "black", linewidth = 0.2, position = "stack"
  ) +
  scale_fill_manual(
    values = hist_colors,
    name   = "",
    breaks = hist_levels,
    drop   = FALSE
  ) +
  labs(
    x = expression(bold(log[10]~"(Protein quantity at 200pg)")),
    y = "Count"
  ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    axis.title   = element_text(size = 14, face = "bold", family = "Arial"),
    axis.text    = element_text(size = 12, color = "black", family = "Arial"),
    legend.text  = element_text(size = 12, family = "Arial"),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

p_hist

# save the Fig. 2 B histogram
ggsave(file.path(plots_dir, "fig2_B_hist.svg"),
       plot = p_hist, device = "svg", width = 9.48, height = 4.68, units = "in")


# Figure 2 rank plots: 0cell and 1x_intact_st ----
# For each condition we average the per-replicate PG.Quantity values for every protein,
# rank the proteins by that average (descending), and plot rank vs log10(average quantity).
# Contaminant proteins (cont_) are grey, E. coli proteins green. Any accession listed in
# `proteins_to_annotate` is labelled with its gene name (PG.Genes) using ggrepel.

rank_point_colors <- c("Contaminant" = "#8F8989", "Ecoli" = "#148509")

make_rank_plot <- function(cond_name) {
  cols <- col_condition_df %>% filter(condition == cond_name) %>% pull(protein_col) %>% unique()
  if(length(cols) == 0) stop(sprintf("No columns found for condition '%s'.", cond_name))

  # Numeric matrix of the replicate quantities for this condition
  qmat <- coli_data %>%
    select(all_of(cols)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(gsub(",", ".", as.character(.))))))

  avg_quantity <- rowMeans(as.matrix(qmat), na.rm = TRUE)  # NaN when all replicates are NA

  rank_df <- tibble(
    PG.ProteinGroups = coli_data$PG.ProteinGroups,
    PG.Genes         = coli_data$PG.Genes,
    avg_quantity     = avg_quantity
  ) %>%
    filter(!is.na(avg_quantity) & !is.nan(avg_quantity) & avg_quantity > 0) %>%
    mutate(protein_type = if_else(str_starts(PG.ProteinGroups, "cont_"), "Contaminant", "Ecoli")) %>%
    arrange(desc(avg_quantity)) %>%
    mutate(rank = row_number())

  # Subset of proteins requested for annotation
  annotate_df <- rank_df %>% filter(PG.ProteinGroups %in% proteins_to_annotate)
  if(length(proteins_to_annotate) > 0) {
    missing <- setdiff(proteins_to_annotate, rank_df$PG.ProteinGroups)
    if(length(missing) > 0)
      warning(sprintf("Condition '%s': these proteins_to_annotate were not found / not quantified: %s",
                      cond_name, paste(missing, collapse = ", ")))
  }

  p_rank <- ggplot(rank_df, aes(x = rank, y = log10(avg_quantity))) +
    geom_point(aes(color = protein_type), size = 1.6, alpha = 0.85) +
    scale_color_manual(
      values = rank_point_colors,
      name   = "Protein type",
      labels = c("Contaminant" = "Contaminant", "Ecoli" = expression(italic("E. coli")))
    ) +
    labs(
      x = "Protein rank",
      y = expression(bold(log[10]~"(average protein quantity)")),
      title = cond_name
    ) +
    theme_classic(base_size = 14, base_family = "Arial") +
    theme(
      plot.title   = element_text(size = 14, face = "bold", family = "Arial", hjust = 0.5),
      axis.title   = element_text(size = 14, face = "bold", family = "Arial"),
      axis.text    = element_text(size = 12, color = "black", family = "Arial"),
      legend.title = element_text(size = 13, face = "bold", family = "Arial"),
      legend.text  = element_text(size = 12, family = "Arial"),
      legend.position = "top"
    )

  # Add gene-name labels only when proteins were requested AND matched
  if(nrow(annotate_df) > 0) {
    p_rank <- p_rank +
      geom_point(
        data = annotate_df,
        aes(x = rank, y = log10(avg_quantity)),
        shape = 21, size = 2.6, stroke = 1, colour = "black", fill = "white"
      ) +
      geom_text_repel(
        data = annotate_df,
        aes(x = rank, y = log10(avg_quantity), label = PG.Genes),
        family = "Arial", size = 4, fontface = "bold",
        box.padding = 0.6, point.padding = 0.3,
        min.segment.length = 0, max.overlaps = Inf,
        segment.color = "black", segment.size = 0.4
      )
  }

  p_rank
}

# The "0cell" rank plot is built with the simple helper (all points grey/green by type).
rank_0cell <- make_rank_plot("0cell")


# ---- Figure 2 rank plot: 1x_intact_st single cell overlaid on the 200pg reference ----
# This plot uses the 200pg sample as a fixed reference "ladder":
#   * BASE LAYER  : every protein quantified at 200pg, ranked by its 200pg average
#                   quantity (descending) and drawn in grey -> this defines the x (rank)
#                   and y (log10 average quantity) coordinate system.
#   * OVERLAY     : the proteins that were also detected in the 1x_intact_st single cell
#                   condition are re-plotted AT THEIR 200pg POSITION (same rank / quantity)
#                   and coloured by how many precursors supported them IN THE SINGLE CELL
#                   data (1, 2 or >=3), using the same three green shades as the histogram.
# This shows which part of the 200pg proteome is recovered from a single cell and how
# well (precursor support) each recovered protein was quantified.
make_overlay_rank_plot <- function(sc_condition = "1x_intact_st",
                                    ref_cols, sc_cols,
                                    min_precursors, max_qvalue) {

  # --- 1. 200pg reference ranking (grey background + coordinate system) ---
  ref_qmat <- coli_data %>%
    select(all_of(ref_cols)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(gsub(",", ".", as.character(.))))))
  ref_avg <- rowMeans(as.matrix(ref_qmat), na.rm = TRUE)  # NaN when all replicates are NA

  ref_rank_df <- tibble(
    PG.ProteinGroups = coli_data$PG.ProteinGroups,
    PG.Genes         = coli_data$PG.Genes,
    avg_quantity     = ref_avg
  ) %>%
    filter(!is.na(avg_quantity) & !is.nan(avg_quantity) & avg_quantity > 0) %>%
    arrange(desc(avg_quantity)) %>%
    mutate(rank = row_number())

  if(nrow(ref_rank_df) == 0) stop("No 200pg reference proteins available for the rank plot.")

  # --- 2. Per-protein precursor category in the single cell condition ---
  # Gather every (protein x replicate) observation in the single cell condition that
  # passes the precursor / Q value thresholds, then average the precursor count per
  # protein across its detected replicates and bin it into 1 / 2 / 3+.
  sc_obs <- gather_quant_with_precursors(
    coli_data,
    quant_cols     = sc_cols,
    min_precursors = min_precursors,
    max_qvalue     = max_qvalue
  )

  sc_protein <- sc_obs %>%
    group_by(PG.ProteinGroups) %>%
    summarise(mean_precursors = mean(precursors, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      prec_cat = case_when(
        mean_precursors < 1.5 ~ "1",
        mean_precursors < 2.5 ~ "2",
        TRUE                  ~ "3+"
      ),
      prec_cat = factor(prec_cat, levels = c("1", "2", "3+"))
    )

  # --- 3. Overlay: single cell proteins positioned at their 200pg rank / quantity ---
  overlay_df <- ref_rank_df %>%
    inner_join(sc_protein %>% select(PG.ProteinGroups, prec_cat), by = "PG.ProteinGroups")

  # Combined colour scale so the grey reference AND the three green categories share one
  # legend. Both layers map `grp` to the same factor levels.
  grp_levels <- c("200pg reference", "1", "2", "3+")
  combined_cols <- c(
    "200pg reference" = "#C8C8C8",
    "1"               = "#A1D99B",
    "2"               = "#41AB5D",
    "3+"              = "#006D2C"
  )
  ref_plot_df <- ref_rank_df %>%
    mutate(grp = factor("200pg reference", levels = grp_levels))
  overlay_plot_df <- overlay_df %>%
    mutate(grp = factor(as.character(prec_cat), levels = grp_levels))

  # --- 4. Annotation (works for both layers; positioned on the shared 200pg coords) ---
  annotate_df <- ref_rank_df %>% filter(PG.ProteinGroups %in% proteins_to_annotate)
  if(length(proteins_to_annotate) > 0) {
    missing <- setdiff(proteins_to_annotate, ref_rank_df$PG.ProteinGroups)
    if(length(missing) > 0)
      warning(sprintf("Overlay rank plot ('%s'): these proteins_to_annotate were not found in the 200pg reference: %s",
                      sc_condition, paste(missing, collapse = ", ")))
  }

  p_rank <- ggplot() +
    # base grey reference
    geom_point(
      data = ref_plot_df,
      aes(x = rank, y = log10(avg_quantity), color = grp),
      size = 1.6, alpha = 0.7
    ) +
    # single cell overlay (drawn on top), coloured by precursor support
    geom_point(
      data = overlay_plot_df,
      aes(x = rank, y = log10(avg_quantity), color = grp),
      size = 2.0, alpha = 0.95
    ) +
    scale_color_manual(
      values = combined_cols,
      name   = "Precursors\n(1x_intact_st)",
      breaks = grp_levels,
      labels = c(
        "200pg reference" = "200pg reference",
        "1"               = "1x cell: n = 1",
        "2"               = "1x cell: n = 2",
        "3+"              = "1x cell: n \u2265 3"
      ),
      drop = FALSE
    ) +
    labs(
      x = "Protein rank (200pg reference)",
      y = expression(bold(log[10]~"(average protein quantity, 200pg)")),
      title = sc_condition
    ) +
    theme_classic(base_size = 14, base_family = "Arial") +
    theme(
      plot.title   = element_text(size = 14, face = "bold", family = "Arial", hjust = 0.5),
      axis.title   = element_text(size = 14, face = "bold", family = "Arial"),
      axis.text    = element_text(size = 12, color = "black", family = "Arial"),
      legend.title = element_text(size = 13, face = "bold", family = "Arial"),
      legend.text  = element_text(size = 12, family = "Arial"),
      legend.position = "top"
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3)))

  # Add gene-name labels only when proteins were requested AND matched
  if(nrow(annotate_df) > 0) {
    p_rank <- p_rank +
      geom_point(
        data = annotate_df,
        aes(x = rank, y = log10(avg_quantity)),
        shape = 21, size = 2.6, stroke = 1, colour = "black", fill = "white"
      ) +
      geom_text_repel(
        data = annotate_df,
        aes(x = rank, y = log10(avg_quantity), label = PG.Genes),
        family = "Arial", size = 4, fontface = "bold",
        box.padding = 0.6, point.padding = 0.3,
        min.segment.length = 0, max.overlaps = Inf,
        segment.color = "black", segment.size = 0.4
      )
  }

  p_rank
}

cols_1x_intact_st <- col_condition_df %>% filter(condition == "1x_intact_st") %>% pull(protein_col) %>% unique()
if(length(cols_1x_intact_st) == 0) stop("No columns found for 1x_intact_st condition.")

rank_1x_intact_st <- make_overlay_rank_plot(
  sc_condition   = "1x_intact_st",
  ref_cols       = cols_200pg,
  sc_cols        = cols_1x_intact_st,
  min_precursors = min_precursors,
  max_qvalue     = max_qvalue
)

rank_0cell
rank_1x_intact_st

# save the Fig. 2 rank plots
ggsave(file.path(plots_dir, "fig2_rankplot_0cell.svg"),
       plot = rank_0cell, device = "svg", width = 6.5, height = 5.2, units = "in")
ggsave(file.path(plots_dir, "fig2_rankplot_1x_intact_st.svg"),
       plot = rank_1x_intact_st, device = "svg", width = 6.5, height = 5.2, units = "in")
