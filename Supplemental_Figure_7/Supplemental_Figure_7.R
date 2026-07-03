# Add user library path if needed
if(dir.exists('~/R/library')) {
  .libPaths(c('~/R/library', .libPaths()))
}

library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(uwot)
library(svglite)
library(cowplot)
library(gridExtra)
library(tidyverse)

# Define the subdirectory where plots should be saved
plots_dir <- "plots"  

# Read the protein quantities tab sep file, and select for columns wanted
BacSCP_data <- read_delim("data/protein_table_supl_fig_7_PGquan_are_with_in3prec.tsv", delim = "\t")

# Create a mapping table for ProteinGroups to Gene names
protein_to_gene <- BacSCP_data %>%
  select(PG.ProteinGroups, PG.Genes) %>%
  distinct()

# Select only PG.ProteinGroups column and columns containing ".raw.PG.Quantity"
BacSCP_data <- BacSCP_data %>%
  select(PG.ProteinGroups, contains(".raw.PG.Quantity"))

# Transpose the data so that samples are rows and proteins are columns
protein_data_t <- as.data.frame(t(BacSCP_data))
colnames(protein_data_t) <- BacSCP_data$PG.ProteinGroups
colnames(protein_data_t) <- make.names(colnames(protein_data_t), unique = TRUE)
protein_data_t <- rownames_to_column(protein_data_t, var = "Filename")
protein_data_t <- protein_data_t %>% slice(-1)

# Extract conditions
protein_data_t <- protein_data_t %>%
  mutate(condition = case_when(
    grepl("THIDTC006_0cell_", Filename) ~ "0cell",
    grepl("BS_dMcsB_chil_10x_", Filename) ~ "10x_BS_chil",
    grepl("BS_dMcsB_HS_10x_", Filename)   ~ "10x_BS_HS",
    grepl("BS_dMcsB_chil_", Filename) & !grepl("10x_", Filename) ~ "1x_BS_chil",
    grepl("BS_dMcsB_HS_",   Filename) & !grepl("10x_", Filename) ~ "1x_BS_HS",
    TRUE ~ NA_character_
  ))

# Identify contaminant and protein columns
cont_cols <- grep("^cont_", colnames(protein_data_t), value = TRUE)
prot_cols <- setdiff(colnames(protein_data_t), c("Filename", "condition", cont_cols))

# ── Helper: extract short name (2 chars before .raw.PG.Quantity) ─────────────
make_short_name <- function(x) {
  sub(".*(.{2})\\.raw\\.PG\\.Quantity.*", "\\1", x)
}

# Figure 6 E: Rank plot -------------------------------------------------------------------------
condition_colors <- c("1x_BS_chil" = "#E69F00", "1x_BS_HS" = "#7B2D8B")

avg_abundances <- protein_data_t %>%
  filter(condition %in% c("1x_BS_chil", "1x_BS_HS")) %>%
  select(-Filename) %>%
  pivot_longer(cols = -condition, names_to = "ProteinGroup", values_to = "Quantity") %>%
  mutate(Quantity = as.numeric(Quantity)) %>%
  filter(!is.na(Quantity), Quantity > 0) %>%
  group_by(condition, ProteinGroup) %>%
  summarise(avg_quantity = mean(Quantity, na.rm = TRUE), .groups = "drop")

avg_abundances <- avg_abundances %>%
  left_join(protein_to_gene, by = c("ProteinGroup" = "PG.ProteinGroups"))

avg_abundances <- avg_abundances %>%
  group_by(condition) %>%
  arrange(desc(avg_quantity)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

target_genes <- c("trypsin", "clpC", "groEL", "groES", "ctc", "pgk", "tpiA", "dnaK")

top1_per_condition <- avg_abundances %>%
  group_by(condition) %>%
  slice_min(rank, n = 1) %>%
  ungroup()

avg_labels <- avg_abundances %>%
  filter(
    paste(condition, ProteinGroup) %in% paste(top1_per_condition$condition, top1_per_condition$ProteinGroup) |
      sapply(strsplit(as.character(PG.Genes), ";"), function(x) any(target_genes %in% trimws(x)))
  ) %>%
  group_by(condition, PG.Genes) %>%
  slice_max(avg_quantity, n = 1) %>%
  ungroup()

rankplot <- ggplot(avg_abundances, aes(x = rank, y = avg_quantity, color = condition)) +
  geom_line(linewidth = 0.8, alpha = 0.5) +
  geom_point(
    data      = avg_labels,
    aes(color = condition),
    size      = 3,
    shape     = 21,
    fill      = "white",
    stroke    = 1.5
  ) +
  geom_label_repel(
    data          = avg_labels,
    aes(label     = PG.Genes, color = condition),
    size          = 3.5,
    fontface      = "bold",
    box.padding   = 0.5,
    point.padding = 0.4,
    segment.size  = 0.4,
    segment.color = "grey40",
    show.legend   = FALSE,
    max.overlaps  = Inf,
    force         = 8,
    fill          = alpha("white", 0.8)
  ) +
  scale_color_manual(
    values = condition_colors,
    labels = c("1x_BS_chil" = "37°C", "1x_BS_HS" = "50°C")
  ) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 10000, by = 25)) +
  labs(
    x     = "Rank (by abundance)",
    y     = "Average Protein Quantity (log10)",
    color = "Condition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title             = element_blank(),
    plot.subtitle          = element_blank(),
    axis.title             = element_text(size = 14, face = "bold"),
    axis.text              = element_text(size = 12, color = "black"),
    axis.text.x            = element_text(angle = 45, hjust = 1),
    legend.title           = element_text(size = 13, face = "bold"),
    legend.text            = element_text(size = 12),
    legend.position        = "inside",
    legend.position.inside = c(0.85, 0.85),
    legend.background      = element_rect(fill = alpha("white", 0.8), color = NA),
    panel.grid.major       = element_blank(),
    panel.grid.minor       = element_blank()
  )

rankplot

# save the rankplot
ggsave(file.path(plots_dir, "fig6_E_rankplot.svg"),
       plot = rankplot, device = "svg", width = 9.48, height = 4.68, units = "in")

# Figure 6 A: Number of quantified proteins -------------------------------------------------------------------------
protein_data_t <- protein_data_t %>%
  mutate(
    n_entries_proteins     = rowSums(!is.na(select(., all_of(prot_cols))) &
                                       select(., all_of(prot_cols)) != 0),
    n_entries_contaminants = rowSums(!is.na(select(., all_of(cont_cols))) &
                                       select(., all_of(cont_cols)) != 0)
  )

plot_data <- protein_data_t %>%
  select(Filename, condition, n_entries_proteins, n_entries_contaminants) %>%
  pivot_longer(
    cols      = starts_with("n_entries"),
    names_to  = "type",
    values_to = "n_entries"
  ) %>%
  mutate(
    type = recode(type,
                  n_entries_proteins     = "B. Subtilis",
                  n_entries_contaminants = "Contaminants"),
    condition = case_when(
      condition == "1x_BS_chil" ~ "37°C",
      condition == "1x_BS_HS"   ~ "50°C",
      TRUE ~ condition
    )
  ) %>%
  mutate(type = factor(type, levels = c("Contaminants", "B. Subtilis")))

summary_data <- plot_data %>%
  group_by(condition, type) %>%
  summarise(mean_entries = mean(n_entries, na.rm = TRUE), .groups = "drop") %>%
  group_by(condition) %>%
  arrange(condition, desc(type)) %>%
  mutate(
    cumulative = cumsum(mean_entries),
    label_y    = cumulative - (mean_entries / 2)
  ) %>%
  ungroup()

dots_data <- protein_data_t %>%
  mutate(
    total_entries = n_entries_proteins + n_entries_contaminants,
    condition = case_when(
      condition == "1x_BS_chil" ~ "37°C",
      condition == "1x_BS_HS"   ~ "50°C",
      TRUE ~ condition
    )
  )

desired_order <- c("0cell", "37°C", "50°C", "10x_BS_chil", "10x_BS_HS")
summary_data$condition <- factor(summary_data$condition, levels = desired_order)
dots_data$condition    <- factor(dots_data$condition,    levels = desired_order)

p <- ggplot() +
  geom_bar(
    data     = summary_data,
    aes(x = condition, y = mean_entries, fill = type),
    stat     = "identity",
    position = "stack",
    width    = 0.7
  ) +
  geom_text(
    data     = summary_data,
    aes(x = condition, y = label_y, label = round(mean_entries), color = type),
    size     = 4,
    fontface = "bold"
  ) +
  geom_jitter(
    data   = dots_data,
    aes(x = condition, y = total_entries),
    shape  = 21,
    fill   = "white",
    color  = "black",
    width  = 0.2,
    height = 0,
    size   = 2.5,
    stroke = 0.5
  ) +
  scale_fill_manual(
    values = c("B. Subtilis" = "#56B4E9", "Contaminants" = "#8F8989"),
    labels = c(
      "B. Subtilis"  = expression(italic("B. Subtilis")),
      "Contaminants" = "Contaminants"
    )
  ) +
  scale_color_manual(
    values = c("B. Subtilis" = "white", "Contaminants" = "black"),
    guide  = "none"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title             = element_blank(),
    axis.title             = element_text(size = 14, face = "bold"),
    axis.text              = element_text(size = 12, color = "black"),
    axis.text.x            = element_text(angle = 45, hjust = 1),
    legend.title           = element_blank(),
    legend.text            = element_text(size = 12),
    legend.position        = "inside",
    legend.position.inside = c(0.15, 0.85),
    legend.background      = element_rect(fill = alpha("white", 0.6), color = NA),
    panel.grid.major       = element_blank(),
    panel.grid.minor       = element_blank(),
    panel.border           = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(
    y    = "Number of quantified proteins",
    x    = "",
    fill = "Entry Type"
  )


p

# save the barplot p
ggsave(file.path(plots_dir, "fig6_A_barplot.svg"),
       plot = p, device = "svg", width = 9.48, height = 4.68, units = "in")


# prepare data for PCA and heatmap ----------------------------------------
protein_data_num <- protein_data_t %>%
  mutate(across(all_of(prot_cols), ~ as.numeric(.)))

mat <- protein_data_num[, prot_cols]
mat <- apply(mat, 2, as.numeric)

# Impute missing values (1st percentile per sample)
impute_min_col <- function(x) {
  # Preserve the vector names even on the all-NA branch, otherwise apply()
  # drops the protein names from the resulting matrix (see fix below).
  if (all(is.na(x))) return(setNames(rep(NA_real_, length(x)), names(x)))
  q <- quantile(x, 0.01, na.rm = TRUE)
  x[is.na(x)] <- q
  return(x)
}

mat_imp <- t(apply(mat, 1, impute_min_col))
# IMPORTANT: re-attach the protein column names. When any row is entirely NA,
# apply() cannot build consistent dimnames and silently drops the column names,
# so as.data.frame() would relabel them V1, V2, ... . The column order is always
# preserved by apply(), so restoring colnames from `mat` (== prot_cols) is safe
# and guarantees pivot_longer(all_of(prot_cols)) can find the columns.
colnames(mat_imp) <- colnames(mat)
mat_imp <- as.data.frame(mat_imp)
mat_imp$condition <- protein_data_t$condition
mat_imp$Filename  <- protein_data_t$Filename

data_long <- mat_imp %>%
  filter(condition %in% c("1x_BS_chil", "1x_BS_HS")) %>%
  pivot_longer(
    cols      = all_of(prot_cols),
    names_to  = "Protein",
    values_to = "Intensity"
  ) %>%
  mutate(log2Int = log2(as.numeric(Intensity) + 1))

data_stats <- data_long %>%
  group_by(Protein) %>%
  summarise(
    log2FC = mean(log2Int[condition == "1x_BS_HS"],   na.rm = TRUE) -
      mean(log2Int[condition == "1x_BS_chil"], na.rm = TRUE),
    pvalue = tryCatch(
      t.test(
        log2Int[condition == "1x_BS_HS"],
        log2Int[condition == "1x_BS_chil"]
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    negLog10P      = -log10(pvalue),
    is_contaminant = FALSE
  )

data_stats <- data_stats %>%
  left_join(protein_to_gene, by = c("Protein" = "PG.ProteinGroups"))

fc_threshold   <- 1.0
pval_threshold <- 0.05

data_stats <- data_stats %>%
  mutate(
    is_significant = !is.na(pvalue) &
      pvalue < pval_threshold &
      abs(log2FC) > fc_threshold,
    gene_label = if_else(is_significant & !is.na(PG.Genes), PG.Genes, NA_character_),
    point_color = factor(
      case_when(
        is_significant ~ "significant",
        TRUE           ~ "protein"
      ),
      levels = c("protein", "significant")
    )
  )

# Figure 6 F: PCA plot: 1x_BS_chil vs 1x_BS_HS -------------------------------------------------------------------------

# 1. Prepare data for PCA — EXCLUDE contaminants
pca_data <- mat_imp %>%
  filter(condition %in% c("1x_BS_chil", "1x_BS_HS"))

# 2. Extract only NON-contaminant protein columns
pca_matrix <- pca_data %>%
  select(all_of(prot_cols)) %>%
  as.matrix()

# 3. Remove columns with zero variance or all NA
col_var    <- apply(pca_matrix, 2, var, na.rm = TRUE)
pca_matrix <- pca_matrix[, !is.na(col_var) & col_var > 0]

# 3b. Remove rows with too many NA values (>50% of columns)
row_na_pct <- rowSums(is.na(pca_matrix)) / ncol(pca_matrix)
valid_rows <- row_na_pct < 0.5
pca_matrix <- pca_matrix[valid_rows, ]
pca_data   <- pca_data[valid_rows, ]

cat("Removed", sum(!valid_rows), "sample(s) with >50% missing values for PCA\n")
cat("PCA will be performed on", nrow(pca_matrix), "samples\n")

# 4. Perform PCA
pca_result <- prcomp(pca_matrix, scale. = TRUE, center = TRUE)

# 5. Extract PC scores and add condition + short name labels
pca_scores           <- as.data.frame(pca_result$x)
pca_scores$condition <- pca_data$condition
pca_scores$Filename  <- pca_data$Filename
pca_scores$short_name <- make_short_name(pca_scores$Filename)

# 6. Calculate variance explained
var_explained <- summary(pca_result)$importance[2, ] * 100

# 7. Plot PCA with sample labels
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text_repel(
    aes(label = short_name),
    size               = 3.5,
    fontface           = "bold",
    box.padding        = 0.5,
    point.padding      = 0.4,
    segment.size       = 0.4,
    segment.color      = "grey40",
    show.legend        = FALSE,
    max.overlaps       = Inf,
    force              = 5,
    fill               = alpha("white", 0.7)
  ) +
  scale_color_manual(
    name   = "Condition",
    values = condition_colors,
    labels = c("1x_BS_chil" = expression(italic("B. Subtilis") ~ "37°C"),
               "1x_BS_HS"   = expression(italic("B. Subtilis") ~ "50°C"))
  ) +
  stat_ellipse(aes(color = condition), level = 0.95, linetype = 2, linewidth = 1.2) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title             = element_blank(),
    axis.title             = element_text(size = 14, face = "bold"),
    axis.text              = element_text(size = 12, color = "black"),
    legend.title           = element_text(size = 13, face = "bold"),
    legend.text            = element_text(size = 12),
    legend.position        = "inside",
    legend.position.inside = c(0.1, 0.1),
    legend.background      = element_rect(fill = alpha("white", 0.8), color = NA),
    panel.grid.major       = element_blank(),
    panel.grid.minor       = element_blank(),
    panel.border           = element_rect(color = "black", fill = NA, linewidth = 1)
  )

pca_plot

# save the PCA pca_plot
ggsave(file.path(plots_dir, "fig6_F_PCA.svg"),
       plot = pca_plot, device = "svg", width = 9.48, height = 4.68, units = "in")

# prepare Heatmap data -------------------------------------------------------------------------
heatmap_data <- mat_imp %>%
  filter(condition %in% c("1x_BS_chil", "1x_BS_HS"))

heatmap_matrix <- heatmap_data %>%
  select(all_of(prot_cols)) %>%
  as.matrix()

# Remove columns with zero variance or all NA
col_var        <- apply(heatmap_matrix, 2, var, na.rm = TRUE)
heatmap_matrix <- heatmap_matrix[, !is.na(col_var) & col_var > 0]

# Remove rows (samples) with too many NA values (>50% of columns)
row_na_pct_hm <- rowSums(is.na(heatmap_matrix)) / ncol(heatmap_matrix)
valid_rows_hm <- row_na_pct_hm < 0.5
heatmap_matrix <- heatmap_matrix[valid_rows_hm, ]
heatmap_data   <- heatmap_data[valid_rows_hm, ]

cat("Removed", sum(!valid_rows_hm), "sample(s) with >50% missing values for heatmap\n")
cat("Heatmap will be created with", nrow(heatmap_matrix), "samples\n")

# Transpose so proteins are rows and samples are columns
heatmap_matrix_t <- t(heatmap_matrix)

# Scale by row (z-score normalization per protein)
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix_t)))

# Short names: always the 2 chars before ".raw.PG.Quantity"
short_names <- make_short_name(heatmap_data$Filename)

# Ensure ALL heatmap matrices carry these sample names
colnames(heatmap_matrix_scaled) <- short_names

# Annotation
annotation_col <- data.frame(
  Condition = heatmap_data$condition,
  row.names = short_names
)

# Sorted annotation (for grouped columns in the "regulated" heatmap)
annotation_col_sorted <- annotation_col %>%
  rownames_to_column("Sample") %>%
  arrange(Condition) %>%
  column_to_rownames("Sample")

# Display annotation (legend shows temperatures)
display_annotation <- annotation_col_sorted %>%
  mutate(Condition = case_when(
    Condition == "1x_BS_chil" ~ "37°C",
    Condition == "1x_BS_HS"   ~ "50°C",
    TRUE ~ Condition
  ))

display_annotation_colors <- list(
  Condition = c("37°C" = "#E69F00", "50°C" = "#7B2D8B")
)

# Figure 6 C: Heatmap: Top 10 most abundant B. Subtilis -------------------------------------------------------------------------

top10_abundant_proteins <- heatmap_data %>%
  select(all_of(prot_cols)) %>%
  summarise(across(everything(), ~ mean(as.numeric(.), na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Protein", values_to = "mean_abundance") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Protein)

heatmap_top10_abund <- heatmap_matrix_scaled[
  rownames(heatmap_matrix_scaled) %in% top10_abundant_proteins, , drop = FALSE
]

# keep same sample labels explicitly
colnames(heatmap_top10_abund) <- colnames(heatmap_matrix_scaled)

top10_gene_labels <- protein_to_gene %>%
  filter(PG.ProteinGroups %in% top10_abundant_proteins) %>%
  mutate(clean_col = make.names(PG.ProteinGroups, unique = TRUE)) %>%
  select(clean_col, PG.Genes)

row_labels_abund <- setNames(
  top10_gene_labels$PG.Genes[match(rownames(heatmap_top10_abund), top10_gene_labels$clean_col)],
  rownames(heatmap_top10_abund)
)
row_labels_abund[is.na(row_labels_abund)] <- rownames(heatmap_top10_abund)[is.na(row_labels_abund)]

heatmap_top10abundant <- pheatmap(
  heatmap_top10_abund,
  annotation_col           = display_annotation[ colnames(heatmap_top10_abund), , drop = FALSE ],
  annotation_colors        = display_annotation_colors,
  labels_col               = colnames(heatmap_top10_abund),     # << explicit
  color                    = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  scale                    = "none",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  show_rownames            = TRUE,
  labels_row               = row_labels_abund,
  show_colnames            = TRUE,                              # << show them here too
  fontsize_col             = 8,
  fontsize_row             = 10,
  main                     = "Heatmap – Top 10 most abundant B. Subtilis proteins",
  border_color             = NA
)

heatmap_top10abundant

# save the heatmap top 10 abundant proteins
ggsave(file.path(plots_dir, "fig6_C_heatmap_top10abundant.svg"),
       plot = heatmap_top10abundant, device = "svg", width = 9.48, height = 4.68, units = "in")

# Figure 6 D: Heatmap: Top 10 regulated proteins + clpC check -------------------------------------------------------------------------
clpc_id <- protein_to_gene %>%
  filter(grepl("clpC", PG.Genes, ignore.case = TRUE)) %>%
  pull(PG.ProteinGroups) %>%
  unique()

top_others <- data_stats %>%
  filter(!is_contaminant, !is.na(log2FC), !is.na(pvalue)) %>%
  filter(!(Protein %in% clpc_id)) %>%
  arrange(desc(is_significant), desc(abs(log2FC))) %>%
  slice_head(n = 9) %>%
  pull(Protein)

top10_regulated_proteins <- c(clpc_id, top_others)

# reorder columns to match sorted annotation
heatmap_top10_reg <- heatmap_matrix_scaled[
  rownames(heatmap_matrix_scaled) %in% top10_regulated_proteins,
  rownames(annotation_col_sorted),
  drop = FALSE
]

# after reordering, labels must follow that order
reg_labels_col <- colnames(heatmap_top10_reg)

top10_reg_gene_labels <- protein_to_gene %>%
  filter(PG.ProteinGroups %in% top10_regulated_proteins) %>%
  mutate(clean_col = make.names(PG.ProteinGroups, unique = TRUE)) %>%
  select(clean_col, PG.Genes)

row_labels_reg <- setNames(
  top10_reg_gene_labels$PG.Genes[match(rownames(heatmap_top10_reg), top10_reg_gene_labels$clean_col)],
  rownames(heatmap_top10_reg)
)
row_labels_reg[is.na(row_labels_reg)] <- rownames(heatmap_top10_reg)[is.na(row_labels_reg)]

heatmap_top10regulated <- pheatmap(
  heatmap_top10_reg,
  annotation_col           = display_annotation[ reg_labels_col, , drop = FALSE ],
  annotation_colors        = display_annotation_colors,
  labels_col               = reg_labels_col,                    # << explicit & ordered
  color                    = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  scale                    = "none",
  cluster_cols             = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method        = "complete",
  show_rownames            = TRUE,
  labels_row               = row_labels_reg,
  show_colnames            = TRUE,                              # << show them here too
  fontsize_col             = 8,
  fontsize_row             = 10,
  main                     = "Heatmap – top 10 regulated B. Subtilis proteins",
  border_color             = NA
)

heatmap_top10regulated

# save the heatmap top 10 regulated proteins
ggsave(file.path(plots_dir, "fig6_D_heatmap_top10regulated.svg"),
       plot = heatmap_top10regulated, device = "svg", width = 9.48, height = 4.68, units = "in")


# Figure 6 B: Volcanoplot -----------------------------------------------------------
data_dir_volcano <- file.path("data/volcano/")

# ── 1. List & categorize files ────────────────────────────────────────────────  
all_files <- list.files(path = data_dir_volcano, pattern = "\\.txt$", full.names = TRUE)  

# volcano files are the files exported from Perseus from the volcano data
volcano_files <- all_files[grepl("volcanodata", all_files)]  

# curve files are the files exported from Perseus from the curves of the respective volcano data
curve_files   <- all_files[grepl("curvedata",   all_files)]

# ── 2. Helper: extract sample ID + condition from filename ────────────────────  
# e.g. "TC006_volcanodata_1pg.txt" → sample = "TC004", condition = "1pg"  
parse_filename <- function(path) {  
  fname <- tools::file_path_sans_ext(basename(path))  
  parts <- strsplit(fname, "_")[[1]]  
  list(  
    sample    = parts[1],                        # e.g. "TC004"  
    condition = parts[length(parts)]             # e.g. "1pg" or "200pg"  
  )  
}


# ── 3. Read all volcano data ──────────────────────────────────────────────────  
volcano_data <- map_dfr(volcano_files, function(f) {  
  meta <- parse_filename(f)  
  read_tsv(f, show_col_types = FALSE) %>%  
    mutate(  
      sample    = meta$sample,  
      condition = meta$condition,  
      source    = basename(f)  
    )  
})

# ── 4. Read all curve data ────────────────────────────────────────────────────  
curve_data <- map_dfr(curve_files, function(f) {  
  meta <- parse_filename(f)  
  read_tsv(f, show_col_types = FALSE) %>%  
    mutate(  
      sample    = meta$sample,  
      condition = meta$condition,  
      source    = basename(f)  
    )  
})

# ── 5. Preview ────────────────────────────────────────────────────────────────  
glimpse(volcano_data)  
glimpse(curve_data)  

# Check which sample/condition combos were loaded  
volcano_data %>% distinct(sample, condition)  
curve_data   %>% distinct(sample, condition)

# ── 6. Plot volcano ────────────────────────────────────────────────────────────────  

# Define genes of interest
genes_of_interest <- c("clpC", "groES", "groEL", "ctc", "pgk", "tpiA", "dnaK")
grey_genes <- c("ctc", "pgk", "groEL")

# Filter data and add a specific color column for labels
v_sc <- volcano_data %>% 
  filter(condition == "sc") %>%
  mutate(
    is_goi = PG.Genes %in% genes_of_interest,
    label_color = if_else(PG.Genes %in% grey_genes, "grey50", "#0072B2")
  )

c_sc <- curve_data %>% filter(condition == "sc")

p_sc <- ggplot(v_sc, aes(x = Difference, y = `-Log(P-value)`)) +
  
  # 1. All points (Grey background)
  geom_point(color = "grey", size = 1.8, alpha = 0.7) +
  
  # 2. Genes of interest (colored by group)
  geom_point(data = v_sc %>% filter(is_goi),
             aes(color = label_color), size = 2, alpha = 0.9) +
  
  # 3. Perseus Curve
  geom_line(data = c_sc, aes(x = x, y = y),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # 4. Labels with conditional coloring
  geom_text_repel(
    data = v_sc %>% filter(is_goi),
    aes(label = PG.Genes, color = label_color, segment.color = label_color),
    size = 3.5,
    fontface = "bold",
    segment.size = 0.3,
    segment.alpha = 0.7,
    min.segment.length = 0,
    box.padding = 0.6,
    point.padding = 0.4,
    force = 2,
    force_pull = 1,
    max.overlaps = Inf,
    direction = "both"
  ) +
  # Map the colors literally from the data column
  scale_color_identity() +
  
  # Formatting & Styling
  coord_cartesian(ylim = c(0, 11), xlim = c(-0.5, 2)) +
  scale_y_continuous(breaks = seq(0, 11, 2)) +
  scale_x_continuous(breaks = seq(-0.5, 2, 0.5)) +
  labs(x = expression("Difference (50 °C " * Delta * "mcsB / 37 °C " * Delta * "mcsB)"),
       y = expression("-log"[10]*" p-value")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

p_sc

# save the volcano plot
ggsave(file.path(plots_dir, "fig6_B_volcanoplot.svg"),
       plot = p_sc, device = "svg", width = 12, height = 6, units = "cm")

