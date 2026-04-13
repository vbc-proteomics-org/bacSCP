library(tidyverse)
library(extrafont)
library(patchwork)

loadfonts(device = "win", quiet = TRUE)

# define directories, read in data and prepare data ------------------------------------
# Define the subdirectory where plots should be saved 
plots_dir <- "plots"  

# Read input files
bacillus_data <- read_delim("data/protein_table_fig3.tsv", delim = "\t")
file_mapping   <- read_delim("data/file_mapping_fig3.tsv", delim = "\t")

# Columns containing PG.Quantity
protein_cols <- grep("PG.Quantity", colnames(bacillus_data), value = TRUE)

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
  base_name   = base_names
) %>%
  left_join(file_mapping %>% select(FileName, condition),
            by = c("base_name" = "FileName"))

if(any(is.na(col_condition_df$condition))) {
  warning("Some filenames could not be matched to conditions.")
  print(col_condition_df[is.na(col_condition_df$condition), ])
}

# Mark contaminant vs B. subtilis
bacillus_data <- bacillus_data %>%
  mutate(protein_type = if_else(str_starts(PG.ProteinGroups, "cont_"), "Contaminant", "B.subtilis"))

# Count quantified proteins per replicate column
plot_data <- map_dfr(protein_cols, function(col) {
  condition <- col_condition_df$condition[col_condition_df$protein_col == col]
  base_name <- col_condition_df$base_name[col_condition_df$protein_col == col]
  data_col  <- bacillus_data[[col]]
  contaminant_count <- sum(!is.na(data_col[bacillus_data$protein_type == "Contaminant"]) &
                             !is.nan(data_col[bacillus_data$protein_type == "Contaminant"]))
  bsubt_count <- sum(!is.na(data_col[bacillus_data$protein_type == "B.subtilis"]) &
                       !is.nan(data_col[bacillus_data$protein_type == "B.subtilis"]))
  total_count <- contaminant_count + bsubt_count
  tibble(protein_col = col, base_name = base_name, condition = condition,
         contaminant_count = contaminant_count, bsubt_count = bsubt_count, total_count = total_count)
})

plot_data <- plot_data %>% filter(!is.na(condition))
exclude_conditions <- c("none")
plot_data <- plot_data %>% filter(!condition %in% exclude_conditions)

# Bar summary with explicit factor levels
bar_data <- plot_data %>%
  group_by(condition) %>%
  summarise(contaminant_count = mean(contaminant_count), bsubt_count = mean(bsubt_count), .groups = "drop") %>%
  pivot_longer(cols = c(contaminant_count, bsubt_count), names_to = "protein_type", values_to = "count") %>%
  mutate(protein_type = recode(protein_type, contaminant_count = "Contaminant", bsubt_count = "B.subtilis")) %>%
  mutate(protein_type = factor(protein_type, levels = c("Contaminant", "B.subtilis")))

# Factor ordering for conditions
all_conditions_present <- unique(plot_data$condition)
desired_order <- c(
  "0cell",
  sort(all_conditions_present[str_starts(all_conditions_present, "1x")]),
  sort(all_conditions_present[str_starts(all_conditions_present, "10x")]),
  "1pg", "5pg", "10pg", "20pg", "50pg", "100pg", "200pg"
)
desired_order_present  <- intersect(desired_order, all_conditions_present)
final_condition_levels <- c(desired_order_present, sort(setdiff(all_conditions_present, desired_order_present)))

plot_data <- plot_data %>% mutate(condition = factor(condition, levels = final_condition_levels))
bar_data  <- bar_data  %>% mutate(condition = factor(condition, levels = final_condition_levels))

# ---- Build annotation table ----
table_row_names <- c("#cells", "intact", "protoplast", "stained")
annotation_cols <- c("condition", "#cells", "intact", "protoplast", "stained")

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

# ---- Shared annotation table function ----
make_table_plot <- function() {
  ggplot(table_long, aes(x = x_pos, y = row_label)) +
    geom_hline(yintercept = seq(0.5, length(table_row_names) + 0.5, by = 1),
               colour = "grey70", linewidth = 0.4) +
    geom_text(aes(label = cell_value), size = 5, family = "Arial", hjust = 0.5, vjust = 0.5) +
    scale_x_continuous(limits = c(0.5, n_conditions + 0.5),
                       breaks = seq_len(n_conditions), expand = c(0, 0)) +
    theme_minimal(base_family = "Arial") +
    theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
          axis.title   = element_blank(),
          axis.text.y  = element_text(size = 15, family = "Arial", hjust = 1),
          panel.grid   = element_blank(),
          plot.margin  = margin(t = 0, r = 5, b = 5, l = 5))
}

# ---- Shared theme ----
shared_theme <- function() {
  list(
    theme_classic(base_size = 12, base_family = "Arial"),
    theme(axis.title.y  = element_text(size = 13, face = "bold", family = "Arial"),
          axis.title.x  = element_blank(),
          axis.text     = element_text(size = 12, color = "black", family = "Arial"),
          axis.text.x   = element_blank(),
          axis.ticks.x  = element_blank(),
          legend.title  = element_text(size = 13, face = "bold", family = "Arial"),
          legend.text   = element_text(size = 12, family = "Arial"),
          axis.line     = element_line(linewidth = 0.8),
          axis.ticks    = element_line(linewidth = 0.8),
          plot.margin   = margin(t = 5, r = 5, b = 0, l = 5))
  )
}

# Colors: Contaminant = grey, B.subtilis = blue
bsubt_colors <- c("Contaminant" = "#8F8989", "B.subtilis" = "#56B4E9")


# Figure 3 A: Stacked bar plot with replicate dots -------------------------------------------------------------------------

bar_annotation <- bar_data %>%
  group_by(condition) %>%
  summarise(total = sum(count), bsubt = sum(count[protein_type == "B.subtilis"]), .groups = "drop")

p1 <- ggplot() +
  geom_bar(data = bar_data, aes(x = condition, y = count, fill = protein_type),
           stat = "identity", position = "stack", width = 0.75, alpha = 0.85) +
  geom_point(data = plot_data, aes(x = condition, y = total_count),
             position = position_jitter(width = 0.15, seed = 42),
             shape = 21, size = 3.2, stroke = 1, colour = "black", fill = "white") +
  geom_text(data = bar_annotation,
            aes(x = condition, y = total + max(total) * 0.05, label = round(bsubt)),
            family = "Arial", size = 4, vjust = 0) +
  scale_fill_manual(values = bsubt_colors, name = "Protein type",
                    labels = c("Contaminant", expression(italic("B. subtilis")))) +
  labs(y = "#quantified proteins") +
  shared_theme()

combined1 <- p1 / make_table_plot() + plot_layout(heights = c(4, 1))
combined1

# save the Fig. 3 A stacked barplot
ggsave(file.path(plots_dir, "fig3_A_barplot.svg"),
       plot = combined1, device = "svg", width = 9.48, height = 4.68, units = "in")

# Figure 3 B: Stacked bar chart of total protein quantity (sum)-------------------------------------------------------------------------

sum_data <- map_dfr(protein_cols, function(col) {
  condition <- col_condition_df$condition[col_condition_df$protein_col == col]
  base_name <- col_condition_df$base_name[col_condition_df$protein_col == col]
  
  if (is.na(condition) || condition %in% exclude_conditions) return(NULL)
  
  data_col <- bacillus_data[[col]]
  
  contaminant_quantity <- sum(
    data_col[bacillus_data$protein_type == "Contaminant"],
    na.rm = TRUE
  )
  
  bsubt_quantity <- sum(
    data_col[bacillus_data$protein_type == "B.subtilis"],
    na.rm = TRUE
  )
  
  total_quantity <- contaminant_quantity + bsubt_quantity
  
  tibble(
    protein_col = col,
    base_name = base_name,
    condition = condition,
    contaminant_quantity = contaminant_quantity,
    bsubt_quantity = bsubt_quantity,
    total_quantity = total_quantity
  )
})

sum_data <- sum_data %>%
  mutate(condition = factor(condition, levels = final_condition_levels))

sum_summary <- sum_data %>%
  group_by(condition) %>%
  summarise(
    mean_contaminant = mean(contaminant_quantity, na.rm = TRUE),
    sd_contaminant   = sd(contaminant_quantity, na.rm = TRUE),
    mean_bsubt       = mean(bsubt_quantity, na.rm = TRUE),
    sd_bsubt         = sd(bsubt_quantity, na.rm = TRUE),
    mean_total       = mean(total_quantity, na.rm = TRUE),
    sd_total         = sd(total_quantity, na.rm = TRUE),
    n                = n(),
    .groups = "drop"
  )

sum_summary_long <- sum_summary %>%
  select(condition, mean_contaminant, mean_bsubt) %>%
  pivot_longer(
    cols = c(mean_contaminant, mean_bsubt),
    names_to = "protein_type",
    values_to = "mean_quantity"
  ) %>%
  mutate(
    protein_type = recode(
      protein_type,
      mean_contaminant = "Contaminant",
      mean_bsubt = "B.subtilis"
    ),
    protein_type = factor(protein_type, levels = c("Contaminant", "B.subtilis"))
  )


# plot the barplot
# --- Lower Panel (0 to 1e4) ---
p3_low <- ggplot() +
  geom_bar(data = sum_summary_long, aes(x = condition, y = mean_quantity, fill = protein_type),
           stat = "identity", position = "stack", width = 0.75, alpha = 0.85) +
  geom_errorbar(data = sum_summary,
                aes(x = condition, ymin = mean_total - sd_total, ymax = mean_total + sd_total),
                width = 0.3, linewidth = 0.6) +
  geom_point(data = sum_data, aes(x = condition, y = total_quantity),
             position = position_jitter(width = 0.15, seed = 42),
             shape = 21, size = 3.2, stroke = 1, colour = "black", fill = "white") +
  scale_fill_manual(values = bsubt_colors, name = "Protein type",
                    labels = c("Contaminant", expression(italic("B. subtilis")))) +
  coord_cartesian(ylim = c(0, 1e4)) +
  scale_y_continuous(breaks = c(0, 5000, 10000),
                     labels = scales::label_scientific()) +
  labs(y = NULL) +
  shared_theme() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, family = "Arial"), # Force fixed size
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 0, l = 5))

# --- Upper Panel (1e4 to 3e5) ---
p3_high <- ggplot() +
  geom_bar(data = sum_summary_long, aes(x = condition, y = mean_quantity, fill = protein_type),
           stat = "identity", position = "stack", width = 0.75, alpha = 0.85) +
  geom_errorbar(data = sum_summary,
                aes(x = condition, ymin = mean_total - sd_total, ymax = mean_total + sd_total),
                width = 0.3, linewidth = 0.6) +
  geom_point(data = sum_data, aes(x = condition, y = total_quantity),
             position = position_jitter(width = 0.15, seed = 42),
             shape = 21, size = 3.2, stroke = 1, colour = "black", fill = "white") +
  scale_fill_manual(values = bsubt_colors, name = "Protein type",
                    labels = c("Contaminant", expression(italic("B. subtilis")))) +
  coord_cartesian(ylim = c(1e4, 3e5)) +
  scale_y_continuous(breaks = c(1e5, 2e5, 3e5),
                     labels = scales::label_scientific()) +
  labs(y = NULL) +
  shared_theme() +
  theme(axis.text.y = element_text(size = 12, family = "Arial"), # Force fixed size
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5))

# --- Combine with Table ---
# We use & theme() to apply the same text scaling rules to the whole assembly
p3_combined <- (p3_high / p3_low / make_table_plot()) +
  plot_layout(heights = c(2, 2, 1)) &
  theme(text = element_text(family = "Arial"))

p3_combined

# save the Fig. 3 B barplot
ggsave(file.path(plots_dir, "fig3_B_barplot.svg"),
       plot = p3_combined, device = "svg", width = 9.48, height = 4.68, units = "in")


# Figure 3 C: CV distribution per condition (proteins in >= 3 replicates, B. subtilis only) -------------------------------------------------------------------------

cv_data <- map_dfr(final_condition_levels, function(cond) {
  
  cond_cols <- col_condition_df %>%
    filter(condition == cond) %>%
    pull(protein_col)
  
  # If < 3 replicates exist → return a placeholder (keeps empty box on x‑axis)
  if(length(cond_cols) < 3) {
    return(tibble(
      condition = cond,
      cv = NA_real_
    ))
  }
  
  qty_matrix <- bacillus_data %>%
    select(PG.ProteinGroups, protein_type, all_of(cond_cols))
  
  qty_long <- qty_matrix %>%
    pivot_longer(cols = all_of(cond_cols),
                 names_to = "replicate",
                 values_to = "quantity") %>%
    filter(!is.na(quantity) & !is.nan(quantity) & quantity > 0) %>%
    filter(protein_type == "B.subtilis")
  
  # Identify proteins quantified in >= 3 replicates
  valid_proteins <- qty_long %>%
    group_by(PG.ProteinGroups) %>%
    summarise(n_valid = n(), .groups = "drop") %>%
    filter(n_valid >= 3) %>%
    pull(PG.ProteinGroups)
  
  # If no valid proteins → return empty placeholder
  if(length(valid_proteins) == 0) {
    return(tibble(
      condition = cond,
      cv = NA_real_
    ))
  }
  
  # Compute CV
  qty_long %>%
    filter(PG.ProteinGroups %in% valid_proteins) %>%
    group_by(PG.ProteinGroups) %>%
    summarise(cv = sd(quantity) / mean(quantity) * 100, .groups = "drop") %>%
    mutate(condition = cond)
  
}) %>%
  mutate(condition = factor(condition, levels = final_condition_levels))

# CV summary for printing 
cv_medians <- cv_data %>%
  group_by(condition) %>%
  summarise(median_cv = median(cv, na.rm = TRUE),
            n_proteins = sum(!is.na(cv)),
            .groups = "drop")

# plot the Boxplot
p5 <- ggplot(cv_data, aes(x = condition, y = cv)) +
  geom_boxplot(
    fill = "#56B4E9", alpha = 0.85,
    outlier.size  = 0.6,
    outlier.alpha = 0.4,
    width = 0.65,
    na.rm = TRUE
  ) +
  geom_text(
    data = cv_medians %>% filter(!is.na(median_cv)),
    # Added paste0 to include the % symbol
    aes(x = condition, y = median_cv, label = paste0(round(median_cv, 1), "%")),
    inherit.aes = FALSE,
    family = "Arial",
    size = 3.5,
    color = "black",
    vjust = 0.5
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = function(x) seq(0, ceiling(max(x, na.rm = TRUE) / 25) * 25, by = 25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(y = expression(bold(italic("B. subtilis")~"protein quantity CV (%)"))) +
  shared_theme() +
  theme(legend.position = "none")

combined3 <- p5 / make_table_plot() + plot_layout(heights = c(4, 1))
combined3

# save the Fig. 3 C boxplot
ggsave(file.path(plots_dir, "fig3_C_boxplot.svg"),
       plot = combined3, device = "svg", width = 9.48, height = 4.68, units = "in")

# Figure 3 D: PCA-------------------------------------------------------------------------

pca_conditions <- c("1x_protop_st", "10x_protop_st")

pca_cols <- col_condition_df %>%
  filter(condition %in% pca_conditions) %>%
  pull(protein_col)

if(length(pca_cols) == 0) stop("No columns found for PCA conditions")

pca_matrix <- bacillus_data %>% select(PG.ProteinGroups, all_of(pca_cols))

pca_matrix_clean <- pca_matrix %>%
  filter(if_all(all_of(pca_cols), ~ !is.na(.) & !is.nan(.) & . > 0))

if(nrow(pca_matrix_clean) < 2) stop("Not enough proteins for PCA")

pca_data <- pca_matrix_clean %>%
  select(-PG.ProteinGroups) %>%
  as.matrix() %>% t()

pca_result <- prcomp(log10(pca_data), center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$sample <- rownames(pca_scores)

pca_scores <- pca_scores %>%
  left_join(col_condition_df %>% select(protein_col, condition, base_name),
            by = c("sample" = "protein_col"))

variance_explained <- summary(pca_result)$importance[2, ] * 100

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = condition, label = base_name)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = condition), level = 0.95, linewidth = 1, linetype = "dashed") +
  scale_color_manual(
    values = c("1x_protop_st" = "#56B4E9", "10x_protop_st" = "#8F8989"),
    name = "Condition"
  ) +
  labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(axis.title   = element_text(size = 13, face = "bold", family = "Arial"),
        axis.text    = element_text(size = 12, color = "black", family = "Arial"),
        legend.title = element_text(size = 13, face = "bold", family = "Arial"),
        legend.text  = element_text(size = 12, family = "Arial"),
        axis.line    = element_line(linewidth = 0.8),
        axis.ticks   = element_line(linewidth = 0.8),
        legend.position = "right")

p_pca

# save the Fig. 3 D PCA
ggsave(file.path(plots_dir, "fig3_D_PCA.svg"),
       plot = p_pca, device = "svg", width = 9.48, height = 4.68, units = "in")
