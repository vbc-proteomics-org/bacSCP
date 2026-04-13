library(tidyverse)
library(extrafont)
library(patchwork)
library(gridExtra)
library(grid)

loadfonts(device = "win", quiet = TRUE)

# define directories, read in data and prepare data ------------------------------------
# Define the subdirectory where plots should be saved
plots_dir <- "plots"  # Subdirectory with plots

# Read input files
coli_data <- read_delim("data/protein_table_fig2.tsv", delim = "\t")
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


# Figure 2 A: Stacked bar plot with replicate dots -------------------------------------------------------------------------

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



# Figure 2 C: Stacked bar chart of total protein quantity (sum) -------------------------------------------------------------------------

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


# Figure 2 zoom in from C: E. coli protein quantity only with reference line -------------------------------------------------------------------------

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


# Figure 2 B: Quantity histogram (NO annotation table) -------------------------------------------------------------------------

coli_ecoli <- coli_data %>% filter(protein_type == "Ecoli")

cols_200pg <- col_condition_df %>% filter(condition == "200pg")        %>% pull(protein_col) %>% unique()
cols_1x    <- col_condition_df %>% filter(condition == "1x_intact_st") %>% pull(protein_col) %>% unique()

if(length(cols_200pg) == 0) stop("No columns found for 200pg condition.")
if(length(cols_1x)    == 0) stop("No columns found for 1x_intact_st condition.")

quant_200pg_long <- coli_ecoli %>%
  select(PG.ProteinGroups, all_of(cols_200pg)) %>%
  pivot_longer(-PG.ProteinGroups, names_to = "protein_col", values_to = "quantity") %>%
  filter(!is.na(quantity) & !is.nan(quantity) & quantity > 0) %>%
  mutate(source = "200pg_all")

quant_1x_long <- coli_ecoli %>%
  select(PG.ProteinGroups, all_of(cols_1x)) %>%
  pivot_longer(-PG.ProteinGroups, names_to = "protein_col", values_to = "quantity") %>%
  filter(!is.na(quantity) & !is.nan(quantity) & quantity > 0)

proteins_in_1x <- unique(quant_1x_long$PG.ProteinGroups)

matched_200pg <- quant_200pg_long %>%
  filter(PG.ProteinGroups %in% proteins_in_1x) %>%
  mutate(source = "matched_from_1x")

if(length(unique(matched_200pg$PG.ProteinGroups)) == 0)
  warning("No proteins quantified in 1x_intact_st were found in 200pg.")

plot_df <- bind_rows(
  quant_200pg_long %>% select(PG.ProteinGroups, quantity, source),
  matched_200pg    %>% select(PG.ProteinGroups, quantity, source)
) %>%
  mutate(source = factor(source,
                         levels = c("200pg_all", "matched_from_1x"),
                         labels = c("200pg (all non contaminant proteins)", "Matched proteins (from 1x_intact_st)")
  ))

p_hist <- ggplot(plot_df, aes(x = log10(quantity), fill = source)) +
  geom_histogram(
    data = subset(plot_df, source == "200pg (all non contaminant proteins)"),
    bins = 50, alpha = 0.6, color = "black", boundary = 0, closed = "left"
  ) +
  geom_histogram(
    data = subset(plot_df, source == "Matched proteins (from 1x_intact_st)"),
    bins = 50, alpha = 0.9, color = "black", boundary = 0, closed = "left"
  ) +
  scale_fill_manual(values = c(
    "200pg (all non contaminant proteins)" = "#808080",
    "Matched proteins (from 1x_intact_st)" = "#148509"
  )) +
  labs(
    x = expression(bold(log[10]~"(Protein quantity at 200pg)")),
    y = "Count", fill = ""
  ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    axis.title  = element_text(size = 14, face = "bold", family = "Arial"),
    axis.text   = element_text(size = 12, color = "black", family = "Arial"),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.position = "top"
  )

p_hist

# save the Fig. 2 B histogram
ggsave(file.path(plots_dir, "fig2_B_hist.svg"),
       plot = p_hist, device = "svg", width = 9.48, height = 4.68, units = "in")