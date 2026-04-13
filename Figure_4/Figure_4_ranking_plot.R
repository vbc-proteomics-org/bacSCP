library(tidyverse)
library(ggrepel)
library(extrafont)

loadfonts(device = "win", quiet = TRUE)

# Define the subdirectory where plots should be safed
plots_dir <- "plots"  # Subdirectory with plots

# Read input files
PG_data <- read_delim("data/rank/protein_table_fig4.tsv", delim = "\t")
file_mapping <- read_delim("data/rank/file_mapping_fig4.tsv", delim = "\t")

# ── 1. Map PG.Quantity columns to Conditions ──────────────────────────────────
pg_qty_cols <- grep("PG\\.Quantity$", colnames(PG_data), value = TRUE)

col_to_cond <- tibble(col_name = pg_qty_cols) %>%
  rowwise() %>%
  mutate(Condition = {
    match_idx <- which(sapply(file_mapping$FileName, function(x) grepl(x, col_name, fixed = TRUE)))
    if(length(match_idx) > 0) file_mapping$Condition[match_idx[1]] else NA_character_
  }) %>%
  ungroup() %>%
  filter(!is.na(Condition))

# ── 2. Process Data: Log2 Transform, Average, and Rank ────────────────────────
pg_ranked <- PG_data %>%
  select(PG.Genes, all_of(col_to_cond$col_name)) %>%
  pivot_longer(-PG.Genes, names_to = "col_name", values_to = "Quantity") %>%
  left_join(col_to_cond, by = "col_name") %>%
  filter(!is.na(Quantity), Quantity > 0) %>%
  group_by(PG.Genes, Condition) %>%
  summarise(mean_log2_qty = mean(log2(Quantity), na.rm = TRUE), .groups = "drop") %>%
  group_by(Condition) %>%
  arrange(desc(mean_log2_qty)) %>%
  mutate(Rank = row_number()) %>%
  ungroup() %>%
  # Extract the numeric value from the condition name (after the first _)
  mutate(InputAmount = as.numeric(str_extract(Condition, "(?<=_)\\d+"))) %>%
  mutate(InputAmount = factor(InputAmount, levels = sort(unique(InputAmount))))

# ── 3. Strict Labeling Logic ──────────────────────────────────────────────────
target_genes <- c("clpC", "groEL", "groES")

make_labels <- function(data) {
  data %>%
    group_by(Condition) %>%
    filter(
      Rank == 1 |
        sapply(strsplit(as.character(PG.Genes), ";"), function(x) any(target_genes %in% x))
    ) %>%
    group_by(Condition, PG.Genes) %>%
    slice_max(mean_log2_qty, n = 1) %>%
    ungroup()
}

# ── 4. Split into two datasets ────────────────────────────────────────────────
pg_50  <- pg_ranked %>% filter(str_detect(Condition, regex("^dY_", ignore_case = TRUE)))
pg_37  <- pg_ranked %>% filter(str_starts(Condition, "37_"))

labels_50 <- make_labels(pg_50)
labels_37 <- make_labels(pg_37)

# ── 5. Shared plot function ───────────────────────────────────────────────────
make_rank_plot <- function(data, labels, title_str) {
  
  # Get the sorted unique input amounts for the legend breaks
  amounts <- sort(unique(data$InputAmount))
  
  ggplot(data, aes(x = Rank, y = mean_log2_qty, color = InputAmount, group = Condition)) +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    geom_point(data = labels, size = 2.5, shape = 21, fill = "white", stroke = 1.2,
               aes(color = InputAmount)) +
    geom_label_repel(
      data = labels,
      aes(label = PG.Genes),
      size = 2.8,
      fontface = "bold",
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey30",
      show.legend = FALSE,
      max.overlaps = 15
    ) +
    # Blue gradient: light blue for low amounts, dark blue for high amounts
       scale_color_manual(
      values = setNames(
        colorRampPalette(c("#c6dbef", "#08306b"))(length(unique(data$InputAmount))),
        levels(data$InputAmount)
      ),
      name   = "Input Amount",
      labels = paste0(levels(data$InputAmount), " pg")
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    labs(
      title    = title_str,
      x        = "Rank",
      y        = "Average log2(Quantity)"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major = element_blank(),
      legend.position   = "right",
      legend.text       = element_text(size = 9),
      legend.title      = element_text(face = "bold", size = 9),
      plot.title        = element_text(face = "bold"),
      plot.subtitle     = element_text(size = 8, color = "grey40")
    )
}

# ── 6. Generate and save both plots ──────────────────────────────────────────
p_50 <- make_rank_plot(pg_50, labels_50, "Protein Rank Plot – 50 °C conditions")
p_37 <- make_rank_plot(pg_37, labels_37, "Protein Rank Plot – 37 °C conditions")

p_50
p_37

ggsave(file.path(plots_dir, "fig4_G_50degC.svg"),
       plot = p_50, device = "svg", width = 16, height = 8, units = "cm")

ggsave(file.path(plots_dir, "fig4_H_37degC.svg"),
       plot = p_37, device = "svg", width = 16, height = 8, units = "cm")
