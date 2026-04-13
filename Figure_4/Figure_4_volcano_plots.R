library(tidyverse)
library(svglite)
library(ggrepel)

#  read in SN data -----------------------------------------------------------

data_dir <- file.path("data/volcano")

# Define the subdirectory where plots should be saved
plots_dir <- "plots"

# ── 1. List & categorize files ────────────────────────────────────────────────  
all_files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)  

# volcano files are the files exported from Perseus from the volcano data
volcano_files <- all_files[grepl("volcanodata", all_files)]  

# curve files are the files exported from Perseus from the curves of the respective volcano data
curve_files   <- all_files[grepl("curvedata",   all_files)]

# ── 2. Helper: extract sample ID + condition from filename ────────────────────  
# e.g. "fig4_volcanodata_1pg.txt" → sample = "fig4", condition = "1pg"  
parse_filename <- function(path) {  
  fname <- tools::file_path_sans_ext(basename(path))  
  parts <- strsplit(fname, "_")[[1]]  
  list(  
    sample    = parts[1],                        # e.g. "fig4"  
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


# ── 6. Plot volcanos ────────────────────────────────────────────────────────────────  

# Define genes of interest which should be labeled in the volcano plots
# add genes here if you want to label/colour other genes 
genes_of_interest <- c("clpC", "groES", "groEL")



# Figure 4 F: volcano plot for 1 pg-------------------------------------------------------------------------
# Filter data for 1pg

v_1pg <- volcano_data %>% filter(condition == "1pg") %>%
  mutate(is_goi = PG.Genes %in% genes_of_interest)
c_1pg <- curve_data %>% filter(condition == "1pg")

p_1pg <- ggplot(v_1pg, aes(x = Difference, y = `-Log(P-value)`)) +
  
  # 1. All points (Grey background)
  geom_point(color = "grey", size = 1.2, alpha = 0.7) +
  
  # 2. Genes of interest (Blue, on top)
  geom_point(data = v_1pg %>% filter(is_goi),
             color = "#0072B2", size = 1.4, alpha = 0.9) +
  
  # 3. Perseus Curve
  geom_line(data = c_1pg, aes(x = x, y = y),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # 4. Labels ONLY for genes of interest
  geom_text_repel(
    data = v_1pg %>% filter(is_goi),
    aes(label = PG.Genes),
    color = "#0072B2",
    size = 2.5,
    fontface = "bold",
    segment.color = "#0072B2",
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
  
  # Formatting & Styling
  coord_cartesian(ylim = c(0, 10), xlim = c(-10, 8)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  scale_x_continuous(breaks = seq(-10, 8, 2)) +
  labs(x = expression("Difference (50 °C 1 pg - 37 °C 1 pg)"),
       y = expression("-log"[10]*" p-value")) +
  annotate("text", x = Inf, y = Inf, label = "1 pg\n 50 °C - 37 °C",  
           color = "black", size = 2.5,  
           hjust = 1.1, vjust = 1.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

p_1pg

ggsave(file.path(plots_dir, "fig4_volcano_1pg.svg"),
       plot = p_1pg, device = "svg", width = 12.7, height = 8, units = "cm")


# Figure 4 E: volcano plot for 5 pg-------------------------------------------------------------------------

# Filter data for 5pg

v_5pg <- volcano_data %>% filter(condition == "5pg") %>%
  mutate(is_goi = PG.Genes %in% genes_of_interest)
c_5pg <- curve_data %>% filter(condition == "5pg")

p_5pg <- ggplot(v_5pg, aes(x = Difference, y = `-Log(P-value)`)) +
  
  # 1. All points (Grey background)
  geom_point(color = "grey", size = 1.2, alpha = 0.7) +
  
  # 2. Genes of interest (Blue, on top)
  geom_point(data = v_5pg %>% filter(is_goi),
             color = "#0072B2", size = 1.4, alpha = 0.9) +
  
  # 3. Perseus Curve
  geom_line(data = c_5pg, aes(x = x, y = y),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # 4. Labels ONLY for genes of interest
  geom_text_repel(
    data = v_5pg %>% filter(is_goi),
    aes(label = PG.Genes),
    color = "#0072B2",
    size = 2.5,
    fontface = "bold",
    segment.color = "#0072B2",
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
  
  # Formatting & Styling
  coord_cartesian(ylim = c(0, 10), xlim = c(-10, 8)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  scale_x_continuous(breaks = seq(-10, 8, 2)) +
  labs(x = expression("Difference (50 °C 5 pg - 37 °C 5 pg)"),
       y = expression("-log"[10]*" p-value")) +
  annotate("text", x = Inf, y = Inf, label = "5 pg\n 50 °C - 37 °C",  
           color = "black", size = 2.5,  
           hjust = 1.1, vjust = 1.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

p_5pg

ggsave(file.path(plots_dir, "fig4_volcano_5pg.svg"),
       plot = p_5pg, device = "svg", width = 12.7, height = 8, units = "cm")


# Figure 4 D: volcano plot for 10 pg-------------------------------------------------------------------------
# Filter data for 10pg

v_10pg <- volcano_data %>% filter(condition == "10pg") %>%
  mutate(is_goi = PG.Genes %in% genes_of_interest)
c_10pg <- curve_data %>% filter(condition == "10pg")

p_10pg <- ggplot(v_10pg, aes(x = Difference, y = `-Log(P-value)`)) +
  
  # 1. All points (Grey background)
  geom_point(color = "grey", size = 1.2, alpha = 0.7) +
  
  # 2. Genes of interest (Blue, on top)
  geom_point(data = v_10pg %>% filter(is_goi),
             color = "#0072B2", size = 1.4, alpha = 0.9) +
  
  # 3. Perseus Curve
  geom_line(data = c_10pg, aes(x = x, y = y),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # 4. Labels ONLY for genes of interest
  geom_text_repel(
    data = v_10pg %>% filter(is_goi),
    aes(label = PG.Genes),
    color = "#0072B2",
    size = 2.5,
    fontface = "bold",
    segment.color = "#0072B2",
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
  
  # Formatting & Styling
  coord_cartesian(ylim = c(0, 10), xlim = c(-10, 8)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  scale_x_continuous(breaks = seq(-10, 8, 2)) +
  labs(x = expression("Difference (50 °C 10 pg - 37 °C 10 pg)"),
       y = expression("-log"[10]*" p-value")) +
  annotate("text", x = Inf, y = Inf, label = "10 pg\n 50 °C - 37 °C",  
           color = "black", size = 2.5,  
           hjust = 1.1, vjust = 1.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

p_10pg

ggsave(file.path(plots_dir, "fig4_volcano_10pg.svg"),
       plot = p_10pg, device = "svg", width = 12.7, height = 8, units = "cm")

# Figure 4 C: volcano plot for 50 pg-------------------------------------------------------------------------
# Filter data for 50pg

v_50pg <- volcano_data %>% filter(condition == "50pg") %>%
  mutate(is_goi = PG.Genes %in% genes_of_interest)
c_50pg <- curve_data %>% filter(condition == "50pg")

p_50pg <- ggplot(v_50pg, aes(x = Difference, y = `-Log(P-value)`)) +
  
  # 1. All points (Grey background)
  geom_point(color = "grey", size = 1.2, alpha = 0.7) +
  
  # 2. Genes of interest (Blue, on top)
  geom_point(data = v_50pg %>% filter(is_goi),
             color = "#0072B2", size = 1.4, alpha = 0.9) +
  
  # 3. Perseus Curve
  geom_line(data = c_50pg, aes(x = x, y = y),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # 4. Labels ONLY for genes of interest
  geom_text_repel(
    data = v_50pg %>% filter(is_goi),
    aes(label = PG.Genes),
    color = "#0072B2",
    size = 2.5,
    fontface = "bold",
    segment.color = "#0072B2",
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
  
  # Formatting & Styling
  coord_cartesian(ylim = c(0, 10), xlim = c(-10, 8)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  scale_x_continuous(breaks = seq(-10, 8, 2)) +
  labs(x = expression("Difference (50 °C 50 pg - 37 °C 50 pg)"),
       y = expression("-log"[10]*" p-value")) +
  annotate("text", x = Inf, y = Inf, label = "50 pg\n 50 °C - 37 °C",  
           color = "black", size = 2.5,  
           hjust = 1.1, vjust = 1.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

p_50pg

ggsave(file.path(plots_dir, "fig4_volcano_50pg.svg"),
       plot = p_50pg, device = "svg", width = 12.7, height = 8, units = "cm")


# Figure 4 B: volcano plot for 200 pg-------------------------------------------------------------------------
# Filter data for 200pg

v_200pg <- volcano_data %>% filter(condition == "200pg") %>%
  mutate(is_goi = PG.Genes %in% genes_of_interest)
c_200pg <- curve_data %>% filter(condition == "200pg")

p_200pg <- ggplot(v_200pg, aes(x = Difference, y = `-Log(P-value)`)) +
  
  # 1. All points (Grey background)
  geom_point(color = "grey", size = 1.2, alpha = 0.7) +
  
  # 2. Genes of interest (Blue, on top)
  geom_point(data = v_200pg %>% filter(is_goi),
             color = "#0072B2", size = 1.4, alpha = 0.9) +
  
  # 3. Perseus Curve
  geom_line(data = c_200pg, aes(x = x, y = y),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # 4. Labels ONLY for genes of interest
  geom_text_repel(
    data = v_200pg %>% filter(is_goi),
    aes(label = PG.Genes),
    color = "#0072B2",
    size = 2.5,
    fontface = "bold",
    segment.color = "#0072B2",
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
  
  # Formatting & Styling
  coord_cartesian(ylim = c(0, 10), xlim = c(-10, 8)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  scale_x_continuous(breaks = seq(-10, 8, 2)) +
  labs(x = expression("Difference (50 °C 200 pg - 37 °C 200 pg)"),
       y = expression("-log"[10]*" p-value")) +
  annotate("text", x = Inf, y = Inf, label = "200 pg\n 50 °C - 37 °C",  
           color = "black", size = 2.5,  
           hjust = 1.1, vjust = 1.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

p_200pg

ggsave(file.path(plots_dir, "fig4_volcano_200pg.svg"),
       plot = p_200pg, device = "svg", width = 12.7, height = 8, units = "cm")


# Figure 4 A: volcano plot for 1000 pg-------------------------------------------------------------------------
# Filter data for 1000pg

v_1000pg <- volcano_data %>% filter(condition == "1000pg") %>%
  mutate(is_goi = PG.Genes %in% genes_of_interest)
c_1000pg <- curve_data %>% filter(condition == "1000pg")

p_1000pg <- ggplot(v_1000pg, aes(x = Difference, y = `-Log(P-value)`)) +
  
  # 1. All points (Grey background)
  geom_point(color = "grey", size = 1.2, alpha = 0.7) +
  
  # 2. Genes of interest (Blue, on top)
  geom_point(data = v_1000pg %>% filter(is_goi),
             color = "#0072B2", size = 1.4, alpha = 0.9) +
  
  # 3. Perseus Curve
  geom_line(data = c_1000pg, aes(x = x, y = y),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # 4. Labels ONLY for genes of interest
  geom_text_repel(
    data = v_1000pg %>% filter(is_goi),
    aes(label = PG.Genes),
    color = "#0072B2",
    size = 2.5,
    fontface = "bold",
    segment.color = "#0072B2",
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
  
  # Formatting & Styling
  coord_cartesian(ylim = c(0, 10), xlim = c(-10, 8)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  scale_x_continuous(breaks = seq(-10, 8, 2)) +
  labs(x = expression("Difference (50 °C 1000 pg - 37 °C 1000 pg)"),
       y = expression("-log"[10]*" p-value")) +
  annotate("text", x = Inf, y = Inf, label = "1000 pg\n 50 °C - 37 °C",  
           color = "black", size = 2.5,  
           hjust = 1.1, vjust = 1.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

p_1000pg

ggsave(file.path(plots_dir, "fig4_volcano_1000pg.svg"),
       plot = p_1000pg, device = "svg", width = 12.7, height = 8, units = "cm")

