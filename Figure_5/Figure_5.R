library(tidyverse)
library(ggrepel)

#  read in SN data -----------------------------------------------------------

data_dir <- file.path("data/volcano/")

# Define the subdirectory where plots should be saved
plots_dir <- "plots"  

# ── 1. List & categorize files ────────────────────────────────────────────────  
all_files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)  

# volcano files are the files exported from Perseus from the volcano data
volcano_files <- all_files[grepl("volcanodata", all_files)]  

# curve files are the files exported from Perseus from the curves of the respective volcano data
curve_files   <- all_files[grepl("curvedata",   all_files)]

# ── 2. Helper: extract sample ID + condition from filename ────────────────────  
# e.g. "fig5_volcanodata_1pg.txt" → sample = "fig5", condition = "1pg"  
parse_filename <- function(path) {  
  fname <- tools::file_path_sans_ext(basename(path))  
  parts <- strsplit(fname, "_")[[1]]  
  list(  
    sample    = parts[1],                        # e.g. "fig5"  
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
# if you want to label more genes, add them here
genes_of_interest <- c("clpC", "groES", "groEL")




# Figure 5 B: dMcsB -------------------------------------------------------------------------

# Filter data for dMcsB
v_dMcsB <- volcano_data %>% filter(condition == "dMcsB") %>%
  mutate(is_goi = PG.Genes %in% genes_of_interest)
c_dMcsB <- curve_data %>% filter(condition == "dMcsB")

p_dMcsB <- ggplot(v_dMcsB, aes(x = Difference, y = `-Log(P-value)`)) +
  
  # 1. All points (Grey background)
  geom_point(color = "grey", size = 1.5, alpha = 0.7) +
  
  # 2. Genes of interest (Blue, on top)
  geom_point(data = v_dMcsB %>% filter(is_goi),
             color = "#0072B2", size = 1.7, alpha = 0.9) +
  
  # 3. Perseus Curve
  geom_line(data = c_dMcsB, aes(x = x, y = y),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # 4. Labels ONLY for genes of interest
  geom_text_repel(
    data = v_dMcsB %>% filter(is_goi),
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
  coord_cartesian(ylim = c(0, 8), xlim = c(-10, 10)) +
  scale_y_continuous(breaks = seq(0, 8, 2)) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  labs(x = expression("Difference (50 °C " * Delta * "mcsB - 37 °C " * Delta * "mcsB)"),
       y = expression("-log"[10]*" p-value")) +
  annotate("text", x = -9.5, y = 6.7,
           label = paste0("\u0394", "mcsB\n 50 °C - 37 °C"),
           color = "black", size = 3.5,
           hjust = 0, vjust = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

p_dMcsB

ggsave(file.path(plots_dir, "fig5_volcano_p_dMcsB.svg"),
       plot = p_dMcsB, device = "svg", width = 12.7, height = 8, units = "cm")



# Figure 5 C: dYwIE -------------------------------------------------------------------------

# Filter data for dYwIE
v_dYwIE <- volcano_data %>% filter(condition == "dYwIE") %>%
  mutate(is_goi = PG.Genes %in% genes_of_interest)
c_dYwIE <- curve_data %>% filter(condition == "dYwIE")

p_dYwIE <- ggplot(v_dYwIE, aes(x = Difference, y = `-Log(P-value)`)) +
  
  # 1. All points (Grey background)
  geom_point(color = "grey", size = 1.5, alpha = 0.7) +
  
  # 2. Genes of interest (Blue, on top)
  geom_point(data = v_dYwIE %>% filter(is_goi),
             color = "#0072B2", size = 1.7, alpha = 0.9) +
  
  # 3. Perseus Curve
  geom_line(data = c_dYwIE, aes(x = x, y = y),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # 4. Labels ONLY for genes of interest
  geom_text_repel(
    data = v_dYwIE %>% filter(is_goi),
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
  coord_cartesian(ylim = c(0, 8), xlim = c(-10, 10)) +
  scale_y_continuous(breaks = seq(0, 8, 2)) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  labs(x = expression("Difference (50 °C " * Delta * "ywIE - 37 °C " * Delta * "ywIE)"),
       y = expression("-log"[10]*" p-value")) +
  annotate("text", x = -9.5, y = 6.7,
           label = paste0("\u0394", "ywIE\n 50 °C - 37 °C"),
           color = "black", size = 3.5,
           hjust = 0, vjust = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

p_dYwIE

ggsave(file.path(plots_dir, "fig5_volcano_p_dYwIE.svg"),
       plot = p_dYwIE, device = "svg", width = 12.7, height = 8, units = "cm")


# Figure 5 A: wt168 -------------------------------------------------------------------------

# Filter data for wt168

v_wt168 <- volcano_data %>% filter(condition == "wt168") %>%
  mutate(is_goi = PG.Genes %in% genes_of_interest)
c_wt168 <- curve_data %>% filter(condition == "wt168")

p_wt168 <- ggplot(v_wt168, aes(x = Difference, y = `-Log(P-value)`)) +
  
  # 1. All points (Grey background)
  geom_point(color = "grey", size = 1.5, alpha = 0.7) +
  
  # 2. Genes of interest (Blue, on top)
  geom_point(data = v_wt168 %>% filter(is_goi),
             color = "#0072B2", size = 1.7, alpha = 0.9) +
  
  # 3. Perseus Curve
  geom_line(data = c_wt168, aes(x = x, y = y),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # 4. Labels ONLY for genes of interest
  geom_text_repel(
    data = v_wt168 %>% filter(is_goi),
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
  coord_cartesian(ylim = c(0, 8), xlim = c(-10, 10)) +
  scale_y_continuous(breaks = seq(0, 8, 2)) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  labs(x = expression("Difference (50 °C wt168 - 37 °C wt168)"),
       y = expression("-log"[10]*" p-value")) +
  annotate("text", x = -9.5, y = 6.7,
           label = "wt168\n 50 °C - 37 °C",
           color = "black", size = 3.5,
           hjust = 0, vjust = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

p_wt168

ggsave(file.path(plots_dir, "fig5_volcano_p_wt168.svg"),
       plot = p_wt168, device = "svg", width = 12.7, height = 8, units = "cm")



# boxplot of genes of interest  (Plot the data: Figure 5 D) --------------------------------------------

# read in data
data_dir_boxplot <- file.path("data/boxplot/")
boxplot_data <- read.delim(file.path(data_dir_boxplot, "protein_table_fig5.txt"), sep = "\t", header = TRUE)

# 1. Pivot the wide matrix to long format  
protein_data_long <- boxplot_data %>%  
  pivot_longer(  
    cols = contains("PG.Quantity"),  
    names_to = "condition_raw",  
    values_to = "Intensity"  
  )

# 2. Process the data for the boxplot and define genes_of_interest_boxplot for this specific plot
genes_of_interest_boxplot <- c("clpC", "groEL", "groES")  

profile_data_box <- protein_data_long %>%  
  rename(Genes = T..PG.Genes) %>%  
  filter(Genes %in% genes_of_interest_boxplot) %>%  
  mutate(condition = case_when(  
    grepl("dMscB_BS_hs", condition_raw) ~ paste0("\u0394", "mcsB 50°C"),  
    grepl("dMscB_BS_st", condition_raw) ~ paste0("\u0394", "mcsB 37°C"),  
    grepl("dYwIE_BS_hs", condition_raw) ~ paste0("\u0394", "ywiE 50°C"),  
    grepl("dYwIE_BS_st", condition_raw) ~ paste0("\u0394", "ywiE 37°C"),  
    grepl("wt168_BS_hs", condition_raw) ~ "wt168 50°C",  
    grepl("wt168_BS_st", condition_raw) ~ "wt168 37°C",  
    TRUE ~ NA_character_  
  )) %>%  
  filter(!is.na(condition))

# 3. Define order and factorize  
condition_order <- c(
  "wt168 37°C", "wt168 50°C",  
  paste0("\u0394", "mcsB 37°C"), paste0("\u0394", "mcsB 50°C"),  
  paste0("\u0394", "ywiE 37°C"), paste0("\u0394", "ywiE 50°C")
)

profile_data_box <- profile_data_box %>%  
  mutate(condition = factor(condition, levels = condition_order))  



# 4. Plot the data: Figure 5 D 

p_profile_boxplot <- profile_data_box %>%  
  ggplot(aes(x = condition,  
             y = Intensity,  
             fill = Genes,  
             color = Genes)) +  
  
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.6,  
               position = position_dodge(width = 0.75)) +  
  
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75,  
                                              jitter.width = 0.15),  
              size = 1.2, alpha = 0.8) +  
  
  geom_vline(xintercept = c(2.5, 4.5), color = "grey40",   
             linetype = "dashed", linewidth = 0.4) +  
  
  scale_fill_manual(values = c(  
    "clpC"  = "#E69F00",  
    "groEL" = "#009E73",  
    "groES" = "#0072B2"  
  ), name = NULL) +  
  scale_color_manual(values = c(  
    "clpC"  = "#E69F00",  
    "groEL" = "#009E73",  
    "groES" = "#0072B2"  
  ), name = NULL) +  
  
  labs(y = "Relative abundance", x = NULL) +  
  
  theme_minimal() +  
  theme(  
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8, color = "black"),  
    axis.text.y  = element_text(size = 9, color = "black"),  
    axis.title.y = element_text(size = 9, color = "black"),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line  = element_line(color = "black"),  
    axis.ticks = element_line(linewidth = 0.3, color = "black"),  
    legend.text = element_text(size = 8),  
    plot.title = element_blank()  
  )  

p_profile_boxplot  

ggsave(file.path(plots_dir, "p_profile_boxplot_fig5.svg"),  
       plot = p_profile_boxplot,  
       device = "svg",  
       width = 12.7, height = 8, units = "cm")