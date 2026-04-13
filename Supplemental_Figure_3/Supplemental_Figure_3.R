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

# HELPER FUNCTION: histogram plot
make_quantity_hist <- function(bulk_long,
                               bacscp_long,
                               bulk_condition,
                               bacscp_condition,
                               highlight_color,
                               bulk_label,
                               bacscp_label
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
    mutate(source = "matched")
  
  if (nrow(matched_bulk) == 0)
    warning(paste("No matched proteins found for BacSCP condition:", bacscp_condition))
  
  plot_df <- bind_rows(
    quant_bulk   %>% select(PG.ProteinGroups, quantity, source),
    matched_bulk %>% select(PG.ProteinGroups, quantity, source)
  ) %>%
    mutate(source = factor(source,
                           levels = c("bulk_all", "matched"),
                           labels = c(bulk_label, bacscp_label)))
  
  ggplot(plot_df, aes(x = log10(quantity), fill = source)) +
    geom_histogram(
      data  = subset(plot_df, source == bulk_label),
      bins  = 50, alpha = 0.6, color = "black", boundary = 0, closed = "left"
    ) +
    geom_histogram(
      data  = subset(plot_df, source == bacscp_label),
      bins  = 50, alpha = 0.9, color = "black", boundary = 0, closed = "left"
    ) +
    scale_fill_manual(values = setNames(c("#808080", highlight_color),
                                        c(bulk_label, bacscp_label))) +
    labs(
      x     = expression(bold(log[10]~"(Protein quantity)")),
      y     = "Count",
      fill  = ""
    ) +
    theme_classic(base_size = 14, base_family = "Arial") +
    theme(
      axis.title      = element_text(size = 14, face = "bold", family = "Arial"),
      axis.text       = element_text(size = 12, color = "black", family = "Arial"),
      legend.text     = element_text(size = 12, family = "Arial"),
      legend.position = "top"
    )
}

# Supplemental Figure 3 B: 50°C — WITHOUT contaminants -------------------------------------------------------------------------

p1 <- make_quantity_hist(
  bulk_long        = bulk200ng_long_nocont,
  bacscp_long      = bacscp_long_nocont,
  bulk_condition   = "200ng dMscB 50°C",
  bacscp_condition = "1x_BS_HS",
  highlight_color  = "#7B2D8B",
  bulk_label       = "200ng dMcsB 50 °C (contaminants removed)",
  bacscp_label     = "Matched proteins (single cell dMcsB 50 °C, contaminants removed, n=20)"
)
p1

# save the histogram
ggsave(file.path(plots_dir, "sup_fig3_B_hist_50°C.svg"),
       plot = p1, device = "svg", width = 10.5, height = 5, units = "in")


# Supplemental Figure 3 A: 37 °C — WITHOUT contaminants -------------------------------------------------------------------------

p2 <- make_quantity_hist(
  bulk_long        = bulk200ng_long_nocont,
  bacscp_long      = bacscp_long_nocont,
  bulk_condition   = "200ng dMscB 37°C",
  bacscp_condition = "1x_BS_chil",
  highlight_color  = "#E69F00",
  bulk_label       = "200ng dMcsB 37 °C (contaminants removed)",
  bacscp_label     = "Matched proteins (single cell dMcsB 37 °C, contaminants removed, n=20)"
)
p2

# save the histogram
ggsave(file.path(plots_dir, "sup_fig3_A_hist_37°C.svg"),
       plot = p2, device = "svg", width = 11, height = 5.5, units = "in")