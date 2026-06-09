library(tidyverse)
library(patchwork)

# === USER VARIABLES ===
input_file <- "ancestry_tracts.txt"
chrom_file <- "chrom_lengths.txt"
output_dir <- "ancestry_plots_freq"

# Resolution of the heatmap and data output in base pairs. 
# Lowered to 1000 (1kb) for finer-scale mapping. 
step_size <- 1000 

# The frequency threshold to define a "high frequency" block for gene mapping
freq_threshold <- 0.1

dir.create(output_dir, showWarnings = FALSE)

# === READ DATA AND RENAME OmyY TO Omy29 ===
chrom_lengths <- read_tsv(chrom_file, col_types = cols()) %>%
  rename(Chr_length = length) %>%
  mutate(LG = str_replace(LG, "OmyY", "Omy29"))

blocks <- read_tsv(input_file, col_types = cols()) %>%
  mutate(Ancestry = str_trim(Ancestry)) %>%
  mutate(LG = str_replace(LG, "OmyY", "Omy29"))

# === DEFINE INTERSPECIFIC INVERSIONS ===
inversions <- tibble(
  LG = c("Omy05", "Omy17", "Omy20", "Omy22", "Omy29"),
  Start_Mb = c(30.102, 36.050, 4.356, 7.308, 2.314),
  End_Mb   = c(87.502, 43.906, 21.129, 22.170, 14.206)
)

# === SORT CHROMOSOMES NUMERICALLY ===
lg_levels <- str_sort(unique(chrom_lengths$LG), numeric = TRUE)

chrom_lengths <- chrom_lengths %>%
  mutate(LG = factor(LG, levels = lg_levels))

inversions <- inversions %>%
  mutate(LG = factor(LG, levels = lg_levels))

# === SORT HETEROZYGOTE COMBINATIONS ALPHABETICALLY ===
# This ensures WCT/IRT becomes IRT/WCT
blocks <- blocks %>%
  mutate(Ancestry = ifelse(grepl("/", Ancestry),
                           sapply(strsplit(Ancestry, "/"), function(x) paste(sort(x), collapse="/")),
                           Ancestry))

# === COLLAPSE ANCESTRY CATEGORIES ===
blocks <- blocks %>%
  mutate(Ancestry = case_when(
    Ancestry %in% c("CRT_AA", "CRT_RR", "CRT_AA/CRT_RR") ~ "CRT",
    Ancestry %in% c("CRT_AA/IRT", "CRT_RR/IRT") ~ "CRT/IRT",
    Ancestry %in% c("CRT_AA/WCT", "CRT_RR/WCT") ~ "CRT/WCT",
    Ancestry %in% c("CRT_AA/YCT", "CRT_RR/YCT") ~ "CRT/YCT",
    TRUE ~ Ancestry
  ))

# =========================================================================
# === PREPARE GENOMIC SAMPLING POINTS ===
# =========================================================================
# Create a sequence of points every 'step_size' along each chromosome
sample_points <- chrom_lengths %>%
  mutate(Position = map(Chr_length, ~ seq(0, .x, by = step_size))) %>%
  unnest(Position) %>%
  dplyr::select(LG, Position)

# =========================================================================
# === FUNCTION: GENERATE HEATMAP PLOT OBJECT & EXTRACT BREAKPOINTS ===
# =========================================================================
plot_frequency_heatmap <- function(target_individuals, output_name, target_ancestries = c("WCT", "IRT/WCT"), plot_tag = NULL, threshold = freq_threshold) {
  
  # 1. Filter population
  pop_blocks <- blocks %>% filter(Individual %in% target_individuals)
  n_inds <- length(target_individuals)
  
  if (nrow(pop_blocks) == 0) {
    warning(paste("Skipping", output_name, "- no matching individuals found.\n"))
    return(NULL)
  }
  
  # 2. Isolate only the blocks matching our target ancestries
  target_blocks <- pop_blocks %>%
    filter(Ancestry %in% target_ancestries)
  
  # 3. Calculate frequency at each sampling point
  freq_data <- sample_points %>%
    mutate(
      Count = map2_dbl(LG, Position, function(chr, pos) {
        sum(target_blocks$LG == chr & target_blocks$Start_bp <= pos & target_blocks$End_bp >= pos)
      }),
      Frequency = Count / n_inds
    )
  
  # =========================================================================
  # === NEW: COLLAPSE POINTS INTO CONTIGUOUS START/STOP BLOCKS ===
  # =========================================================================
  high_freq_blocks <- freq_data %>%
    filter(Frequency >= threshold) %>%
    arrange(LG, Position) %>%
    group_by(LG) %>%
    # Create a grouping variable that increments whenever there is a gap 
    # larger than step_size, meaning it's a new continuous block
    mutate(block_id = cumsum(c(1, diff(Position) > step_size))) %>%
    group_by(LG, block_id) %>%
    summarize(
      Start_bp = min(Position),
      End_bp = max(Position) + step_size, # Add step_size so the block covers the whole final bin
      Mean_Frequency = mean(Frequency),
      Max_Frequency = max(Frequency),
      .groups = "drop"
    ) %>%
    mutate(Population = output_name) %>%
    dplyr::select(Population, LG, Start_bp, End_bp, Mean_Frequency, Max_Frequency)
  
  # Save the BED-like block data to a TSV file
  safe_name <- gsub(" ", "_", output_name)
  data_output_path <- file.path(output_dir, paste0(safe_name, "_HighFreq_Breakpoints.txt"))
  write_tsv(high_freq_blocks, data_output_path)
  cat("Saved contiguous breakpoint data for", output_name, "to:", data_output_path, "\n")
  # =========================================================================
  
  # Prepare for plotting by converting 0s to NA so they don't plot
  plot_data <- freq_data %>%
    mutate(Frequency = ifelse(Frequency == 0, NA, Frequency))
  
  # 4. Create Plot
  p <- ggplot() +
    # Background solid line for chromosomes
    geom_segment(data = chrom_lengths, aes(
      x = 0, xend = Chr_length/1e6,
      y = LG, yend = LG
    ), color = "gray85", linewidth = 3, lineend = "round") +
    
    # Heatmap tiles mapping frequency to color
    geom_tile(data = plot_data, aes(
      x = Position/1e6, 
      y = LG, 
      fill = Frequency
    ), width = step_size/1e6, height = 0.62) +
    
    # Inversion boundaries (empty boxes overlayed on top) - Reduced linewidth to 0.3
    geom_rect(data = inversions, aes(
      xmin = Start_Mb, xmax = End_Mb,
      ymin = as.numeric(LG) - 0.45, ymax = as.numeric(LG) + 0.45,
      color = "Inversion"
    ), fill = NA, linewidth = 0.3) +
    
    scale_y_discrete() + 
    
    # Colorblind-friendly Blue to Orange gradient in 5% increments
    scale_fill_gradient(
      low = "#4575b4",               # Deep Blue
      high = "#f46d43",              # Vibrant Orange
      limits = c(0, 0.4),            
      breaks = seq(0, 1, by = 0.05), # 5% bins
      na.value = "transparent", 
      labels = scales::percent_format(accuracy = 1), 
      name = "WCT\nAncestry"
    ) +
    scale_color_manual(name = NULL, values = c("Inversion" = "black")) +
    
    # Apply dynamic title WITH italicized N, and removed subtitle
    labs(
      title = bquote(paste(.(output_name), " Consensus (", italic(N) == .(n_inds), ")")),
      subtitle = NULL,
      tag = plot_tag,
      x = "Position (Mb)", y = "Chromosome"
    ) +
    
    # Molecular Ecology Publication Theme
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, color = "black", hjust = 0.5),
      plot.tag = element_text(size = 16, face = "bold"),
      plot.tag.position = "topleft",
      legend.position = "right", 
      legend.box = "vertical"
    )
  
  return(p)
}

# =========================================================================
# === RUN PLOTS & GENERATE MULTI-PANEL FIGURE ===
# =========================================================================

# 1. Read the explicit sample lists
bear_samples <- readLines("upper_bear_samples.txt")
efyaak_samples <- readLines("efyaak_samples.txt")

# 2. Generate the heatmap plot objects focusing ONLY on WCT and IRT/WCT
p_yaak_heat <- plot_frequency_heatmap(
  target_individuals = efyaak_samples,
  output_name = "EF Yaak",
  target_ancestries = c("WCT", "IRT/WCT"),
  plot_tag = "a."
)

p_bear_upper_heat <- plot_frequency_heatmap(
  target_individuals = bear_samples,
  output_name = "Upper Bear",
  target_ancestries = c("WCT", "IRT/WCT"),
  plot_tag = "b."
)

# 3. Save INDIVIDUAL plots as PNGs
if(!is.null(p_yaak_heat)) {
  yaak_png_path <- file.path(output_dir, "EF_Yaak_WCT_Heatmap.png")
  ggsave(filename = yaak_png_path, plot = p_yaak_heat, width = 12, height = 6, dpi = 300)
  cat("Saved individual PNG:", yaak_png_path, "\n")
}

if(!is.null(p_bear_upper_heat)) {
  bear_png_path <- file.path(output_dir, "Upper_Bear_WCT_Heatmap.png")
  ggsave(filename = bear_png_path, plot = p_bear_upper_heat, width = 12, height = 6, dpi = 300)
  cat("Saved individual PNG:", bear_png_path, "\n")
}

# 4. Combine EF Yaak and Upper Bear into a 2-panel figure using patchwork & save as PDF
if(!is.null(p_yaak_heat) & !is.null(p_bear_upper_heat)) {
  combined_figure <- (p_yaak_heat / p_bear_upper_heat) + 
    plot_layout(guides = "collect")
  
  # Save the combined figure as a PDF
  combined_pdf_path <- file.path(output_dir, "Figure_WCT_Heatmap_Combined.pdf")
  ggsave(filename = combined_pdf_path, plot = combined_figure, width = 12, height = 12, dpi = 300)
  cat("Saved combined heatmap figure to:", combined_pdf_path, "\n")
}