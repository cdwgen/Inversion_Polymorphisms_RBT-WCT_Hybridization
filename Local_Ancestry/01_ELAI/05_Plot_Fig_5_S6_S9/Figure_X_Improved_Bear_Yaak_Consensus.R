library(tidyverse)
library(patchwork) # Added for multi-panel plot combination

# === USER VARIABLES ===
input_file <- "ancestry_tracts.txt"
chrom_file <- "chrom_lengths.txt"
output_dir <- "ancestry_plots"

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

# === DEFINE FACTORS AND COLORS ===
ancestry_levels <- c(
  "Uncertain", "CRT", "IRT", "WCT", "YCT",
  "CRT/IRT", "CRT/WCT", "CRT/YCT",
  "IRT/WCT", "IRT/YCT", "WCT/YCT"
)

ancestry_colors <- c(
  "Uncertain" = "gray70",
  "CRT"       = "#66c2a5",
  "IRT"       = "#fc8d62",
  "WCT"       = "#8da0cb",
  "YCT"       = "#e5c494",
  "CRT/IRT"   = "#1b9e77",
  "CRT/WCT"   = "#e78ac3",
  "CRT/YCT"   = "#ffd92f",
  "IRT/WCT"   = "#a6d854",
  "IRT/YCT"   = "#e7298a",
  "WCT/YCT"   = "#d95f02"
)

blocks <- blocks %>%
  mutate(
    Ancestry = factor(Ancestry, levels = ancestry_levels),
    LG = factor(LG, levels = lg_levels) 
  )

# =========================================================================
# === FUNCTION: GENERATE CONSENSUS PLOT OBJECT ===
# =========================================================================
# Modified to return the ggplot object instead of immediately saving it
plot_consensus <- function(target_individuals, output_name, baseline_ancestry, filter_majority = FALSE, plot_tag = NULL) {
  
  # 1. Filter by EXACT matches only using the provided list
  pop_blocks <- blocks %>% filter(Individual %in% target_individuals)
  
  # Exit function early if no data matches
  if (nrow(pop_blocks) == 0) {
    warning(paste("Skipping", output_name, "- no matching individuals found.\n"))
    return(NULL)
  }
  
  # Calculate number of unique individuals making up this consensus based on the explicitly provided list
  n_inds <- length(target_individuals)
  
  # 2. Filter blocks based on logic
  if (filter_majority) {
    majority_ancestry <- pop_blocks %>%
      group_by(Individual) %>%
      count(Ancestry, wt = End_bp - Start_bp + 1, name = "total_bp") %>%
      slice_max(total_bp, n = 1, with_ties = FALSE) %>%
      dplyr::select(Individual, Majority = Ancestry) %>%
      ungroup()
    
    plot_blocks <- pop_blocks %>%
      left_join(majority_ancestry, by = "Individual") %>%
      filter(Ancestry != Majority) %>%
      filter(Ancestry != baseline_ancestry)
  } else {
    plot_blocks <- pop_blocks %>%
      filter(Ancestry != baseline_ancestry)
  }
  
  # -------------------------------------------------------------------------
  # EDITED HERE: Sort blocks so "Uncertain" is on the very bottom, followed 
  # by "CRT/IRT", and everything else is plotted on top.
  # -------------------------------------------------------------------------
  plot_blocks <- plot_blocks %>%
    arrange(Ancestry != "Uncertain", Ancestry != "CRT/IRT")
  
  # 3. Use output_name directly for the title
  plot_title <- output_name
  
  # 4. Create Plot
  p <- ggplot() +
    # Background dotted line for chromosomes
    geom_segment(data = chrom_lengths, aes(
      x = 0, xend = Chr_length/1e6,
      y = LG, yend = LG
    ), color = "black", linetype = "dotted", linewidth = 0.5) +
    
    # Inversion boundaries
    geom_rect(data = inversions, aes(
      xmin = Start_Mb, xmax = End_Mb,
      ymin = as.numeric(LG) - 0.35, ymax = as.numeric(LG) + 0.35,
      fill = "Inversion"
    ), color = "black", linewidth = 0.3) +
    
    # Ancestry blocks
    geom_segment(data = plot_blocks, aes(
      x = Start_bp/1e6, xend = End_bp/1e6,
      y = LG, yend = LG,
      color = Ancestry
    ), linewidth = 2, alpha = 0.7) +
    
    scale_y_discrete() + 
    scale_color_manual(values = ancestry_colors, drop = TRUE) +
    scale_fill_manual(name = NULL, values = c("Inversion" = "transparent")) +
    
    # Apply dynamic title, subtitle (N = X), and optional tag
    labs(
      title = plot_title,
      subtitle = bquote(italic(N) == .(n_inds)),
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
      strip.background = element_rect(fill = "gray90", color = "black"),
      panel.spacing = unit(0.3, "lines"),
      legend.position = "right", 
      legend.box = "vertical",
      legend.text = element_text(size = 10),
      legend.title = element_blank()
    )
  
  return(p)
}

# =========================================================================
# === RUN PLOTS & GENERATE MULTI-PANEL FIGURE ===
# =========================================================================

# 1. Read the explicit sample lists
bear_samples <- readLines("upper_bear_samples.txt")
efyaak_samples <- readLines("efyaak_samples.txt")
bear20_samples <- readLines("BearCreek_20_samples.txt")
bear21_samples <- readLines("Bear_21_samples.txt")
lower_bear_samples <- c(bear20_samples, bear21_samples)

# 2. Generate the individual plot objects
p_yaak <- plot_consensus(
  target_individuals = efyaak_samples,
  output_name = "Combined EF Yaak Consensus (EFYaak17_275, EForkYaak_20, EFYaak_21, EFYaak_22, EFYaak_24)",
  baseline_ancestry = "IRT",
  filter_majority = FALSE,
  plot_tag = "a."
)

p_bear_upper <- plot_consensus(
  target_individuals = bear_samples,
  output_name = "Combined Upper Bear Consensus (Bear_17, BearUp_22, Bear_24)",
  baseline_ancestry = "IRT",
  filter_majority = TRUE,
  plot_tag = "b."
)

# 3. Combine EF Yaak and Upper Bear into a 2-panel figure using patchwork
if(!is.null(p_yaak) & !is.null(p_bear_upper)) {
  
  # NEW: Save individual plots with their "a." and "b." labels before combining
  ggsave(file.path(output_dir, "EF_Yaak_Consensus_Individual.pdf"), p_yaak, width=12, height=7, dpi=300)
  ggsave(file.path(output_dir, "EF_Yaak_Consensus_Individual.png"), p_yaak, width=12, height=7, dpi=300)
  
  ggsave(file.path(output_dir, "Upper_Bear_Consensus_Individual.pdf"), p_bear_upper, width=12, height=7, dpi=300)
  ggsave(file.path(output_dir, "Upper_Bear_Consensus_Individual.png"), p_bear_upper, width=12, height=7, dpi=300)
  
  # Combine them
  combined_figure <- (p_yaak / p_bear_upper) + 
    plot_layout(guides = "collect")
  
  # Save the combined figure as a PDF for the main text
  combined_pdf_path <- file.path(output_dir, "Figure_EFYaak_UpperBear_Combined.pdf")
  combined_png_path <- file.path(output_dir, "Figure_EFYaak_UpperBear_Combined.png")
  
  ggsave(filename = combined_pdf_path, plot = combined_figure, width = 12, height = 12, dpi = 300)
  ggsave(filename = combined_png_path, plot = combined_figure, width = 12, height = 12, dpi = 300)
  cat("Saved combined 2-panel figure to:", combined_pdf_path, "and .png\n")
}

# 4. (Optional) Run & Save Individual Lower Bear Plots as separate PDFs and PNGs
p_lower <- plot_consensus(lower_bear_samples, "Combined Lower Bear Consensus", "IRT", TRUE)
if(!is.null(p_lower)) {
  ggsave(file.path(output_dir, "Combined_Lower_Bear.pdf"), p_lower, width=12, height=7, dpi=300)
  ggsave(file.path(output_dir, "Combined_Lower_Bear.png"), p_lower, width=12, height=7, dpi=300)
}

# Added plot_tag = "a."
p_bear20 <- plot_consensus(bear20_samples, "Bear_20 Consensus", "IRT", TRUE, plot_tag = "a.")
if(!is.null(p_bear20)) {
  ggsave(file.path(output_dir, "Bear_20_Consensus.pdf"), p_bear20, width=12, height=7, dpi=300)
  ggsave(file.path(output_dir, "Bear_20_Consensus.png"), p_bear20, width=12, height=7, dpi=300)
}

# Added plot_tag = "b."
p_bear21 <- plot_consensus(bear21_samples, "Bear_21 Consensus", "IRT", TRUE, plot_tag = "b.")
if(!is.null(p_bear21)) {
  ggsave(file.path(output_dir, "Bear_21_Consensus.pdf"), p_bear21, width=12, height=7, dpi=300)
  ggsave(file.path(output_dir, "Bear_21_Consensus.png"), p_bear21, width=12, height=7, dpi=300)
}

cat("Script complete!\n")