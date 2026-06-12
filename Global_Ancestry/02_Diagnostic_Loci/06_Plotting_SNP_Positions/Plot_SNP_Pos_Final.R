library(tidyverse)

# ==============================================================================
# 1. SETUP & CONFIGURATION
# ==============================================================================

# Input Files
snp_file   <- "snps_pos.tsv" 
chrom_file <- "chrom_lengths.txt"

# Diagnostic Files Dictionary
diag_files <- list(
  "RBT" = "RBT.txt",
  "WCT" = "WCT.txt",
  "YCT" = "YCT.txt",
  "IRT" = "IRT.txt",
  "CRT" = "CRT.txt"
)

# Colors
diag_colors <- c(
  "YCT" = "#E3A63D", "WCT" = "#72428A", "RBT" = "#2C5AA0",
  "IRT" = "#E88634", "CRT" = "#62A2A3"
)
default_diag_color <- "#D55E00"

# Output settings
img_width  <- 12
img_height <- 8
img_dpi    <- 300

# Initialize Counter (starts at 1 for 'a')
plot_counter <- 1 

# ==============================================================================
# 2. COMMON DATA PROCESSING
# ==============================================================================
message(">>> Loading and cleaning data...")

# Load & Clean Chromosomes
chrom_lengths <- read_tsv(chrom_file, col_types = cols()) %>%
  rename(CHROM = LG, CHR_LEN = length) %>%
  mutate(CHROM = as.character(CHROM), CHR_LEN_MB = CHR_LEN / 1e6) %>%
  mutate(CHROM = ifelse(CHROM == "OmyY", "Omy29", CHROM))

# Load & Clean Background SNPs
snps_bg <- read_tsv(snp_file, col_names = c("CHROM", "POS"), col_types = "cd") %>%
  mutate(CHROM = ifelse(CHROM == "OmyY", "Omy29", CHROM)) %>%
  filter(CHROM %in% chrom_lengths$CHROM) %>%
  distinct() # <-- Added to remove duplicate SNPs

# Filter Chromosomes to match SNPs
chrom_lengths <- chrom_lengths %>% filter(CHROM %in% snps_bg$CHROM)

# Establish Order
chrom_order <- chrom_lengths %>%
  mutate(chr_num = suppressWarnings(as.numeric(gsub("Omy", "", CHROM)))) %>%
  arrange(chr_num) %>%
  pull(CHROM)

# Apply Factors & Create Backbone
chrom_lengths  <- chrom_lengths %>% mutate(CHROM = factor(CHROM, levels = chrom_order))
snps_bg        <- snps_bg       %>% mutate(CHROM = factor(CHROM, levels = chrom_order))
chrom_backbone <- chrom_lengths %>% transmute(CHROM, x0 = 0, x1 = CHR_LEN_MB)

# ==============================================================================
# 3. PLOT A: GENOME-WIDE DISTRIBUTION
# ==============================================================================
current_label <- paste0(letters[plot_counter], ".") # Creates "a."
message(paste0(">>> Generating Plot ", current_label))

p1 <- ggplot() +
  geom_segment(data = chrom_backbone, aes(x = x0, xend = x1, y = CHROM, yend = CHROM),
               linewidth = 5, color = "#E5E7EB", lineend = "round") +
  geom_linerange(data = snps_bg,
                 aes(x = POS / 1e6, ymin = as.numeric(CHROM) - 0.3, ymax = as.numeric(CHROM) + 0.3),
                 color = "#2C3E50", alpha = 0.3, linewidth = 0.15) +
  labs(
    title = "Genome-wide SNP Distribution",
    subtitle = paste0("Total SNPs: ", format(nrow(snps_bg), big.mark=",")),
    x = "Genomic Position (Mb)", y = NULL,
    tag = current_label 
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  theme_minimal(base_size = 14) +
  theme(
    # --- UPDATED SIZE TO MATCH TITLE ---
    plot.tag = element_text(face = "bold", size = 18, hjust = 0), 
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    
    plot.subtitle = element_text(color = "grey40", margin = margin(b = 20)),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    axis.text.x = element_text(color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linetype = "dashed")
  )

ggsave(paste0(letters[plot_counter], "_Genome_Wide.png"), p1, width = img_width, height = img_height, dpi = img_dpi, bg = "white")
plot_counter <- plot_counter + 1

# ==============================================================================
# 4. PLOTS B-Z: DIAGNOSTIC OVERLAYS
# ==============================================================================

for (set_name in names(diag_files)) {
  
  file_path  <- diag_files[[set_name]]
  this_color <- if(set_name %in% names(diag_colors)) diag_colors[[set_name]] else default_diag_color
  
  current_label <- paste0(letters[plot_counter], ".") # Creates "b.", "c.", etc.
  message(paste0(">>> Generating Plot ", current_label, " : ", set_name))
  
  # Load specific diagnostic file
  diag_snps <- read_tsv(file_path, col_names = c("CHROM", "POS"), col_types = "cd") %>%
    mutate(CHROM = ifelse(CHROM == "OmyY", "Omy29", CHROM)) %>%
    filter(CHROM %in% chrom_lengths$CHROM) %>%
    mutate(CHROM = factor(CHROM, levels = chrom_order))
  
  if(nrow(diag_snps) == 0) { warning("No SNPs for ", set_name); next }
  
  p_diag <- ggplot() +
    geom_segment(data = chrom_backbone, aes(x = x0, xend = x1, y = CHROM, yend = CHROM),
                 linewidth = 5, color = "#E5E7EB", lineend = "round") +
    geom_linerange(data = snps_bg,
                   aes(x = POS/1e6, ymin = as.numeric(CHROM) - 0.3, ymax = as.numeric(CHROM) + 0.3),
                   color = "#999999", alpha = 0.1, linewidth = 0.1) +
    geom_linerange(data = diag_snps,
                   aes(x = POS/1e6, ymin = as.numeric(CHROM) - 0.4, ymax = as.numeric(CHROM) + 0.4),
                   color = this_color, alpha = 1.0, linewidth = 0.8) +
    labs(
      title = paste0("Diagnostic Set: ", set_name),
      subtitle = paste0(format(nrow(diag_snps), big.mark=","), " markers"),
      x = "Genomic Position (Mb)", y = NULL,
      tag = current_label
    ) +
    scale_x_continuous(expand = c(0.01, 0)) +
    theme_minimal(base_size = 14) +
    theme(
      # --- UPDATED SIZE TO MATCH TITLE ---
      plot.tag = element_text(face = "bold", size = 18, hjust = 0),
      plot.title = element_text(face = "bold", size = 18, color = "black"),
      
      plot.subtitle = element_text(color = this_color, face = "bold", margin = margin(b=20)),
      axis.text.y = element_text(color = "black", size = 10, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92", linetype = "dashed")
    )
  
  # Save
  ggsave(paste0(letters[plot_counter], "_Diagnostic_", set_name, ".png"), p_diag, width = img_width, height = img_height, dpi = img_dpi, bg = "white")
  
  plot_counter <- plot_counter + 1
}

message(">>> All plots saved.")