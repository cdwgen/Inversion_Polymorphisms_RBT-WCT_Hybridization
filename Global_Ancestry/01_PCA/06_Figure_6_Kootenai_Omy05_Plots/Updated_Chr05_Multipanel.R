#!/usr/bin/env Rscript

## PCA for Chr05 Inversions - Fixed Layout (Named Mapping)
## Top Row: Panel A (60%) | Spacer (25%) | Legend (15%)
## Bottom Row: Panel B (50%) | Panel C (50%)

library(vcfR)
library(adegenet)
library(ade4)
library(ggplot2)
library(patchwork) 
library(dplyr)

# -------------------------------------------------------------------------
# 1. Configuration & Themes
# -------------------------------------------------------------------------

theme_molecol <- function() {
  theme_classic(base_size = 12, base_family = "sans") +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "plain", color = "black"),
      legend.text = element_text(size = 10, color = "black"),
      legend.title = element_text(face = "bold", size = 11, color = "black"),
      legend.title.align = 0.5,
      legend.position = "right", 
      # Important: Snap legend to top-left of its container box
      legend.justification = c(0, 1), 
      plot.margin = margin(5, 5, 5, 5)
    )
}

# Mapping and Colors
display_map <- c(
  "WCT_AA" = "WCT[AA]",
  "WCT_AxYCT_A" = "WCT[A] %*% YCT[A]",
  "RBT_RR" = "RBT[RR]",
  "RBT_AR" = "RBT[AR]",
  "WCT_A/RBT_R" = "WCT[A]/RBT[R]",
  "RBT_AA" = "RBT[AA]",
  "WCT_A/RBT_A" = "WCT[A]/RBT[A]",
  "WCT_A/RBT_A;RBT_AA" = "WCT[A]/RBT[A] * \";\" * RBT[AA]",
  "RBT_AR;WCT_AxRBT_R" = "RBT[AR] * \";\" * WCT[A] %*% RBT[R]",
  "YCT_AA" = "YCT[AA]",
  "RBT_AA;YCT_AxRBT_A" = "RBT[AA] * \";\" * YCT[A] %*% RBT[A]",
  "WCT_A/RBT_A_Inv1" = "'WCT[A]/RBT[A] (Inv2)'",
  "YCT_A/RBT_A_Inv2_Rec" = "'YCT[A]/RBT[A] (Inv2 Recombinant)'",
  "Recombinant" = "'Recombinant'"
)

custom_colors <- c(
  "RBT_AA"             = "#d07562",
  "RBT_AA;YCT_AxRBT_A" = "#f88dc5",
  "RBT_AR"             = "#65ae63",
  "RBT_AR;WCT_AxRBT_R" = "#cda13d",
  "RBT_RR"             = "#529c96",
  "WCT_AxYCT_A"        = "#797095",
  "WCT_A/RBT_A"        = "#ffa222",
  "WCT_A/RBT_A;RBT_AA" = "#fff342",
  "WCT_A/RBT_R"        = "#95709c",
  "WCT_AA"             = "#e73032",
  "YCT_AA"             = "#c47166",
  "WCT_A/RBT_A_Inv1"     = "grey50",
  "YCT_A/RBT_A_Inv2_Rec" = "grey50",
  "Recombinant"          = "grey50"
)

legend_order <- c(
  "RBT_AA", "RBT_AR", "RBT_RR", "WCT_AA", "YCT_AA",
  "WCT_A/RBT_R", "WCT_A/RBT_A", "WCT_A/RBT_A;RBT_AA",
  "RBT_AR;WCT_AxRBT_R", "RBT_AA;YCT_AxRBT_A", "WCT_AxYCT_A"
)

# -------------------------------------------------------------------------
# 2. PCA Function
# -------------------------------------------------------------------------

run_pca_internal <- function(vcf_path, popmap_path, filter_yct = FALSE) {
  message("Processing: ", vcf_path)
  
  vcf <- read.vcfR(vcf_path, verbose = FALSE)
  
  # Fix duplicate IDs
  dup_ids <- duplicated(vcf@fix[, "ID"])
  if (any(dup_ids)) {
    vcf@fix[dup_ids, "ID"] <- paste0(vcf@fix[dup_ids, "ID"], "_", seq_len(sum(dup_ids)))
  }
  genind_data <- vcfR2genind(vcf)
  
  popmap <- read.table(popmap_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(popmap) <- c("sample_id", "Omy05_Genotype", "group")
  
  if (filter_yct) {
    popmap <- popmap[!grepl("YCT", popmap$Omy05_Genotype), ]
    common_inds <- intersect(indNames(genind_data), popmap$sample_id)
    genind_data <- genind_data[common_inds, ]
    popmap <- popmap[match(common_inds, popmap$sample_id), ]
  } else {
    ind_names <- indNames(genind_data)
    popmap <- popmap[match(ind_names, popmap$sample_id), ]
  }
  
  pop(genind_data) <- as.factor(popmap$Omy05_Genotype)
  
  genind_tab <- tab(genind_data)
  genind_tab[is.na(genind_tab)] <- colMeans(genind_tab, na.rm = TRUE)[col(genind_tab)][is.na(genind_tab)]
  genind_data@tab <- genind_tab
  
  pca_result <- dudi.pca(genind_tab, center = TRUE, scale = FALSE, scannf = FALSE, nf = 3)
  pct_var <- round((pca_result$eig / sum(pca_result$eig)) * 100, 2)
  
  df <- data.frame(
    Axis1 = pca_result$li[,1],
    Axis2 = pca_result$li[,2],
    Omy05_Genotype = popmap$Omy05_Genotype,
    stringsAsFactors = FALSE
  )
  return(list(df = df, var = pct_var))
}

# -------------------------------------------------------------------------
# 3. Generate Data
# -------------------------------------------------------------------------

# A: Full Inversion (No YCT)
data_a <- run_pca_internal("02_Subset_VCF/Omy05_RBT_Full_subset.recode.vcf", "popmap_Omy05.txt", filter_yct = TRUE)
# B: Inv1
data_b <- run_pca_internal("02_Subset_VCF/Omy05_RBT_Inv1_subset.recode.vcf", "popmap_Omy05.txt")
# C: Inv2
data_c <- run_pca_internal("02_Subset_VCF/Omy05_RBT_Inv2_subset.recode.vcf", "popmap_Omy05.txt")

# -------------------------------------------------------------------------
# 4. Construct Plots
# -------------------------------------------------------------------------

make_gg <- function(pca_data, show_legend = TRUE) {
  ggplot(pca_data$df, aes(Axis1, Axis2, color = Omy05_Genotype)) +
    geom_point(size = 3, alpha = 0.9, show.legend = show_legend) +
    scale_color_manual(
      values = custom_colors,
      labels = function(breaks) { parse(text = display_map[breaks]) },
      breaks = legend_order,
      limits = legend_order,
      name = "Chr05 Genotype"
    ) +
    theme_molecol() +
    labs(
      title = NULL,
      x = paste0("PC1 (", pca_data$var[1], "%)"),
      y = paste0("PC2 (", pca_data$var[2], "%)")
    )
}

# Create Plot Objects
p_top       <- make_gg(data_a, show_legend = FALSE) # Panel A
p_bot_left  <- make_gg(data_b, show_legend = TRUE)  # Panel B (Legend Source)
p_bot_right <- make_gg(data_c, show_legend = TRUE)  # Panel C (Legend Source)

# -------------------------------------------------------------------------
# 5. Assemble Layout (NAMED MAPPING)
# -------------------------------------------------------------------------

message("Assembling final figure with strict NAMED mapping...")

# STRICT 20-COLUMN GRID
# Row 1: A (12u = 60%) | S (5u = 25%) | L (3u = 15%)
# Row 2: B (10u = 50%) | C (10u = 50%)
# The letters in 'design' match the arguments in wrap_plots() below.

layout_design <- "
AAAAAAAAAAAASSSSSLLL
BBBBBBBBBBCCCCCCCCCC
"

# We pass the plots as NAMED arguments. This prevents Panel B from accidentally jumping to 'S' or 'L'.
final_plot <- wrap_plots(
  A = p_top,
  B = p_bot_left,
  C = p_bot_right,
  L = guide_area(),  # This catches the legend
  S = plot_spacer(), # This creates the empty space
  design = layout_design
) + 
  plot_layout(guides = "collect", heights = c(1.3, 1)) + 
  plot_annotation(tag_levels = 'a', tag_suffix = '.') &
  theme(plot.tag = element_text(face = 'bold', size = 16))

# -------------------------------------------------------------------------
# 6. Save
# -------------------------------------------------------------------------

ggsave("Figure_Chr05_Balanced_Multipanel.pdf", final_plot, width = 12, height = 10, dpi = 500)
ggsave("Figure_Chr05_Balanced_Multipanel.png", final_plot, width = 12, height = 10, dpi = 500)

message("Done! Saved Figure_Chr05_Balanced_Multipanel.")
print(final_plot)