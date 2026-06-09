#!/usr/bin/env Rscript

# ============================================================
# Arm-to-Arm Linkage (Fusion/Fission Test)
# Test linkage between "p-arm" (Seg1) and "q-arm" (Seg2).
# ============================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(broom)

# --- Directories and files ---
input_dir <- "./counts_output"
output_dir <- "fusion_tests"
dir.create(output_dir, showWarnings = FALSE)

set_name <- "RBT"
counts_file <- file.path(input_dir, paste0(set_name, "_counts.tsv"))
if(!file.exists(counts_file)) stop("Counts file not found: ", counts_file)

hybrid_table_file <- "ancestry_table.tsv"
if(!file.exists(hybrid_table_file)) stop("Missing hybrid table: ", hybrid_table_file)

strict_classes_file <- "triangle_hybrid_classifications.tsv"
if(!file.exists(strict_classes_file)) stop("Run triangle script first! Missing: ", strict_classes_file)

yct_metrics_file <- file.path(input_dir, "YCT_diagnostic_metrics.tsv")
if(!file.exists(yct_metrics_file)) stop("Missing YCT diagnostic metrics: ", yct_metrics_file)

duplicate_file <- "duplicate_samples.txt"
whitelist_file  <- "whitelist_samples.txt"

# --- Known F1 samples ---
true_f1s <- c("GOC1_03_18","LSM_14_11","ParmenterCr_20_007",
              "W_14_015","KRBY_22_06","KRBY_22_02","W_14_008","SML_20_157")

# --- Fusion Targets (Centromere Breakpoints) ---
fusion_list <- list(
  "Omy08" = 41651000,
  "Omy11" = 39014000,
  "Omy15" = 39379000,
  "Omy01" = 45000000 # Control
)

# --- MASTER MANUSCRIPT PALETTE (9 Colors, Okabe-Ito + Grey) ---
color_dict <- c(
  "Unadmixed RBT"      = "#D55E00", # Vermilion
  "Unadmixed WCT"      = "#0072B2", # Blue
  "Unadmixed YCT"      = "#009E73", # Teal
  "F1"                 = "#F0E442", # Yellow
  "F2"                 = "#CC79A7", # Pink/Rose
  "Complex Admixture"  = "#999999", # Grey
  "Backcross (to RBT)" = "#E69F00", # Gold
  "Backcross (to WCT)" = "#56B4E9", # Sky Blue
  "Backcross (to YCT)" = "#000000"  # Black
)

# --- LEGEND LABELS ---
label_dict <- c(
  "Unadmixed RBT"      = "Unadmixed RBT",
  "Unadmixed WCT"      = "Unadmixed WCT",
  "Unadmixed YCT"      = "Unadmixed YCT",
  "F1"                 = "F1",
  "F2"                 = "F2",
  "Complex Admixture"  = "Complex Admixture",
  "Backcross (to RBT)" = "Backcross (to RBT)",
  "Backcross (to WCT)" = "Backcross (to WCT)",
  "Backcross (to YCT)" = "Backcross (to YCT)"
)

# --- Unified Theme ---
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "right"
    )
}

# --- Read genotype counts ---
df <- read_tsv(counts_file, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
  mutate(
    CHROM = sapply(strsplit(Locus,"_"), `[`, 1),
    POS   = as.integer(sapply(strsplit(Locus,"_"), `[`, 2))
  )
sample_cols <- setdiff(colnames(df), c("Locus","CHROM","POS"))
df[sample_cols] <- map_df(df[sample_cols], ~ as.numeric(.x))

# --- Read ancestry, strict classifications, and YCT metrics ---
hybrid_table <- read_tsv(hybrid_table_file, show_col_types = FALSE)
strict_classes <- read_tsv(strict_classes_file, show_col_types = FALSE)
yct_data <- read_tsv(yct_metrics_file, show_col_types = FALSE) %>%
  rename(YCT_Het = Heterozygosity, YCT_Raw = Raw_Ancestry)

duplicate_ids <- if(file.exists(duplicate_file)) readLines(duplicate_file) %>% .[. != ""] else character()
whitelist_ids <- if(file.exists(whitelist_file)) readLines(whitelist_file) %>% .[. != ""] else hybrid_table$Sample_ID

hybrid_table <- hybrid_table %>%
  filter(!Sample_ID %in% duplicate_ids & Sample_ID %in% whitelist_ids)

# --- Hierarchical Hybrid Classification Logic ---
hybrid_table <- hybrid_table %>%
  left_join(strict_classes, by = "Sample_ID") %>%
  left_join(yct_data, by = "Sample_ID") %>%
  rowwise() %>%
  mutate(
    hybrid_class = case_when(
      # 1. Biological override for known F1s
      Sample_ID %in% true_f1s ~ "F1",
      # 2. Inherit strict class from triangle plot if present
      !is.na(Hybrid_Class) ~ as.character(Hybrid_Class),
      # 3. Strict 99% check for YCTs
      pYCT >= 0.99 | (!is.na(YCT_Raw) & YCT_Raw >= 0.99) ~ "Unadmixed YCT",
      # 4. Strict YCT Backcrosses (Geometric bounds)
      !is.na(YCT_Het) & !is.na(YCT_Raw) & YCT_Raw >= 0.50 & 
        (YCT_Het - 2 * (1 - YCT_Raw)) >= -0.1 & (YCT_Het - 2 * (1 - YCT_Raw)) <= 0.1 ~ "Backcross (to YCT)",
      # 5. Fallback
      TRUE ~ "Complex Admixture"
    )
  ) %>%
  ungroup() %>%
  mutate(hybrid_class = factor(hybrid_class, levels = names(color_dict)))

# --- Heterozygosity function ---
het_fraction <- function(x) {
  x <- as.numeric(x)
  valid <- x %in% c(0,1,2)
  if(sum(valid) == 0) return(NA_real_)
  sum(x==1, na.rm=TRUE)/sum(valid)
}

save_plot <- function(plot, name, width=10, height=8){
  ggsave(file.path(output_dir, paste0(name,".png")), plot, width=width, height=height, dpi=300)
  ggsave(file.path(output_dir, paste0(name,".pdf")), plot, width=width, height=height, device="pdf", dpi=300)
}

# --- Main Analysis Loop ---
combined_plot_data <- list()

for(chrom in names(fusion_list)){
  breakpoint <- fusion_list[[chrom]]
  message("Processing ", chrom, " breakpoint ", breakpoint)
  
  seg1 <- df %>% filter(CHROM==chrom & POS < breakpoint) %>% select(all_of(sample_cols))
  seg2 <- df %>% filter(CHROM==chrom & POS >= breakpoint) %>% select(all_of(sample_cols))
  
  if(nrow(seg1) == 0 | nrow(seg2) == 0){
    message("  Skipping (insufficient loci)")
    next
  }
  
  valid_samples <- intersect(colnames(seg1), hybrid_table$Sample_ID)
  seg1_het <- sapply(seg1[valid_samples], het_fraction)
  seg2_het <- sapply(seg2[valid_samples], het_fraction)
  
  plot_df <- tibble(
    individual = valid_samples,
    seg1_het = seg1_het[valid_samples],
    seg2_het = seg2_het[valid_samples],
    chromosome = chrom
  ) %>%
    left_join(hybrid_table %>% select(Sample_ID, hybrid_class), by=c("individual"="Sample_ID")) %>%
    filter(!is.na(seg1_het) & !is.na(seg2_het) & !is.na(hybrid_class))
  
  if(nrow(plot_df) < 3) next
  combined_plot_data[[chrom]] <- plot_df
  
  # Boxplot
  plot_df_box <- plot_df %>%
    mutate(
      seg1_status = case_when(
        is.na(seg1_het) ~ NA_character_,
        between(seg1_het, 0.80, 1.0) ~ "het",
        TRUE ~ "non_het"
      ),
      seg1_status = factor(seg1_status, levels=c("non_het","het"))
    ) %>% filter(!is.na(seg1_status))
  
  p_box <- ggplot(plot_df_box, aes(x=seg1_status, y=seg2_het, fill=hybrid_class)) +
    geom_boxplot(outlier.size=0.5, alpha=0.8) +
    facet_wrap(~hybrid_class, scales="free_x") +
    scale_fill_manual(values=color_dict, labels=label_dict, drop=TRUE) + 
    labs(
      title = paste0(chrom, " Arm Coupling: p-arm vs q-arm"),
      x = "p-arm Genotype", 
      y = "q-arm Heterozygosity"
    ) +
    theme_pub() + theme(legend.position="none")
  
  save_plot(p_box, paste0("Boxplot_", chrom, "_Fusion_Test"), width=12, height=6)
}

# --- MANUSCRIPT FIGURE: Faceted Grid ---
if(length(combined_plot_data) > 0){
  message("Generating Manuscript Grid Figure...")
  combined_df <- bind_rows(combined_plot_data)
  combined_df <- combined_df %>%
    mutate(chromosome = ifelse(chromosome == "Omy01", "Chr01 (Control)", chromosome))
  
  correlations <- combined_df %>%
    group_by(chromosome) %>%
    summarise(r = cor(seg1_het, seg2_het, method="pearson", use="complete.obs"), .groups="drop")
  
  combined_df <- combined_df %>%
    left_join(correlations, by="chromosome") %>%
    mutate(chrom_label = sprintf('atop(bold("%s"), bold("(") * italic(r) ~ bold("= %s)"))', chromosome, sprintf("%.2f", r)))
  
  p_grid <- ggplot(combined_df, aes(x=seg1_het, y=seg2_het, fill=hybrid_class)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey80") + 
    geom_point(shape=21, color="black", size=2.5, stroke=0.3, alpha=0.85) +
    facet_wrap(~chrom_label, ncol=2, scales = "free", labeller = label_parsed) + 
    scale_fill_manual(values=color_dict, labels=label_dict, drop=TRUE) + 
    scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) +
    scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) +
    labs(
      x = "Arm 1 Heterozygosity", 
      y = "Arm 2 Heterozygosity", 
      fill = "Inferred Hybrid Class"
    ) +
    theme_pub() +
    theme(
      aspect.ratio = 1,
      panel.spacing = unit(1, "lines"),
      legend.position = "right", 
      legend.background = element_rect(fill = "white", color = NA) 
    )
  
  save_plot(p_grid, "Manuscript_Fig_Fusion_Summary_2x2", width=10, height=8) 
  message("Saved Manuscript Figure.")
}

cat("Done.\n")