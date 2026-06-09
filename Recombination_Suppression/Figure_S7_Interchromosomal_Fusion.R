#!/usr/bin/env Rscript

# ============================================================
# Manuscript Script: Inversion vs Target Chromosome (Fusion Test)
# Comparison: Omy20 Inversion vs Omy28 (Test) & Omy10 (Control)
# Style: Manuscript Ready (Right Legend, No Title, Control First)
# ============================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(broom)

# --- Directories and files ---
input_dir <- "./counts_output"
output_dir <- "inversion_vs_fusion_Omy20_Omy28"
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

# --- CONFIGURATION ---
inv_key_to_test <- "Chr20" 
target_chrom    <- "Omy28"
control_chrom   <- "Omy10"

# --- Inversion definitions ---
inversions <- list(
  "Chr20" = list(chrom="Omy20", start=4356000, end=21129000)
)

# --- COLOR PALETTE (Okabe-Ito, Aligned to Triangle Plot) ---
color_dict <- c(
  "Unadmixed RBT"      = "#D55E00", 
  "Unadmixed WCT"      = "#0072B2", 
  "Unadmixed YCT"      = "#009E73", 
  "F1"                 = "#F0E442", 
  "F2"                 = "#CC79A7", 
  "Complex Admixture"  = "#999999", 
  "Backcross (to RBT)" = "#E69F00", 
  "Backcross (to WCT)" = "#56B4E9", 
  "Backcross (to YCT)" = "#000000"  
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

# --- Theme ---
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black"),
      legend.position = "right",
      legend.justification = "center",
      legend.margin = margin(l = 10) 
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

# --- Read ancestry, classifications, and YCT metrics ---
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
    final_class = case_when(
      !is.na(Hybrid_Class) ~ as.character(Hybrid_Class),
      pRBT >= 0.99 ~ "Unadmixed RBT",
      pWCT >= 0.99 ~ "Unadmixed WCT",
      pYCT >= 0.99 | (!is.na(YCT_Raw) & YCT_Raw >= 0.99) ~ "Unadmixed YCT",
      !is.na(YCT_Het) & !is.na(YCT_Raw) & YCT_Raw >= 0.50 & 
        (YCT_Het - 2 * (1 - YCT_Raw)) >= -0.1 & (YCT_Het - 2 * (1 - YCT_Raw)) <= 0.1 ~ "Backcross (to YCT)",
      TRUE ~ "Complex Admixture" 
    )
  ) %>%
  ungroup() %>%
  mutate(hybrid_class = factor(final_class, levels = names(color_dict)))

# --- 6. Print YCT Backcross Individuals ---
yct_bx <- hybrid_table %>% filter(final_class == "Backcross (to YCT)") %>% pull(Sample_ID)
cat("\n--- Individuals classified as Backcross (to YCT) ---\n")
if(length(yct_bx) > 0) print(yct_bx) else cat("None found.\n")
cat("---------------------------------------------------\n\n")

# --- Heterozygosity function ---
het_fraction <- function(x) {
  x <- as.numeric(x)
  valid <- x %in% c(0,1,2)
  if(sum(valid) == 0) return(NA_real_)
  sum(x==1, na.rm=TRUE)/sum(valid)
}

save_plot <- function(plot, name, width=8, height=7){
  ggsave(file.path(output_dir, paste0(name,".png")), plot, width=width, height=height, dpi=300)
  ggsave(file.path(output_dir, paste0(name,".pdf")), plot, width=width, height=height, device="pdf", dpi=300)
}

# ============================================================
# Analysis & Plotting
# ============================================================

inv <- inversions[[inv_key_to_test]]
inv_loci <- df %>% filter(CHROM==inv$chrom & POS>=inv$start & POS<=inv$end) %>% select(all_of(sample_cols))
valid_samples <- intersect(colnames(inv_loci), hybrid_table$Sample_ID)

inv_het_vector <- sapply(inv_loci[valid_samples], het_fraction)

target_loci <- df %>% filter(CHROM == target_chrom) %>% select(all_of(sample_cols))
target_het_vector <- sapply(target_loci[valid_samples], het_fraction)

control_loci <- df %>% filter(CHROM == control_chrom) %>% select(all_of(sample_cols))
control_het_vector <- sapply(control_loci[valid_samples], het_fraction)

label_test <- paste0(target_chrom, " (Ocl20 Fusion)")
label_ctrl <- "Chr10 (Unfused)" 

df_test <- tibble(
  individual = valid_samples,
  x_het = inv_het_vector,
  y_het = target_het_vector,
  type_raw = label_test
)

df_control <- tibble(
  individual = valid_samples,
  x_het = inv_het_vector,
  y_het = control_het_vector,
  type_raw = label_ctrl
)

combined_df <- bind_rows(df_test, df_control) %>%
  left_join(hybrid_table %>% select(Sample_ID, hybrid_class), by=c("individual"="Sample_ID")) %>%
  filter(!is.na(x_het) & !is.na(y_het) & !is.na(hybrid_class))

stats <- combined_df %>%
  group_by(type_raw) %>%
  summarize(r_val = cor(x_het, y_het, method="pearson", use="complete.obs"), .groups="drop")

combined_df <- combined_df %>%
  left_join(stats, by="type_raw") %>%
  mutate(
    facet_label = sprintf('atop(bold("%s"), bold("(") * italic(r) ~ bold("= %s)"))', type_raw, sprintf("%.2f", r_val))
  )

unique_labels <- unique(combined_df$facet_label)
label_test_final <- unique_labels[grep(target_chrom, unique_labels)]
label_ctrl_final <- unique_labels[grep("Chr10", unique_labels)] 

combined_df$facet_label <- factor(combined_df$facet_label, 
                                  levels = c(label_ctrl_final, label_test_final))

message("Generating ordered side-by-side plot...")

p_facet <- ggplot(combined_df, aes(x=x_het, y=y_het, fill=hybrid_class)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey80") + 
  geom_point(shape=21, color="black", size=2.5, stroke=0.3, alpha=0.85) +
  facet_wrap(~facet_label, ncol=2, scales="free", labeller = label_parsed) +
  scale_fill_manual(values=color_dict, labels=label_dict, drop=TRUE) + 
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) +
  scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1)) +
  labs(
    title = NULL, 
    x = "Chr20 Inversion Heterozygosity",
    y = "Chromosome-wide Heterozygosity",
    fill = "Inferred Hybrid Class"
  ) +
  theme_pub() +
  theme(
    aspect.ratio = 1,
    panel.spacing = unit(1.0, "lines"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA) 
  )

save_plot(p_facet, paste0("Manuscript_Fig_Fusion_Test_Final"), width=11, height=5)
cat("Done. Plot saved as 'Manuscript_Fig_Fusion_Test_Final.pdf'\n")