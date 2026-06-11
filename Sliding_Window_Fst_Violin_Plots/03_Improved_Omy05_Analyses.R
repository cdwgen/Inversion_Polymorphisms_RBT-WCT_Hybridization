library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(broom)
library(effsize)
library(purrr)
library(rstatix)
library(patchwork)

# ==== Theme ====
theme_molecol <- function() {
  theme_classic(base_size = 12, base_family = "sans") +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12, color = "black"),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
}

# ==== Colors ====
# Retaining your original colors, adding a distinct blue for the chromosome background
cb_palette <- c(
  "Genome-wide" = "#004488",
  "Omy05 background" = "#56B4E9", 
  "Inversion" = "#E69F00"
)

# ==== Load ONLY windowed FST ====
fst_files <- list.files(pattern = "^fst_.*_vs_.*\\.windowed\\.weir\\.fst$") %>%
  .[!grepl("_global|_site", .)]

load_fst <- function(file) {
  read_tsv(file, comment = "#", col_types = cols()) %>%
    rename(scaffold = CHROM, start = BIN_START, end = BIN_END, fst = MEAN_FST) %>%
    mutate(pair = str_remove_all(file, "^fst_|\\.windowed\\.weir\\.fst$"))
}

# Keep negatives for accurate stats
fst_data <- bind_rows(lapply(fst_files, load_fst))

# ==== Labels ====
pair_expressions <- c(
  "CRT_AA_vs_CRT_RR" = expression(CRT[AA] ~ "vs" ~ CRT[RR]),
  "CRT_AA_vs_IRT"    = expression(CRT[AA] ~ "vs IRT"),
  "CRT_AA_vs_WCT"    = expression(CRT[AA] ~ "vs WCT"),
  "CRT_AA_vs_YCT"    = expression(CRT[AA] ~ "vs YCT"),
  "CRT_RR_vs_IRT"    = expression(CRT[RR] ~ "vs IRT"),
  "CRT_RR_vs_WCT"    = expression(CRT[RR] ~ "vs WCT"),
  "CRT_RR_vs_YCT"    = expression(CRT[RR] ~ "vs YCT"),
  "IRT_vs_WCT"       = expression("IRT vs WCT"),
  "IRT_vs_YCT"       = expression("IRT vs YCT"),
  "WCT_vs_YCT"       = expression("WCT vs YCT")
)

# ==== Inversion coordinates ====
inv_start <- 30102000
inv_end   <- 87502000

# ==== Define regions (ALL THREE) ====
fst_data <- fst_data %>%
  mutate(region_type = case_when(
    scaffold == "Omy05" & start >= inv_start & end <= inv_end ~ "Inversion",
    scaffold == "Omy05" ~ "Omy05 background",
    TRUE ~ "Genome-wide"
  )) %>%
  mutate(region_type = factor(region_type,
                              levels = c("Genome-wide", "Omy05 background", "Inversion")))

# ==== Statistical tests ====
run_all_tests <- function(df) {
  
  inv <- df$fst[df$region_type == "Inversion"]
  omy <- df$fst[df$region_type == "Omy05 background"]
  gen <- df$fst[df$region_type == "Genome-wide"]
  
  # ---- Test 1: Inversion vs Genome-wide ----
  if(length(inv) >= 5 & length(gen) >= 5) {
    test1 <- wilcox.test(inv, gen, alternative = "greater", exact = FALSE)
    eff1  <- effsize::cliff.delta(inv, gen)$estimate
  } else {
    test1 <- list(p.value = NA)
    eff1  <- NA
  }
  
  # ---- Test 2: Inversion vs Omy05 ----
  if(length(inv) >= 5 & length(omy) >= 5) {
    test2 <- wilcox.test(inv, omy, alternative = "greater", exact = FALSE)
    eff2  <- effsize::cliff.delta(inv, omy)$estimate
  } else {
    test2 <- list(p.value = NA)
    eff2  <- NA
  }
  
  data.frame(
    P_genome = test1$p.value,
    Effect_genome = eff1,
    P_omy05 = test2$p.value,
    Effect_omy05 = eff2
  )
}

stats_results <- fst_data %>%
  group_by(pair) %>%
  reframe(run_all_tests(cur_data()))

# ==== Summary table ====
table_s2 <- fst_data %>%
  group_by(pair, region_type) %>%
  summarize(
    Mean_Fst   = round(mean(fst), 3),
    Median_Fst = round(median(fst), 3),
    StDev      = round(sd(fst), 3),
    N_Windows  = n(),
    .groups = "drop"
  ) %>%
  left_join(stats_results, by = "pair") %>%
  mutate(
    P_genome = case_when(
      P_genome < 0.001 ~ "<0.001",
      TRUE ~ as.character(round(P_genome, 3))
    ),
    P_omy05 = case_when(
      P_omy05 < 0.001 ~ "<0.001",
      TRUE ~ as.character(round(P_omy05, 3))
    ),
    Effect_genome = as.character(round(Effect_genome, 3)),
    Effect_omy05 = as.character(round(Effect_omy05, 3))
  ) %>%
  mutate(
    P_genome = ifelse(region_type == "Inversion", P_genome, "—"),
    Effect_genome = ifelse(region_type == "Inversion", Effect_genome, "—"),
    P_omy05 = ifelse(region_type == "Inversion", P_omy05, "—"),
    Effect_omy05 = ifelse(region_type == "Inversion", Effect_omy05, "—")
  ) %>%
  rename(
    `Taxon Pair` = pair,
    Region = region_type,
    `Mean FST` = Mean_Fst,
    `Median FST` = Median_Fst,
    `# Windows` = N_Windows,
    `P (vs Genome)` = P_genome,
    `Effect (vs Genome)` = Effect_genome,
    `P (vs Omy05)` = P_omy05,
    `Effect (vs Omy05)` = Effect_omy05
  )

write_csv(table_s2, "Table_S5_Fst_Statistics_Cleaned.csv")

# ==== Violin Plot ====
plot_violins_annotated <- function(data_subset) {
  
  # Explicitly define position dodge to perfectly align violins and boxplots
  pd <- position_dodge(width = 0.9)
  
  ggplot(data_subset, aes(x = pair, y = fst, fill = region_type)) +
    geom_violin(position = pd, scale = "width", trim = FALSE, alpha = 0.8, color = "white", linewidth = 0.2) +
    geom_boxplot(position = pd, width = 0.2, outlier.size = 0.1, alpha = 0.9, color = "black") +
    scale_x_discrete(labels = pair_expressions) +
    scale_fill_manual(
      values = cb_palette,
      labels = c("Genome-wide", "Chr05 Background", "Chr05 Complex") 
    ) +
    coord_cartesian(ylim = c(0, 1.05)) +  # Restored from original to clip negative tails cleanly
    labs(
      title = NULL, 
      y = expression(italic(F)[ST]), 
      x = NULL,
      fill = "Region" 
    ) +
    theme_molecol() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color="black"),
      legend.text = element_text(size = 11),
      legend.title = element_text(face = "bold", size = 12)
    )
}

# Generate and save plots
p_violin_all <- plot_violins_annotated(fst_data)

rbt_pairs_ids <- c("CRT_AA_vs_CRT_RR", "CRT_AA_vs_IRT", "CRT_RR_vs_IRT")
p_violin_rbt <- fst_data %>%
  filter(pair %in% rbt_pairs_ids) %>%
  plot_violins_annotated()

ggsave("Figure_Violin_All_Final.pdf", p_violin_all, width = 12, height = 6, dpi=300)
ggsave("Figure_Violin_All_Final.png", p_violin_all, width = 12, height = 6, dpi=600)
ggsave("Figure_Violin_RBT_Final.pdf", p_violin_rbt, width = 8, height = 6, dpi=300)

print(p_violin_all)

# ==== Export Word-Friendly Condensed Table ====
table_s2_word <- table_s2 %>%
  mutate(
    # Combine Mean, SD, and Median
    `FST: Mean ± SD (Median)` = paste0(`Mean FST`, " ± ", StDev, " (", `Median FST`, ")"),
    
    # Combine Effect Size and P-value for Genome comparison
    `vs Genome: Effect (p)` = case_when(
      Region == "Inversion" ~ paste0(`Effect (vs Genome)`, " (", `P (vs Genome)`, ")"),
      TRUE ~ "—"
    ),
    
    # Combine Effect Size and P-value for Omy05 comparison
    `vs Chr05: Effect (p)` = case_when(
      Region == "Inversion" ~ paste0(`Effect (vs Omy05)`, " (", `P (vs Omy05)`, ")"),
      TRUE ~ "—"
    )
  ) %>%
  # Select only the condensed columns
  select(
    `Taxon Pair`, 
    Region, 
    `FST: Mean ± SD (Median)`, 
    `# Windows`, 
    `vs Genome: Effect (p)`, 
    `vs Chr05: Effect (p)`
  )

write_csv(table_s2_word, "Table_S5_WordFriendly.csv")