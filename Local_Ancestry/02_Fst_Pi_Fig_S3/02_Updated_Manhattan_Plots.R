# ==== Load Libraries ====
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)

# ==== Custom Theme ====
theme_molecol <- function() {
  theme_classic(base_size = 12, base_family = "sans") +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12, color = "black"),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
}

# ==== Labels ====
custom_labeller_exprs <- c(
  "CRT_AA_vs_CRT_RR" = 'bold(CRT[AA]~"vs"~CRT[RR])',
  "CRT_AA_vs_IRT"    = 'bold(CRT[AA]~"vs IRT")',
  "CRT_AA_vs_WCT"    = 'bold(CRT[AA]~"vs WCT")',
  "CRT_AA_vs_YCT"    = 'bold(CRT[AA]~"vs YCT")',
  "CRT_RR_vs_IRT"    = 'bold(CRT[RR]~"vs IRT")',
  "CRT_RR_vs_WCT"    = 'bold(CRT[RR]~"vs WCT")',
  "CRT_RR_vs_YCT"    = 'bold(CRT[RR]~"vs YCT")',
  "IRT_vs_WCT"       = 'bold("IRT vs WCT")',
  "IRT_vs_YCT"       = 'bold("IRT vs YCT")',
  "WCT_vs_YCT"       = 'bold("WCT vs YCT")',
  "CRT_AA" = 'bold(CRT[AA])',
  "CRT_RR" = 'bold(CRT[RR])',
  "IRT"    = 'bold("IRT")',
  "WCT"    = 'bold("WCT")',
  "YCT"    = 'bold("YCT")'
)

# ==== Load π and FST files ====

pi_files <- list.files(pattern = "^pi_.*\\.windowed\\.pi$")

# FIX: only load TRUE windowed FST (exclude global + site)
fst_files <- list.files(pattern = "^fst_.*_vs_.*\\.windowed\\.weir\\.fst$") %>%
  .[!grepl("_global|_site", .)]

load_pi <- function(file) {
  read_tsv(file, comment = "#", col_types = cols()) %>%
    rename(scaffold = CHROM, start = BIN_START, end = BIN_END, pi = PI) %>%
    mutate(pop = str_remove_all(file, "^pi_|\\.windowed\\.pi$"))
}

load_fst <- function(file) {
  read_tsv(file, comment = "#", col_types = cols()) %>%
    rename(scaffold = CHROM, start = BIN_START, end = BIN_END, fst = MEAN_FST) %>%
    mutate(pair = str_remove_all(file, "^fst_|\\.windowed\\.weir\\.fst$"))
}

pi_data <- bind_rows(lapply(pi_files, load_pi))

# CRITICAL FIX: Floor negative FST values to 0 strictly for clean visualization
fst_data <- bind_rows(lapply(fst_files, load_fst)) %>%
  mutate(fst = ifelse(fst < 0, 0, fst))

# ==== Set Factor Levels for Plot Ordering ====
# This forces ggplot to render the FST pairs in the exact order of your supplementary table
pair_order <- c(
  "CRT_AA_vs_CRT_RR",
  "CRT_AA_vs_IRT",
  "CRT_AA_vs_WCT",
  "CRT_AA_vs_YCT",
  "CRT_RR_vs_IRT",
  "CRT_RR_vs_WCT",
  "CRT_RR_vs_YCT",
  "IRT_vs_WCT",
  "IRT_vs_YCT",
  "WCT_vs_YCT"
)

fst_data <- fst_data %>%
  mutate(pair = factor(pair, levels = pair_order))

# ==== Chromosome order ====
chrom_order <- c(sprintf("Omy%02d", 1:28), "OmyY", sprintf("Omy%02d", 30:32))

# ==== Offsets ====
build_offsets <- function(df, scaffold_order) {
  df %>%
    group_by(scaffold) %>%
    summarize(max_pos = max(end, na.rm = TRUE), .groups = "drop") %>%
    filter(scaffold %in% scaffold_order) %>%
    mutate(scaffold = factor(scaffold, levels = scaffold_order)) %>%
    arrange(scaffold) %>%
    mutate(offset = lag(cumsum(max_pos + 1e6), default = 0))
}

scaffold_offsets <- build_offsets(bind_rows(
  pi_data %>% select(scaffold, end),
  fst_data %>% select(scaffold, end)
), chrom_order)

scaffold_labels <- scaffold_offsets %>%
  mutate(center = offset + max_pos / 2)

# ==== Add positions ====
add_positions <- function(df, offsets) {
  df %>%
    inner_join(offsets, by = "scaffold") %>%
    mutate(pos = start + offset)
}

pi_data <- add_positions(pi_data, scaffold_offsets)
fst_data <- add_positions(fst_data, scaffold_offsets)

# ==== Plot functions ====

plot_metric_facet <- function(df, yvar, groupvar, ylab_expr, output_file) {
  
  p <- ggplot(df, aes(x = pos, y = .data[[yvar]])) +
    geom_point(aes(color = scaffold), size = 0.3, alpha = 0.8) +
    scale_color_manual(values = rep(c("gray30", "skyblue3"),
                                    length.out = n_distinct(df$scaffold))) +
    geom_vline(xintercept = scaffold_offsets$offset,
               linetype = "dotted", color = "gray60", linewidth = 0.3) +
    scale_x_continuous(breaks = scaffold_labels$center,
                       labels = scaffold_labels$scaffold,
                       expand = expansion(mult = c(0.01, 0.01))) +
    facet_wrap(as.formula(paste("~", groupvar)),
               scales = "free_y", ncol = 1,
               labeller = as_labeller(custom_labeller_exprs,
                                      default = label_parsed)) +
    labs(x = "Chromosome", y = ylab_expr) +
    theme_molecol() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  
  ggsave(output_file, p,
         width = 12,
         height = 1.8 * length(unique(df[[groupvar]])))
}

plot_metric_individual <- function(df, yvar, groupvar, ylab_expr, prefix) {
  for (grp in unique(df[[groupvar]])) {
    
    df_sub <- df %>% filter(.data[[groupvar]] == grp)
    parsed_title <- parse(text = custom_labeller_exprs[grp])
    
    p <- ggplot(df_sub, aes(x = pos, y = .data[[yvar]])) +
      geom_point(aes(color = scaffold), size = 0.4, alpha = 0.8) +
      scale_color_manual(values = rep(c("gray30", "skyblue3"),
                                      length.out = n_distinct(df_sub$scaffold))) +
      geom_vline(xintercept = scaffold_offsets$offset,
                 linetype = "dotted", color = "gray60", linewidth = 0.3) +
      scale_x_continuous(breaks = scaffold_labels$center,
                         labels = scaffold_labels$scaffold,
                         expand = expansion(mult = c(0.01, 0.01))) +
      labs(x = "Chromosome", y = ylab_expr, title = parsed_title) +
      theme_molecol() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank()
      )
    
    ggsave(paste0(prefix, "_", grp, ".png"), p, width = 12, height = 4)
  }
}

plot_metric_Omy05 <- function(df, yvar, groupvar, ylab_expr, prefix) {
  
  df_Omy05 <- df %>% filter(scaffold == "Omy05")
  
  for (grp in unique(df_Omy05[[groupvar]])) {
    
    df_grp <- df_Omy05 %>% filter(.data[[groupvar]] == grp)
    
    title_expr <- paste0('"', prefix, ' "~',
                         custom_labeller_exprs[grp],
                         '~" Omy05"')
    
    p <- ggplot(df_grp, aes(x = start, y = .data[[yvar]])) +
      geom_point(size = 0.6, alpha = 0.9, color = "steelblue") +
      labs(x = "Position (bp)", y = ylab_expr,
           title = parse(text = title_expr)) +
      theme_molecol()
    
    ggsave(paste0(prefix, "_Omy05_", grp, ".pdf"),
           p, width = 8, height = 4)
  }
}

# ==== Run plots ====

plot_metric_facet(pi_data, "pi", "pop",
                  expression(pi~"(nucleotide diversity)"),
                  "pi_manhattan_facet.pdf")

plot_metric_individual(pi_data, "pi", "pop",
                       expression(pi~"(nucleotide diversity)"),
                       "pi_manhattan")

plot_metric_Omy05(pi_data, "pi", "pop",
                  expression(pi~"(nucleotide diversity)"),
                  "pi")

plot_metric_facet(fst_data, "fst", "pair",
                  expression(italic(F)[ST]),
                  "fst_manhattan_facet.png")

plot_metric_individual(fst_data, "fst", "pair",
                       expression(italic(F)[ST]),
                       "fst_manhattan")

plot_metric_Omy05(fst_data, "fst", "pair",
                  expression(italic(F)[ST]),
                  "fst")