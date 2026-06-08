# ==== Load Libraries ====
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(purrr)
library(Metrics) # for rmse

# ==== File paths ====
admix_file <- "ADMIXTURE.tsv"
diag_file <- "input_ancestry.tsv"
norm_file <- "diagnostics_normalized.tsv"
elai_file <- "ELAI.tsv"

# ==== Read data ====
admix <- read_tsv(admix_file, col_types = cols())
diag  <- read_tsv(diag_file, col_types = cols())
norm  <- read_tsv(norm_file, col_types = cols())
elai  <- read_tsv(elai_file, col_types = cols())

# ==== Convert to long format ====

admix_long <- admix %>%
  select(Sample_ID,
         CRT = CRT_admix,
         IRT = IRT_admix,
         RBT = RBT_admix,
         WCT = WCT_admix,
         YCT = YCT_admix) %>%
  pivot_longer(cols = -Sample_ID, names_to = "Ancestry", values_to = "ADMIXTURE")

diag_long <- diag %>%
  select(Sample_ID,
         CRT = pCRT,
         IRT = pIRT,
         RBT = pRBT,
         WCT = pWCT,
         YCT = pYCT) %>%
  pivot_longer(cols = -Sample_ID, names_to = "Ancestry", values_to = "Diagnostics")

norm_long <- norm %>%
  select(Sample_ID,
         CRT = norm_pCRT,
         IRT = norm_pIRT,
         RBT = norm_pRBT,
         WCT = norm_pWCT,
         YCT = norm_pYCT) %>%
  pivot_longer(cols = -Sample_ID, names_to = "Ancestry", values_to = "Normalized")

elai_long <- elai %>%
  select(Sample_ID,
         CRT = CRT,
         IRT = IRT,
         RBT = pRBT,
         WCT = WCT,
         YCT = YCT) %>%
  pivot_longer(cols = -Sample_ID, names_to = "Ancestry", values_to = "ELAI")

# ==== Merge into one long-format table ====
combined <- diag_long %>%
  left_join(norm_long, by = c("Sample_ID", "Ancestry")) %>%
  left_join(admix_long, by = c("Sample_ID", "Ancestry")) %>%
  left_join(elai_long, by = c("Sample_ID", "Ancestry"))

# Save combined table
write_tsv(combined, "ancestry_estimate_comparison.tsv")

# ==== Summary statistics ====
# Compare each method against ADMIXTURE as baseline
methods <- c("Diagnostics", "Normalized", "ELAI")

summary_stats <- combined %>%
  group_by(Ancestry) %>%
  summarise(
    across(all_of(methods),
           list(
             cor = ~ cor(.x, ADMIXTURE, use = "complete.obs"),
             rmse = ~ rmse(ADMIXTURE, .x)
           ),
           .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

write_tsv(summary_stats, "ancestry_summary_statistics.tsv")
print(summary_stats)

# ==== Visualization ====
plot_data <- combined %>%
  pivot_longer(cols = c(Diagnostics, Normalized, ADMIXTURE, ELAI),
               names_to = "Method", values_to = "Estimate")

# Line plot: estimates by sample
ggplot(plot_data, aes(x = Method, y = Estimate, color = Method, group = Sample_ID)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ Ancestry, scales = "free_y") +
  labs(title = "Ancestry Estimates by Method",
       x = "Method", y = "Ancestry Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ancestry_estimate_comparison_plot.pdf", width = 12, height = 7)

# Jitter plot: distribution across methods
ggplot(plot_data, aes(x = Method, y = Estimate, color = Method)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ Ancestry, scales = "free_y") +
  labs(title = "Ancestry Estimates by Method",
       x = "Method", y = "Ancestry Proportion") +
  theme_bw()
ggsave("ancestry_estimate_jitter_plot.pdf", width = 12, height = 7)

# List of methods
methods <- c("Diagnostics", "Normalized", "ADMIXTURE", "ELAI")

# Function to calculate all pairwise correlations + RMSE for a given ancestry
pairwise_stats <- function(df, ancestry) {
  df_sub <- df %>% filter(Ancestry == ancestry) %>% select(all_of(methods))
  
  # Generate all pairwise combinations
  combs <- combn(methods, 2, simplify = FALSE)
  
  # Compute correlations + RMSE
  map_dfr(combs, ~ {
    tibble(
      Ancestry = ancestry,
      Method1 = .x[1],
      Method2 = .x[2],
      Correlation = cor(df_sub[[.x[1]]], df_sub[[.x[2]]], use = "complete.obs"),
      RMSE = rmse(df_sub[[.x[1]]], df_sub[[.x[2]]])
    )
  })
}

# Apply to all ancestries
ancestries <- unique(combined$Ancestry)
pairwise_stats_df <- map_dfr(ancestries, ~ pairwise_stats(combined, .x))

# View results
print(pairwise_stats_df)
write_tsv(pairwise_stats_df, "ancestry_pairwise_stats.tsv")
