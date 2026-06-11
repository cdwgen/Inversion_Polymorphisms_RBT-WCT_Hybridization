#!/usr/bin/env Rscript

# ============================================================
# Habitat-genotype relationships (Final Polish)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(DirichletReg)
})

# ============================================================
# 1. Configuration: Publication Style
# ============================================================
theme_molecol <- function() {
  theme_classic(base_size = 12, base_family = "sans") +
    theme(
      text = element_text(color = "black"),
      
      # Axis text size 11 for balance
      axis.text = element_text(color = "black", size = 11),
      
      # FIX: Changed face from "bold" to "plain"
      axis.title = element_text(face = "plain", color = "black", size = 13),
      
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12, color = "black"),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
}

# Custom Palette
geno_colors <- c("RBT_AA" = "#cb6651", "RBT_AR" = "#54a451", "RBT_RR" = "#3f908a")

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------
df <- read.table("Omy05_table.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)

df$Habitat_Type <- factor(trimws(df$Habitat_Type))
df$RBT_AA <- as.numeric(df$RBT_AA)
df$RBT_AR <- as.numeric(df$RBT_AR)
df$RBT_RR <- as.numeric(df$RBT_RR)
df$N       <- as.numeric(df$N)
df_clean <- df[!is.na(df$Habitat_Type), ]

habitat_lookup <- c(
  "River" = "River",
  "Tributary" = "Tributary",
  "Lake" = "Lake",
  "Valley_Floor" = "Valley Floor",
  "Reservoir" = "Reservoir"
)
df_clean$Habitat_Display <- habitat_lookup[as.character(df_clean$Habitat_Type)]
df_clean$Habitat_Display <- factor(df_clean$Habitat_Display,
                                   levels = c("Tributary", "Valley Floor", "River", "Lake", "Reservoir"))

# ------------------------------------------------------------
# Stats (Summary Tables + Tests)
# ------------------------------------------------------------
df_clean <- df_clean %>%
  mutate(
    AA_count = round(RBT_AA * N),
    AR_count = round(RBT_AR * N),
    RR_count = round(RBT_RR * N)
  )

geno_counts <- df_clean %>%
  group_by(Habitat_Display) %>%
  summarise(
    AA = sum(AA_count),
    AR = sum(AR_count),
    RR = sum(RR_count),
    Total_N = sum(N),
    .groups = "drop"
  )
write.table(geno_counts, "habitat_geno_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

geno_matrix <- as.matrix(geno_counts[, c("AA", "AR", "RR")])
rownames(geno_matrix) <- geno_counts$Habitat_Display
chi_res <- chisq.test(geno_matrix, simulate.p.value = TRUE, B = 10000)

df_dir <- DR_data(df_clean[, c("RBT_AA", "RBT_AR", "RBT_RR")])
fit_dir <- DirichReg(df_dir ~ Habitat_Display, data = df_clean)

# Save Stats
sink("stats_results.txt")
cat("========== Chi-square on counts ==========\n")
print(chi_res)
cat("\n========== Dirichlet regression ==========\n")
print(summary(fit_dir))
sink()

# ------------------------------------------------------------
# Visualization
# ------------------------------------------------------------
df_long <- df_clean %>%
  select(Habitat_Display, RBT_AA, RBT_AR, RBT_RR) %>%
  pivot_longer(cols = starts_with("RBT_"), names_to = "Genotype", values_to = "Frequency")

df_long$Genotype <- factor(df_long$Genotype, levels = c("RBT_AA", "RBT_AR", "RBT_RR"))

# Define Labels (Subscripts)
geno_labels <- c(
  "RBT_AA" = expression(RBT["AA"]), 
  "RBT_AR" = expression(RBT["AR"]), 
  "RBT_RR" = expression(RBT["RR"])
)

# 1. Stacked Bar Plot
p1 <- ggplot(df_long, aes(x = Habitat_Display, y = Frequency, fill = Genotype)) +
  stat_summary(fun = mean, geom = "bar", position = "stack", width = 0.9) +
  
  scale_fill_manual(
    values = geno_colors, 
    labels = geno_labels, 
    name = "Chr05 Genotype" # FIX: Updated Legend Title
  ) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  
  labs(y = "Mean Frequency", x = "Habitat Type") +
  theme_molecol() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    # Ensure legend title is bold (defined in theme, but good to double check visual consistency)
    legend.title = element_text(face = "bold", size = 11)
  )

# 2. Boxplots
p2 <- ggplot(df_long, aes(x = Habitat_Display, y = Frequency, fill = Genotype)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.size = 1.0, linewidth = 0.4) +
  facet_wrap(~Genotype, nrow = 1, labeller = as_labeller(c(
    "RBT_AA" = "RBT[AA]", 
    "RBT_AR" = "RBT[AR]", 
    "RBT_RR" = "RBT[RR]"
  ), default = label_parsed)) +
  scale_fill_manual(values = geno_colors, guide = "none") +
  labs(y = "Frequency", x = "Habitat Type") +
  theme_molecol() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# Save plots
ggsave("Chr05_habitat_stacked_bar.pdf", p1, width = 8, height = 5, dpi = 300)
ggsave("Chr05_habitat_stacked_bar.png", p1, width = 8, height = 5, dpi = 300)

ggsave("Chr05_habitat_boxplots.pdf", p2, width = 10, height = 4, dpi = 300)
ggsave("Chr05_habitat_boxplots.png", p2, width = 10, height = 4, dpi = 300)

message("Plots saved with 'Chr05 Genotype' legend and plain axis titles.")