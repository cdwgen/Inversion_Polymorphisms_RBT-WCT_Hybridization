# Load input
df <- read.table("input_ancestry.tsv", header = TRUE, sep = "\t")

# If pRBT isn't in the table, calculate it as sum of CRT and IRT
if (!"pRBT" %in% colnames(df)) {
  df$pRBT <- df$pCRT + df$pIRT
}

# Avoid division by zero when pCRT + pIRT = 0
sum_CRT_IRT <- df$pCRT + df$pIRT
sum_CRT_IRT[sum_CRT_IRT == 0] <- 1  # avoid divide-by-zero

# Adjust within pRBT
df$adj_pCRT <- df$pRBT * (df$pCRT / sum_CRT_IRT)
df$adj_pIRT <- df$pRBT * (df$pIRT / sum_CRT_IRT)

# Compute total ancestry (excluding original pRBT)
df$total <- df$adj_pCRT + df$adj_pIRT + df$pWCT + df$pYCT

# Normalize
df$norm_pCRT <- df$adj_pCRT / df$total
df$norm_pIRT <- df$adj_pIRT / df$total
df$norm_pWCT <- df$pWCT / df$total
df$norm_pYCT <- df$pYCT / df$total

# Select relevant columns
df_out <- df[, c("Sample_ID", "norm_pCRT", "norm_pIRT", "norm_pWCT", "norm_pYCT")]

# Save output
write.table(df_out, "diagnostics_normalized.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
