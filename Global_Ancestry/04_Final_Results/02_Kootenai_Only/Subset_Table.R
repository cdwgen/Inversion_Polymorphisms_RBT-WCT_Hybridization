# ============================================================
# Subset & Reorder Table by Sample List (SAFE VERSION)
# RStudio-friendly version (no command line args)
# ============================================================

  library(readr)
  library(dplyr)

# ------------------------------------------------------------
# 0. USER: Define paths here
# ------------------------------------------------------------
big_file   <- "../01_all_samples_all_methods/All_Samples_All_Methods.tsv"
order_file <- "desired_order.tsv"
out_file   <- "Kootenai_only_all_methods.tsv"

cat("\n=============================================\n")
cat(" Subset + Order Script (RStudio Version)\n")
cat("=============================================\n")
cat("Input main table:  ", big_file, "\n")
cat("Desired order file:", order_file, "\n")
cat("Output file:       ", out_file, "\n\n")

# ------------------------------------------------------------
# 1. Read input data
# ------------------------------------------------------------
cat("Reading large table...\n")
big <- read_tsv(big_file, col_types = cols(.default = "c"))

if (!"Sample_ID" %in% names(big)) {
  stop("ERROR: The main table does NOT contain a 'Sample_ID' column.")
}

cat("Reading order list...\n")
order_list <- read_tsv(order_file, col_types = cols(.default = "c"))

if (!"Sample_ID" %in% names(order_list)) {
  stop("ERROR: The order list file must contain a column named 'Sample_ID'.")
}

order_vec <- order_list$Sample_ID

# ------------------------------------------------------------
# 2. Safety checks for mismatches
# ------------------------------------------------------------
cat("\nChecking for mismatched sample IDs...\n")

missing_in_big <- setdiff(order_vec, big$Sample_ID)
extra_in_big   <- setdiff(big$Sample_ID, order_vec)

if (length(missing_in_big) > 0) {
  cat("\n❌ ERROR: These sample IDs are in the order list but NOT in the big table:\n")
  print(missing_in_big)
  stop("\nFix your sample names before running again.\n")
}

if (length(extra_in_big) > 0) {
  cat("\n⚠ WARNING: These sample IDs are in the big table but not in the order list:\n")
  print(extra_in_big)
  cat("\n(They will be excluded from the output.)\n")
}

# ------------------------------------------------------------
# 3. Subset + reorder safely
# ------------------------------------------------------------
cat("\nSubsetting + ordering rows...\n")

subset_df <- big %>%
  filter(Sample_ID %in% order_vec) %>%
  slice(match(order_vec, Sample_ID))

if (nrow(subset_df) != length(order_vec)) {
  stop("ERROR: Output row count does NOT match order list length. Something went wrong.")
}

# ------------------------------------------------------------
# 4. Write output
# ------------------------------------------------------------
cat("\nWriting output file:", out_file, "\n")
write_tsv(subset_df, out_file)

cat("\n✅ DONE! Output successfully written.\n")
cat("=============================================\n\n")
