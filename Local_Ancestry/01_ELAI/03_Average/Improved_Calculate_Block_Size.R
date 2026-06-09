#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tibble)

# =========================
# USER CONFIG
# =========================

ancestry_labels <- c("CRT", "IRT", "WCT", "YCT")
sample_list <- readLines("sample_order.txt")

ps21_dir <- "."
S_original <- 5
N <- length(sample_list)

# =========================
# CORE THRESHOLDS (ANCESTRY DOSAGES FOR HET CALLING)
# =========================

homo_thresh <- 1.85        # Strong single ancestry (2 copies, ~2.0)
uncertain_thresh <- 1.30   # Low-confidence floor for TOTAL signal

# Heterozygote boundaries (1 copy each, ~1.0 + ~1.0)
het_low <- 0.85            # Lower bound for a single hybrid allele
het_high <- 1.15           # Upper bound for a single hybrid allele
het_sum_thresh <- 1.85     # Combined top-two sum must be as strong as a homozygote

# =========================
# BLOCK FILTERING
# =========================

merge_max_snps <- 5
merge_max_bp   <- 100000
min_snps       <- 10
min_block_bp   <- 200000   # Set to 200kb to protect older, real backcross tracts

# =========================
# CALLING FUNCTION (FIXED LOGIC)
# =========================

call_states_vectorized <- function(dosage) {

  M <- ncol(dosage)
  t_dosage <- t(dosage)

  # Get #1 ancestry
  max_idx <- max.col(t_dosage, ties.method = "first")
  max_vals <- dosage[cbind(max_idx, seq_len(M))]

  # Get #2 ancestry
  dosage_temp <- dosage
  dosage_temp[cbind(max_idx, seq_len(M))] <- -1

  sec_idx <- max.col(t(dosage_temp), ties.method = "first")
  sec_vals <- dosage[cbind(sec_idx, seq_len(M))]

  # Calculate Total Top Signal
  sum12 <- max_vals + sec_vals

  states <- rep("Uncertain", M)

  # -----------------------------
  # STEP 1: UNCERTAIN (WEAK SIGNAL)
  # -----------------------------
  # FIXED: Now checks if the COMBINED top two signal is too weak
  is_uncertain <- sum12 < uncertain_thresh

  # -----------------------------
  # STEP 2: HOMOZYGOUS / DOMINANT
  # -----------------------------
  is_homo <- (!is_uncertain) & (max_vals >= homo_thresh)

  states[is_homo] <- ancestry_labels[max_idx[is_homo]]

  # -----------------------------
  # STEP 3: TRUE MIXTURE (THE 1.0 + 1.0 TEST)
  # -----------------------------
  # Both top ancestries must be near 1.0, AND their sum must be strong.
  is_het <- (!is_uncertain) &
            (!is_homo) &
            (max_vals >= het_low & max_vals <= het_high) &
            (sec_vals >= het_low & sec_vals <= het_high) &
            (sum12 >= het_sum_thresh)

  if (any(is_het)) {
    l1 <- ancestry_labels[max_idx[is_het]]
    l2 <- ancestry_labels[sec_idx[is_het]]

    states[is_het] <- paste(pmin(l1, l2), pmax(l1, l2), sep = "/")
  }

  return(states)
}

# =========================
# MERGE SHORT TRACTS (PROTECTED HET SEGMENTS)
# =========================

merge_short_blocks <- function(df) {
  df <- as.data.frame(df)

  repeat {
    n <- nrow(df)
    if (n <= 2) break

    merged <- FALSE

    for (i in 2:(n-1)) {

      left  <- df[i-1,]
      mid   <- df[i,]
      right <- df[i+1,]

      # Skip uncertain blocks
      if (mid$Ancestry == "Uncertain") next

      # PROTECTION: do NOT merge confident heterozygous tracts
      is_confident_het <- grepl("/", mid$Ancestry)

      if (mid$Ancestry != left$Ancestry &&
          left$Ancestry == right$Ancestry &&
          left$Ancestry != "Uncertain" &&
          !is_confident_het &&
          (mid$n_SNPs <= merge_max_snps ||
           mid$Block_length_bp <= merge_max_bp)) {

        df[i-1, "SNP_end"] <- right$SNP_end
        df[i-1, "End_bp"]  <- right$End_bp
        df[i-1, "n_SNPs"]  <- left$n_SNPs + mid$n_SNPs + right$n_SNPs
        df[i-1, "Block_length_bp"] <- df[i-1, "End_bp"] - df[i-1, "Start_bp"] + 1

        df <- df[-c(i, i+1), ]
        merged <- TRUE
        break
      }
    }

    if (!merged) break
  }

  return(as_tibble(df))
}

# =========================
# BLOCK CALLING
# =========================

call_blocks <- function(ps21_file, snpinfo_file, LG) {

  yy <- scan(ps21_file, quiet = TRUE)
  snpinfo <- fread(snpinfo_file)

  M <- nrow(snpinfo)

  if (length(yy) != S_original * M * N) {
    stop(paste("Dimension mismatch:", LG))
  }

  dim(yy) <- c(S_original, M, N)
  pos_vec <- snpinfo$pos

  results <- list()

  for (ind in seq_len(N)) {

    raw_dosage <- yy[,,ind]

    # collapse CRT_AA + CRT_RR
    crt_combined <- raw_dosage[1, ] + raw_dosage[2, ]
    dosage <- rbind(crt_combined, raw_dosage[3:5, ])

    states <- call_states_vectorized(dosage)

    rle_states <- rle(states)
    ends <- cumsum(rle_states$lengths)
    starts <- c(1, head(ends, -1) + 1)

    df <- tibble(
      Individual = sample_list[ind],
      LG = LG,
      SNP_start = starts,
      SNP_end = ends,
      Start_bp = pos_vec[starts],
      End_bp = pos_vec[ends],
      Ancestry = rle_states$values,
      n_SNPs = rle_states$lengths
    ) %>%
      mutate(Block_length_bp = End_bp - Start_bp + 1)

    df <- merge_short_blocks(df)

    df <- df %>%
      filter(n_SNPs >= min_snps) %>%
      filter(Block_length_bp >= min_block_bp) %>%
      mutate(
        Block_length_cM = Block_length_bp / 1e6,
        Block_length_M  = Block_length_bp / 1e8
      )

    results[[ind]] <- df
  }

  bind_rows(results)
}

# =========================
# MAIN LOOP
# =========================

ps21_files <- list.files(ps21_dir, pattern="\\.ps21\\.txt$", full.names=TRUE)
chromosomes <- gsub("\\.ps21\\.txt$", "", basename(ps21_files))

all_blocks <- list()

for (i in seq_along(ps21_files)) {

  LG <- chromosomes[i]
  ps21_file <- ps21_files[i]
  snpinfo_file <- file.path(ps21_dir, paste0(LG, ".snpinfo.txt"))

  if (!file.exists(snpinfo_file)) {
    warning("Skipping ", LG)
    next
  }

  cat("Processing", LG, "...\n")
  all_blocks[[LG]] <- call_blocks(ps21_file, snpinfo_file, LG)
}

all_blocks_df <- bind_rows(all_blocks)

# =========================
# OUTPUT
# =========================

write.table(all_blocks_df, "ancestry_tracts.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

block_summary <- all_blocks_df %>%
  filter(Ancestry != "Uncertain") %>%
  group_by(Individual, Ancestry) %>%
  summarise(
    mean_block_bp = mean(Block_length_bp),
    median_block_bp = median(Block_length_bp),
    mean_block_M = mean(Block_length_M),
    n_blocks = n(),
    generations_since_admix = 1 / mean(Block_length_M),
    .groups="drop"
  )

write.table(block_summary, "ancestry_block_summary.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

cat("\nDone. Output files generated.\n")