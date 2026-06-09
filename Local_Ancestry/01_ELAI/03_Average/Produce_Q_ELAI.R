# ==== Load Libraries ====
library(data.table)
library(dplyr)
library(tibble)
library(tidyr)

# ==== User Config ====
ancestry_labels <- c("CRT_AA", "CRT_RR", "IRT", "WCT", "YCT")
sample_list <- readLines("sample_order.txt")  # one sample per line
ps21_dir <- "."  # folder with OmyXX.ps21.txt files
chrom_lengths <- fread("chrom_lengths.txt")   # columns: LG length

S <- length(ancestry_labels)
N <- length(sample_list)

# ==== Helper Function: Compute per-chromosome ancestry ====
compute_chrom_ancestry <- function(ps21_file) {
  yy <- scan(ps21_file, what=double(), quiet=TRUE)
  
  # Sanity check: total length divisible by S*N?
  total_len <- length(yy)
  if (total_len %% (S * N) != 0) {
    warning(ps21_file, ": length ", total_len,
            " is NOT divisible by S*N = ", S*N, 
            " → skipping this file")
    return(NULL)
  }
  
  M <- total_len / (S * N)
  dim(yy) <- c(S, M, N)
  
  mat <- matrix(0, nrow=N, ncol=S,
                dimnames=list(sample_list, ancestry_labels))
  
  for (ind in seq_len(N)) {
    ancestry_sums <- rowSums(yy[,,ind])
    mat[ind, ] <- ancestry_sums / sum(ancestry_sums)
  }
  
  as.data.frame(mat) %>% rownames_to_column("Individual")
}

# ==== Main ====
all_perchrom <- list()
ps21_files <- list.files(ps21_dir, pattern="\\.ps21\\.txt$", full.names=TRUE)
chromosomes <- gsub("\\.ps21\\.txt$", "", basename(ps21_files))

for (i in seq_along(ps21_files)) {
  LG <- chromosomes[i]
  cat("Processing", LG, "...\n")
  
  ps21_file <- ps21_files[i]
  if (!file.exists(ps21_file)) {
    warning("Skipping ", LG, " (missing file)")
    next
  }
  
  df <- compute_chrom_ancestry(ps21_file)
  if (is.null(df)) next  # skip if sample check failed
  
  df$Chromosome <- LG
  df$Chr_length <- chrom_lengths$length[chrom_lengths$LG == LG]
  
  all_perchrom[[LG]] <- df
}

# Combine all chromosomes
if (length(all_perchrom) == 0) stop("No valid chromosomes to process!")
all_data <- bind_rows(all_perchrom)

# ==== Output per-chromosome files ====
# Long format
per_chrom_long <- all_data %>%
  pivot_longer(cols = all_of(ancestry_labels),
               names_to = "Ancestry",
               values_to = "Fraction")
write.table(per_chrom_long, "per_chrom_ancestry_long.txt", sep="\t",
            quote=FALSE, row.names=FALSE)

# Wide format
per_chrom_wide <- all_data %>%
  select(Individual, Chromosome, all_of(ancestry_labels))
write.table(per_chrom_wide, "per_chrom_ancestry_wide.txt", sep="\t",
            quote=FALSE, row.names=FALSE)

# ==== Genome-wide Q values (weighted by chromosome length) ====
Q_genomewide <- all_data %>%
  group_by(Individual) %>%
  summarise(across(all_of(ancestry_labels),
                   ~ sum(.x * Chr_length) / sum(Chr_length)),
            .groups="drop")
write.table(Q_genomewide, "genomewide_Q_weighted.txt", sep="\t",
            quote=FALSE, row.names=FALSE)
