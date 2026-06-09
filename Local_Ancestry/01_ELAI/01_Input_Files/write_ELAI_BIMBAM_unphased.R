### === USER INPUTS =======================================================

# Path to full, unphased or phased, multi-population VCF (must be bgzipped or plain VCF)
vcf_path <- "/mnt/ceph/well0766/UMT_Redband/04_vcf/01_STACKS/01_global_ancestry/02_Concat_VCF/02_remove_CCT/all_RADCap_10x.filtered.noCCT.recode.vcf.gz"

# Directory containing sample list files (e.g., CRT.txt, IRT.txt, etc.)
sample_list_dir <- "/mnt/ceph/well0766/UMT_Redband/06_local_ancestry/01_ELAI/01_all_chrom_unphased/01_Input_Files/sets_rapture"

# Output directory (will be organized by population, chromosome, etc.)
output_dir <- "/mnt/ceph/well0766/UMT_Redband/06_local_ancestry/01_ELAI/01_all_chrom_unphased/01_Input_Files"

### =======================================================================


### === LIBRARIES AND FUNCTIONS ==========================================

# Load required library
library(vcfR)

# Read sample list from text file
read_samples <- function(file) {
  readLines(file)
}

# Subset a vcfR object to a specific set of sample names
subset_vcf <- function(vcf, samples) {
  sample_indices <- which(colnames(vcf@gt)[-1] %in% samples) + 1  # +1 to include FORMAT column
  if (length(sample_indices) == 0) stop("No samples from list found in VCF.")
  vcf@gt <- vcf@gt[, c(1, sample_indices)]
  return(vcf)
}

# Convert a vcfR object into ELAI-compatible BIMBAM format (per chromosome)
writeELAI <- function(vcf, group_name, output_base) {
  chroms <- unique(vcf@fix[, "CHROM"])
  
  for (chrom in chroms) {
    cat("🧬 Processing group:", group_name, "| Chromosome:", chrom, "\n")
    
    idx <- which(vcf@fix[, "CHROM"] == chrom)
    if (length(idx) == 0) next
    
    fix <- vcf@fix[idx, ]
    gt <- vcf@gt[idx, ]
    samples <- colnames(gt)[-1]
    n_samples <- length(samples)
    
    gt_matrix <- extract.gt(vcf[idx, ], element = "GT", as.numeric = FALSE, return.alleles = FALSE)
    
    ref <- fix[, "REF"]
    alt <- fix[, "ALT"]
    ids <- fix[, "ID"]
    pos <- fix[, "POS"]
    marker_ids <- ifelse(ids == ".", paste0("SNP_", chrom, "_", pos), ids)
    
    gmat <- matrix(NA, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
    
    for (i in 1:nrow(gt_matrix)) {
      for (j in 1:ncol(gt_matrix)) {
        g <- gt_matrix[i, j]
        if (is.na(g) || g %in% c(".", "./.", ".|.")) {
          gmat[i, j] <- "??"
        } else {
          a <- unlist(strsplit(g, "[|/]"))
          if (length(a) < 2 || "." %in% a) {
            gmat[i, j] <- "??"
          } else {
            a1 <- ifelse(a[1] == "0", ref[i],
                         ifelse(a[1] == "1", alt[i], "?"))
            a2 <- ifelse(a[2] == "0", ref[i],
                         ifelse(a[2] == "1", alt[i], "?"))
            gmat[i, j] <- paste0(a1, a2)
          }
        }
      }
    }
    
    # Output paths
    geno_dir <- file.path(output_base, group_name, "genos")
    loci_dir <- file.path(output_base, group_name, "loci")
    dir.create(geno_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(loci_dir, recursive = TRUE, showWarnings = FALSE)
    
    geno_file <- file.path(geno_dir, paste0(chrom, ".txt"))
    loci_file <- file.path(loci_dir, paste0(chrom, ".txt"))
    
    # === Write Genotype File (no "=" for unphased) ===
    cat(ncol(gmat), file = geno_file, sep = "\n")
    cat(nrow(gmat), file = geno_file, append = TRUE, sep = "\n")
    cat(paste("IND", paste(samples, collapse = ", "), sep = ", "), file = geno_file, append = TRUE, sep = "\n")
    
    for (i in 1:nrow(gmat)) {
      line <- paste(marker_ids[i], paste(gmat[i, ], collapse = ", "), sep = ", ")
      cat(line, file = geno_file, append = TRUE, sep = "\n")
    }
    
    # === Write Marker Position File ===
    loci_data <- data.frame(marker_ids, pos, chrom)
    write.table(loci_data, loci_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ", ")
  }
  
  cat("Finished group:", group_name, "\n\n")
}

### =======================================================================


### === MAIN SCRIPT =======================================================

cat("Loading full VCF...\n")
vcf_all <- read.vcfR(vcf_path)

sample_files <- list.files(sample_list_dir, pattern = "\\.txt$", full.names = TRUE)

for (sfile in sample_files) {
  group <- gsub("\\.txt$", "", basename(sfile))
  samples <- read_samples(sfile)
  cat("Group:", group, "| N samples:", length(samples), "\n")
  
  vcf_group <- subset_vcf(vcf_all, samples)
  writeELAI(vcf_group, group, output_dir)
}

cat("All groups processed and ELAI inputs written to:", output_dir, "\n")

