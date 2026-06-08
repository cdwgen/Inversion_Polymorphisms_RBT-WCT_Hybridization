#!/usr/bin/env Rscript

# Load necessary libraries
library(vcfR)
library(adegenet)
library(ade4)
library(ggplot2)
library(RColorBrewer)
library(htmlwidgets)
library(plotly)
library(dplyr)

#### PCA from VCF and popmap ####

# Read the VCF file
vcf <- read.vcfR("/mnt/ceph/well0766/UMT_Redband/04_vcf/01_STACKS/01_global_ancestry/02_Concat_VCF/03_remove_dups/all_RADCap_10x.filtered.noCCT.nodups.recode.vcf")

# Fix duplicate SNP IDs (Edit 1)
dup_ids <- duplicated(vcf@fix[, "ID"])
if (any(dup_ids)) {
  message("Found duplicated SNP IDs, renaming them...")
  vcf@fix[dup_ids, "ID"] <- paste0(vcf@fix[dup_ids, "ID"], "_", seq_len(sum(dup_ids)))
}

# Convert VCF data to genind object
genind_data <- vcfR2genind(vcf)

# Read popmap file (three columns: sample_id, population, group)
popmap <- read.table("no_dup_popmap.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(popmap) <- c("sample_id", "population", "group")

# Ensure individual names match order
ind_names <- indNames(genind_data)
popmap <- popmap[match(ind_names, popmap$sample_id), ]

# Check for mismatches
if (any(is.na(popmap$population))) {
  stop("Mismatch between VCF sample IDs and popmap. Please check your files.")
}

# Assign population to genind object
pop(genind_data) <- as.factor(popmap$population)

# Impute missing values
genind_tab <- tab(genind_data)
genind_tab[is.na(genind_tab)] <- colMeans(genind_tab, na.rm = TRUE)[col(genind_tab)][is.na(genind_tab)]
genind_data@tab <- genind_tab

# Run PCA
pca_result <- dudi.pca(genind_tab, center = TRUE, scale = FALSE, scannf = FALSE, nf = 3)

# % variance explained
percent_variation <- round((pca_result$eig / sum(pca_result$eig)) * 100, 2)

# Prepare dataframe
pca_df <- data.frame(pca_result$li)
pca_df$sample_id <- ind_names
pca_df$population <- popmap$population
pca_df$group <- popmap$group

# Set color palette
pop_colors <- setNames(colorRampPalette(brewer.pal(8, "Set1"))(length(unique(pca_df$population))),
                       unique(pca_df$population))

# Static 2D ggplot
p <- ggplot(pca_df, aes(x = Axis1, y = Axis2, color = population)) +
  geom_point(size = 3) +
  scale_color_manual(values = pop_colors) +
  theme_classic() +
  labs(title = "All Samples (RADCap 2025)",
       x = paste0("PC1 (", percent_variation[1], "%)"),
       y = paste0("PC2 (", percent_variation[2], "%)"),
       color = "Population") +
  theme(legend.position = "right")

# Show and save static plot
print(p)
ggsave("PCA_All.png", plot = p, device = "png", width = 10, height = 10, units = "in", dpi = 500)

# Save PCA coordinates
write.csv(pca_df, "PCA_All_coordinates.csv", row.names = FALSE)

### Interactive 2D Plot ###
p_interactive_2D <- plot_ly(
  data = pca_df,
  x = ~Axis1,
  y = ~Axis2,
  type = 'scatter',
  mode = 'markers',
  color = ~population,
  text = ~paste("Sample:", sample_id,
                "<br>Population:", population,
                "<br>Group:", group),
  hoverinfo = 'text',
  marker = list(size = 6)
) %>%
  layout(
    title = "Interactive PCA Plot (2D)",
    xaxis = list(title = paste0("PC1 (", percent_variation[1], "%)")),
    yaxis = list(title = paste0("PC2 (", percent_variation[2], "%)"))
  )

# Save 2D plot
saveWidget(p_interactive_2D, file = "PCA_All_Interactive_2D.html")

### Interactive 3D Plot ###
p_interactive_3D <- plot_ly(
  data = pca_df,
  x = ~Axis1,
  y = ~Axis2,
  z = ~Axis3,
  type = 'scatter3d',
  mode = 'markers',
  color = ~population,
  text = ~paste("Sample:", sample_id,
                "<br>Population:", population,
                "<br>Group:", group),
  hoverinfo = 'text',
  marker = list(size = 4)
) %>%
  layout(
    title = "Interactive PCA Plot (3D)",
    scene = list(
      xaxis = list(title = paste0("PC1 (", percent_variation[1], "%)")),
      yaxis = list(title = paste0("PC2 (", percent_variation[2], "%)")),
      zaxis = list(title = paste0("PC3 (", percent_variation[3], "%)"))
    )
  )

# Save 3D plot
saveWidget(p_interactive_3D, file = "PCA_All_Interactive_3D.html")