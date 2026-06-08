library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)

# Function to create one plot for a given LG
make_lg_plot <- function(bed, pop_map, global_map, lg_to_plot, out_name) {
  # Join to global map
  renamed_map <- inner_join(pop_map, global_map %>% select(name, Global_LG, Global_Cm_order), by = "name")
  
  merged <- inner_join(bed, renamed_map, by = "name") %>%
    filter(Global_LG == lg_to_plot)
  
  if (nrow(merged) == 0) return(NULL)
  
  # Convert physical position to Mb
  merged <- merged %>% mutate(pos_mb = bp / 1e6)
  
  # Build simplified plot
  p <- ggplot(merged, aes(x = pos_mb, y = Global_Cm_order)) +
    geom_point(aes(color = strand), size = 1, alpha = 0.8) +
    scale_x_continuous(name = "Physical position (Mb))") +
    scale_y_continuous(name = "Genetic position (cM from global map)") +
    scale_color_manual(values = c("+" = "blue", "-" = "red")) +
    labs(title = paste0(out_name, ": ", lg_to_plot, " markers aligned to Omy", gsub("LG", "", lg_to_plot))) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      text = element_text(size = 14)
    )
  
  return(p)
}

# Wrapper for a population
plot_all_lgs <- function(bed_file, map_file, global_map, out_prefix) {
  bed <- read_tsv(bed_file, col_names = FALSE)
  colnames(bed) <- c("chr", "start", "end", "name", "score", "strand")
  bed <- bed %>% mutate(bp = (start + end) / 2)
  
  pop_map <- read_tsv(map_file)  # must contain 'name', 'LG', and 'cM' columns
  
  lg_list <- paste0("LG", sprintf("%02d", 1:29))
  
  plots <- lapply(lg_list, function(lg) {
    make_lg_plot(bed, pop_map, global_map, lg, out_prefix)
  })
  
  plots <- Filter(Negate(is.null), plots)
  
  if (length(plots) > 0) {
    pdf(paste0(out_prefix, "_all_LG_alignments.pdf"), width = 10, height = 6)
    for (p in plots) print(p)
    dev.off()
  } else {
    warning(paste("No markers aligned for", out_prefix))
  }
}

# Load global map with correct column names
global_map <- read_tsv("global_map_clean.tsv", col_names = FALSE)
colnames(global_map) <- c("CatID", "name", "Global_LG", "Global_Cm_order", "Arlee_LG", "Arlee_pos")

# Run for each population
plot_all_lgs("EFYaak.bed",   "EFYaak_linkage_map.tsv",   global_map, "EFYaak")
plot_all_lgs("Gerrard.bed",  "Gerrard_linkage_map.tsv",  global_map, "Gerrard")
plot_all_lgs("Wolf.bed",     "Wolf_linkage_map.tsv",     global_map, "Wolf")

# Prepare a linkage-style map for global markers (run this once)
write_tsv(
  global_map %>% select(name, LG = Global_LG, cM = Global_Cm_order),
  "global_linkage_map.tsv"
)

# Then run:
plot_all_lgs("Global.bed", "global_linkage_map.tsv", global_map, "Global")
