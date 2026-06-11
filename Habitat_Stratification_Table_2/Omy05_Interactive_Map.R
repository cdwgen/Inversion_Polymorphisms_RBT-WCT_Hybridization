#!/usr/bin/env Rscript
  
  library(leaflet)
library(ggplot2)
library(htmlwidgets)
library(dplyr)
library(base64enc)

# --------------------------------------------------------------
# Read & clean data
# --------------------------------------------------------------
df <- read.table("Omy05_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df$Lat  <- as.numeric(df$Lat)
df$Long <- as.numeric(df$Long)
df$Long[df$Long > 0] <- -df$Long[df$Long > 0]

# --------------------------------------------------------------
# Make Habitat_Type human-readable
# --------------------------------------------------------------
habitat_lookup <- c(
  "River" = "River",
  "Tributary" = "Tributary",
  "Lake" = "Lake",
  "Valley_Floor" = "Valley Floor"
)

df$Habitat_Display <- habitat_lookup[df$Habitat_Type]

# --------------------------------------------------------------
# Pie chart colors & labels
# --------------------------------------------------------------
pie_colors <- c(
  "RBT_AA" = "#cb6651",
  "RBT_AR" = "#54a451",
  "RBT_RR" = "#3f908a"
)

popup_labels <- c(
  "RBT<sub>AA</sub>",
  "RBT<sub>AR</sub>",
  "RBT<sub>RR</sub>"
)

# --------------------------------------------------------------
# Step 1: Generate pie PNGs as base64 strings
# --------------------------------------------------------------
pie_base64 <- vector("character", nrow(df))

for(i in 1:nrow(df)){
  pie_data <- df[i, c("RBT_AA", "RBT_AR", "RBT_RR")]
  pie_df <- data.frame(
    genotype = names(pie_data),
    proportion = as.numeric(pie_data)
  )
  
  tmpfile <- tempfile(fileext = ".png")
  
  p <- ggplot(pie_df, aes(x = "", y = proportion, fill = genotype)) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = pie_colors) +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(tmpfile, plot = p, width = 2, height = 2, dpi = 144)
  
  pie_base64[i] <- paste0("data:image/png;base64,", base64encode(tmpfile))
}

# --------------------------------------------------------------
# Step 2: Build leaflet map
# --------------------------------------------------------------
icons_list <- icons(
  iconUrl = pie_base64,
  iconWidth = 50,
  iconHeight = 50
)

m <- leaflet(df) %>%
  addTiles() %>%
  addMarkers(
    lng = df$Long,
    lat = df$Lat,
    icon = icons_list,
    popup = paste0(
      "<b>Population:</b> ", df$Population, "<br>",
      "<b>Group:</b> ", df$Group, "<br>",
      "<b>Year:</b> ", df$Year, "<br>",
      "<b>Habitat type:</b> ", df$Habitat_Display, "<br>",
      "<b>Associated life history:</b> ", df$Associated_Life_History, "<br>",
      "<b>N:</b> ", df$N, "<br>",
      "<b>Genotype proportions:</b><br>",
      "RBT<sub>AA</sub>: ", df$RBT_AA, "<br>",
      "RBT<sub>AR</sub>: ", df$RBT_AR, "<br>",
      "RBT<sub>RR</sub>: ", df$RBT_RR
    )
  ) %>%
  addLegend(
    position = "bottomright",
    colors = pie_colors,
    labels = popup_labels,
    title = "Genotype proportions",
    opacity = 1
  )

# --------------------------------------------------------------
# Step 3: Save HTML
# --------------------------------------------------------------
m
htmlwidgets::saveWidget(m, "Omy05_genotype_map_subscript.html", selfcontained = TRUE)