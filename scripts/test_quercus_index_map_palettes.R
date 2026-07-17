library(readxl)
library(tidyverse)
library(maps)
library(colorspace)
library(viridisLite)
library(Polychrome)

setwd("/Users/kuowenhsi/Library/CloudStorage/OneDrive-MissouriBotanicalGarden/General - IMLS National Leadership Grant 2023/Manuscript/Quercus")

metadata <- read_excel("data/metadata/Quercus_metadata_with_index.xlsx", sheet = "metadata")

label_y_offset <- 0.12
index_circle_size_range <- c(6, 10)
lon_padding <- 0.75
lat_padding <- 0.75

plot_data <- metadata %>%
  filter(
    !is.na(index),
    !is.na(index_new_county_state),
    !is.na(Latitude),
    !is.na(Longitude)
  ) %>%
  mutate(
    index = as.character(index),
    county_group = factor(index_new_county_state, levels = unique(index_new_county_state[order(index)]))
  )

county_levels <- levels(plot_data$county_group)

county_index_data <- plot_data %>%
  group_by(index, county_group) %>%
  summarise(
    Longitude = mean(Longitude),
    Latitude = mean(Latitude),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(index) %>%
  mutate(
    n_label = paste0("n = ", n),
    circle_size = index_circle_size_range[1] +
      (n - min(n)) / (max(n) - min(n)) *
        (index_circle_size_range[2] - index_circle_size_range[1]),
    n_label_latitude = Latitude - (label_y_offset + circle_size * 0.025)
  )

state_names <- plot_data %>%
  filter(!is.na(State)) %>%
  pull(State) %>%
  str_to_lower() %>%
  unique()

state_map <- map_data("state") %>%
  filter(region %in% state_names)

county_map <- map_data("county") %>%
  filter(region %in% state_names)

x_limits <- range(plot_data$Longitude, na.rm = TRUE) + c(-lon_padding, lon_padding)
y_limits <- range(plot_data$Latitude, na.rm = TRUE) + c(-lat_padding, lat_padding)

palette_tests <- list(
  colorspace_dark3 = qualitative_hcl(length(county_levels), palette = "Dark 3"),
  viridis = viridis(length(county_levels)),
  turbo = turbo(length(county_levels)),
  polychrome = palette36.colors(length(county_levels))
)

make_index_map <- function(palette_name, palette_values) {
  names(palette_values) <- county_levels

  ggplot() +
    geom_polygon(
      data = state_map,
      aes(x = long, y = lat, group = group),
      fill = "grey96",
      color = "grey45",
      linewidth = 0.45
    ) +
    geom_path(
      data = county_map,
      aes(x = long, y = lat, group = group),
      color = "grey78",
      linewidth = 0.2
    ) +
    geom_point(
      data = county_index_data,
      aes(x = Longitude, y = Latitude, fill = county_group, size = n),
      shape = 21,
      color = "black",
      stroke = 0.45,
      alpha = 0.9
    ) +
    geom_text(
      data = county_index_data,
      aes(x = Longitude, y = Latitude, label = index),
      size = 3.2,
      fontface = "bold",
      color = "black"
    ) +
    geom_text(
      data = county_index_data,
      aes(x = Longitude, y = n_label_latitude, label = n_label),
      size = 2.8,
      fontface = "bold",
      color = "black",
      check_overlap = FALSE
    ) +
    scale_fill_manual(values = palette_values) +
    scale_size(
      range = index_circle_size_range,
      breaks = c(1, 5, 10, 20, 40)
    ) +
    coord_quickmap(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    labs(
      title = paste("Quercus County Sample Size Index -", palette_name),
      x = "Longitude",
      y = "Latitude",
      fill = "County",
      size = "Sample size"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_line(color = "grey90", linewidth = 0.25),
      legend.position = "right",
      legend.key.height = unit(0.42, "cm"),
      legend.text = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 15)
    )
}

walk2(
  names(palette_tests),
  palette_tests,
  \(palette_name, palette_values) {
    output_file <- file.path("figures", "metadata_maps", paste0("Quercus_metadata_county_index_map_", palette_name, ".png"))
    ggsave(
      output_file,
      make_index_map(palette_name, palette_values),
      width = 11,
      height = 5.5,
      dpi = 300
    )
    message("Saved ", output_file)
  }
)
