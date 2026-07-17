library(readxl)
library(tidyverse)
library(maps)
library(pals)

# Editing notes:
# - Change setwd() if this folder moves or if you run the script on another computer.
# - Change metadata_file if the Excel file name or location changes.
# - Change repeated_shapes to use a different sequence of point shapes.
# - county_colors uses pals::glasbey(), a discrete palette for distinguishing many groups.
# - Change label_y_offset to move the "n = ?" labels up or down.
# - Change index_circle_size_range to make the indexed-map circles smaller/larger.
# - Change lon_padding and lat_padding to zoom the map in or out.
# - Change the ggsave width/height/dpi values near the bottom for figure sizes.
# - This script uses index_new_county_state from Quercus_metadata_with_index.xlsx
#   so maps match the generated metadata file.

setwd("/Users/kuowenhsi/Library/CloudStorage/OneDrive-MissouriBotanicalGarden/General - IMLS National Leadership Grant 2023/Manuscript/Quercus")

metadata_file <- "data/metadata/Quercus_metadata_with_index.xlsx"
metadata_sheet <- "metadata"
repeated_shapes <- c(21, 22, 23, 24)
label_y_offset <- 0.12
index_circle_size_range <- c(6, 10)
lon_padding <- 0.75
lat_padding <- 0.75

metadata <- read_excel(metadata_file, sheet = metadata_sheet)

required_columns <- c("index", "index_new_county_state", "State", "Latitude", "Longitude")
missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}

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

contrast_text_color <- function(fill_colors) {
  rgb_values <- col2rgb(fill_colors) / 255
  luminance <- 0.2126 * rgb_values[1, ] + 0.7152 * rgb_values[2, ] + 0.0722 * rgb_values[3, ]
  ifelse(luminance < 0.45, "white", "black")
}

county_colors <- glasbey(length(county_levels)) %>%
  set_names(county_levels)

county_shapes <- repeated_shapes %>%
  rep(length.out = length(county_levels)) %>%
  set_names(county_levels)

county_counts <- plot_data %>%
  group_by(county_group) %>%
  summarise(
    Longitude = mean(Longitude),
    Latitude = mean(Latitude) - label_y_offset,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", n))

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
    n_label_latitude = Latitude - (label_y_offset + circle_size * 0.025),
    index_label_color = contrast_text_color(county_colors[as.character(county_group)])
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

quercus_map <- ggplot() +
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
    data = plot_data,
    aes(x = Longitude, y = Latitude, fill = county_group, shape = county_group),
    color = "black",
    size = 3.2,
    stroke = 0.35,
    alpha = 0.9
  ) +
  geom_text(
    data = county_counts,
    aes(x = Longitude, y = Latitude, label = label),
    size = 2.7,
    fontface = "bold",
    color = "black",
    check_overlap = FALSE
  ) +
  scale_shape_manual(values = county_shapes) +
  scale_fill_manual(values = county_colors) +
  coord_quickmap(xlim = x_limits, ylim = y_limits, expand = FALSE) +
  labs(
    title = "Quercus Sampling Localities",
    x = "Longitude",
    y = "Latitude",
    fill = "County",
    shape = "County"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_line(color = "grey90", linewidth = 0.25),
    legend.position = "right",
    legend.key.height = unit(0.42, "cm"),
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 15)
  )

ggsave("figures/metadata_maps/Quercus_metadata_county_map.png", quercus_map, width = 11, height = 5.5, dpi = 300)
ggsave("figures/metadata_maps/Quercus_metadata_county_map.pdf", quercus_map, width = 11, height = 5.5)

quercus_index_map <- ggplot() +
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
    aes(x = Longitude, y = Latitude, label = index, color = index_label_color),
    size = 3.2,
    fontface = "bold"
  ) +
  geom_text(
    data = county_index_data,
    aes(x = Longitude, y = n_label_latitude, label = n_label),
    size = 2.8,
    fontface = "bold",
    color = "black",
    check_overlap = FALSE
  ) +
  scale_fill_manual(values = county_colors) +
  scale_color_identity(guide = "none") +
  scale_size(
    range = index_circle_size_range,
    breaks = c(1, 5, 10, 20, 40)
  ) +
  coord_quickmap(xlim = x_limits, ylim = y_limits, expand = FALSE) +
  labs(
    title = "Quercus County Sample Size Index",
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

ggsave("figures/metadata_maps/Quercus_metadata_county_index_map.png", quercus_index_map, width = 11, height = 5.5, dpi = 300)
ggsave("figures/metadata_maps/Quercus_metadata_county_index_map.pdf", quercus_index_map, width = 11, height = 5.5)

message("Saved Quercus_metadata_county_map.png and Quercus_metadata_county_map.pdf")
message("Saved Quercus_metadata_county_index_map.png and Quercus_metadata_county_index_map.pdf")
message("Records plotted: ", nrow(plot_data))
message("Records skipped for missing County/Latitude/Longitude: ", nrow(metadata) - nrow(plot_data))
