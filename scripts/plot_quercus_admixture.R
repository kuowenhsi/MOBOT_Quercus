library(tidyverse)
library(readxl)
library(sf)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(scatterpie)
library(pals)

# Plot ADMIXTURE CV error, individual ancestry barplots, and map scatterpies.
# Outputs are written to figures/admixture.

project_dir <- "/Users/kuowenhsi/Library/CloudStorage/OneDrive-MissouriBotanicalGarden/General - IMLS National Leadership Grant 2023/Manuscript/Quercus"
setwd(project_dir)

admixture_dir <- file.path(project_dir, "data", "admixture")
output_dir <- file.path(project_dir, "figures", "admixture")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

metadata_file <- file.path(project_dir, "data", "metadata", "Quercus_metadata_with_index.xlsx")
metadata_sheet <- "metadata"
fam_file <- file.path(admixture_dir, "Quercus_afterPhase_MAF001_mLD.fam")
q_prefix <- file.path(admixture_dir, "Quercus_afterPhase_MAF001_mLD")
log_pattern <- "^Quercus_admix_20260522_K([0-9]+)_cv10[.]log$"

k_cv <- 1:20
k_plot <- 2:10

theme_quercus <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = "grey55"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

read_cv_log <- function(path) {
  lines <- readr::read_lines(path, progress = FALSE)
  cv_line <- lines[str_detect(lines, "CV error")]
  if (length(cv_line) == 0) {
    return(tibble(K = NA_integer_, cv_error = NA_real_, file = basename(path)))
  }

  tibble(
    K = as.integer(str_match(cv_line[1], "K=([0-9]+)")[, 2]),
    cv_error = as.numeric(str_match(cv_line[1], ":\\s*([0-9.]+)")[, 2]),
    file = basename(path)
  )
}

read_q <- function(k) {
  q_file <- sprintf("%s.%d.Q", q_prefix, k)
  if (!file.exists(q_file)) {
    stop("Missing Q file: ", q_file)
  }

  cluster_cols <- paste0("cluster_", seq_len(k))
  readr::read_table(
    q_file,
    col_names = cluster_cols,
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(admixture_row = row_number(), .before = 1)
}

prepare_k_data <- function(k, sample_metadata) {
  q <- read_q(k)
  if (nrow(q) != nrow(sample_metadata)) {
    stop("K = ", k, " Q rows (", nrow(q), ") do not match fam rows (", nrow(sample_metadata), ").")
  }

  bind_cols(sample_metadata, q %>% select(-admixture_row)) %>%
    filter(!is.na(State), !is.na(index), !is.na(index_new_county_state), !is.na(Latitude), !is.na(Longitude)) %>%
    mutate(
      K = k,
      state_order = as.integer(factor(State, levels = unique(State[order(state_abbr, State)]))),
      population_label = index_new_county_state,
      population_label = factor(population_label, levels = unique(population_label[order(index, State)])),
      sample_order = factor(Sample_Name, levels = Sample_Name)
    )
}

plot_barplot <- function(k_data, k) {
  cluster_cols <- paste0("cluster_", seq_len(k))

  bar_data <- k_data %>%
    pivot_longer(
      cols = all_of(cluster_cols),
      names_to = "cluster",
      values_to = "ancestry"
    ) %>%
    mutate(cluster = factor(cluster, levels = cluster_cols))

  ggplot(bar_data, aes(x = sample_order, y = ancestry, fill = cluster)) +
    geom_col(width = 1, linewidth = 0) +
    facet_grid(
      . ~ index,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
    scale_fill_brewer(
      palette = "Set3",
      labels = paste("Cluster", seq_len(k)),
      name = sprintf("K = %d", k)
    ) +
    labs(
      title = sprintf("Quercus ADMIXTURE Ancestry Proportions (K = %d)", k),
      x = "Population index",
      y = "Ancestry proportion"
    ) +
    theme_quercus(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(0.05, "lines"),
      strip.text.x = element_text(size = 7),
      legend.key.width = unit(0.8, "cm")
    )
}

compute_map_extent <- function(point_tbl) {
  point_sf <- st_as_sf(point_tbl, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(3857)
  bbox <- st_bbox(point_sf)
  x_span <- bbox$xmax - bbox$xmin
  y_span <- bbox$ymax - bbox$ymin
  x_pad <- x_span * 0.14
  y_pad <- y_span * 0.32
  st_bbox(
    c(
      xmin = as.numeric(bbox["xmin"]) - as.numeric(x_pad),
      ymin = as.numeric(bbox["ymin"]) - as.numeric(y_pad),
      xmax = as.numeric(bbox["xmax"]) + as.numeric(x_pad),
      ymax = as.numeric(bbox["ymax"]) + as.numeric(y_pad)
    ),
    crs = st_crs(point_sf)
  )
}

plot_scatterpie_map <- function(k_data, k, state_sf, county_sf, map_bbox, nudge_tbl, base_pie_radius) {
  cluster_cols <- paste0("cluster_", seq_len(k))

  pie_data <- k_data %>%
    group_by(index, population_label) %>%
    summarise(
      Longitude = mean(Longitude, na.rm = TRUE),
      Latitude = mean(Latitude, na.rm = TRUE),
      state_abbr = paste(sort(unique(state_abbr)), collapse = "/"),
      n = n(),
      across(all_of(cluster_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(index) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(3857) %>%
    cbind(st_coordinates(.)) %>%
    st_drop_geometry() %>%
    as_tibble() %>%
    left_join(nudge_tbl, by = "index") %>%
    mutate(
      dx = replace_na(dx, 0),
      dy = replace_na(dy, 0),
      X_offset = X + dx,
      Y_offset = Y + dy,
      pie_radius = scales::rescale(sqrt(n), to = c(base_pie_radius * 0.75, base_pie_radius * 1.35)),
      label = paste0(index, "\n", state_abbr),
      n_label = paste0("n = ", n),
      label_y = Y_offset - pie_radius * 1.85,
      n_label_y = Y_offset + pie_radius * 1.35
    )

  line_tbl <- pie_data %>%
    select(index, X_anchor = X, Y_anchor = Y, X_offset, Y_offset) %>%
    pivot_longer(
      cols = c(X_anchor, Y_anchor, X_offset, Y_offset),
      names_to = c(".value", "position"),
      names_pattern = "([XY])_(anchor|offset)"
    ) %>%
    mutate(position = factor(position, levels = c("anchor", "offset"))) %>%
    arrange(index, position)

  ggplot() +
    geom_sf(
      data = state_sf,
      fill = "grey96",
      color = "grey45",
      linewidth = 0.45
    ) +
    geom_sf(
      data = county_sf,
      fill = NA,
      color = "grey80",
      linewidth = 0.18
    ) +
    geom_line(
      data = line_tbl,
      aes(x = X, y = Y, group = index),
      color = "grey30",
      linewidth = 0.3,
      alpha = 0.8
    ) +
    geom_point(
      data = pie_data,
      aes(x = X, y = Y),
      color = "black",
      size = 1.4
    ) +
    scatterpie::geom_scatterpie(
      data = pie_data,
      aes(x = X_offset, y = Y_offset, r = pie_radius),
      cols = cluster_cols,
      color = "black",
      linewidth = 0.25,
      alpha = 1
    ) +
    geom_text(
      data = pie_data,
      aes(x = X_offset, y = n_label_y, label = n_label),
      size = 2.6,
      fontface = "bold"
    ) +
    geom_text(
      data = pie_data,
      aes(x = X_offset, y = label_y, label = label),
      size = 2.8,
      fontface = "bold",
      lineheight = 0.85
    ) +
    scale_fill_brewer(
      palette = "Set3",
      labels = paste("Cluster", seq_len(k)),
      name = sprintf("K = %d", k)
    ) +
    annotation_scale(location = "br", width_hint = 0.22) +
    coord_sf(
      xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
      ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
      expand = FALSE,
      crs = 3857
    ) +
    labs(
      title = sprintf("Quercus ADMIXTURE Population Mean Ancestry (K = %d)", k),
      x = NULL,
      y = NULL
    ) +
    theme_quercus(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      legend.key.width = unit(0.8, "cm")
    )
}

metadata <- read_excel(metadata_file, sheet = metadata_sheet) %>%
  mutate(Sample_Name = as.character(Sample_Name))

required_columns <- c(
  "Sample_Name", "State", "state_abbr", "index", "index_new_county_state",
  "Latitude", "Longitude"
)
missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0) {
  stop("Missing required metadata column(s): ", paste(missing_columns, collapse = ", "))
}

fam <- readr::read_table(
  fam_file,
  col_names = c("FID", "Sample_Name", "PID", "MID", "Sex", "Phenotype"),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  mutate(admixture_row = row_number(), Sample_Name = as.character(Sample_Name))

metadata_for_join <- metadata %>%
  filter(!is.na(Sample_Name)) %>%
  select(all_of(required_columns), County, Species, Accession, Planting_Number)

sample_metadata <- fam %>%
  left_join(metadata_for_join, by = "Sample_Name")

unmatched_samples <- sample_metadata %>%
  filter(is.na(State) | is.na(index) | is.na(Latitude) | is.na(Longitude)) %>%
  pull(Sample_Name)

if (length(unmatched_samples) > 0) {
  warning(
    "The following fam sample(s) lack complete metadata and will be omitted from plots: ",
    paste(unmatched_samples, collapse = ", ")
  )
}

sample_metadata %>%
  readr::write_csv(file.path(output_dir, "Quercus_admixture_sample_metadata_join.csv"))

log_files <- list.files(admixture_dir, pattern = log_pattern, full.names = TRUE)
cv_tbl <- purrr::map_dfr(log_files, read_cv_log) %>%
  filter(K %in% k_cv) %>%
  right_join(tibble(K = k_cv), by = "K") %>%
  arrange(K)

missing_cv <- setdiff(k_cv, cv_tbl$K)
missing_cv <- cv_tbl %>%
  filter(is.na(cv_error)) %>%
  pull(K)
if (length(missing_cv) > 0) {
  warning("Missing CV error values for K = ", paste(missing_cv, collapse = ", "))
}

readr::write_csv(cv_tbl, file.path(output_dir, "Quercus_ADMIXTURE_CV_error_K01_K20.csv"))

cv_plot_data <- cv_tbl %>%
  filter(!is.na(cv_error))

cv_caption <- if (length(missing_cv) > 0) {
  paste("CV error not found in log file for K =", paste(missing_cv, collapse = ", "))
} else {
  NULL
}

cv_plot <- ggplot(cv_plot_data, aes(x = K, y = cv_error)) +
  geom_line(color = "grey25", linewidth = 0.5) +
  geom_point(size = 2.4, color = "#0072B2") +
  scale_x_continuous(breaks = k_cv) +
  labs(
    title = "Quercus ADMIXTURE Cross-validation Error",
    x = "K",
    y = "CV error",
    caption = cv_caption
  ) +
  theme_quercus(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "Quercus_ADMIXTURE_CV_error_K01_K20.png"), cv_plot, width = 7.5, height = 4.5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Quercus_ADMIXTURE_CV_error_K01_K20.pdf"), cv_plot, width = 7.5, height = 4.5, bg = "white")

map_metadata <- sample_metadata %>%
  filter(!is.na(State), !is.na(Latitude), !is.na(Longitude))

state_names <- map_metadata %>%
  pull(state_abbr) %>%
  unique()

population_points <- map_metadata %>%
  group_by(index, index_new_county_state) %>%
  summarise(
    Longitude = mean(Longitude, na.rm = TRUE),
    Latitude = mean(Latitude, na.rm = TRUE),
    .groups = "drop"
  )

map_bbox <- compute_map_extent(population_points)
map_bbox_sfc <- st_as_sfc(map_bbox)

state_sf <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(postal %in% state_names) %>%
  st_transform(3857) %>%
  st_make_valid() %>%
  st_intersection(map_bbox_sfc)

county_sf <- maps::map("county", fill = TRUE, plot = FALSE) %>%
  st_as_sf() %>%
  mutate(
    ID = as.character(ID),
    state = str_extract(ID, "^[^,]+")
  ) %>%
  filter(state %in% str_to_lower(unique(map_metadata$State))) %>%
  st_set_crs(4326) %>%
  st_transform(3857) %>%
  st_make_valid() %>%
  st_intersection(map_bbox_sfc)

# Hard-coded nudges for scatterpie positions, in EPSG:3857 meters.
# Edit dx/dy values to manually move pies while keeping black dots at true population means.
nudge_tbl <- tribble(
  ~index, ~dx, ~dy,
  "01", -45000, -35000,
  "02", -65000,  75000,
  "03", -35000, -70000,
  "04",  70000,  70000,
  "05",  85000, -15000,
  "06", -25000, -45000,
  "07", -35000,  30000,
  "08",  45000,  30000,
  "09", -25000, -65000,
  "10",  45000, -35000,
  "11",  95000,  35000,
  "12",  45000,  35000,
  "13",  45000,  25000
)

base_pie_radius <- min(
  as.numeric(map_bbox["xmax"] - map_bbox["xmin"]),
  as.numeric(map_bbox["ymax"] - map_bbox["ymin"])
) * 0.032

all_k_long <- vector("list", length(k_plot))
names(all_k_long) <- paste0("K", k_plot)

for (k in k_plot) {
  message("Plotting K = ", k)
  k_data <- prepare_k_data(k, sample_metadata)
  cluster_cols <- paste0("cluster_", seq_len(k))

  all_k_long[[paste0("K", k)]] <- k_data %>%
    select(K, Sample_Name, State, state_abbr, index, population_label, Longitude, Latitude, all_of(cluster_cols)) %>%
    pivot_longer(all_of(cluster_cols), names_to = "cluster", values_to = "ancestry")

  bar_plot <- plot_barplot(k_data, k)
  map_plot <- plot_scatterpie_map(k_data, k, state_sf, county_sf, map_bbox, nudge_tbl, base_pie_radius)

  ggsave(
    file.path(output_dir, sprintf("Quercus_ADMIXTURE_barplot_K%02d.png", k)),
    bar_plot,
    width = 13,
    height = 5.2,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(output_dir, sprintf("Quercus_ADMIXTURE_barplot_K%02d.pdf", k)),
    bar_plot,
    width = 13,
    height = 5.2,
    bg = "white"
  )

  ggsave(
    file.path(output_dir, sprintf("Quercus_ADMIXTURE_scatterpie_map_K%02d.png", k)),
    map_plot,
    width = 10.5,
    height = 6.2,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(output_dir, sprintf("Quercus_ADMIXTURE_scatterpie_map_K%02d.pdf", k)),
    map_plot,
    width = 10.5,
    height = 6.2,
    bg = "white"
  )
}

bind_rows(all_k_long) %>%
  readr::write_csv(file.path(output_dir, "Quercus_ADMIXTURE_individual_ancestry_long_K02_K10.csv"))

message("Done. Outputs saved to: ", output_dir)
