library(readxl)
library(tidyverse)
library(writexl)

# Editing notes:
# - Change setwd() if this folder moves or if you run this script on another computer.
# - Change metadata_file if the original Excel file name changes.
# - Change output_file if you want a different output file name.
# - Change combined_counties or combined_county_name to merge different county groups.
# - The new index is ordered west-to-east by mean Longitude of each new_county.
# - index_new_county_state is formatted like "01-Caddo (LA)".
# - The original metadata already has an index column, so this script renames it
#   to original_index before adding the new county-order index column.

setwd("/Users/kuowenhsi/Library/CloudStorage/OneDrive-MissouriBotanicalGarden/General - IMLS National Leadership Grant 2023/Manuscript/Quercus")

metadata_file <- "data/metadata/Quercus_metadata.xlsx"
output_file <- "data/metadata/Quercus_metadata_with_index.xlsx"

combined_counties <- c("Nevada", "Ouachita", "Poinsett")
combined_county_name <- "Nevada-Ouachita-Poinsett"

metadata <- read_excel(metadata_file)

required_columns <- c("County", "State", "Latitude", "Longitude")
missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}

state_lookup <- tibble(
  State = state.name,
  state_abbr = state.abb
)

metadata_with_new_county <- metadata %>%
  rename(original_index = index) %>%
  mutate(
    new_county = case_when(
      County %in% combined_counties ~ combined_county_name,
      TRUE ~ County
    )
  ) %>%
  left_join(state_lookup, by = "State")

county_index_lookup <- metadata_with_new_county %>%
  filter(
    !is.na(new_county),
    !is.na(Latitude),
    !is.na(Longitude)
  ) %>%
  group_by(new_county) %>%
  summarise(
    mean_longitude = mean(Longitude),
    mean_latitude = mean(Latitude),
    state_abbr = paste(sort(unique(state_abbr[!is.na(state_abbr)])), collapse = "/"),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(mean_longitude) %>%
  mutate(
    index = str_pad(row_number(), width = 2, pad = "0"),
    index_new_county_state = paste0(index, "-", new_county, " (", state_abbr, ")")
  ) %>%
  select(index, new_county, index_new_county_state, state_abbr, mean_longitude, mean_latitude, n)

metadata_indexed <- metadata_with_new_county %>%
  left_join(
    county_index_lookup %>% select(index, new_county, index_new_county_state),
    by = "new_county"
  ) %>%
  relocate(index, new_county, index_new_county_state, .after = County)

write_xlsx(
  list(
    metadata = metadata_indexed,
    county_index_lookup = county_index_lookup
  ),
  output_file
)

message("Saved ", output_file)
message("Rows written: ", nrow(metadata_indexed))
message("County groups indexed: ", nrow(county_index_lookup))
