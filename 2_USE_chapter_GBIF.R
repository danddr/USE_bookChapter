# ================================================
#     Downloading and Cleaning Data from GBIF
# ================================================
setwd("/home/daniele/Documents/working_files/USE_bookChapter/")
# GBIF requires user credentials for data downloads. Set your GBIF username and password here:
# Details: https://docs.ropensci.org/rgbif/articles/gbif_credentials.html
usethis::edit_r_environ()
# For more details on using CoordinateCleaner, refer to:
# https://ropensci.github.io/CoordinateCleaner/articles/Cleaning_GBIF_data_with_CoordinateCleaner.html

# Load Required Libraries
library(dplyr)
library(tidyr)
library(rgbif)
library(maps)
library(CoordinateCleaner)
library(countrycode)

# Define Target Species
# To obtain occurrence records, it is preferable to use a taxonKey.
myspecies <- c("Acer pseudoplatanus")
taxonKey <- name_backbone(myspecies)$usageKey
taxonKey

# Download Occurrence Data
# Specify download filters such as country, year, and geospatial presence.
res <- occ_download(
  pred_in("taxonKey", taxonKey),
  pred_in("country", "IT"), # set the country
  pred_gte("year", 1990),
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus", "PRESENT"),
  format = "SIMPLE_CSV"
)

# Check download status and retrieve metadata
occ_download_wait(res)
occ_download_meta(res) # download key: 0010041-241024112534372

# GBIF citation DOI
rgbif::gbif_citation("0010041-241024112534372") # in our case https://doi.org/10.15468/dl.dy7rn8

# Import Data
# Once downloaded, import the data
GBIFpre <- occ_download_get(key = res, overwrite = TRUE) %>% 
  occ_download_import()
str(GBIFpre)
write.csv(GBIFpre, "outputs/acer_pseudoplatanus.csv")

# Select Relevant Columns for Mapping and Cleaning
myspecies_coords <- GBIFpre %>% 
  dplyr::select(
    "decimalLongitude", "decimalLatitude", "individualCount",
    "coordinateUncertaintyInMeters", "institutionCode", "elevation",
    "basisOfRecord", "issue", countryCode, species
  )

# Quick Map of Occurrence Data
map("world", xlim = range(myspecies_coords$decimalLongitude), # if the map doesn't appear right at first, run this command again
    ylim = range(myspecies_coords$decimalLatitude))
points(myspecies_coords[, c("decimalLongitude", "decimalLatitude")], pch = ".")

# Clean the Dataset
# Remove records with absence, zero abundance, and large coordinate uncertainty (> 5 km).
myspecies_coords <- myspecies_coords %>%
  drop_na(decimalLatitude, decimalLongitude) %>%
  filter(basisOfRecord %in% c("OCCURRENCE", "HUMAN_OBSERVATION")) %>%
  filter(coordinateUncertaintyInMeters < 5000) %>%
  filter(!is.na(individualCount) | individualCount > 0)

# Convert country codes to ISO3 format for cleaner processing
myspecies_coords$iso3 <- countrycode::countrycode(myspecies_coords$countryCode, 
                                                  origin = 'iso2c', 
                                                  destination = 'iso3c')

# Flag potentially problematic records using CoordinateCleaner
flags <- CoordinateCleaner::clean_coordinates(
  x = myspecies_coords,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  countries = "iso3",
  species = "species",
  tests = c("capitals", "centroids", "equal", "gbif", "zeros", "seas", "duplicates")
)

summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

# Exclude problematic records
dat_cl <- myspecies_coords[flags$.summary, ]

# Map Cleaned Occurrence Data
map("world", xlim = range(dat_cl$decimalLongitude), 
    ylim = range(dat_cl$decimalLatitude))
points(dat_cl[, c("decimalLongitude", "decimalLatitude")], pch = ".")

# Save the Cleaned Data
saveRDS(dat_cl, "outputs/gbifClean.RDS")
