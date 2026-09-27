source("config/paths.R")
check_six_inputs()

# ============================================================
# GPW POPULATION DENSITY EXTRACTION IN R
# Keeps ALL rows.
# Pre-2000 rows retained with population outputs = NA.
# Temporal matching:
# 2000–2004 -> 2000
# 2005–2009 -> 2005
# 2010–2014 -> 2010
# 2015–2019 -> 2015
# 2020+     -> 2020
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(stringr)
})

# -----------------------------
# 1. PATHS
# -----------------------------

gpw_dir <- file.path(raster_root, "Population Density")

lon_col <- "Long"
lat_col <- "Lat"
year_col <- "Year"

# -----------------------------
# 2. FIND INPUT CSV FILES
# -----------------------------

input_files <- datasets$Input

cat("Found CSV files:\n")
print(basename(input_files))

# -----------------------------
# 3. READ GPW RASTERS
# -----------------------------
# Expected filenames must contain 2000, 2005, 2010, 2015, 2020

gpw_files <- list.files(
  gpw_dir,
  pattern = "\\.(tif|tiff|img|asc|grd|nc)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(gpw_files) == 0) {
  stop("No GPW raster files found. Check gpw_dir and file extensions.")
}

get_year_from_filename <- function(x) {
  yr <- str_extract(basename(x), "2000|2005|2010|2015|2020")
  as.integer(yr)
}

gpw_years <- sapply(gpw_files, get_year_from_filename)

if (any(is.na(gpw_years))) {
  stop(
    "Could not detect year in some GPW files:\n",
    paste(basename(gpw_files[is.na(gpw_years)]), collapse = "\n")
  )
}

gpw_list <- setNames(
  lapply(gpw_files, rast),
  as.character(gpw_years)
)

required_years <- c("2000", "2005", "2010", "2015", "2020")

missing_years <- setdiff(required_years, names(gpw_list))
if (length(missing_years) > 0) {
  stop("Missing GPW rasters for years: ", paste(missing_years, collapse = ", "))
}

cat("\nLoaded GPW rasters:\n")
print(names(gpw_list))

# -----------------------------
# 4. TEMPORAL MATCHING FUNCTION
# -----------------------------

assign_pop_year <- function(y) {
  dplyr::case_when(
    is.na(y) ~ NA_real_,
    y < 2000 ~ NA_real_,
    y >= 2000 & y <= 2004 ~ 2000,
    y >= 2005 & y <= 2009 ~ 2005,
    y >= 2010 & y <= 2014 ~ 2010,
    y >= 2015 & y <= 2019 ~ 2015,
    y >= 2020 ~ 2020
  )
}

# -----------------------------
# 5. PROCESS ONE CSV
# -----------------------------

process_one_file <- function(csv_path) {

  cat("\nProcessing:", basename(csv_path), "\n")

  dat <- read_extraction_points(csv_path)

  if (!all(c(lon_col, lat_col, year_col) %in% names(dat))) {
    stop("Missing one of Long, Lat, Year in: ", basename(csv_path))
  }

  dat$row_id_pop <- seq_len(nrow(dat))

  dat[[lon_col]] <- suppressWarnings(as.numeric(dat[[lon_col]]))
  dat[[lat_col]] <- suppressWarnings(as.numeric(dat[[lat_col]]))
  dat[[year_col]] <- suppressWarnings(as.numeric(substr(as.character(dat[[year_col]]), 1, 4)))

  dat$pop_density_year_used_Final <- assign_pop_year(dat[[year_col]])
  dat$pop_density_people_per_km2_Final <- NA_real_

  valid_coord <- which(
    !is.na(dat[[lon_col]]) &
      !is.na(dat[[lat_col]]) &
      dat[[lon_col]] >= -180 &
      dat[[lon_col]] <= 180 &
      dat[[lat_col]] >= -90 &
      dat[[lat_col]] <= 90
  )

  valid_extract <- valid_coord[
    !is.na(dat$pop_density_year_used_Final[valid_coord])
  ]

  if (length(valid_extract) == 0) {
    warning("No valid post-2000 rows with coordinates in ", basename(csv_path))
  } else {

    pts <- vect(
      dat[valid_extract, ],
      geom = c(lon_col, lat_col),
      crs = "EPSG:4326"
    )

    for (yr in c(2000, 2005, 2010, 2015, 2020)) {

      idx <- valid_extract[dat$pop_density_year_used_Final[valid_extract] == yr]

      if (length(idx) == 0) next

      cat("  Extracting GPW year:", yr, "| rows:", length(idx), "\n")

      r <- gpw_list[[as.character(yr)]]

      pts_year <- vect(
        dat[idx, ],
        geom = c(lon_col, lat_col),
        crs = "EPSG:4326"
      )

      pts_year <- project(pts_year, crs(r))

      vals <- terra::extract(r, pts_year, ID = FALSE)

      # If raster has multiple layers, use first layer
      dat$pop_density_people_per_km2_Final[idx] <- vals[[1]]
    }
  }

  dat$row_id_pop <- NULL

  out_path <- file.path(
    file.path(dirname(dirname(csv_path)), "Extracted"),
    paste0(
      tools::file_path_sans_ext(basename(csv_path)),
      "_POP_DENSITY_FINAL.csv"
    )
  )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(dat, out_path)

  cat("Saved:", out_path, "\n")
  cat("Rows kept:", nrow(dat), "\n")
  cat("Pre-2000 / missing pop year rows:",
      sum(is.na(dat$pop_density_year_used_Final)), "\n")
  cat("Missing extracted population values:",
      sum(is.na(dat$pop_density_people_per_km2_Final)), "\n")
}

# -----------------------------
# 6. RUN ALL FILES
# -----------------------------

for (f in input_files) {
  process_one_file(f)
}
