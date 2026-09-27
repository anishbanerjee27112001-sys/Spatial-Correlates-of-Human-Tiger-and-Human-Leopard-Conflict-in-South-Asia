source("config/paths.R")
check_six_inputs()

# ============================================================
# GLAD 30 m FRAGMENTATION EXTRACTION
# ============================================================
#
# Calculates FRAGSTATS-style forest fragmentation metrics
# from GLAD GLCLU2020 30 m categorical land-cover rasters.
#
# Forest definition:
#   GLAD classes 5–13
#
# Metrics:
#   1. Forest patch density              patches / 100 ha
#   2. Mean forest patch size            ha
#   3. Forest edge density               m / ha
#   4. Forest aggregation index          %
#
# Spatial scales:
#   1 km
#   2 km
#   5 km
#
# GLAD temporal matching:
#   < 2000      -> excluded
#   2000–2004   -> 2000
#   2005–2009   -> 2005
#   2010–2014   -> 2010
#   2015–2019   -> 2015
#   >= 2020     -> 2020
#
# ============================================================

# ------------------------------------------------------------
# 0. PACKAGES
# ------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(landscapemetrics)
library(stringr)

# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

# Folder containing GLAD 30 m categorical rasters.
#
# Filenames need to contain:
# 2000, 2005, 2010, 2015 or 2020
#
# Example:
# GLAD_LCLUC_2000.tif

# Input HLC / HTC datasets

# Output folder

# ------------------------------------------------------------
# 2. METRIC CRS
# ------------------------------------------------------------

# Metric CRS used for South Asia.
#
# Buffers will therefore be interpreted in metres.
#
# IMPORTANT:
# fragmentation metrics should not be calculated directly
# on geographic latitude/longitude rasters.

metric_crs <- "EPSG:7755"

# ------------------------------------------------------------
# 3. FIND GLAD RASTERS
# ------------------------------------------------------------

glad_files <- list.files(
  glad_dir,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

# Extract GLAD reference year from filename

get_glad_year <- function(x) {

  as.numeric(
    str_extract(
      basename(x),
      "(2000|2005|2010|2015|2020)"
    )
  )

}

glad_lookup <- data.frame(
  year = sapply(
    glad_files,
    get_glad_year
  ),
  file = glad_files
) |>
  filter(!is.na(year)) |>
  arrange(year)

print(glad_lookup)

# ------------------------------------------------------------
# Check required GLAD years
# ------------------------------------------------------------

required_years <- c(
  2000,
  2005,
  2010,
  2015,
  2020
)

missing_years <- setdiff(
  required_years,
  glad_lookup$year
)

if (length(missing_years) > 0) {

  warning(
    "Missing GLAD rasters for: ",
    paste(
      missing_years,
      collapse = ", "
    )
  )

}

# ------------------------------------------------------------
# Check for duplicate rasters per year
# ------------------------------------------------------------

duplicate_years <- glad_lookup |>
  count(year) |>
  filter(n > 1)

if (nrow(duplicate_years) > 0) {

  warning(
    "More than one GLAD raster detected for some years. ",
    "The first matching raster will be used."
  )

  print(duplicate_years)

}

# ------------------------------------------------------------
# 4. TEMPORAL MATCHING
# ------------------------------------------------------------

# Correct GLAD temporal bins:
#
# < 2000      -> NA / excluded
# 2000–2004   -> 2000
# 2005–2009   -> 2005
# 2010–2014   -> 2010
# 2015–2019   -> 2015
# >= 2020     -> 2020

assign_glad_year <- function(year) {

  case_when(

    is.na(year) ~ NA_real_,

    year < 2000 ~ NA_real_,

    year >= 2000 & year <= 2004 ~ 2000,

    year >= 2005 & year <= 2009 ~ 2005,

    year >= 2010 & year <= 2014 ~ 2010,

    year >= 2015 & year <= 2019 ~ 2015,

    year >= 2020 ~ 2020,

    TRUE ~ NA_real_

  )

}

# ------------------------------------------------------------
# 5. GLAD FOREST MASK
# ------------------------------------------------------------

# GLAD forest classes:
#
# 5–13 = tree / forest-related land-cover classes
#
# Keep:
#
# 1 = forest
# 0 = non-forest
#
# Non-forest must NOT be converted to NA because fragmentation
# metrics require forest/non-forest boundaries.

make_forest_binary <- function(lc) {

  forest <- ifel(
    lc >= 5 & lc <= 13,
    1,
    0
  )

  names(forest) <- "forest"

  return(forest)

}

# ------------------------------------------------------------
# 6. FRAGMENTATION METRIC FUNCTION
# ------------------------------------------------------------

calc_frag_metrics <- function(forest_crop) {

  # ----------------------------------------------------------
  # Extract available raster values
  # ----------------------------------------------------------

  vals <- values(
    forest_crop,
    na.rm = TRUE
  )

  # ----------------------------------------------------------
  # Empty raster / outside coverage
  # ----------------------------------------------------------

  if (length(vals) == 0) {

    return(
      data.frame(

        forest_patch_density = NA_real_,

        forest_mean_patch_size = NA_real_,

        forest_edge_density = NA_real_,

        forest_aggregation = NA_real_

      )
    )

  }

  # ----------------------------------------------------------
  # No forest in buffer
  # ----------------------------------------------------------
  #
  # PD = 0
  # AREA_MN = 0
  # ED = 0
  #
  # AI is undefined where there is no forest class.

  if (!any(vals == 1, na.rm = TRUE)) {

    return(
      data.frame(

        forest_patch_density = 0,

        forest_mean_patch_size = 0,

        forest_edge_density = 0,

        forest_aggregation = NA_real_

      )
    )

  }

  # ----------------------------------------------------------
  # Convert to categorical raster
  # ----------------------------------------------------------

  forest_crop <- as.factor(
    forest_crop
  )

  # ----------------------------------------------------------
  # Patch density
  #
  # Units:
  # patches / 100 ha
  # ----------------------------------------------------------

  pd <- tryCatch(

    lsm_c_pd(
      forest_crop,
      directions = 8
    ) |>
      filter(class == 1) |>
      pull(value),

    error = function(e) NA_real_

  )

  # ----------------------------------------------------------
  # Mean patch area
  #
  # Units:
  # hectares
  # ----------------------------------------------------------

  area_mn <- tryCatch(

    lsm_c_area_mn(
      forest_crop,
      directions = 8
    ) |>
      filter(class == 1) |>
      pull(value),

    error = function(e) NA_real_

  )

  # ----------------------------------------------------------
  # Edge density
  #
  # Units:
  # metres / hectare
  # ----------------------------------------------------------

  ed <- tryCatch(

    lsm_c_ed(
      forest_crop,
      directions = 8
    ) |>
      filter(class == 1) |>
      pull(value),

    error = function(e) NA_real_

  )

  # ----------------------------------------------------------
  # Aggregation index
  #
  # Units:
  # percent
  # ----------------------------------------------------------

  ai <- tryCatch(

    lsm_c_ai(
      forest_crop,
      directions = 8
    ) |>
      filter(class == 1) |>
      pull(value),

    error = function(e) NA_real_

  )

  # ----------------------------------------------------------
  # Return one-row dataframe
  # ----------------------------------------------------------

  data.frame(

    forest_patch_density =
      ifelse(
        length(pd) == 0,
        NA_real_,
        pd[1]
      ),

    forest_mean_patch_size =
      ifelse(
        length(area_mn) == 0,
        NA_real_,
        area_mn[1]
      ),

    forest_edge_density =
      ifelse(
        length(ed) == 0,
        NA_real_,
        ed[1]
      ),

    forest_aggregation =
      ifelse(
        length(ai) == 0,
        NA_real_,
        ai[1]
      )

  )

}

# ------------------------------------------------------------
# 7. YEAR-COLUMN DETECTION
# ------------------------------------------------------------

get_year_column <- function(df) {

  if ("Year" %in% names(df)) {

    return("Year")

  }

  if ("YearMidpoint" %in% names(df)) {

    return("YearMidpoint")

  }

  if ("Year.Midpoint" %in% names(df)) {

    return("Year.Midpoint")

  }

  if ("Year Mid-Point" %in% names(df)) {

    return("Year Mid-Point")

  }

  if ("Year_bg" %in% names(df)) {

    return("Year_bg")

  }

  stop(
    "No usable year column found."
  )

}

# ------------------------------------------------------------
# 8. COORDINATE CHECK
# ------------------------------------------------------------

check_coordinates <- function(df) {

  if (!all(c("Long", "Lat") %in% names(df))) {

    stop(
      "Dataset must contain Long and Lat columns."
    )

  }

  df$Long <- as.numeric(
    df$Long
  )

  df$Lat <- as.numeric(
    df$Lat
  )

  bad_coords <- which(
    is.na(df$Long) |
      is.na(df$Lat) |
      df$Long < -180 |
      df$Long > 180 |
      df$Lat < -90 |
      df$Lat > 90
  )

  if (length(bad_coords) > 0) {

    warning(
      length(bad_coords),
      " rows contain invalid coordinates and will be excluded."
    )

  }

  df <- df |>
    filter(
      !is.na(Long),
      !is.na(Lat),
      Long >= -180,
      Long <= 180,
      Lat >= -90,
      Lat <= 90
    )

  return(df)

}

# ------------------------------------------------------------
# 9. MAIN EXTRACTION FUNCTION
# ------------------------------------------------------------

extract_glad_frag_for_file <- function(
  csv_path,
  species_label,
  background_label
) {

  message(
    "\n============================================"
  )

  message(
    "Processing: ",
    basename(csv_path)
  )

  message(
    "============================================"
  )

  # ----------------------------------------------------------
  # Read input dataset
  # ----------------------------------------------------------

  df_original <- read_extraction_points(csv_path)

  # Permanent unique ID from original dataset

  df_original$row_id_frag <- seq_len(
    nrow(df_original)
  )

  # Working copy

  df <- df_original

  # ----------------------------------------------------------
  # Check coordinates
  # ----------------------------------------------------------

  df <- check_coordinates(
    df
  )

  # ----------------------------------------------------------
  # Detect year column
  # ----------------------------------------------------------

  year_col <- get_year_column(
    df
  )

  df$Year_frag <- suppressWarnings(
    as.numeric(
      df[[year_col]]
    )
  )

  # ----------------------------------------------------------
  # Assign GLAD reference year
  # ----------------------------------------------------------

  df$GLAD_year_frag <- assign_glad_year(
    df$Year_frag
  )

  # ----------------------------------------------------------
  # Report observations outside temporal coverage
  # ----------------------------------------------------------

  n_invalid_year <- sum(
    is.na(df$GLAD_year_frag)
  )

  n_pre2000 <- sum(
    !is.na(df$Year_frag) &
      df$Year_frag < 2000
  )

  if (n_pre2000 > 0) {

    message(
      "Excluding ",
      n_pre2000,
      " observations dated before 2000."
    )

  }

  if (n_invalid_year > n_pre2000) {

    message(
      "Additional rows with missing/invalid year: ",
      n_invalid_year - n_pre2000
    )

  }

  # ----------------------------------------------------------
  # Remove observations not covered by GLAD
  # ----------------------------------------------------------

  df <- df |>
    filter(
      !is.na(GLAD_year_frag)
    )

  if (nrow(df) == 0) {

    stop(
      "No observations remain after temporal filtering."
    )

  }

  # ----------------------------------------------------------
  # Print temporal distribution
  # ----------------------------------------------------------

  message(
    "\nGLAD temporal distribution:"
  )

  print(
    table(
      df$GLAD_year_frag
    )
  )

  # ----------------------------------------------------------
  # Convert points to sf
  # ----------------------------------------------------------

  pts_sf <- st_as_sf(

    df,

    coords = c(
      "Long",
      "Lat"
    ),

    crs = 4326,

    remove = FALSE

  )

  # ----------------------------------------------------------
  # Convert to terra vector
  # ----------------------------------------------------------

  pts_vect <- vect(
    pts_sf
  )

  # ----------------------------------------------------------
  # Project points to metric CRS
  # ----------------------------------------------------------

  pts_metric <- project(
    pts_vect,
    metric_crs
  )

  # ----------------------------------------------------------
  # Buffer sizes
  # ----------------------------------------------------------

  buffer_sizes <- c(

    "1km" = 1000,

    "2km" = 2000,

    "5km" = 5000

  )

  all_results <- list()

  # ----------------------------------------------------------
  # Process each GLAD year separately
  # ----------------------------------------------------------

  for (
    glad_year in
    sort(
      unique(
        df$GLAD_year_frag
      )
    )
  ) {

    # --------------------------------------------------------
    # Locate raster
    # --------------------------------------------------------

    glad_matches <- glad_lookup$file[
      glad_lookup$year == glad_year
    ]

    if (length(glad_matches) == 0) {

      warning(
        "No GLAD raster found for year: ",
        glad_year
      )

      next

    }

    glad_file <- glad_matches[1]

    message(
      "\nGLAD year: ",
      glad_year
    )

    message(
      "Raster: ",
      basename(glad_file)
    )

    # --------------------------------------------------------
    # Load categorical raster
    # --------------------------------------------------------

    lc <- rast(
      glad_file
    )

    # Use first layer if raster unexpectedly has multiple bands

    if (nlyr(lc) > 1) {

      warning(
        basename(glad_file),
        " contains multiple layers. Using first layer."
      )

      lc <- lc[[1]]

    }

    # --------------------------------------------------------
    # Check raster CRS
    # --------------------------------------------------------

    if (
      is.na(crs(lc)) ||
      crs(lc) == ""
    ) {

      stop(
        "GLAD raster has no CRS: ",
        glad_file
      )

    }

    # --------------------------------------------------------
    # Create binary forest raster
    # --------------------------------------------------------

    forest <- make_forest_binary(
      lc
    )

    # --------------------------------------------------------
    # Select points belonging to this temporal bin
    # --------------------------------------------------------

    pts_year <- pts_metric[
      pts_metric$GLAD_year_frag == glad_year,
    ]

    message(
      "Observations: ",
      nrow(pts_year)
    )

    # --------------------------------------------------------
    # Crop source raster before projection
    # --------------------------------------------------------
    #
    # This avoids projecting the entire regional/global
    # 30 m GLAD raster.
    #
    # Use a 5.5 km envelope around all observations so that
    # every 5 km extraction buffer is fully covered.
    # --------------------------------------------------------

    points_union <- aggregate(
      pts_year
    )

    year_extent_metric <- buffer(
      points_union,
      width = 5500
    )

    # Transform crop area into original raster CRS

    year_extent_src <- project(
      year_extent_metric,
      crs(forest)
    )

    # --------------------------------------------------------
    # Check whether raster and points overlap
    # --------------------------------------------------------

    overlap_test <- relate(
      year_extent_src,
      as.polygons(
        ext(forest),
        crs = crs(forest)
      ),
      "intersects"
    )

    if (!any(overlap_test)) {

      warning(
        "Points for ",
        glad_year,
        " do not overlap GLAD raster extent."
      )

      next

    }

    # --------------------------------------------------------
    # Crop raster
    # --------------------------------------------------------

    forest_subset <- crop(
      forest,
      year_extent_src
    )

    # --------------------------------------------------------
    # Project cropped categorical raster to metric CRS
    # --------------------------------------------------------
    #
    # Nearest-neighbour resampling is mandatory for
    # categorical land-cover data.
    #
    # res = 30 retains approximately 30 m spatial resolution.
    # --------------------------------------------------------

    forest_metric <- project(

      forest_subset,

      metric_crs,

      method = "near",

      res = 30

    )

    # Ensure binary integer classes after reprojection

    forest_metric <- ifel(
      forest_metric >= 0.5,
      1,
      0
    )

    names(
      forest_metric
    ) <- "forest"

    # --------------------------------------------------------
    # Process each buffer size
    # --------------------------------------------------------

    for (
      buf_name in names(
        buffer_sizes
      )
    ) {

      radius_m <- buffer_sizes[
        [buf_name]
      ]

      message(
        "  Buffer: ",
        buf_name
      )

      res_list <- vector(
        "list",
        nrow(pts_year)
      )

      # ------------------------------------------------------
      # Point-by-point extraction
      # ------------------------------------------------------

      for (
        i in seq_len(
          nrow(pts_year)
        )
      ) {

        # ----------------------------------------------------
        # Circular buffer
        # ----------------------------------------------------

        buf_i <- buffer(
          pts_year[i, ],
          width = radius_m
        )

        # ----------------------------------------------------
        # Crop and mask
        # ----------------------------------------------------

        forest_crop <- tryCatch(

          {

            temp <- crop(
              forest_metric,
              buf_i
            )

            mask(
              temp,
              buf_i
            )

          },

          error = function(e) {

            message(
              "    Crop/mask error at row ",
              pts_year$row_id_frag[i],
              ": ",
              e$message
            )

            NULL

          }

        )

        # ----------------------------------------------------
        # Calculate metrics
        # ----------------------------------------------------

        if (is.null(forest_crop)) {

          metrics <- data.frame(

            forest_patch_density = NA_real_,

            forest_mean_patch_size = NA_real_,

            forest_edge_density = NA_real_,

            forest_aggregation = NA_real_

          )

        } else {

          metrics <- calc_frag_metrics(
            forest_crop
          )

        }

        # ----------------------------------------------------
        # Attach identifiers
        # ----------------------------------------------------

        metrics$row_id_frag <-
          pts_year$row_id_frag[i]

        metrics$buffer <-
          buf_name

        metrics$GLAD_year <-
          glad_year

        res_list[[i]] <-
          metrics

        # Optional progress message every 1000 points

        if (
          i %% 1000 == 0
        ) {

          message(
            "    Completed ",
            i,
            " / ",
            nrow(pts_year)
          )

        }

      }

      # ------------------------------------------------------
      # Save year × buffer results
      # ------------------------------------------------------

      all_results[
        [paste(
          glad_year,
          buf_name,
          sep = "_"
        )]
      ] <- bind_rows(
        res_list
      )

    }

    # --------------------------------------------------------
    # Release large raster objects
    # --------------------------------------------------------

    rm(
      forest,
      forest_subset,
      forest_metric,
      lc
    )

    gc()

  }

  # ----------------------------------------------------------
  # Combine all extraction results
  # ----------------------------------------------------------

  long_res <- bind_rows(
    all_results
  )

  if (nrow(long_res) == 0) {

    stop(
      "No fragmentation metrics were successfully extracted."
    )

  }

  # ----------------------------------------------------------
  # Convert long -> wide
  # ----------------------------------------------------------

  wide_res <- long_res |>

    select(

      row_id_frag,

      buffer,

      forest_patch_density,

      forest_mean_patch_size,

      forest_edge_density,

      forest_aggregation

    ) |>

    pivot_wider(

      id_cols =
        row_id_frag,

      names_from =
        buffer,

      values_from = c(

        forest_patch_density,

        forest_mean_patch_size,

        forest_edge_density,

        forest_aggregation

      )

    )

  # ----------------------------------------------------------
  # Join metrics back to ORIGINAL dataset
  # ----------------------------------------------------------
  #
  # Important:
  #
  # Pre-2000 observations remain in the output but their
  # GLAD fragmentation variables are NA.
  #
  # They are excluded from extraction, not silently removed
  # from the original dataset.
  # ----------------------------------------------------------

  temporal_info <- df |>
    select(
      row_id_frag,
      Year_frag,
      GLAD_year_frag
    )

  out <- df_original |>

    left_join(
      temporal_info,
      by = "row_id_frag"
    ) |>

    left_join(
      wide_res,
      by = "row_id_frag"
    )

  # ----------------------------------------------------------
  # Rename variables
  # ----------------------------------------------------------

  names(out) <- names(out) |>

    str_replace(
      "^forest_patch_density_",
      "GLAD_Forest_Patch_Density_patches_per_100ha_"
    ) |>

    str_replace(
      "^forest_mean_patch_size_",
      "GLAD_Forest_Mean_Patch_Size_ha_"
    ) |>

    str_replace(
      "^forest_edge_density_",
      "GLAD_Forest_Edge_Density_m_per_ha_"
    ) |>

    str_replace(
      "^forest_aggregation_",
      "GLAD_Forest_Aggregation_Index_"
    )

  # ----------------------------------------------------------
  # Output filename
  # ----------------------------------------------------------

  out_dir <- file.path(dirname(dirname(csv_path)), "Extracted")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (prefix in c("GLAD_Forest_Patch_Density_patches_per_100ha_",
                   "GLAD_Forest_Mean_Patch_Size_ha_", "GLAD_Forest_Edge_Density_m_per_ha_")) {
    idx <- startsWith(names(out), prefix) & !grepl("_FINAL$", names(out))
    names(out)[idx] <- paste0(names(out)[idx], "_FINAL")
  }
  out_file <- file.path(

    out_dir,

    paste0(

      species_label,
      "_",
      background_label,
      "_Merged_Final_GLAD_Fragmentation.csv"

    )

  )

  # ----------------------------------------------------------
  # Save
  # ----------------------------------------------------------

  write.csv(

    out,

    out_file,

    row.names = FALSE

  )

  message(
    "\nSaved: ",
    out_file
  )

  return(
    out_file
  )

}

# ============================================================
# 10. FILES TO PROCESS
# ============================================================

files_to_process <- data.frame(csv_path = datasets$Input,
  species = datasets$Species, background = datasets$Method)
outputs <- mapply(extract_glad_frag_for_file,
  csv_path = files_to_process$csv_path,
  species_label = files_to_process$species,
  background_label = files_to_process$background,
  SIMPLIFY = TRUE, USE.NAMES = FALSE)
print(outputs)
