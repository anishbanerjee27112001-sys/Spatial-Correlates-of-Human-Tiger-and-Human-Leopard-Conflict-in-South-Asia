source("config/paths.R")
check_six_inputs()

# ============================================================
# TEMPORALLY MATCHED TOTAL LIVESTOCK DENSITY EXTRACTION
# Combines buffalo, cattle, chicken, goats, pigs, sheep
# Uses DA rasters for 2010 and 2015 + combined 2020 raster
# Extracts mean livestock density within 1 km buffer
# ============================================================

library(terra)
library(dplyr)
library(stringr)

# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

rasters_2010_dir <- file.path(raster_root, "2010")
rasters_2015_dir <- file.path(raster_root, "2015")

raster_2020_file <- file.path(raster_root, "GLW4-2020.D-DA.GLEAM3-ALL-LU.tif")

# ------------------------------------------------------------
# 2. FUNCTION TO LOAD AND SUM DA SPECIES RASTERS
# ------------------------------------------------------------

make_total_livestock_raster <- function(folder_path, year_label) {
  
  files <- list.files(
    folder_path,
    pattern = "\\.tif$|\\.tiff$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  # Keep only DA density rasters
  files <- files[str_detect(basename(files), regex("Da|DA", ignore_case = TRUE))]
  
  if (length(files) == 0) {
    stop("No DA tif files found for ", year_label, " in: ", folder_path)
  }
  
  message("\n", year_label, " DA rasters found:")
  print(basename(files))
  
  r_list <- lapply(files, rast)
  
  # Align all rasters to the first raster
  template <- r_list[[1]]
  
  r_list_aligned <- lapply(r_list, function(r) {
    if (!compareGeom(template, r, stopOnError = FALSE)) {
      r <- resample(r, template, method = "bilinear")
    }
    return(r)
  })
  
  r_stack <- rast(r_list_aligned)
  
  # Sum all species densities
  total <- app(r_stack, fun = sum, na.rm = TRUE)
  names(total) <- paste0("livestock_total_", year_label)
  
  return(total)
}

# ------------------------------------------------------------
# 3. CREATE TOTAL LIVESTOCK RASTERS
# ------------------------------------------------------------

ls_2010_total <- make_total_livestock_raster(rasters_2010_dir, "2010")
ls_2015_total <- make_total_livestock_raster(rasters_2015_dir, "2015")
ls_2020_total <- rast(raster_2020_file)

if (nlyr(ls_2020_total) > 1) {
  ls_2020_total <- ls_2020_total[[1]]
}

names(ls_2020_total) <- "livestock_total_2020"

# ------------------------------------------------------------
# 4. TEMPORAL MATCHING FUNCTION
# ------------------------------------------------------------

assign_livestock_year <- function(year) {
  ifelse(year <= 2010, 2010,
  ifelse(year >= 2011 & year <= 2019, 2015,
  ifelse(year >= 2020, 2020, NA)))
}

# ------------------------------------------------------------
# 5. YEAR COLUMN DETECTION
# ------------------------------------------------------------

get_year_column <- function(df) {
  possible_cols <- c(
    "Year",
    "YearMidpoint",
    "Year.Midpoint",
    "Year Mid-Point",
    "Year_bg"
  )
  
  year_col <- possible_cols[possible_cols %in% names(df)][1]
  
  if (is.na(year_col)) {
    stop("No recognised year column found.")
  }
  
  return(year_col)
}

# ------------------------------------------------------------
# 6. EXTRACTION FUNCTION
# ------------------------------------------------------------

extract_total_livestock <- function(csv_path, output_name) {
  
  message("\nProcessing: ", basename(csv_path))
  
  df <- read_extraction_points(csv_path)
  
  year_col <- get_year_column(df)
  df$livestock_input_year <- as.numeric(substr(as.character(df[[year_col]]), 1, 4))
  df$livestock_raster_year <- assign_livestock_year(df$livestock_input_year)
  
  # Convert points to spatial vector
  pts <- vect(
    df,
    geom = c("Long", "Lat"),
    crs = "EPSG:4326",
    keepgeom = TRUE
  )
  
  df$livestock_density_1km <- NA_real_
  
  # -------------------------
  # 2010 extraction
  # -------------------------
  idx_2010 <- which(df$livestock_raster_year == 2010)
  
  if (length(idx_2010) > 0) {
    pts_2010 <- project(pts[idx_2010, ], crs(ls_2010_total))
    buf_2010 <- buffer(pts_2010, width = 1000)
    
    vals <- extract(
      ls_2010_total,
      buf_2010,
      fun = mean,
      na.rm = TRUE
    )
    
    df$livestock_density_1km[idx_2010] <- vals[, 2]
  }
  
  # -------------------------
  # 2015 extraction
  # -------------------------
  idx_2015 <- which(df$livestock_raster_year == 2015)
  
  if (length(idx_2015) > 0) {
    pts_2015 <- project(pts[idx_2015, ], crs(ls_2015_total))
    buf_2015 <- buffer(pts_2015, width = 1000)
    
    vals <- extract(
      ls_2015_total,
      buf_2015,
      fun = mean,
      na.rm = TRUE
    )
    
    df$livestock_density_1km[idx_2015] <- vals[, 2]
  }
  
  # -------------------------
  # 2020 extraction
  # -------------------------
  idx_2020 <- which(df$livestock_raster_year == 2020)
  
  if (length(idx_2020) > 0) {
    pts_2020 <- project(pts[idx_2020, ], crs(ls_2020_total))
    buf_2020 <- buffer(pts_2020, width = 1000)
    
    vals <- extract(
      ls_2020_total,
      buf_2020,
      fun = mean,
      na.rm = TRUE
    )
    
    df$livestock_density_1km[idx_2020] <- vals[, 2]
  }
  
  out_dir <- file.path(dirname(dirname(csv_path)), "Extracted")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, output_name)
  write.csv(df, out_file, row.names = FALSE)
  
  message("Saved: ", out_file)
  return(out_file)
}

# ------------------------------------------------------------
# 7. FILES TO PROCESS
# ------------------------------------------------------------

files_to_process <- data.frame(csv_path = datasets$Input,
  output_name = paste0(datasets$Label, "_Livestock.csv"))
outputs <- mapply(extract_total_livestock, files_to_process$csv_path,
                  files_to_process$output_name)
print(outputs)
