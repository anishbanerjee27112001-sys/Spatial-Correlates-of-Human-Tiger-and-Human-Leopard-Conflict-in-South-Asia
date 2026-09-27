# Project paths. Run R scripts from the repository root.
# Set HCC_DATA_DIR to use a different computer or data location.
base_dir <- Sys.getenv("HCC_DATA_DIR", unset = "C:/Users/Anish Banerjee/Desktop/HCC/Spatial Determinants of HCC/Final Datasets")
methods <- c("AccessibleArea", "BiasCorrected", "UniformRandom")
datasets <- expand.grid(Method = methods, Species = c("HLC", "HTC"),
                        stringsAsFactors = FALSE)
datasets$Label <- paste(datasets$Species, datasets$Method, "Merged_Final", sep = "_")
datasets$Input <- file.path(base_dir, paste(datasets$Species, "Data"), "Background",
                            paste0(datasets$Label, ".csv"))
datasets$Extracted <- file.path(base_dir, paste(datasets$Species, "Data"), "Extracted")
conflict_files <- setNames(file.path(base_dir, paste(c("HLC", "HTC"), "Data"),
                          paste(c("HLC", "HTC"), "Final Conflict Data.xlsx")), c("HLC", "HTC"))
# Raster locations supplied in the original document; change here if relocated.
raster_root <- Sys.getenv("HCC_RASTER_DIR", unset = "C:/Users/Anish Banerjee/Desktop/HCC/Rasters")
glad_dir <- Sys.getenv("HCC_GLAD_DIR", unset = "C:/Users/Anish Banerjee/Desktop/HWC Paper/GLAD_30m")
worldclim_file <- Sys.getenv("HCC_WORLDCLIM_FILE", unset = file.path(raster_root, "wc2.1_30s_bio.tif"))
check_six_inputs <- function() {
  absent <- datasets$Input[!file.exists(datasets$Input)]
  if (length(absent)) stop("Missing input files:\n", paste(absent, collapse = "\n"))
  invisible(TRUE)
}
# Use the same stable row identifiers in R and GEE extraction outputs.
# Existing predictor columns are omitted to prevent duplicate/stale outputs.
read_extraction_points <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("Long", "Lat", "Year", "Conflict")
  absent <- setdiff(required, names(x))
  if (length(absent)) stop(basename(path), " missing: ", paste(absent, collapse = ", "))
  if (!"Record_ID" %in% names(x)) {
    x$Record_ID <- sprintf("%s_%07d", tools::file_path_sans_ext(basename(path)), seq_len(nrow(x)))
  }
  x$Record_ID <- as.character(x$Record_ID)
  if (anyNA(x$Record_ID) || any(!nzchar(x$Record_ID)) || anyDuplicated(x$Record_ID)) stop("Invalid Record_ID in ", path)
  for (nm in required) x[[nm]] <- suppressWarnings(as.numeric(as.character(x[[nm]])))
  if (any(!is.finite(x$Long) | !is.finite(x$Lat) | abs(x$Long) > 180 | abs(x$Lat) > 90)) stop("Invalid coordinates in ", path)
  if (any(!is.finite(x$Year))) stop("Year must already be assigned in ", path)
  if (any(!x$Conflict %in% c(0, 1))) stop("Conflict must be 0 or 1 in ", path)
  x[c("Record_ID", required)]
}
