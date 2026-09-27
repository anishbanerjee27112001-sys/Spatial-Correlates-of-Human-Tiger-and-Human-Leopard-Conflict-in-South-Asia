source("config/paths.R")
check_six_inputs()

library(terra)
if (!file.exists(worldclim_file)) stop("WorldClim raster not found: ", worldclim_file)
bio <- rast(worldclim_file)
ids <- c(3, 4, 5, 6, 12, 16, 17)
# Match explicit layer names; never infer band order from an unnamed stack.
layer_ids <- suppressWarnings(as.integer(sub("^.*bio_?", "", names(bio))))
idx <- match(ids, layer_ids)
if (anyNA(idx)) stop("WorldClim layers must have names ending bio3/bio_3, etc. Found: ", paste(names(bio), collapse = ", "))
bio <- bio[[idx]]
names(bio) <- paste0("bio", ids)
for (i in seq_len(nrow(datasets))) {
  dat <- read_extraction_points(datasets$Input[i])
  pts <- vect(dat, geom = c("Long", "Lat"), crs = "EPSG:4326")
  pts <- project(pts, crs(bio))
  values <- terra::extract(bio, pts, ID = FALSE)
  out <- cbind(dat, values)
  dir.create(datasets$Extracted[i], recursive = TRUE, showWarnings = FALSE)
  path <- file.path(datasets$Extracted[i], paste0(datasets$Label[i], "_Bioclim.csv"))
  write.csv(out, path, row.names = FALSE)
  message("Saved: ", path)
}
