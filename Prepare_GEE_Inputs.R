source("config/paths.R")
library(readxl)
for (sp in names(conflict_files)) {
  if (!file.exists(conflict_files[[sp]])) stop("Missing conflict workbook: ", conflict_files[[sp]])
  x <- as.data.frame(read_excel(conflict_files[[sp]], sheet = 1))
  required <- c("Long", "Lat", "Year")
  if (!all(required %in% names(x))) stop(sp, " workbook requires Long, Lat and Year columns.")
  if ("Conflict" %in% names(x) && any(is.na(x$Conflict) | x$Conflict != 1)) stop(sp, " conflict workbook contains non-presence rows.")
  for (nm in required) x[[nm]] <- suppressWarnings(as.numeric(as.character(x[[nm]])))
  if (any(!is.finite(x$Long) | !is.finite(x$Lat) | abs(x$Long) > 180 | abs(x$Lat) > 90 | !is.finite(x$Year))) stop(sp, " workbook has invalid coordinates or years.")
  x$Conflict <- 1L
  x$Record_ID <- sprintf("%s_Conflict_%07d", sp, seq_len(nrow(x)))
  folder <- file.path(base_dir, paste(sp, "Data"), "GEE_Input")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  write.csv(x[c("Record_ID", "Long", "Lat", "Year", "Conflict")],
            file.path(folder, paste0(sp, "_Final_Conflict_Data.csv")), row.names = FALSE)
}
message("Conflict CSVs prepared. Preparing the six extraction-input CSVs next.")
check_six_inputs()
for (i in seq_len(nrow(datasets))) {
  x <- read_extraction_points(datasets$Input[i])
  folder <- file.path(base_dir, paste(datasets$Species[i], "Data"), "GEE_Input")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file.path(folder, paste0(datasets$Label[i], ".csv")), row.names = FALSE)
}
message("Upload the GEE_Input CSVs to Earth Engine using the filename without .csv as asset name.")
