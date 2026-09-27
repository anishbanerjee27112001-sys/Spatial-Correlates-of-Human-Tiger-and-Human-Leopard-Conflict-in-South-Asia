packages <- c("readr", "readxl", "dplyr", "tidyr", "purrr", "caret", "pROC",
  "car", "broom", "sf", "blockCV", "tibble", "terra", "stringr", "landscapemetrics")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
