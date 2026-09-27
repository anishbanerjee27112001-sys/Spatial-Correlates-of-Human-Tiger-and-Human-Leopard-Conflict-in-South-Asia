# ============================================================
# HLC GLM
# Human-leopard conflict: binomial generalized linear models
#
# Compares three predictor scales, five spatial cross-validation
# block sizes and three background-sampling methods.
# Preprocessing is fitted within training folds; PR-AUC ranks models.
# Only full models are fitted and exported for the manuscript.
# ============================================================

# ============================================================
# 0. PACKAGES
# ============================================================

required_packages <- c(
  "readr", "dplyr", "tidyr", "purrr", "caret", "pROC",
  "car", "broom", "sf", "blockCV", "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Run Setup_Packages.R first. Missing: ",
       paste(missing_packages, collapse = ", "))
}

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(caret)
library(pROC)
library(car)
library(broom)
library(sf)
library(blockCV)
library(tibble)

# ============================================================
# 1. GLOBAL SETTINGS
# ============================================================

source("config/paths.R")

buffer_options <- c(1, 2, 5)

block_sizes_to_test <- c(10000, 20000, 30000, 40000, 50000)

k_folds_main <- 5
seed_main <- 123

cor_threshold_main <- 0.70
vif_threshold_main <- 5

# ============================================================
# 2. A PRIORI PREDICTOR SET
# ============================================================

base_vars <- c(
  "pop_density_people_per_km2_Final",
  "livestock_density_1km",
  "bio3", "bio4", "bio5", "bio6",
  "bio12", "bio16", "bio17",
  "dist_to_road_km_FINAL",
  "dist_to_settlement_km_FINAL",
  "GLAD_dist_to_forest_km_Final",
  "JRC_dist_to_water_km_Final"
)

get_landcover_vars <- function(buffer_km) {

  c(
    paste0("GLAD_pct_bare_", buffer_km, "km_Final"),
    paste0("GLAD_pct_built_", buffer_km, "km_Final"),
    paste0("GLAD_pct_cropland_", buffer_km, "km_Final"),
    paste0("GLAD_pct_forest_", buffer_km, "km_Final"),
    paste0("GLAD_pct_grassland_", buffer_km, "km_Final"),
    paste0("GLAD_pct_shrubland_", buffer_km, "km_Final"),
    paste0("GLAD_SHDI_Unitless_", buffer_km, "km_FINAL")
  )
}

get_buffer_vars <- function(buffer_km) {

  suffix_m <- ifelse(
    buffer_km == 1,
    "1000",
    ifelse(buffer_km == 2, "2000", "5000")
  )

  c(
    paste0("NDVI_mean_", buffer_km, "km"),
    paste0(
      "road_density_",
      suffix_m,
      "m_buf_km_per_km2"
    ),
    paste0("mean_water_occurrence_", buffer_km, "km"),
    paste0("mean_water_seasonality_", buffer_km, "km"),
    paste0("slope_mean_", buffer_km, "km_buf"),
    paste0("elev_mean_", buffer_km, "km_buf"),
    paste0("aspect_circmean_", buffer_km, "km_buf"),
    paste0(
      "pct_hansen_treecover_loss_area_",
      buffer_km,
      "km"
    ),
    paste0(
      "pct_hansen_baseline_treecover_lost_",
      buffer_km,
      "km"
    )
  )
}

get_fragmentation_vars <- function(buffer_km) {

  c(
    paste0(
      "GLAD_Forest_Patch_Density_patches_per_100ha_",
      buffer_km,
      "km_FINAL"
    ),
    paste0(
      "GLAD_Forest_Mean_Patch_Size_ha_",
      buffer_km,
      "km_FINAL"
    ),
    paste0(
      "GLAD_Forest_Edge_Density_m_per_ha_",
      buffer_km,
      "km_FINAL"
    )
  )
}

make_candidate_predictors <- function(buffer_km) {

  unique(c(
    base_vars,
    get_landcover_vars(buffer_km),
    get_buffer_vars(buffer_km),
    get_fragmentation_vars(buffer_km)
  ))
}

# ============================================================
# 4. GENERAL HELPERS
# ============================================================

bt <- function(x) {
  paste0("`", gsub("`", "``", x), "`")
}

make_formula <- function(response, predictors) {

  if (length(predictors) == 0) {
    stop("No predictors available for model.")
  }

  as.formula(
    paste(
      bt(response),
      "~",
      paste(bt(predictors), collapse = " + ")
    )
  )
}

calc_pr_auc <- function(obs, prob) {

  obs_num <- as.integer(as.character(obs))
  thresholds <- sort(unique(prob), decreasing = TRUE)

  pr <- purrr::map_dfr(thresholds, function(th) {

    pred <- ifelse(prob >= th, 1, 0)

    tp <- sum(pred == 1 & obs_num == 1, na.rm = TRUE)
    fp <- sum(pred == 1 & obs_num == 0, na.rm = TRUE)
    fn <- sum(pred == 0 & obs_num == 1, na.rm = TRUE)

    tibble(
      precision = ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp)),
      recall = ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
    )
  }) %>%
    filter(!is.na(precision), !is.na(recall)) %>%
    arrange(recall)

  if (nrow(pr) < 2) return(NA_real_)

  sum(
    diff(pr$recall) *
      (head(pr$precision, -1) + tail(pr$precision, -1)) / 2,
    na.rm = TRUE
  )
}

# ============================================================
# 5. SPATIALLY WITHHELD PERMUTATION IMPORTANCE
# ============================================================

calculate_permutation_importance <- function(
    model,
    test_data,
    response = "Conflict",
    predictors,
    n_repeats = 30,
    seed = 123
) {

  if (length(predictors) == 0) return(tibble())

  obs <- test_data[[response]]

  # Baseline discrimination on the unchanged held-out spatial fold.
  baseline_prob <- predict(
    model,
    newdata = test_data,
    type = "response"
  )

  baseline_auc <- as.numeric(
    pROC::auc(
      pROC::roc(
        response = obs,
        predictor = baseline_prob,
        levels = c(0, 1),
        direction = "<",
        quiet = TRUE
      )
    )
  )

  purrr::map_dfr(
    seq_along(predictors),
    function(j) {

      predictor_j <- predictors[j]

      purrr::map_dfr(
        seq_len(n_repeats),
        function(r) {

          # Reproducible, predictor- and repeat-specific permutation.
          set.seed(seed + j * 10000 + r)

          permuted_test <- test_data
          permuted_test[[predictor_j]] <- sample(
            permuted_test[[predictor_j]],
            size = nrow(permuted_test),
            replace = FALSE
          )

          perm_prob <- predict(
            model,
            newdata = permuted_test,
            type = "response"
          )

          perm_auc <- as.numeric(
            pROC::auc(
              pROC::roc(
                response = obs,
                predictor = perm_prob,
                levels = c(0, 1),
                direction = "<",
                quiet = TRUE
              )
            )
          )

          tibble(
            predictor = predictor_j,
            permutation_repeat = r,
            baseline_auc = baseline_auc,
            permuted_auc = perm_auc,
            delta_auc = baseline_auc - perm_auc
          )
        }
      )
    }
  )
}

# ============================================================
# 6. RAW MODEL DATA
# ============================================================

prepare_raw_model_data <- function(df, predictors) {

  required_core <- c("Conflict", "Lat", "Long", "Year")

  missing_core <- setdiff(required_core, names(df))

  if (length(missing_core) > 0) {
    stop(
      "Missing core columns: ",
      paste(missing_core, collapse = ", ")
    )
  }

  missing_predictors <- setdiff(predictors, names(df))

  if (length(missing_predictors) > 0) {
    stop(
      "Missing required predictor columns: ",
      paste(missing_predictors, collapse = ", ")
    )
  }

  dat <- df %>%
    select(all_of(c(required_core, predictors))) %>%
    mutate(
      Year = suppressWarnings(
        as.integer(substr(as.character(Year), 1, 4))
      )
    ) %>%
    filter(
      !is.na(Lat),
      !is.na(Long),
      is.finite(Lat),
      is.finite(Long),
      !is.na(Conflict),

      # Restrict BOTH conflict and background records to the
      # study period beginning in 2000.
      !is.na(Year),
      Year >= 2000
    )

  dat$Conflict <- factor(dat$Conflict, levels = c(0, 1))

  # Safety check: no pre-2000 observation can enter a model.
  if (any(dat$Year < 2000, na.rm = TRUE)) {
    stop("Pre-2000 observations remain after Year >= 2000 filtering.")
  }

  dat
}

# ============================================================
# 7. CIRCULAR ASPECT
# ============================================================

transform_aspect <- function(dat, aspect_col) {

  if (!(aspect_col %in% names(dat))) {
    stop("Aspect column missing: ", aspect_col)
  }

  aspect_rad <- dat[[aspect_col]] * pi / 180

  dat[[paste0(aspect_col, "_northness")]] <- cos(aspect_rad)
  dat[[paste0(aspect_col, "_eastness")]] <- sin(aspect_rad)

  dat %>%
    select(-all_of(aspect_col))
}

# ============================================================
# 8. TRAINING-FOLD IMPUTATION
# ============================================================

fit_imputation_params <- function(train, response = "Conflict") {

  pred_names <- setdiff(
    names(train),
    c(response, "Lat", "Long", "Year")
  )

  num_names <- pred_names[
    vapply(train[, pred_names, drop = FALSE], is.numeric, logical(1))
  ]

  if (length(num_names) == 0) return(numeric())

  vapply(
    train[, num_names, drop = FALSE],
    function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) return(NA_real_)
      median(x, na.rm = TRUE)
    },
    numeric(1)
  )
}

apply_imputation <- function(dat, medians) {

  out <- dat

  for (v in names(medians)) {

    if (!(v %in% names(out))) next

    bad <- is.na(out[[v]]) | !is.finite(out[[v]])
    out[[v]][bad] <- medians[[v]]
  }

  out
}

# ============================================================
# 9. NEAR-ZERO VARIANCE
# ============================================================

find_nzv_drops <- function(train, response = "Conflict") {

  x <- train %>%
    select(-any_of(c(response, "Lat", "Long", "Year")))

  if (ncol(x) == 0) return(character())

  nzv <- caret::nearZeroVar(
    x,
    saveMetrics = TRUE
  )

  rownames(
    nzv[nzv$nzv, , drop = FALSE]
  )
}

# ============================================================
# 10. CORRELATION FILTER
# ============================================================

find_correlation_drops <- function(
    train,
    response = "Conflict",
    cutoff = 0.70
) {

  x <- train %>%
    select(-any_of(c(response, "Lat", "Long", "Year")))

  numeric_names <- names(x)[
    vapply(x, is.numeric, logical(1))
  ]

  if (length(numeric_names) <= 1) {
    return(character())
  }

  current <- numeric_names
  dropped <- character()

  repeat {

    if (length(current) <= 1) break

    cor_mat <- suppressWarnings(
      cor(
        x[, current, drop = FALSE],
        use = "pairwise.complete.obs"
      )
    )

    diag(cor_mat) <- 0

    cor_abs <- abs(cor_mat)
    cor_abs[!is.finite(cor_abs)] <- 0

    max_cor <- max(cor_abs, na.rm = TRUE)

    if (!is.finite(max_cor) || max_cor < cutoff) break

    # Identify the most highly correlated pair currently retained.
    idx <- which(
      cor_abs == max_cor,
      arr.ind = TRUE
    )[1, ]

    v1 <- rownames(cor_abs)[idx[1]]
    v2 <- colnames(cor_abs)[idx[2]]

    others_v1 <- setdiff(current, v1)
    others_v2 <- setdiff(current, v2)

    mean_cor_v1 <- if (length(others_v1) > 0) {
      mean(cor_abs[v1, others_v1], na.rm = TRUE)
    } else {
      0
    }

    mean_cor_v2 <- if (length(others_v2) > 0) {
      mean(cor_abs[v2, others_v2], na.rm = TRUE)
    } else {
      0
    }

    if (mean_cor_v1 > mean_cor_v2) {

      drop_var <- v1

    } else if (mean_cor_v2 > mean_cor_v1) {

      drop_var <- v2

    } else {

      # Deterministic final tie-break for exact equality.
      drop_var <- sort(c(v1, v2))[2]
    }

    dropped <- c(dropped, drop_var)
    current <- setdiff(current, drop_var)
  }

  unique(dropped)
}

# ============================================================
# 11. TRAINING-FOLD SCALING
# ============================================================

fit_scaling_params <- function(
    train,
    response = "Conflict"
) {

  pred_names <- setdiff(
    names(train),
    c(response, "Lat", "Long", "Year")
  )

  num_names <- pred_names[
    vapply(train[, pred_names, drop = FALSE], is.numeric, logical(1))
  ]

  if (length(num_names) == 0) {
    return(
      tibble(
        variable = character(),
        mean = numeric(),
        sd = numeric()
      )
    )
  }

  tibble(
    variable = num_names,
    mean = vapply(
      train[, num_names, drop = FALSE],
      mean,
      numeric(1),
      na.rm = TRUE
    ),
    sd = vapply(
      train[, num_names, drop = FALSE],
      sd,
      numeric(1),
      na.rm = TRUE
    )
  ) %>%
    mutate(
      sd = ifelse(!is.finite(sd) | sd == 0, 1, sd)
    )
}

apply_scaling <- function(dat, scaling_params) {

  out <- dat

  if (nrow(scaling_params) == 0) return(out)

  for (i in seq_len(nrow(scaling_params))) {

    v <- scaling_params$variable[i]

    if (!(v %in% names(out))) next

    out[[v]] <- (
      out[[v]] - scaling_params$mean[i]
    ) / scaling_params$sd[i]
  }

  out
}

# ============================================================
# 12. VIF FILTER
# ============================================================

reduce_vif <- function(
    df,
    response = "Conflict",
    threshold = 5
) {

  current_df <- df
  dropped <- character()

  repeat {

    predictors <- setdiff(names(current_df), response)

    if (length(predictors) <= 1) break

    model <- tryCatch(
      suppressWarnings(
        glm(
          make_formula(response, predictors),
          data = current_df,
          family = binomial,
          control = glm.control(maxit = 100)
        )
      ),
      error = function(e) NULL
    )

    if (is.null(model)) break

    vif_vals <- tryCatch(
      car::vif(model),
      error = function(e) NULL
    )

    if (is.null(vif_vals)) break

    if (is.matrix(vif_vals)) {
      vif_vals <- vif_vals[, ncol(vif_vals)]
    }

    vif_vals <- vif_vals[is.finite(vif_vals)]

    if (length(vif_vals) == 0) break

    max_vif <- max(vif_vals, na.rm = TRUE)

    if (!is.finite(max_vif) || max_vif < threshold) break

    # Remove the predictor with the greatest remaining VIF.
    drop_var <- names(which.max(vif_vals))

    dropped <- c(dropped, drop_var)

    current_df <- current_df %>%
      select(-all_of(drop_var))
  }

  list(
    data = current_df,
    dropped = unique(dropped)
  )
}

# ============================================================
# 13. FIT TRAINING PREPROCESSOR
# ============================================================

fit_fold_preprocessor <- function(
    train_raw,
    buffer_km,
    response = "Conflict",
    cor_cutoff = 0.70,
    vif_threshold = 5
) {

  aspect_col <- paste0(
    "aspect_circmean_",
    buffer_km,
    "km_buf"
  )

  # Circular aspect transformation
  train <- transform_aspect(
    train_raw,
    aspect_col
  )

  # Training-only imputation
  imputation_params <- fit_imputation_params(
    train,
    response = response
  )

  train <- apply_imputation(
    train,
    imputation_params
  )

  # Drop variables completely missing in the training fold
  all_missing_vars <- names(imputation_params)[
    is.na(imputation_params)
  ]

  if (length(all_missing_vars) > 0) {
    train <- train %>%
      select(-any_of(all_missing_vars))
  }

  # Training-only near-zero variance
  nzv_drops <- find_nzv_drops(
    train,
    response = response
  )

  if (length(nzv_drops) > 0) {
    train <- train %>%
      select(-any_of(nzv_drops))
  }

  # Training-only correlation filtering
  cor_drops <- find_correlation_drops(
    train,
    response = response,
    cutoff = cor_cutoff
  )

  if (length(cor_drops) > 0) {
    train <- train %>%
      select(-any_of(cor_drops))
  }

  # Training-only scaling
  scaling_params <- fit_scaling_params(
    train,
    response = response
  )

  train <- apply_scaling(
    train,
    scaling_params
  )

  # Remove non-predictive metadata before VIF
  train_model <- train %>%
    select(-any_of(c("Lat", "Long", "Year")))

  # Training-only VIF
  vif_result <- reduce_vif(
    train_model,
    response = response,
    threshold = vif_threshold
  )

  train_model <- vif_result$data
  vif_drops <- vif_result$dropped

  final_predictors <- setdiff(
    names(train_model),
    response
  )

  list(
    train_model = train_model,
    aspect_col = aspect_col,
    imputation_params = imputation_params,
    all_missing_vars = all_missing_vars,
    nzv_drops = nzv_drops,
    cor_drops = cor_drops,
    scaling_params = scaling_params,
    vif_drops = vif_drops,
    final_predictors = final_predictors
  )
}

apply_fold_preprocessor <- function(
    test_raw,
    prep,
    response = "Conflict"
) {

  test <- transform_aspect(
    test_raw,
    prep$aspect_col
  )

  test <- apply_imputation(
    test,
    prep$imputation_params
  )

  test <- test %>%
    select(-any_of(prep$all_missing_vars)) %>%
    select(-any_of(prep$nzv_drops)) %>%
    select(-any_of(prep$cor_drops)) %>%
    select(-any_of(prep$vif_drops))

  test <- apply_scaling(
    test,
    prep$scaling_params
  )

  test_model <- test %>%
    select(
      any_of(c(response, prep$final_predictors))
    )

  missing_after <- setdiff(
    prep$final_predictors,
    names(test_model)
  )

  if (length(missing_after) > 0) {
    stop(
      "Test fold missing retained predictors: ",
      paste(missing_after, collapse = ", ")
    )
  }

  test_model
}

# ============================================================
# 14. METRICS
# ============================================================

calculate_fold_metrics <- function(
    obs,
    pred_prob
) {

  roc_obj <- pROC::roc(
    response = obs,
    predictor = pred_prob,
    quiet = TRUE
  )

  auc_val <- as.numeric(
    pROC::auc(roc_obj)
  )

  pr_auc <- calc_pr_auc(
    obs,
    pred_prob
  )

  obs_num <- as.integer(
    as.character(obs)
  )

  brier <- mean(
    (pred_prob - obs_num)^2,
    na.rm = TRUE
  )

  # Clamp probabilities for numerically stable log-loss.
  eps <- 1e-15
  p_clip <- pmin(
    pmax(pred_prob, eps),
    1 - eps
  )

  logloss <- -mean(
    obs_num * log(p_clip) +
      (1 - obs_num) * log(1 - p_clip),
    na.rm = TRUE
  )

  tibble(
    auc = auc_val,
    pr_auc = pr_auc,
    brier = brier,
    logloss = logloss
  )
}

# ============================================================
# 15. SPATIAL CV GLM
# ============================================================

fit_block_cv_glm <- function(
    df,
    predictors,
    buffer_km,
    block_size,
    k = 5,
    seed = 123,
    response = "Conflict",
    cor_cutoff = 0.70,
    vif_threshold = 5,
    calculate_importance = FALSE,
    permutation_repeats = 30
) {

  set.seed(seed)

  raw_df <- prepare_raw_model_data(
    df,
    predictors
  )

  sf_df <- st_as_sf(
    raw_df,
    coords = c("Long", "Lat"),
    crs = 4326,
    remove = FALSE
  )

  sb <- blockCV::cv_spatial(
    x = sf_df,
    column = response,
    k = k,
    size = block_size,
    hexagon = FALSE,
    selection = "random",
    iteration = 100,  # Candidate fold assignments; not repeated model CV.
    progress = FALSE
  )

  fold_metrics <- list()
  fold_predictors <- list()
  fold_coefficients <- list()
  fold_preprocessing <- list()
  fold_predictions <- list()
  fold_permutation_importance <- list()

  for (i in seq_along(sb$folds_list)) {

    train_ids <- sb$folds_list[[i]][[1]]
    test_ids <- sb$folds_list[[i]][[2]]

    train_raw <- sf_df[train_ids, ] %>%
      st_drop_geometry()

    test_raw <- sf_df[test_ids, ] %>%
      st_drop_geometry()

    # AUC is undefined if one class is absent
    if (
      length(unique(train_raw[[response]])) < 2 ||
      length(unique(test_raw[[response]])) < 2
    ) {
      warning(
        "Skipping fold ",
        i,
        ": one response class is absent."
      )
      next
    }

    prep <- fit_fold_preprocessor(
      train_raw = train_raw,
      buffer_km = buffer_km,
      response = response,
      cor_cutoff = cor_cutoff,
      vif_threshold = vif_threshold
    )

    train <- prep$train_model

    test <- apply_fold_preprocessor(
      test_raw = test_raw,
      prep = prep,
      response = response
    )

    predictors_final <- prep$final_predictors

    if (length(predictors_final) == 0) {
      warning(
        "Skipping fold ",
        i,
        ": no predictors retained."
      )
      next
    }

    glm_fit <- suppressWarnings(
      glm(
        make_formula(
          response,
          predictors_final
        ),
        data = train,
        family = binomial,
        control = glm.control(maxit = 100)
      )
    )

    pred_prob <- predict(
      glm_fit,
      newdata = test,
      type = "response"
    )

    # Optional out-of-fold permutation importance. Permutation is performed
    # only in the held-out spatial fold, after training-derived preprocessing.
    if (calculate_importance) {
      fold_permutation_importance[[i]] <-
        calculate_permutation_importance(
          model = glm_fit,
          test_data = test,
          response = response,
          predictors = predictors_final,
          n_repeats = permutation_repeats,
          seed = seed + i * 100000
        ) %>%
        mutate(
          fold = i,
          .before = 1
        )
    }

    metric_row <- calculate_fold_metrics(
      obs = test[[response]],
      pred_prob = pred_prob
    ) %>%
      mutate(
        fold = i,
        n_train = nrow(train),
        n_test = nrow(test),
        n_predictors = length(predictors_final),
        .before = 1
      )

    fold_metrics[[i]] <- metric_row

    # Save out-of-fold predictions for ROC/PR/calibration/figure creation.
    fold_predictions[[i]] <- tibble(
      fold = i,
      observed = as.integer(as.character(test[[response]])),
      predicted_probability = pred_prob
    )

    fold_predictors[[i]] <- tibble(
      fold = i,
      predictor = predictors_final
    )

    fold_coefficients[[i]] <- broom::tidy(
      glm_fit
    ) %>%
      filter(term != "(Intercept)") %>%
      transmute(
        fold = i,
        predictor = term,
        estimate = estimate,
        std_error = std.error,
        odds_ratio = exp(estimate),
        sign = case_when(
          estimate > 0 ~ 1L,
          estimate < 0 ~ -1L,
          TRUE ~ 0L
        )
      )

    fold_preprocessing[[i]] <- bind_rows(
      tibble(
        fold = i,
        stage = "NZV",
        variable = prep$nzv_drops
      ),
      tibble(
        fold = i,
        stage = "Correlation",
        variable = prep$cor_drops
      ),
      tibble(
        fold = i,
        stage = "VIF",
        variable = prep$vif_drops
      )
    )
  }

  metrics_df <- bind_rows(
    fold_metrics
  )

  if (nrow(metrics_df) == 0) {
    stop(
      "No valid spatial CV folds were fitted."
    )
  }

  predictor_stability <- bind_rows(
    fold_predictors
  ) %>%
    count(
      predictor,
      name = "folds_retained"
    ) %>%
    mutate(
      retention_proportion =
        folds_retained / nrow(metrics_df)
    ) %>%
    arrange(
      desc(retention_proportion),
      predictor
    )

  coefficient_stability <- bind_rows(
    fold_coefficients
  ) %>%
    group_by(predictor) %>%
    summarise(
      folds_estimated = n(),
      mean_beta = mean(estimate, na.rm = TRUE),
      sd_beta = sd(estimate, na.rm = TRUE),
      median_beta = median(estimate, na.rm = TRUE),
      min_beta = min(estimate, na.rm = TRUE),
      max_beta = max(estimate, na.rm = TRUE),
      mean_odds_ratio = mean(odds_ratio, na.rm = TRUE),
      positive_folds = sum(sign > 0, na.rm = TRUE),
      negative_folds = sum(sign < 0, na.rm = TRUE),
      dominant_sign_proportion = max(
        positive_folds,
        negative_folds
      ) / folds_estimated,
      .groups = "drop"
    ) %>%
    mutate(
      coefficient_stable =
        folds_estimated >= 3 &
        dominant_sign_proportion >= 0.80
    ) %>%
    arrange(
      desc(folds_estimated),
      desc(dominant_sign_proportion),
      predictor
    )

  preprocessing_log <- bind_rows(
    fold_preprocessing
  )

  permutation_importance_raw <- bind_rows(
    fold_permutation_importance
  )

  if (nrow(permutation_importance_raw) > 0) {

    permutation_importance_summary <-
      permutation_importance_raw %>%
      group_by(predictor) %>%
      summarise(
        folds_evaluated = n_distinct(fold),
        n_permutations = n(),
        mean_delta_auc = mean(delta_auc, na.rm = TRUE),
        sd_delta_auc = sd(delta_auc, na.rm = TRUE),
        median_delta_auc = median(delta_auc, na.rm = TRUE),
        q025_delta_auc = quantile(delta_auc, 0.025, na.rm = TRUE),
        q975_delta_auc = quantile(delta_auc, 0.975, na.rm = TRUE),
        mean_baseline_auc = mean(baseline_auc, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_delta_auc))

  } else {

    permutation_importance_summary <- tibble()
  }

  # ----------------------------------------------------------
  # Final full-data model for coefficient reporting only
  # ----------------------------------------------------------

  full_prep <- fit_fold_preprocessor(
    train_raw = raw_df,
    buffer_km = buffer_km,
    response = response,
    cor_cutoff = cor_cutoff,
    vif_threshold = vif_threshold
  )

  final_df <- full_prep$train_model
  final_predictors <- full_prep$final_predictors

  final_model <- glm(
    make_formula(
      response,
      final_predictors
    ),
    data = final_df,
    family = binomial,
    control = glm.control(maxit = 100)
  )

  list(
    fold_metrics = metrics_df,
    fold_predictions = bind_rows(fold_predictions),
    predictor_stability = predictor_stability,
    coefficient_stability = coefficient_stability,
    preprocessing_log = preprocessing_log,

    final_model = final_model,
    final_predictors = final_predictors,

    # Retained for later model-predicted/partial-response figures.
    # Predictors are on the same processed/scaled scale used by the fitted GLM.
    final_model_data = final_df,

    permutation_importance_raw = permutation_importance_raw,
    permutation_importance_summary = permutation_importance_summary,

    final_scaling_params = full_prep$scaling_params,
    final_imputation_params = full_prep$imputation_params,
    final_nzv_drops = full_prep$nzv_drops,
    final_cor_drops = full_prep$cor_drops,
    final_vif_drops = full_prep$vif_drops,

    summary = metrics_df %>%
      summarise(
        mean_auc = mean(auc, na.rm = TRUE),
        sd_auc = sd(auc, na.rm = TRUE),
        mean_pr_auc = mean(pr_auc, na.rm = TRUE),
        sd_pr_auc = sd(pr_auc, na.rm = TRUE),
        mean_brier = mean(brier, na.rm = TRUE),
        sd_brier = sd(brier, na.rm = TRUE),
        mean_logloss = mean(logloss, na.rm = TRUE),
        sd_logloss = sd(logloss, na.rm = TRUE),
        successful_folds = n()
      )
  )
}

# ============================================================
# 16. RUN ONE SPECIES
# ============================================================

run_species_analysis <- function(
    analysis_name,
    data_dir,
    output_dir
) {

  if (!dir.exists(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE
    )
  }

  background_files <- c(
    AccessibleArea = file.path(
      data_dir, "Background",
      paste0(
        analysis_name,
        "_AccessibleArea_Merged_Final.csv"
      )
    ),
    BiasCorrected = file.path(
      data_dir, "Background",
      paste0(
        analysis_name,
        "_BiasCorrected_Merged_Final.csv"
      )
    ),
    UniformRandom = file.path(
      data_dir, "Background",
      paste0(
        analysis_name,
        "_UniformRandom_Merged_Final.csv"
      )
    )
  )

  missing_files <- background_files[
    !file.exists(background_files)
  ]

  if (length(missing_files) > 0) {
    stop(
      "Missing files:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  }

  data_list <- lapply(
    background_files,
    function(f) {
      read_csv(
        f,
        show_col_types = FALSE,
        na = c(
          "",
          "NA",
          "NaN",
          "NULL",
          "null"
        )
      )
    }
  )

  selection_results <- tibble()

  coefficient_tables <- list()
  fold_metric_tables <- list()
  predictor_stability_tables <- list()
  fold_prediction_tables <- list()
  final_model_objects <- list()

  for (background_name in names(data_list)) {

    for (buffer_km in buffer_options) {

      candidate_predictors <- make_candidate_predictors(
        buffer_km = buffer_km
      )

      # Strict consistency check
      missing_predictors <- setdiff(
        candidate_predictors,
        names(data_list[[background_name]])
      )

      if (length(missing_predictors) > 0) {
        stop(
          analysis_name,
          " / ",
          background_name,
          " / ",
          buffer_km,
          " km is missing required columns:\n",
          paste(
            missing_predictors,
            collapse = "\n"
          )
        )
      }

      for (block_size in block_sizes_to_test) {

          model_type <- "Full"

          model_id <- paste(
            analysis_name,
            background_name,
            paste0(buffer_km, "km"),
            paste0(block_size / 1000, "kmBlock"),
            model_type,
            sep = "_"
          )

          cat(
            "\nRunning:",
            model_id,
            "\n"
          )

          res <- fit_block_cv_glm(
            df = data_list[[background_name]],
            predictors = candidate_predictors,
            buffer_km = buffer_km,
            block_size = block_size,
            k = k_folds_main,
            seed = seed_main,
            cor_cutoff = cor_threshold_main,
            vif_threshold = vif_threshold_main
          )

          fold_metric_tables[[model_id]] <-
            res$fold_metrics %>%
            mutate(
              Model_ID = model_id,
              Background = background_name,
              Buffer_km = buffer_km,
              Block_Size_m = block_size,
              Model_Type = model_type
            )

          fold_prediction_tables[[model_id]] <-
            res$fold_predictions %>%
            mutate(
              Model_ID = model_id,
              Background = background_name,
              Buffer_km = buffer_km,
              Block_Size_m = block_size,
              Model_Type = model_type
            )

          predictor_stability_tables[[model_id]] <-
            res$predictor_stability %>%
            mutate(
              Model_ID = model_id,
              Background = background_name,
              Buffer_km = buffer_km,
              Block_Size_m = block_size,
              Model_Type = model_type
            )

          coef_tbl <- broom::tidy(
            res$final_model
          ) %>%
            filter(
              term != "(Intercept)"
            ) %>%
            mutate(
              Model_ID = model_id,
              Background = background_name,
              Buffer_km = buffer_km,
              Block_Size_m = block_size,
              Model_Type = model_type,

              # Continuous predictors were standardized.
              # For linear terms OR is per +1 SD.
              odds_ratio = exp(estimate),
              conf_low = exp(
                estimate - 1.96 * std.error
              ),
              conf_high = exp(
                estimate + 1.96 * std.error
              ),

              effect_direction = case_when(
                estimate > 0 ~ "Positive",
                estimate < 0 ~ "Negative",
                TRUE ~ "Zero"
              )
            )

          coefficient_tables[[model_id]] <-
            coef_tbl

          # Save model object + preprocessing information for later figures
          # or alternative prediction grids.
          final_model_objects[[model_id]] <- list(
            model = res$final_model,
            final_predictors = res$final_predictors,
            final_model_data = res$final_model_data,
            scaling_params = res$final_scaling_params,
            imputation_params = res$final_imputation_params,
            nzv_drops = res$final_nzv_drops,
            correlation_drops = res$final_cor_drops,
            vif_drops = res$final_vif_drops,
            buffer_km = buffer_km,
            block_size_m = block_size,
            background = background_name,
            model_type = model_type
          )

          selection_results <- bind_rows(
            selection_results,
            tibble(
              Model_ID = model_id,
              Background = background_name,
              Buffer_km = buffer_km,
              Block_Size_m = block_size,
              Model_Type = model_type,

              N_Candidate_Predictors =
                length(candidate_predictors),

              N_Final_Predictors =
                length(res$final_predictors),

              Final_Predictors = paste(
                res$final_predictors,
                collapse = "; "
              ),

              Mean_AUC =
                res$summary$mean_auc,

              SD_AUC =
                res$summary$sd_auc,

              Mean_PR_AUC =
                res$summary$mean_pr_auc,

              SD_PR_AUC =
                res$summary$sd_pr_auc,

              Mean_Brier =
                res$summary$mean_brier,

              SD_Brier =
                res$summary$sd_brier,

              Mean_LogLoss =
                res$summary$mean_logloss,

              SD_LogLoss =
                res$summary$sd_logloss,

              Successful_Folds =
                res$summary$successful_folds
            )
          )
      }
    }
  }

  all_model_results <- selection_results

  scale_block_summary <- all_model_results %>%
    group_by(Buffer_km, Block_Size_m) %>%
    summarise(
      Mean_PR_AUC_Across_Backgrounds =
        mean(Mean_PR_AUC, na.rm = TRUE),
      Mean_AUC_Across_Backgrounds =
        mean(Mean_AUC, na.rm = TRUE),
      Mean_Brier_Across_Backgrounds =
        mean(Mean_Brier, na.rm = TRUE),
      Mean_LogLoss_Across_Backgrounds =
        mean(Mean_LogLoss, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Block_km = Block_Size_m / 1000
    ) %>%
    select(
      Buffer_km,
      Block_km,
      Mean_PR_AUC_Across_Backgrounds,
      Mean_AUC_Across_Backgrounds,
      Mean_Brier_Across_Backgrounds,
      Mean_LogLoss_Across_Backgrounds
    ) %>%
    arrange(
      desc(Mean_PR_AUC_Across_Backgrounds),
      desc(Mean_AUC_Across_Backgrounds),
      Mean_Brier_Across_Backgrounds
    )

  optimal_setting <- scale_block_summary %>%
    slice(1)

  optimal_buffer_km <- optimal_setting$Buffer_km[[1]]
  optimal_block_m <- optimal_setting$Block_km[[1]] * 1000

  # Direct background comparison at a COMMON optimal buffer + block.
  background_comparison_optimal <- all_model_results %>%
    filter(
      Buffer_km == optimal_buffer_km,
      Block_Size_m == optimal_block_m
    ) %>%
    arrange(
      # Keep the same ranking hierarchy used for scale/block selection:
      # PR-AUC primary, AUC secondary, Brier score tertiary.
      desc(Mean_PR_AUC),
      desc(Mean_AUC),
      Mean_Brier
    )

  best_background_optimal <- background_comparison_optimal %>%
    slice(1)

  best_model_id <- best_background_optimal$Model_ID[[1]]
  best_background_name <- best_background_optimal$Background[[1]]

  backgrounds_opt <- background_comparison_optimal$Background

  permutation_importance_optimal <- list()

  for (background_name in backgrounds_opt) {

    candidate_predictors_opt <- make_candidate_predictors(
      optimal_buffer_km
    )

    importance_result <- fit_block_cv_glm(
      df = data_list[[background_name]],
      predictors = candidate_predictors_opt,
      buffer_km = optimal_buffer_km,
      block_size = optimal_block_m,
      k = k_folds_main,
      seed = seed_main,
      cor_cutoff = cor_threshold_main,
      vif_threshold = vif_threshold_main,
      calculate_importance = TRUE,
      permutation_repeats = 30
    )

    permutation_importance_optimal[[background_name]] <-
      importance_result$permutation_importance_raw %>%
      mutate(
        Background = background_name,
        Buffer_km = optimal_buffer_km,
        Block_km = optimal_block_m / 1000,
        .before = 1
      )
  }

  figure3_permutation_raw <- bind_rows(
    permutation_importance_optimal
  )

  # First average repeated permutations within each spatial fold. This prevents
  # repeated permutations from being mistaken for independent CV replicates.
  figure3_permutation_fold <- figure3_permutation_raw %>%
    group_by(
      Background,
      Buffer_km,
      Block_km,
      fold,
      predictor
    ) %>%
    summarise(
      baseline_auc = first(.data$baseline_auc),
      mean_permuted_auc = mean(.data$permuted_auc, na.rm = TRUE),
      fold_mean_delta_auc = mean(.data$delta_auc, na.rm = TRUE),
      sd_delta_auc_within_fold = stats::sd(.data$delta_auc, na.rm = TRUE),
      permutation_repeats = dplyr::n_distinct(.data$permutation_repeat),
      .groups = "drop"
    )

  # Background-specific importance summarised across spatial folds.
  figure3_permutation_by_background <- figure3_permutation_fold %>%
    group_by(
      Background,
      Buffer_km,
      Block_km,
      predictor
    ) %>%
    summarise(
      folds_evaluated = n_distinct(.data$fold),
      sd_delta_auc_across_folds = if (
        n_distinct(.data$fold) > 1
      ) {
        stats::sd(.data$fold_mean_delta_auc, na.rm = TRUE)
      } else {
        NA_real_
      },
      median_delta_auc = median(.data$fold_mean_delta_auc, na.rm = TRUE),
      min_delta_auc = min(.data$fold_mean_delta_auc, na.rm = TRUE),
      max_delta_auc = max(.data$fold_mean_delta_auc, na.rm = TRUE),
      mean_baseline_auc = mean(.data$baseline_auc, na.rm = TRUE),
      mean_delta_auc = mean(.data$fold_mean_delta_auc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(
      Background,
      desc(mean_delta_auc)
    )

  figure3_permutation_consensus <- figure3_permutation_by_background %>%
    group_by(predictor) %>%
    summarise(
      backgrounds_evaluated = n_distinct(.data$Background),
      total_folds_evaluated = sum(.data$folds_evaluated),
      sd_delta_auc_across_backgrounds = if (
        n_distinct(.data$Background) > 1
      ) {
        stats::sd(.data$mean_delta_auc, na.rm = TRUE)
      } else {
        NA_real_
      },
      min_background_delta_auc = min(.data$mean_delta_auc, na.rm = TRUE),
      max_background_delta_auc = max(.data$mean_delta_auc, na.rm = TRUE),
      mean_delta_auc = mean(.data$mean_delta_auc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_delta_auc))

  optimal_full_ids <- background_comparison_optimal$Model_ID

  optimal_coef_full <- bind_rows(
    coefficient_tables[optimal_full_ids]
  ) %>%
    transmute(
      Model_ID,
      Background,
      variable = term,
      Estimate = estimate,
      SE = std.error,
      p = p.value,
      OR = odds_ratio,
      Direction = case_when(
        estimate > 0 ~ "Positive",
        estimate < 0 ~ "Negative",
        TRUE ~ "Zero"
      ),
      Significant = as.integer(p.value < 0.05)
    )

  aspect_raw_opt <- paste0(
    "aspect_circmean_",
    optimal_buffer_km,
    "km_buf"
  )

  candidate_universe <- make_candidate_predictors(
    optimal_buffer_km
  )

  candidate_universe <- setdiff(
    candidate_universe,
    aspect_raw_opt
  )

  candidate_universe <- unique(c(
    candidate_universe,
    paste0(aspect_raw_opt, "_northness"),
    paste0(aspect_raw_opt, "_eastness")
  ))

  backgrounds_opt <- background_comparison_optimal$Background

  consensus_long <- tidyr::expand_grid(
    variable = candidate_universe,
    Background = backgrounds_opt
  ) %>%
    left_join(
      optimal_coef_full,
      by = c("variable", "Background")
    ) %>%
    mutate(
      Retained = as.integer(!is.na(Estimate)),
      Positive = as.integer(Retained == 1 & Estimate > 0),
      Negative = as.integer(Retained == 1 & Estimate < 0),
      Zero = as.integer(Retained == 1 & Estimate == 0),
      Significant = replace_na(Significant, 0L)
    )

  variable_consensus_base <- consensus_long %>%
    group_by(variable) %>%
    summarise(
      Total_models = n(),
      Retained_n = sum(Retained),
      Retained_prop = mean(Retained),

      Positive_n = sum(Positive),
      Negative_n = sum(Negative),
      Zero_n = sum(Zero),

      Positive_prop_when_retained =
        ifelse(Retained_n > 0, Positive_n / Retained_n, NA_real_),

      Negative_prop_when_retained =
        ifelse(Retained_n > 0, Negative_n / Retained_n, NA_real_),

      Direction_consistency =
        ifelse(
          Retained_n > 0,
          pmax(Positive_n, Negative_n, Zero_n) / Retained_n,
          NA_real_
        ),

      Dominant_direction = case_when(
        Retained_n == 0 ~ NA_character_,
        Positive_n == pmax(Positive_n, Negative_n, Zero_n) ~ "Positive",
        Negative_n == pmax(Positive_n, Negative_n, Zero_n) ~ "Negative",
        TRUE ~ "Zero"
      ),

      Significant_n = sum(Significant),
      Significant_prop_all_models =
        Significant_n / Total_models,

      Significant_prop_when_retained =
        ifelse(
          Retained_n > 0,
          Significant_n / Retained_n,
          NA_real_
        ),

      Mean_Estimate =
        ifelse(Retained_n > 0, mean(Estimate, na.rm = TRUE), NA_real_),

      SD_Estimate =
        ifelse(Retained_n > 1, sd(Estimate, na.rm = TRUE), NA_real_),

      Median_Estimate =
        ifelse(Retained_n > 0, median(Estimate, na.rm = TRUE), NA_real_),

      Mean_OR =
        ifelse(Retained_n > 0, mean(OR, na.rm = TRUE), NA_real_),

      SD_OR =
        ifelse(Retained_n > 1, sd(OR, na.rm = TRUE), NA_real_),

      Minimum_p =
        ifelse(Retained_n > 0, min(p, na.rm = TRUE), NA_real_),

      Maximum_p =
        ifelse(Retained_n > 0, max(p, na.rm = TRUE), NA_real_),

      Stability = case_when(
        Retained_prop == 1 &
          Direction_consistency == 1 ~ "Very stable",

        Retained_prop >= 2/3 &
          Direction_consistency >= 0.80 ~ "Stable",

        Retained_prop >= 1/3 ~ "Variable",

        TRUE ~ "Not retained"
      ),

      .groups = "drop"
    )

  # Background-specific wide columns matching the reference file style.
  consensus_wide <- consensus_long %>%
    select(
      variable,
      Background,
      Retained,
      Estimate,
      SE,
      p,
      OR,
      Direction,
      Significant
    ) %>%
    pivot_wider(
      names_from = Background,
      values_from = c(
        Retained,
        Estimate,
        SE,
        p,
        OR,
        Direction,
        Significant
      ),
      names_glue = "{.value}_{Background}"
    )

  variable_consensus <- variable_consensus_base %>%
    left_join(
      consensus_wide,
      by = "variable"
    ) %>%
    arrange(
      desc(Retained_prop),
      desc(Direction_consistency),
      desc(Significant_prop_when_retained),
      variable
    )

  # EXPORTS — CORE MODEL FILES + FIGURE 3 IMPORTANCE FILES

  all_fold_performance <- bind_rows(fold_metric_tables)

  all_oof_predictions <- bind_rows(fold_prediction_tables)

  optimal_full_coefficients <- bind_rows(
    coefficient_tables[optimal_full_ids]
  )

  final_coefficients <- optimal_full_coefficients

  all_retention <- bind_rows(predictor_stability_tables)

  model_objects_to_save <- list(
    full_models = final_model_objects,
    optimal_buffer_km = optimal_buffer_km,
    optimal_block_km = optimal_block_m / 1000,
    best_background = best_background_name
  )

  # Manuscript outputs contain full models only.
  write_csv(selection_results,
            file.path(output_dir, "01_AllModelComparison.csv"))
  write_csv(scale_block_summary,
            file.path(output_dir, "02_ScaleBlockRanking.csv"))
  write_csv(optimal_setting,
            file.path(output_dir, "03_OptimalScaleBlock.csv"))
  write_csv(background_comparison_optimal,
            file.path(output_dir, "04_BackgroundComparison_Optimal.csv"))
  write_csv(all_fold_performance,
            file.path(output_dir, "05_FoldPerformance.csv"))
  write_csv(final_coefficients,
            file.path(output_dir, "06_FinalModelCoefficients.csv"))
  write_csv(all_retention,
            file.path(output_dir, "07_PredictorRetention.csv"))
  write_csv(variable_consensus,
            file.path(output_dir, "08_VariableConsensus.csv"))
  write_csv(all_oof_predictions,
            file.path(output_dir, "11_OutOfFoldPredictions.csv"))
  saveRDS(model_objects_to_save,
          file.path(output_dir, "12_FinalModelObjects.rds"))
  write_csv(figure3_permutation_raw,
            file.path(output_dir, "13_PermutationImportance_Raw.csv"))
  write_csv(figure3_permutation_by_background,
            file.path(output_dir, "14_PermutationImportance_ByBackground_Corrected.csv"))
  write_csv(figure3_permutation_consensus,
            file.path(output_dir, "15_PermutationImportance_Consensus_Corrected.csv"))

  cat(
    "\n============================================================\n"
  )

  cat(
    analysis_name,
    " GLM ANALYSIS COMPLETE\n"
  )

  cat(
    paste0(
      "Optimal linear buffer: ", optimal_buffer_km, " km\n",
      "Optimal linear spatial block: ", optimal_block_m / 1000, " km\n",
      "Best background at optimal setting: ", best_background_name, "\n"
    )
  )

  cat(
    "Compared buffers: 1, 2 and 5 km; spatial blocks: 10, 20, 30, 40 and 50 km.\n"
  )

  cat(
    "============================================================\n"
  )

  writeLines(capture.output(sessionInfo()),
             file.path(output_dir, "16_SessionInfo.txt"))

  print(optimal_setting)
  print(background_comparison_optimal)
  print(variable_consensus %>% filter(Stability %in% c("Very stable", "Stable")))

  invisible(
    list(
      all_results = selection_results,
      scale_block_summary = scale_block_summary,
      optimal_setting = optimal_setting,
      background_comparison_optimal = background_comparison_optimal,
      best_background_optimal = best_background_optimal,
      variable_consensus = variable_consensus,
      permutation_importance_raw = figure3_permutation_raw,
      permutation_importance_by_background = figure3_permutation_by_background,
      permutation_importance_consensus = figure3_permutation_consensus
    )
  )
}

# ============================================================
# 17. RUN HLC ONLY
# ============================================================

hlc_data_dir <- file.path(
  base_dir,
  "HLC Data"
)

hlc_output_dir <- file.path(
  hlc_data_dir,
  "HLC_GLM_Final_Outputs"
)

hlc_results <- run_species_analysis(
  analysis_name = "HLC",
  data_dir = hlc_data_dir,
  output_dir = hlc_output_dir
)

  # END — HLC
