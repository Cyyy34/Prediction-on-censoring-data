# =============================================================================
# CASE STUDY: Censored Data Imputation Comparison (FIXED v2)
# Methods: MICE Post-hoc, MICE Truncated Normal, Tobit
# Fix: LOD bounds now use known assay detection limits (hardcoded lookup table)
#      instead of inferring from data min/max
# =============================================================================

library(mice)
library(truncnorm)
library(AER)
library(dplyr)
library(readxl)
library(glmnet)
library(pROC)
library(caret)

# =============================================================================
# 0. PATHS & SETTINGS
# =============================================================================
data_path  <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Case_Study_Results/cleaned_dataset_with_indicators.xlsx"
output_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Case_Study_Results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

TARGET      <- "TYPE_CA_A09"
N_BOOT      <- 100
SEL_THRESH  <- 0.50
N_FOLDS     <- 5
M_IMP       <- 5
MAXIT_MICE  <- 5
MAXIT_TOBIT <- 5
set.seed(123)

# =============================================================================
# 1. VARIABLE DEFINITIONS
# =============================================================================
BIOMARKER_VARS <- c(
  "ADAM12_V1", "ENDOGLIN_V1", "SFLT1_V1", "VEGF_V1", "AFP_V1",
  "fbHCG_V1", "INHIBINA_V1", "PAPPA_V1", "PLGF_V1",
  "ADAM12_V2", "ENDOGLIN_V2", "SFLT1_V2", "VEGF_V2", "AFP_V2",
  "fbHCG_V2", "INHIBINA_V2", "PAPPA_V2", "PLGF_V2"
)

INDICATOR_VARS <- c(
  "ADAM12c_V1", "ADAM12c_V2", "ENDOGLINc_V1", "ENDOGLINc_V2",
  "SFLT1c_V1", "SFLT1c_V2", "VEGFc_V1", "VEGFc_V2",
  "AFPc_V1", "AFPc_V2", "fbHCGc_V1", "fbHCGc_V2",
  "INHIBINAc_V1", "INHIBINAc_V2", "PAPPAc_V1", "PAPPAc_V2",
  "PLGFc_V1", "PLGFc_V2"
)

CATEGORICAL_COLS <- c(
  "CRace", "BMI_Cat", "PreGestDM", "S02C01", "U3DA03", "AgeCat_V1",
  "V1AD12a", "V1AD12b", "V1AD12c", "V1AD12e", "V1AE01", "V1AF04",
  "V2AE06c", "V2AE06d", "V2AE06e", "V2AE09", "poverty",
  "Education", "V1AD08", "V1AD08a", "V1AD08b", "V1AG04", "V1AG05",
  "V1AG06", "V1AG07", "V1AG10", "V1AH03", "VXXB01ae_V1a", "V2AD03",
  "V2AD03a", "V2AH02", "V2AH05", "V2IA07", "VXXB01ae_V2a", "VXXB01aa_V2a",
  "V3AD03", "V3AD03a", "V3AD07d1", "V3AF02", "V3AF05", "V3AG01", "V3AG02",
  "V3AG03", "V3AG04", "V3AG05", "V3AG06", "V3AG07", "V3AG08", "VXXB01ae_V3a",
  "VXXB01aa_V3a", "V3AG09", "V3AG10"
)

NUMERIC_COLS <- c(
  "BMI", "Age_at_V1", "PctFedPoverty",
  "STRESSTABSAMOUNT", "STRESSTABSYEARS", "VITD_MCG", "SUP_VITD",
  "VITAMINDAMOUNT", "VITAMINDYEARS", "U1CC07", "U1CD07",
  "Height", "VXXB01aa_V1a", "U2BB02", "U2BC03b", "U2CC07", "U2CD07",
  "V2BA01_KG", "VXXB01aa_V2a", "BMI_V2", "U3BB02", "U3BC03b", "U3CC07",
  "U3CD07", "U3DB02_weeks", "U3DB02_days", "U3DC01",
  "U3DD02a1", "U3DD02a2", "U3DD02a3", "U3DD02b1", "U3DD02b2", "U3DD02b3",
  "U3DD03a1", "U3DD03a2", "U3DD03a3", "U3DD03b1", "U3DD03b2", "U3DD03b3",
  "U3DD04a1", "U3DD04a2", "U3DD04a3", "U3DD04b1", "U3DD04b2", "U3DD04b3",
  "U3DF02a1", "U3DF02a2", "U3DF02a3", "U3DF02b1", "U3DF02b2", "U3DF02b3",
  "U3DF03a1", "U3DF03a2", "U3DF03a3", "U3DF03b1", "U3DF03b2", "U3DF03b3",
  "U3DF04a1", "U3DF04a2", "U3DF04a3", "U3DF04b1", "U3DF04b2", "U3DF04b3",
  "V3BA01_KG", "VXXB01aa_V3a", "BMI_V3",
  "ADAM12_V1", "ENDOGLIN_V1", "SFLT1_V1", "VEGF_V1", "AFP_V1",
  "fbHCG_V1", "INHIBINA_V1", "PAPPA_V1", "PLGF_V1",
  "ADAM12_V2", "ENDOGLIN_V2", "SFLT1_V2", "VEGF_V2", "AFP_V2",
  "fbHCG_V2", "INHIBINA_V2", "PAPPA_V2", "PLGF_V2"
)

# =============================================================================
# 2. KNOWN ASSAY DETECTION LIMITS (replaces infer_lod)
# Source: columns_to_process_up / columns_to_process_lower
# lower = -Inf means no left censoring for that variable
# upper = Inf  means no right censoring for that variable
# =============================================================================
build_lod_info <- function(indicator_vars) {
  # Each entry: biomarker -> list(lower, upper, ind_col)
  # -Inf / Inf = no censoring on that side
  lod_table <- list(
    ADAM12_V1   = list(lower = 1e-9,  upper = 20,     ind_col = "ADAM12c_V1"),
    ADAM12_V2   = list(lower = 1e-9,  upper = 20,     ind_col = "ADAM12c_V2"),
    ENDOGLIN_V1 = list(lower = 1e-9,  upper = 10,     ind_col = "ENDOGLINc_V1"),
    ENDOGLIN_V2 = list(lower = 1e-9,  upper = 10,     ind_col = "ENDOGLINc_V2"),
    SFLT1_V1    = list(lower = -Inf,  upper = Inf,    ind_col = "SFLT1c_V1"),   # commented out in original
    SFLT1_V2    = list(lower = -Inf,  upper = Inf,    ind_col = "SFLT1c_V2"),   # commented out in original
    VEGF_V1     = list(lower = 0.48,  upper = Inf,    ind_col = "VEGFc_V1"),    # lower only
    VEGF_V2     = list(lower = 0.48,  upper = Inf,    ind_col = "VEGFc_V2"),    # lower only
    AFP_V1      = list(lower = 0.1,   upper = 225,    ind_col = "AFPc_V1"),
    AFP_V2      = list(lower = 0.1,   upper = 225,    ind_col = "AFPc_V2"),
    fbHCG_V1    = list(lower = -Inf,  upper = Inf,    ind_col = "fbHCGc_V1"),   # commented out in original
    fbHCG_V2    = list(lower = 0.24,  upper = Inf,    ind_col = "fbHCGc_V2"),   # lower only
    INHIBINA_V1 = list(lower = 6.25,  upper = 1492,   ind_col = "INHIBINAc_V1"),
    INHIBINA_V2 = list(lower = 6.25,  upper = 1492,   ind_col = "INHIBINAc_V2"),
    PAPPA_V1    = list(lower = 1,     upper = 2300,   ind_col = "PAPPAc_V1"),
    PAPPA_V2    = list(lower = -Inf,  upper = Inf,    ind_col = "PAPPAc_V2"),   # commented out in original
    PLGF_V1     = list(lower = 0.5,   upper = 1341.3, ind_col = "PLGFc_V1"),
    PLGF_V2     = list(lower = 0.5,   upper = Inf,    ind_col = "PLGFc_V2")     # lower only in original
  )
  
  # Sanity check: warn if ind_col not in indicator_vars
  for (bv in names(lod_table)) {
    ic <- lod_table[[bv]]$ind_col
    if (!is.na(ic) && !ic %in% indicator_vars)
      warning("Indicator column ", ic, " for ", bv, " not found in INDICATOR_VARS")
  }
  
  # Print summary
  cat("\nLOD lookup table:\n")
  cat(sprintf("  %-20s  lower=%10s  upper=%10s\n", "Biomarker", "lower", "upper"))
  for (bv in names(lod_table)) {
    info <- lod_table[[bv]]
    cat(sprintf("  %-20s  lower=%10s  upper=%10s\n",
                bv,
                ifelse(is.infinite(info$lower), "-Inf", info$lower),
                ifelse(is.infinite(info$upper),  "Inf", info$upper)))
  }
  lod_table
}

# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================
clean_column_names <- function(nms) {
  nms <- gsub("[^[:alnum:]_]", "_", nms)
  nms <- gsub("^([0-9])", "X\\1", nms)
  nms
}

# Build _ind columns (Below_Lower / Above_Upper / NA) from raw c_ columns
build_indicator_df <- function(df, lod_info) {
  df_ind <- df
  for (bv in names(lod_info)) {
    info    <- lod_info[[bv]]
    ind_col <- info$ind_col
    new_ind <- paste0(bv, "_ind")
    if (is.na(ind_col) || !ind_col %in% colnames(df)) {
      df_ind[[new_ind]] <- NA_character_
      next
    }
    raw_ind <- df[[ind_col]]
    status  <- rep(NA_character_, nrow(df))
    # Only flag as censored if that side actually has a finite limit
    if (is.finite(info$lower))
      status[!is.na(raw_ind) & raw_ind == "<"] <- "Below_Lower"
    if (is.finite(info$upper))
      status[!is.na(raw_ind) & raw_ind == ">"] <- "Above_Upper"
    df_ind[[new_ind]] <- status
  }
  df_ind
}

# Set censored biomarker values to NA
set_censored_to_na <- function(df, lod_info) {
  for (bv in names(lod_info)) {
    info    <- lod_info[[bv]]
    ind_col <- info$ind_col
    if (is.na(ind_col) || !ind_col %in% colnames(df) || !bv %in% colnames(df)) next
    raw_ind  <- df[[ind_col]]
    censored <- (!is.na(raw_ind)) &
      ((raw_ind == "<" & is.finite(info$lower)) |
         (raw_ind == ">" & is.finite(info$upper)))
    df[[bv]][censored] <- NA_real_
    df[[bv]] <- suppressWarnings(as.numeric(df[[bv]]))
  }
  df
}

prepare_feature_matrix <- function(df_imp, target_var, cat_cols, num_cols,
                                   fill_na = FALSE) {
  drop_cols <- c(target_var,
                 grep("c_V[12]$|_ind$", names(df_imp), value = TRUE))
  X <- df_imp[, !names(df_imp) %in% drop_cols, drop = FALSE]
  
  for (col in intersect(cat_cols, colnames(X)))
    X[[col]] <- as.integer(as.factor(X[[col]]))
  for (col in intersect(num_cols, colnames(X)))
    X[[col]] <- suppressWarnings(as.numeric(X[[col]]))
  for (col in colnames(X))
    if (is.character(X[[col]]))
      X[[col]] <- suppressWarnings(as.numeric(X[[col]]))
  
  if (fill_na) {
    for (col in colnames(X))
      if (is.numeric(X[[col]]) && any(is.na(X[[col]])))
        X[[col]][is.na(X[[col]])] <- mean(X[[col]], na.rm = TRUE)
  }
  
  single_lv <- sapply(X, function(v) length(unique(v[!is.na(v)])) <= 1)
  if (any(single_lv)) X <- X[, !single_lv, drop = FALSE]
  colnames(X) <- clean_column_names(colnames(X))
  X
}

# =============================================================================
# 4. DATA LOADING & PREPROCESSING
# =============================================================================
cat("Loading data...\n")
df_raw <- read_excel(data_path)
cat("Raw dimensions:", nrow(df_raw), "x", ncol(df_raw), "\n")

df_raw <- df_raw[, !colnames(df_raw) %in% c("check column", "PublicID")]
df_raw <- df_raw[!is.na(df_raw[[TARGET]]), ]
cat("After removing missing outcome:", nrow(df_raw), "rows\n")

keep_vars <- unique(c(CATEGORICAL_COLS, NUMERIC_COLS, INDICATOR_VARS, TARGET))
keep_vars <- intersect(keep_vars, colnames(df_raw))
df_raw    <- df_raw[, keep_vars]

non_biomarker_vars <- setdiff(colnames(df_raw), c(BIOMARKER_VARS, INDICATOR_VARS, TARGET))
miss_rate <- colMeans(is.na(df_raw[, non_biomarker_vars, drop = FALSE]))
drop_vars <- names(miss_rate[miss_rate > 0.30])
if (length(drop_vars) > 0) {
  cat("Dropping", length(drop_vars), "non-biomarker variables with >30% missing:\n")
  cat(paste(drop_vars, collapse = ", "), "\n")
  df_raw <- df_raw[, !colnames(df_raw) %in% drop_vars]
}
CATEGORICAL_COLS <- intersect(CATEGORICAL_COLS, colnames(df_raw))
NUMERIC_COLS     <- intersect(NUMERIC_COLS,     colnames(df_raw))
cat("Working dimensions:", nrow(df_raw), "x", ncol(df_raw), "\n")

# Build LOD info from known assay limits
lod_info <- build_lod_info(INDICATOR_VARS)

# Set censored values to NA
df_raw <- set_censored_to_na(df_raw, lod_info)

# Build _ind columns
df_work <- build_indicator_df(df_raw, lod_info)

# Sanity check
cat("\nSanity check - censored counts per biomarker:\n")
for (bv in names(lod_info)) {
  if (!bv %in% colnames(df_work)) next
  ind_col <- paste0(bv, "_ind")
  n_na    <- sum(is.na(df_work[[bv]]))
  n_below <- if (ind_col %in% names(df_work)) sum(df_work[[ind_col]] == "Below_Lower", na.rm=TRUE) else 0
  n_above <- if (ind_col %in% names(df_work)) sum(df_work[[ind_col]] == "Above_Upper", na.rm=TRUE) else 0
  if (n_na > 0 || n_below > 0 || n_above > 0)
    cat(sprintf("  %-20s NA=%d  Below=%d  Above=%d\n", bv, n_na, n_below, n_above))
}

y_all <- as.numeric(as.factor(df_work[[TARGET]])) - 1
cat("\nEvent rate:", round(mean(y_all), 4), "\n")

# =============================================================================
# 5. IMPUTATION FUNCTIONS
# =============================================================================
CURRENT_VAR_NAME <- NULL
global_df_ref    <- NULL
global_lod_info  <- NULL

get_bounds_case <- function(var_name, n_rows, data_df) {
  info <- global_lod_info[[var_name]]
  if (is.null(info)) return(list(lower = rep(-Inf, n_rows), upper = rep(Inf, n_rows)))
  ind_col   <- paste0(var_name, "_ind")
  lower_vec <- rep(-Inf, n_rows)
  upper_vec <- rep( Inf, n_rows)
  if (ind_col %in% names(data_df)) {
    indicators <- data_df[[ind_col]]
    if (is.finite(info$upper))
      lower_vec[!is.na(indicators) & indicators == "Above_Upper"] <- info$upper
    if (is.finite(info$lower))
      upper_vec[!is.na(indicators) & indicators == "Below_Lower"] <- info$lower
  }
  list(lower = lower_vec, upper = upper_vec)
}

mice.impute.posthoc_case <- function(y, ry, x, wy = NULL, ...) {
  imputed <- mice.impute.norm(y, ry, x, wy = wy, ...)
  bounds  <- get_bounds_case(CURRENT_VAR_NAME, length(y), global_df_ref)
  if (is.null(wy)) wy <- !ry
  imputed <- pmax(imputed, bounds$lower[wy])
  imputed <- pmin(imputed, bounds$upper[wy])
  imputed
}

mice.impute.truncnorm_case <- function(y, ry, x, wy = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  x      <- cbind(1, as.matrix(x))
  parm   <- mice:::.norm.draw(y, ry, x)
  x_miss <- x[wy, , drop = FALSE]
  mu     <- x_miss %*% parm$beta
  sigma  <- parm$sigma
  bounds <- get_bounds_case(CURRENT_VAR_NAME, length(y), global_df_ref)
  rtruncnorm(sum(wy),
             a    = bounds$lower[wy],
             b    = bounds$upper[wy],
             mean = mu,
             sd   = sigma)
}

make_mice_imputer_case <- function(method_type) {
  function(y, ry, x, wy = NULL, ...) {
    var_name <- NA
    for (col in names(global_df_ref)) {
      if (is.numeric(global_df_ref[[col]]) &&
          length(y) == nrow(global_df_ref) &&
          all(abs(y[ry] - global_df_ref[[col]][ry]) < 1e-10, na.rm = TRUE)) {
        var_name <- col; break
      }
    }
    assign("CURRENT_VAR_NAME", var_name, envir = .GlobalEnv)
    if (method_type == "posthoc") {
      mice.impute.posthoc_case(y, ry, x, wy, ...)
    } else {
      mice.impute.truncnorm_case(y, ry, x, wy, ...)
    }
  }
}

pool_mice_to_single <- function(mice_obj, orig_df, target_var) {
  long_dat <- complete(mice_obj, action = "long", include = FALSE)
  pooled   <- long_dat %>%
    group_by(.id) %>%
    summarise(across(everything(), function(x) {
      if (is.numeric(x)) mean(x) else dplyr::first(x)
    })) %>%
    select(-.id, -.imp) %>%
    as.data.frame()
  pooled[[target_var]] <- orig_df[[target_var]]
  pooled
}

run_mice_imputation <- function(df_input, lod_info_list, target_var,
                                method_type = "posthoc",
                                m = M_IMP, maxit = MAXIT_MICE) {
  global_df_ref   <<- df_input
  global_lod_info <<- lod_info_list
  
  # Variables that actually have finite censoring limits
  censored_vars <- names(lod_info_list)[sapply(lod_info_list, function(x)
    is.finite(x$lower) || is.finite(x$upper))]
  
  ind_cols_internal <- grep("_ind$", names(df_input), value = TRUE)
  raw_c_cols        <- grep("c_V[12]$", names(df_input), value = TRUE)
  df_mice           <- df_input[, !names(df_input) %in% raw_c_cols]
  
  ini  <- mice(df_mice, maxit = 0, print = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  
  exclude_pred <- c(ind_cols_internal, target_var)
  pred[, intersect(exclude_pred, colnames(pred))] <- 0
  meth[c(ind_cols_internal, target_var)] <- ""
  
  wrapper_name <- paste0(method_type, "_case_wrapper")
  assign(paste0("mice.impute.", wrapper_name),
         make_mice_imputer_case(method_type),
         envir = .GlobalEnv)
  meth[intersect(censored_vars, names(meth))] <- wrapper_name
  
  imp <- mice(df_mice, method = meth, predictorMatrix = pred,
              m = m, maxit = maxit, printFlag = FALSE)
  pool_mice_to_single(imp, df_mice, target_var)
}

run_tobit_imputation <- function(df_input, lod_info_list, target_var,
                                 m = M_IMP, max_iter = MAXIT_TOBIT) {
  censored_vars <- names(lod_info_list)[sapply(lod_info_list, function(x)
    is.finite(x$lower) || is.finite(x$upper))]
  
  ind_cols_internal <- grep("_ind$", names(df_input), value = TRUE)
  raw_c_cols        <- grep("c_V[12]$", names(df_input), value = TRUE)
  df_base           <- df_input[, !names(df_input) %in% raw_c_cols]
  
  impute_once <- function(seed_val) {
    set.seed(seed_val)
    df_f <- df_base
    for (v in names(df_f))
      if (is.numeric(df_f[[v]]) && any(is.na(df_f[[v]])))
        df_f[[v]][is.na(df_f[[v]])] <- mean(df_f[[v]], na.rm = TRUE)
    
    for (iter in seq_len(max_iter)) {
      for (bv in censored_vars) {
        info    <- lod_info_list[[bv]]
        ind_col <- paste0(bv, "_ind")
        if (!bv %in% names(df_f) || !ind_col %in% names(df_f)) next
        
        mask_lower <- !is.na(df_f[[ind_col]]) & df_f[[ind_col]] == "Below_Lower"
        mask_upper <- !is.na(df_f[[ind_col]]) & df_f[[ind_col]] == "Above_Upper"
        if (!any(mask_lower) && !any(mask_upper)) next
        
        y_temp             <- df_f[[bv]]
        if (any(mask_lower)) y_temp[mask_lower] <- info$lower
        if (any(mask_upper)) y_temp[mask_upper] <- info$upper
        
        exclude_cols <- c(bv, ind_cols_internal, target_var)
        predictors   <- setdiff(names(df_f), exclude_cols)
        predictors   <- predictors[sapply(df_f[predictors], is.numeric)]
        if (length(predictors) == 0) next
        
        data_fit        <- df_f
        data_fit[[bv]]  <- y_temp
        
        # Use -Inf/Inf for sides with no censoring
        tobit_left  <- if (is.finite(info$lower)) info$lower else -Inf
        tobit_right <- if (is.finite(info$upper)) info$upper else  Inf
        
        tryCatch({
          fit       <- tobit(as.formula(paste(bv, "~", paste(predictors, collapse = "+"))),
                             left = tobit_left, right = tobit_right,
                             data = data_fit)
          sigma_hat <- fit$scale
          mu_hat    <- predict(fit, newdata = data_fit)
          
          if (any(mask_lower))
            df_f[[bv]][mask_lower] <- rtruncnorm(
              sum(mask_lower), a = -Inf, b = info$lower,
              mean = mu_hat[mask_lower], sd = sigma_hat)
          if (any(mask_upper))
            df_f[[bv]][mask_upper] <- rtruncnorm(
              sum(mask_upper), a = info$upper, b = Inf,
              mean = mu_hat[mask_upper], sd = sigma_hat)
        }, error = function(e)
          message("Tobit failed for ", bv, " iter ", iter, ": ", e$message))
      }
    }
    df_f
  }
  
  imp_list  <- lapply(seq_len(m), function(k) impute_once(1000 + k))
  df_pooled <- imp_list[[1]]
  avg_cols  <- setdiff(names(df_pooled)[sapply(df_pooled, is.numeric)], target_var)
  for (col in avg_cols) {
    mat              <- sapply(imp_list, function(d) d[[col]])
    df_pooled[[col]] <- rowMeans(mat)
  }
  df_pooled[[target_var]] <- df_base[[target_var]]
  df_pooled
}

# =============================================================================
# 6. RUN IMPUTATION & SAVE
# =============================================================================
cat("\n=== Running MICE Post-hoc imputation ===\n")
df_mice_posthoc <- run_mice_imputation(df_work, lod_info, TARGET, method_type = "posthoc")
cat("NA in biomarkers after MICE_Posthoc:\n")
for (v in names(lod_info)) if (v %in% names(df_mice_posthoc))
  cat(sprintf("  %-20s NA=%d\n", v, sum(is.na(df_mice_posthoc[[v]]))))
write.csv(df_mice_posthoc, file.path(output_dir, "imputed_mice_posthoc.csv"), row.names = FALSE)

cat("\n=== Running MICE Truncated Normal imputation ===\n")
df_mice_trunc <- run_mice_imputation(df_work, lod_info, TARGET, method_type = "truncnorm")
cat("NA in biomarkers after MICE_Trunc:\n")
for (v in names(lod_info)) if (v %in% names(df_mice_trunc))
  cat(sprintf("  %-20s NA=%d\n", v, sum(is.na(df_mice_trunc[[v]]))))
write.csv(df_mice_trunc, file.path(output_dir, "imputed_mice_trunc.csv"), row.names = FALSE)

cat("\n=== Running Tobit imputation ===\n")
df_tobit <- run_tobit_imputation(df_work, lod_info, TARGET)
cat("NA in biomarkers after Tobit:\n")
for (v in names(lod_info)) if (v %in% names(df_tobit))
  cat(sprintf("  %-20s NA=%d\n", v, sum(is.na(df_tobit[[v]]))))
write.csv(df_tobit, file.path(output_dir, "imputed_tobit.csv"), row.names = FALSE)

cat("\nAll imputed datasets saved to:", output_dir, "\n")
cat("Next step: run case_study_pipeline.R for each method.\n")




########################################visual
