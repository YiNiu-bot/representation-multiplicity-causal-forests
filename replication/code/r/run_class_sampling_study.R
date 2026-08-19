#!/usr/bin/env Rscript

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_flag)) {
  normalizePath(sub("^--file=", "", script_flag[[1L]]))
} else {
  normalizePath("replication/code/r/run_class_sampling_study.R")
}
repository_dir <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
replication_dir <- file.path(repository_dir, "replication")
run_scope <- match.arg(
  tolower(Sys.getenv("RI_RUN_SCOPE", "all")),
  c("all", "simulation")
)
rlib <- Sys.getenv("RI_RLIB", "")
if (nzchar(rlib) && dir.exists(rlib)) .libPaths(c(rlib, .libPaths()))

suppressPackageStartupMessages({
  library(grf)
  library(haven)
})
expected_grf_version <- "2.4.0"
allow_unpinned_grf <- Sys.getenv("RI_ALLOW_UNPINNED_GRF", "0") == "1"
if (as.character(packageVersion("grf")) != expected_grf_version &&
    !allow_unpinned_grf) {
  stop(
    sprintf(
      "Expected grf %s, found %s. Set RI_ALLOW_UNPINNED_GRF=1 only for a non-replication run.",
      expected_grf_version,
      packageVersion("grf")
    ),
    call. = FALSE
  )
}
source(file.path(replication_dir, "code", "r", "class_sampled_forest.R"))

output_dir <- Sys.getenv(
  "RI_OUTPUT_DIR",
  file.path(repository_dir, "rerun_outputs", "class_sampling")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source_dir <- Sys.getenv(
  "RI_SOURCE_DIR",
  file.path(replication_dir, "data", "kenya")
)

physical_cores <- parallel::detectCores(logical = FALSE)
if (is.na(physical_cores)) physical_cores <- 28L
threads <- as.integer(Sys.getenv("RI_THREADS", as.character(physical_cores - 1L)))
real_trees <- as.integer(Sys.getenv("RI_REAL_TREES", "2000"))
real_seed_count <- as.integer(Sys.getenv("RI_REAL_SEEDS", "3"))
mc_trees <- as.integer(Sys.getenv("RI_MC_TREES", "800"))
mc_reps <- as.integer(Sys.getenv("RI_MC_REPS", "20"))
mc_n_train <- as.integer(Sys.getenv("RI_MC_N_TRAIN", "3000"))
mc_n_test <- as.integer(Sys.getenv("RI_MC_N_TEST", "30000"))
base_seed <- 20260808L

required_data <- c(
  "GE_HH-Analysis_AllHHs.dta",
  "GE_HH-Survey-BL_Analysis_AllHHs.dta"
)
missing_data <- required_data[!file.exists(file.path(source_dir, required_data))]
if (run_scope == "all" && length(missing_data)) {
  stop(
    sprintf(
      "Missing required public data under %s: %s",
      source_dir,
      paste(missing_data, collapse = ", ")
    ),
    call. = FALSE
  )
}
expected_data_md5 <- c(
  "GE_HH-Analysis_AllHHs.dta" = "6841329ee31e3fee6413651bf8a1c5ff",
  "GE_HH-Survey-BL_Analysis_AllHHs.dta" = "35b0b1384ea3cf82e38bd77c73442f3c"
)
allow_data_mismatch <- Sys.getenv("RI_ALLOW_DATA_MISMATCH", "0") == "1"
if (run_scope == "all") {
  observed_data_md5 <- unname(tools::md5sum(file.path(source_dir, required_data)))
  names(observed_data_md5) <- required_data
  bad_md5 <- names(expected_data_md5)[
    tolower(observed_data_md5[names(expected_data_md5)]) != expected_data_md5
  ]
  if (length(bad_md5) && !allow_data_mismatch) {
    stop(
      sprintf(
        "Input checksum mismatch for: %s. Set RI_ALLOW_DATA_MISMATCH=1 only for a non-replication run.",
        paste(bad_md5, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

weighted_mean <- function(x, w) sum(w * x) / sum(w)

cluster_mean <- function(x, w, cluster) {
  estimate <- weighted_mean(x, w)
  influence <- w * (x - estimate) / sum(w)
  cluster_influence <- rowsum(influence, cluster, reorder = FALSE)[, 1]
  g <- length(cluster_influence)
  list(
    estimate = estimate,
    se = sqrt(g / (g - 1) * sum(cluster_influence^2))
  )
}

top_share <- function(score, share, id = seq_along(score)) {
  k <- floor(length(score) * share)
  selected <- rep(FALSE, length(score))
  if (k > 0L) selected[order(-score, id)[seq_len(k)]] <- TRUE
  selected
}

comparison_metrics <- function(tau, selected, reference_tau, reference_selected) {
  intersection <- sum(selected & reference_selected)
  union <- sum(selected | reference_selected)
  c(
    maximum_absolute_prediction_difference = max(abs(tau - reference_tau)),
    cate_correlation = cor(tau, reference_tau),
    cate_rank_correlation = cor(tau, reference_tau, method = "spearman"),
    top_set_overlap = if (union == 0L) 1 else intersection / union,
    assignment_switching = mean(selected != reference_selected)
  )
}

median_impute <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

zero_impute <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x
}

append_feature <- function(design, name, train, target = NULL) {
  design$train <- cbind(design$train, setNames(data.frame(train), name))
  if (!is.null(design$target)) {
    stopifnot(!is.null(target))
    design$target <- cbind(design$target, setNames(data.frame(target), name))
  }
  design$train <- as.matrix(design$train)
  if (!is.null(design$target)) design$target <- as.matrix(design$target)
  design
}

pooled_values <- function(X.train, X.target, feature) {
  if (is.null(X.target)) X.train[, feature] else c(X.train[, feature], X.target[, feature])
}

split_pooled <- function(values, n.train, has.target) {
  list(
    train = values[seq_len(n.train)],
    target = if (has.target) values[-seq_len(n.train)] else NULL
  )
}

add_exact_aliases <- function(design, features) {
  for (feature in features) {
    design <- append_feature(
      design,
      paste0(feature, "__admin_alias"),
      design$train[, feature],
      if (is.null(design$target)) NULL else design$target[, feature]
    )
  }
  design
}

add_unit_recodings <- function(design, features) {
  n.train <- nrow(design$train)
  has.target <- !is.null(design$target)
  for (feature in features) {
    all_x <- pooled_values(design$train, design$target, feature)
    scale_x <- sd(all_x)
    if (!is.finite(scale_x) || scale_x == 0) next
    standardized <- (all_x - mean(all_x)) / scale_x
    split_x <- split_pooled(standardized, n.train, has.target)
    design <- append_feature(
      design,
      paste0(feature, "__standardized"),
      split_x$train,
      split_x$target
    )
  }
  design
}

add_complements <- function(design, features) {
  n.train <- nrow(design$train)
  has.target <- !is.null(design$target)
  for (feature in features) {
    all_x <- pooled_values(design$train, design$target, feature)
    complemented <- min(all_x) + max(all_x) - all_x
    split_x <- split_pooled(complemented, n.train, has.target)
    design <- append_feature(
      design,
      paste0(feature, "__complement"),
      split_x$train,
      split_x$target
    )
  }
  design
}

add_monotone_bundle <- function(design, features) {
  n.train <- nrow(design$train)
  has.target <- !is.null(design$target)
  for (feature in features) {
    all_x <- pooled_values(design$train, design$target, feature)
    shifted <- all_x - min(all_x)
    logged <- log1p(shifted)
    ranked <- rank(all_x, ties.method = "average")
    log_split <- split_pooled(logged, n.train, has.target)
    rank_split <- split_pooled(ranked, n.train, has.target)
    design <- append_feature(
      design,
      paste0(feature, "__log"),
      log_split$train,
      log_split$target
    )
    design <- append_feature(
      design,
      paste0(feature, "__rank"),
      rank_split$train,
      rank_split$target
    )
  }
  design
}

make_representation <- function(
    X.train,
    X.target = NULL,
    scenario,
    aliases,
    units,
    complements,
    monotone
) {
  design <- list(train = as.matrix(X.train), target = if (is.null(X.target)) NULL else as.matrix(X.target))
  if (scenario %in% c("administrative_aliases", "routine_bundle")) {
    design <- add_exact_aliases(design, aliases)
  }
  if (scenario %in% c("unit_recodings", "routine_bundle")) {
    design <- add_unit_recodings(design, units)
  }
  if (scenario %in% c("complements", "routine_bundle")) {
    design <- add_complements(design, complements)
  }
  if (scenario %in% c("raw_log_rank", "routine_bundle")) {
    design <- add_monotone_bundle(design, monotone)
  }
  design
}

fit_forest <- function(
    X.train,
    Y,
    W,
    X.target,
    mtry,
    seed,
    num.trees,
    Y.hat,
    W.hat,
    sample.weights = NULL,
    clusters = NULL
) {
  mtry <- as.numeric(mtry)
  if (length(mtry) != 1L || !is.finite(mtry) || mtry <= 0) {
    stop("mtry must be one positive finite number.", call. = FALSE)
  }
  arguments <- list(
    X = X.train,
    Y = Y,
    W = W,
    Y.hat = Y.hat,
    W.hat = W.hat,
    mtry = mtry,
    num.trees = as.integer(num.trees),
    seed = as.integer(seed),
    num.threads = threads
  )
  if (!is.null(sample.weights)) arguments$sample.weights <- sample.weights
  if (!is.null(clusters)) arguments$clusters <- clusters
  fit <- do.call(grf::causal_forest, arguments)
  predictions <- if (is.null(X.target)) {
    as.numeric(predict(fit)$predictions)
  } else {
    as.numeric(predict(fit, X.target)$predictions)
  }
  list(
    fit = fit,
    predictions = predictions,
    input_dimension = ncol(X.train),
    raw_dimension = ncol(X.train),
    fitted_dimension = ncol(X.train),
    semantic_dimension = NA_integer_,
    mtry = mtry,
    class_map = NULL
  )
}

fit_method <- function(
    design,
    method,
    Y,
    W,
    base_dimension,
    base_mtry,
    seed,
    num.trees,
    Y.hat,
    W.hat,
    sample.weights = NULL,
    clusters = NULL
) {
  if (method == "ordinary") {
    return(fit_forest(
      design$train, Y, W, design$target, base_mtry, seed, num.trees,
      Y.hat, W.hat, sample.weights, clusters
    ))
  }

  if (method == "exact_deduplication") {
    reduced <- drop_exact_columns(design$train, design$target)
    result <- fit_forest(
      reduced$X.train, Y, W, reduced$X.target, base_mtry, seed, num.trees,
      Y.hat, W.hat, sample.weights, clusters
    )
    result$input_dimension <- ncol(design$train)
    return(result)
  }

  if (method == "dimension_adjusted_mtry") {
    adjusted_mtry <- ceiling(base_mtry * ncol(design$train) / base_dimension)
    return(fit_forest(
      design$train, Y, W, design$target, adjusted_mtry, seed, num.trees,
      Y.hat, W.hat, sample.weights, clusters
    ))
  }

  if (method == "class_sampled") {
    arguments <- list(
      X.train = design$train,
      Y = Y,
      W = W,
      X.target = design$target,
      mtry = base_mtry,
      Y.hat = Y.hat,
      W.hat = W.hat,
      num.trees = as.integer(num.trees),
      seed = as.integer(seed),
      num.threads = threads
    )
    if (!is.null(sample.weights)) arguments$sample.weights <- sample.weights
    if (!is.null(clusters)) arguments$clusters <- clusters
    result <- do.call(class_sampled_causal_forest, arguments)
    result$input_dimension <- ncol(design$train)
    return(result)
  }

  stop(sprintf("Unknown method: %s", method), call. = FALSE)
}

write_class_map <- function(class_map, sample_name, scenario, seed) {
  class_map$sample <- sample_name
  class_map$scenario <- scenario
  class_map$seed <- seed
  class_map[, c(
    "sample", "scenario", "seed", "raw_column", "raw_name", "class_id",
    "representative", "representative_name"
  )]
}

summarize_simulation <- function(results) {
  groups <- unique(results[, c("scenario", "method")])
  metrics <- c(
    "cate_rmse", "cate_correlation", "cate_rank_correlation",
    "representation_cate_correlation", "representation_rank_correlation",
    "assignment_switching", "top_set_overlap", "policy_value",
    "policy_regret", "policy_value_difference", "regret_difference"
  )
  rows <- vector("list", nrow(groups) * length(metrics))
  cursor <- 0L
  for (i in seq_len(nrow(groups))) {
    keep <- results$scenario == groups$scenario[i] & results$method == groups$method[i]
    for (metric in metrics) {
      values <- results[[metric]][keep]
      n.values <- length(values)
      standard_deviation <- if (n.values > 1L) sd(values) else NA_real_
      standard_error <- if (n.values > 1L) {
        standard_deviation / sqrt(n.values)
      } else {
        NA_real_
      }
      critical <- if (n.values > 1L) qt(0.975, df = n.values - 1L) else NA_real_
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        scenario = groups$scenario[i],
        method = groups$method[i],
        metric = metric,
        replications = n.values,
        mean = mean(values),
        sd = standard_deviation,
        se_mean = standard_error,
        ci95_low = mean(values) - critical * standard_error,
        ci95_high = mean(values) + critical * standard_error,
        q025 = unname(quantile(values, 0.025)),
        median = median(values),
        q975 = unname(quantile(values, 0.975)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

cat(sprintf(
  "Using %d threads and %d Monte Carlo replications with %d trees.\n",
  threads, mc_reps, mc_trees
))

scenarios <- c(
  "administrative_aliases",
  "unit_recodings",
  "complements",
  "raw_log_rank",
  "routine_bundle"
)
methods <- c(
  "ordinary",
  "exact_deduplication",
  "dimension_adjusted_mtry",
  "class_sampled"
)

real_results <- NULL
if (run_scope == "all") {
cat("Preparing the public Kenya cash-transfer audit.\n")
analysis <- read_dta(
  file.path(source_dir, "GE_HH-Analysis_AllHHs.dta"),
  col_select = c(
    "hhid_key", "eligible", "baselined", "treat", "village_code",
    "hhweight_EL", "p2_consumption_pc_wins_PPP"
  )
)
baseline <- read_dta(
  file.path(source_dir, "GE_HH-Survey-BL_Analysis_AllHHs.dta"),
  col_select = c(
    "hhid_key", "hhsize1_BL", "female_BL", "haschildhh_BL",
    "has_child_sch_BL", "numchild3_BL", "numchild6_BL", "widowed_BL",
    "age60_BL", "emp_BL", "selfemp_BL", "any_livestock",
    "land_ownland_BL", "own_land_acres_BL", "h1_4_radiotv_BL",
    "num_meals_yest", "num_meals_yest_protein"
  )
)

idx <- match(analysis$hhid_key, baseline$hhid_key)
d <- cbind(analysis, baseline[idx, setdiff(names(baseline), "hhid_key")])
keep <- d$eligible == 1 & d$baselined == 1 & !is.na(d$treat) &
  !is.na(d$village_code) & !is.na(d$hhweight_EL) & d$hhweight_EL > 0 &
  !is.na(d$p2_consumption_pc_wins_PPP)
d <- d[keep, ]

X_real <- cbind(
  hh_size = median_impute(d$hhsize1_BL),
  female_head = zero_impute(d$female_BL),
  has_children = zero_impute(d$haschildhh_BL),
  child_in_school = zero_impute(d$has_child_sch_BL),
  child_under_3 = as.numeric(zero_impute(d$numchild3_BL) > 0),
  child_under_6 = as.numeric(zero_impute(d$numchild6_BL) > 0),
  widow = zero_impute(d$widowed_BL),
  elder = zero_impute(d$age60_BL),
  employed = zero_impute(d$emp_BL),
  self_employed = zero_impute(d$selfemp_BL),
  livestock = zero_impute(d$any_livestock),
  owns_land = zero_impute(d$land_ownland_BL),
  owns_quarter_acre = as.numeric(zero_impute(d$own_land_acres_BL) >= 0.25),
  owns_radio_tv = as.numeric(zero_impute(d$h1_4_radiotv_BL) > 0),
  meals = median_impute(d$num_meals_yest),
  protein_meals = median_impute(d$num_meals_yest_protein)
)
Y_real <- as.numeric(d$p2_consumption_pc_wins_PPP)
W_real <- as.numeric(d$treat)
weights_real <- as.numeric(d$hhweight_EL)
clusters_real <- as.integer(factor(d$village_code))
base_dimension_real <- ncol(X_real)
base_mtry_real <- base_dimension_real
budget_share_real <- 0.5

real_aliases <- c("hh_size", "employed", "owns_land", "owns_radio_tv")
real_units <- c("hh_size", "meals", "protein_meals")
real_complements <- c("female_head", "employed", "owns_land")
real_monotone <- c("hh_size", "meals", "protein_meals")

e_hat_real <- rep(weighted_mean(W_real, weights_real), length(W_real))
y_forest_real <- regression_forest(
  X_real,
  Y_real,
  sample.weights = weights_real,
  clusters = clusters_real,
  num.trees = real_trees,
  mtry = base_mtry_real,
  seed = base_seed - 1L,
  num.threads = threads
)
y_hat_real <- as.numeric(predict(y_forest_real)$predictions)
gamma_real <- (W_real - e_hat_real) * (Y_real - y_hat_real) /
  (e_hat_real * (1 - e_hat_real))

real_rows <- list()
real_class_maps <- list()
real_seeds <- base_seed + seq_len(real_seed_count) - 1L
for (seed in real_seeds) {
  canonical_design <- make_representation(
    X_real,
    scenario = "canonical",
    aliases = real_aliases,
    units = real_units,
    complements = real_complements,
    monotone = real_monotone
  )
  references <- list()
  for (method in methods) {
    cat(sprintf("Real canonical fit: method=%s, seed=%d.\n", method, seed))
    reference <- fit_method(
      canonical_design, method, Y_real, W_real,
      base_dimension_real, base_mtry_real, seed, real_trees,
      y_hat_real, e_hat_real, weights_real, clusters_real
    )
    reference_selected <- top_share(
      reference$predictions,
      budget_share_real,
      d$hhid_key
    )
    reference_value <- cluster_mean(
      as.numeric(reference_selected) * gamma_real,
      weights_real,
      clusters_real
    )
    reference$selected <- reference_selected
    reference$value <- reference_value
    references[[method]] <- reference
    if (!is.null(reference$class_map)) {
      real_class_maps[[length(real_class_maps) + 1L]] <- write_class_map(
        reference$class_map, "real_data", "canonical", seed
      )
    }
    real_rows[[length(real_rows) + 1L]] <- data.frame(
      scenario = "canonical",
      method = method,
      seed = seed,
      n = nrow(X_real),
      clusters = length(unique(clusters_real)),
      raw_dimension = reference$input_dimension,
      fitted_dimension = reference$fitted_dimension,
      semantic_dimension = reference$semantic_dimension,
      mtry = reference$mtry,
      budget_share = budget_share_real,
      maximum_absolute_prediction_difference = 0,
      cate_correlation = 1,
      cate_rank_correlation = 1,
      top_set_overlap = 1,
      assignment_switching = 0,
      pseudo_policy_value = reference_value$estimate,
      pseudo_policy_value_se = reference_value$se,
      pseudo_policy_value_difference = 0,
      pseudo_policy_value_difference_se = 0,
      stringsAsFactors = FALSE
    )
    gc(verbose = FALSE)
  }

  for (scenario in scenarios) {
    design <- make_representation(
      X_real,
      scenario = scenario,
      aliases = real_aliases,
      units = real_units,
      complements = real_complements,
      monotone = real_monotone
    )
    for (method in methods) {
      cat(sprintf(
        "Real recoding fit: scenario=%s, method=%s, seed=%d, p=%d.\n",
        scenario, method, seed, ncol(design$train)
      ))
      fitted <- fit_method(
        design, method, Y_real, W_real,
        base_dimension_real, base_mtry_real, seed, real_trees,
        y_hat_real, e_hat_real, weights_real, clusters_real
      )
      selected <- top_share(fitted$predictions, budget_share_real, d$hhid_key)
      metrics <- comparison_metrics(
        fitted$predictions,
        selected,
        references[[method]]$predictions,
        references[[method]]$selected
      )
      value <- cluster_mean(
        as.numeric(selected) * gamma_real,
        weights_real,
        clusters_real
      )
      paired <- cluster_mean(
        as.numeric(selected - references[[method]]$selected) * gamma_real,
        weights_real,
        clusters_real
      )
      if (!is.null(fitted$class_map)) {
        real_class_maps[[length(real_class_maps) + 1L]] <- write_class_map(
          fitted$class_map, "real_data", scenario, seed
        )
      }
      real_rows[[length(real_rows) + 1L]] <- data.frame(
        scenario = scenario,
        method = method,
        seed = seed,
        n = nrow(X_real),
        clusters = length(unique(clusters_real)),
        raw_dimension = fitted$input_dimension,
        fitted_dimension = fitted$fitted_dimension,
        semantic_dimension = fitted$semantic_dimension,
        mtry = fitted$mtry,
        budget_share = budget_share_real,
        maximum_absolute_prediction_difference = metrics[[
          "maximum_absolute_prediction_difference"
        ]],
        cate_correlation = metrics[["cate_correlation"]],
        cate_rank_correlation = metrics[["cate_rank_correlation"]],
        top_set_overlap = metrics[["top_set_overlap"]],
        assignment_switching = metrics[["assignment_switching"]],
        pseudo_policy_value = value$estimate,
        pseudo_policy_value_se = value$se,
        pseudo_policy_value_difference = paired$estimate,
        pseudo_policy_value_difference_se = paired$se,
        stringsAsFactors = FALSE
      )
      gc(verbose = FALSE)
    }
  }
}

real_results <- do.call(rbind, real_rows)
write.csv(
  real_results,
  file.path(output_dir, "real_natural_recodings.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, real_class_maps),
  file.path(output_dir, "real_class_maps.csv"),
  row.names = FALSE
)
}

cat("Running the repeated-sample randomized-trial design.\n")
base_dimension_mc <- 20L
base_mtry_mc <- 6L
budget_share_mc <- 0.20
mc_aliases <- c("x01", "x08", "x20")
mc_units <- c("x02", "x05", "x06")
mc_complements <- c("x03", "x07", "x09")
mc_monotone <- c("x01", "x02", "x03")

response_functions <- function(X) {
  group_a <- X[, "x01"] > 0.70 & X[, "x03"] > 0.35
  group_b <- X[, "x02"] > 0.70 & X[, "x03"] <= 0.35
  tau <- 0.05 + 1.40 * group_a + 1.05 * group_b +
    0.25 * (X[, "x04"] - 0.5)
  mu <- 1 + X[, "x05"] + 0.5 * X[, "x06"]^2 - 0.5 * X[, "x07"]
  list(mu = mu, tau = tau, group_a = group_a, group_b = group_b)
}

mc_rows <- list()
mc_class_maps <- list()
for (replication in seq_len(mc_reps)) {
  data_seed <- base_seed + 10000L + replication
  forest_seed <- base_seed + 20000L + replication
  set.seed(data_seed)
  X_train <- matrix(
    runif(mc_n_train * base_dimension_mc),
    nrow = mc_n_train,
    ncol = base_dimension_mc
  )
  X_test <- matrix(
    runif(mc_n_test * base_dimension_mc),
    nrow = mc_n_test,
    ncol = base_dimension_mc
  )
  colnames(X_train) <- colnames(X_test) <- sprintf(
    "x%02d",
    seq_len(base_dimension_mc)
  )

  train_signal <- response_functions(X_train)
  test_signal <- response_functions(X_test)
  W_train <- rbinom(mc_n_train, 1, 0.5)
  Y_train <- train_signal$mu + W_train * train_signal$tau + rnorm(mc_n_train)
  Y_hat_train <- train_signal$mu + 0.5 * train_signal$tau
  W_hat_train <- rep(0.5, mc_n_train)

  canonical_design <- make_representation(
    X_train,
    X_test,
    "canonical",
    mc_aliases,
    mc_units,
    mc_complements,
    mc_monotone
  )
  optimal_selected <- top_share(
    test_signal$tau,
    budget_share_mc,
    seq_len(mc_n_test)
  )
  optimal_value <- mean(optimal_selected * test_signal$tau)
  canonical_classes <- certify_split_classes(
    canonical_design$train,
    canonical_design$target
  )
  references <- list()
  for (method in methods) {
    cat(sprintf(
      "Monte Carlo canonical fit %02d/%02d: method=%s.\n",
      replication, mc_reps, method
    ))
    reference <- fit_method(
      canonical_design, method, Y_train, W_train,
      base_dimension_mc, base_mtry_mc, forest_seed, mc_trees,
      Y_hat_train, W_hat_train
    )
    reference_selected <- top_share(
      reference$predictions,
      budget_share_mc,
      seq_len(mc_n_test)
    )
    reference_value <- mean(reference_selected * test_signal$tau)
    reference_regret <- optimal_value - reference_value
    reference$selected <- reference_selected
    reference$value <- reference_value
    reference$regret <- reference_regret
    references[[method]] <- reference
    if (!is.null(reference$class_map)) {
      mc_class_maps[[length(mc_class_maps) + 1L]] <- write_class_map(
        reference$class_map,
        "simulation",
        "canonical",
        replication
      )
    }
    mc_rows[[length(mc_rows) + 1L]] <- data.frame(
      replication = replication,
      data_seed = data_seed,
      forest_seed = forest_seed,
      scenario = "canonical",
      method = method,
      n_train = mc_n_train,
      n_test = mc_n_test,
      raw_dimension = reference$input_dimension,
      fitted_dimension = reference$fitted_dimension,
      semantic_dimension = canonical_classes$semantic_dimension,
      mtry = reference$mtry,
      budget_share = budget_share_mc,
      maximum_absolute_prediction_difference = 0,
      cate_rmse = sqrt(mean((reference$predictions - test_signal$tau)^2)),
      cate_correlation = cor(reference$predictions, test_signal$tau),
      cate_rank_correlation = cor(
        reference$predictions,
        test_signal$tau,
        method = "spearman"
      ),
      representation_cate_correlation = 1,
      representation_rank_correlation = 1,
      assignment_switching = 0,
      top_set_overlap = 1,
      policy_value = reference_value,
      policy_regret = reference_regret,
      policy_value_difference = 0,
      regret_difference = 0,
      share_group_a_treated = mean(test_signal$group_a[reference_selected]),
      share_group_b_treated = mean(test_signal$group_b[reference_selected]),
      stringsAsFactors = FALSE
    )
    gc(verbose = FALSE)
  }

  for (scenario in scenarios) {
    design <- make_representation(
      X_train,
      X_test,
      scenario,
      mc_aliases,
      mc_units,
      mc_complements,
      mc_monotone
    )
    declared_classes <- certify_split_classes(design$train, design$target)
    for (method in methods) {
      cat(sprintf(
        "Monte Carlo fit %02d/%02d: scenario=%s, method=%s, p=%d, G=%d.\n",
        replication, mc_reps, scenario, method, ncol(design$train),
        declared_classes$semantic_dimension
      ))
      fitted <- fit_method(
        design, method, Y_train, W_train,
        base_dimension_mc, base_mtry_mc, forest_seed, mc_trees,
        Y_hat_train, W_hat_train
      )
      selected <- top_share(
        fitted$predictions,
        budget_share_mc,
        seq_len(mc_n_test)
      )
      metrics <- comparison_metrics(
        fitted$predictions,
        selected,
        references[[method]]$predictions,
        references[[method]]$selected
      )
      value <- mean(selected * test_signal$tau)
      regret <- optimal_value - value
      if (!is.null(fitted$class_map)) {
        mc_class_maps[[length(mc_class_maps) + 1L]] <- write_class_map(
          fitted$class_map,
          "simulation",
          scenario,
          replication
        )
      }
      mc_rows[[length(mc_rows) + 1L]] <- data.frame(
        replication = replication,
        data_seed = data_seed,
        forest_seed = forest_seed,
        scenario = scenario,
        method = method,
        n_train = mc_n_train,
        n_test = mc_n_test,
        raw_dimension = fitted$input_dimension,
        fitted_dimension = fitted$fitted_dimension,
        semantic_dimension = declared_classes$semantic_dimension,
        mtry = fitted$mtry,
        budget_share = budget_share_mc,
        maximum_absolute_prediction_difference = metrics[[
          "maximum_absolute_prediction_difference"
        ]],
        cate_rmse = sqrt(mean((fitted$predictions - test_signal$tau)^2)),
        cate_correlation = cor(fitted$predictions, test_signal$tau),
        cate_rank_correlation = cor(
          fitted$predictions,
          test_signal$tau,
          method = "spearman"
        ),
        representation_cate_correlation = metrics[["cate_correlation"]],
        representation_rank_correlation = metrics[["cate_rank_correlation"]],
        assignment_switching = metrics[["assignment_switching"]],
        top_set_overlap = metrics[["top_set_overlap"]],
        policy_value = value,
        policy_regret = regret,
        policy_value_difference = value - references[[method]]$value,
        regret_difference = regret - references[[method]]$regret,
        share_group_a_treated = mean(test_signal$group_a[selected]),
        share_group_b_treated = mean(test_signal$group_b[selected]),
        stringsAsFactors = FALSE
      )
      gc(verbose = FALSE)
    }
  }
}

mc_results <- do.call(rbind, mc_rows)
mc_summary <- summarize_simulation(mc_results)
write.csv(
  mc_results,
  file.path(output_dir, "simulation_natural_recodings.csv"),
  row.names = FALSE
)
write.csv(
  mc_summary,
  file.path(output_dir, "simulation_natural_recodings_summary.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, mc_class_maps),
  file.path(output_dir, "simulation_class_maps.csv"),
  row.names = FALSE
)

invariance_check <- subset(mc_results, method == "class_sampled")
if (max(abs(invariance_check$assignment_switching)) > 0 ||
    max(abs(invariance_check$policy_value_difference)) > 0 ||
    max(invariance_check$maximum_absolute_prediction_difference) > 0 ||
    min(invariance_check$representation_cate_correlation) < 1 - 1e-14) {
  stop("Class-sampled invariance check failed.", call. = FALSE)
}

real_invariance_check <- NULL
if (run_scope == "all") {
  real_invariance_check <- subset(real_results, method == "class_sampled")
  real_invariance_check <- subset(real_invariance_check, scenario != "canonical")
  if (max(abs(real_invariance_check$assignment_switching)) > 0 ||
      max(real_invariance_check$maximum_absolute_prediction_difference) > 0 ||
      min(real_invariance_check$cate_correlation) < 1 - 1e-14) {
    stop("Real-data class-sampled invariance check failed.", call. = FALSE)
  }
}

metadata <- c(
  sprintf("completed_utc: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("run_scope: %s", run_scope),
  sprintf("grf_version: %s", packageVersion("grf")),
  sprintf("haven_version: %s", packageVersion("haven")),
  sprintf("allow_unpinned_grf: %s", allow_unpinned_grf),
  sprintf("allow_data_mismatch: %s", allow_data_mismatch),
  sprintf("threads: %d", threads),
  sprintf("real_trees: %d", real_trees),
  sprintf("real_seed_count: %d", real_seed_count),
  sprintf("mc_trees: %d", mc_trees),
  sprintf("mc_reps: %d", mc_reps),
  sprintf("mc_n_train: %d", mc_n_train),
  sprintf("mc_n_test: %d", mc_n_test),
  sprintf("mc_class_sampled_max_switching: %.17g", max(invariance_check$assignment_switching)),
  sprintf(
    "mc_class_sampled_max_prediction_difference: %.17g",
    max(invariance_check$maximum_absolute_prediction_difference)
  ),
  sprintf("mc_class_sampled_max_value_difference: %.17g", max(abs(invariance_check$policy_value_difference)))
)
if (run_scope == "all") {
  metadata <- c(
    metadata,
    sprintf("real_n: %d", nrow(X_real)),
    sprintf("real_clusters: %d", length(unique(clusters_real))),
    sprintf(
      "real_class_sampled_max_switching: %.17g",
      max(real_invariance_check$assignment_switching)
    ),
    sprintf(
      "real_class_sampled_max_prediction_difference: %.17g",
      max(real_invariance_check$maximum_absolute_prediction_difference)
    )
  )
}
writeLines(metadata, file.path(output_dir, "run_metadata.txt"))
if (run_scope == "all") {
  data_md5 <- tools::md5sum(file.path(source_dir, required_data))
  writeLines(
    sprintf("%s  %s", unname(data_md5), basename(names(data_md5))),
    file.path(output_dir, "data_md5.txt")
  )
}
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

cat("Completed class-sampled GRF extension.\n")
print(subset(
  mc_summary,
  metric %in% c("assignment_switching", "policy_regret")
))
