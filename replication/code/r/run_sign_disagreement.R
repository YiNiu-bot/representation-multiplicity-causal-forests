#!/usr/bin/env Rscript

started_at <- Sys.time()
script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_flag)) {
  normalizePath(sub("^--file=", "", script_flag[[1L]]))
} else {
  normalizePath("replication/code/r/run_sign_disagreement.R")
}
repository_dir <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
replication_dir <- file.path(repository_dir, "replication")
rlib <- Sys.getenv("RI_RLIB", "")
if (nzchar(rlib) && dir.exists(rlib)) .libPaths(c(rlib, .libPaths()))
suppressPackageStartupMessages(library(grf))
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

threads <- as.integer(Sys.getenv("RI_THREADS", "15"))
replications <- as.integer(Sys.getenv("RI_ROBUST_REPS", "20"))
n_test <- as.integer(Sys.getenv("RI_ROBUST_N_TEST", "30000"))
base_seed <- as.integer(Sys.getenv("RI_ROBUST_SEED", "20260810"))
output_dir <- Sys.getenv(
  "RI_OUTPUT_DIR",
  file.path(repository_dir, "rerun_outputs", "sign_disagreement")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

append_feature <- function(design, name, train, target) {
  design$train <- cbind(design$train, setNames(data.frame(train), name))
  design$target <- cbind(design$target, setNames(data.frame(target), name))
  design$train <- as.matrix(design$train)
  design$target <- as.matrix(design$target)
  design
}

split_pooled <- function(values, n.train) {
  list(train = values[seq_len(n.train)], target = values[-seq_len(n.train)])
}

monotone_transformations <- function(values) {
  shifted <- values - min(values)
  standardized <- (values - mean(values)) / sd(values)
  list(
    standardized = standardized,
    complement = min(values) + max(values) - values,
    units = 1000 * values + 250,
    log = log1p(shifted),
    rank = rank(values, ties.method = "average"),
    square = shifted^2,
    exponential = exp(values),
    cube = shifted^3,
    square_root = sqrt(shifted),
    logistic = plogis(values),
    arctangent = atan(values),
    inverse_hyperbolic_sine = asinh(values),
    reciprocal = 1 / (1 + values),
    normal_cdf = pnorm(standardized),
    fourth_power = shifted^4,
    hyperbolic_sine = sinh(values)
  )
}

enrich_feature <- function(design, feature, multiplicity) {
  if (multiplicity == 0L) return(design)
  n.train <- nrow(design$train)
  values <- c(design$train[, feature], design$target[, feature])
  transformations <- monotone_transformations(values)
  stopifnot(multiplicity <= length(transformations))
  for (name in names(transformations)[seq_len(multiplicity)]) {
    transformed <- split_pooled(transformations[[name]], n.train)
    design <- append_feature(
      design,
      paste0(feature, "__", name),
      transformed$train,
      transformed$target
    )
  }
  design
}

make_representation <- function(X.train, X.test, feature, multiplicity) {
  enrich_feature(
    list(train = X.train, target = X.test),
    feature,
    multiplicity
  )
}

default_mtry <- function(dimension) {
  min(ceiling(sqrt(dimension) + 20), dimension)
}

fit_ordinary <- function(
    design,
    Y,
    W,
    Y.hat,
    W.hat,
    trees,
    mtry_rule,
    seed
) {
  arguments <- list(
    X = design$train,
    Y = Y,
    W = W,
    Y.hat = Y.hat,
    W.hat = W.hat,
    num.trees = trees,
    seed = seed,
    num.threads = threads
  )
  if (mtry_rule == "fixed_6") arguments$mtry <- 6L
  fit <- do.call(causal_forest, arguments)
  as.numeric(predict(fit, design$target)$predictions)
}

fit_class_sampled <- function(
    design,
    Y,
    W,
    Y.hat,
    W.hat,
    trees,
    seed
) {
  fit <- class_sampled_causal_forest(
    design$train,
    Y,
    W,
    X.target = design$target,
    mtry = 6L,
    Y.hat = Y.hat,
    W.hat = W.hat,
    num.trees = trees,
    seed = seed,
    num.threads = threads
  )
  fit$predictions
}

conditional_mean <- function(values, keep) {
  if (!any(keep)) NA_real_ else mean(values[keep])
}

comparison_metrics <- function(x2_prediction, x1_prediction, truth) {
  reversal <- x2_prediction * x1_prediction < 0
  estimated_margin_010 <- pmin(abs(x1_prediction), abs(x2_prediction)) >= 0.10
  estimated_margin_025 <- pmin(abs(x1_prediction), abs(x2_prediction)) >= 0.25
  true_margin_010 <- abs(truth) >= 0.10
  true_margin_025 <- abs(truth) >= 0.25
  c(
    cate_correlation = cor(x2_prediction, x1_prediction),
    cate_rank_correlation = cor(x2_prediction, x1_prediction, method = "spearman"),
    sign_reversal = mean(reversal),
    positive_to_negative = mean(x1_prediction > 0 & x2_prediction < 0),
    negative_to_positive = mean(x1_prediction < 0 & x2_prediction > 0),
    strong_reversal_010 = mean(reversal & estimated_margin_010),
    strong_reversal_025 = mean(reversal & estimated_margin_025),
    reversal_given_estimated_margin_010 = conditional_mean(
      reversal,
      estimated_margin_010
    ),
    reversal_given_estimated_margin_025 = conditional_mean(
      reversal,
      estimated_margin_025
    ),
    reversal_given_true_margin_010 = conditional_mean(reversal, true_margin_010),
    reversal_given_true_margin_025 = conditional_mean(reversal, true_margin_025),
    x1_wrong_sign_true_margin_010 = conditional_mean(
      x1_prediction * truth < 0,
      true_margin_010
    ),
    x2_wrong_sign_true_margin_010 = conditional_mean(
      x2_prediction * truth < 0,
      true_margin_010
    ),
    x1_wrong_sign_true_margin_025 = conditional_mean(
      x1_prediction * truth < 0,
      true_margin_025
    ),
    x2_wrong_sign_true_margin_025 = conditional_mean(
      x2_prediction * truth < 0,
      true_margin_025
    ),
    mean_absolute_prediction_change = mean(abs(x2_prediction - x1_prediction))
  )
}

summarize_long <- function(data, keys) {
  groups <- unique(data[, keys, drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    keep <- rep(TRUE, nrow(data))
    for (key in keys) keep <- keep & data[[key]] == groups[[key]][i]
    values <- data$value[keep]
    standard_error <- sd(values) / sqrt(length(values))
    critical <- qt(0.975, df = length(values) - 1L)
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        replications = length(values),
        mean = mean(values),
        sd = sd(values),
        se_mean = standard_error,
        ci95_low = mean(values) - critical * standard_error,
        ci95_high = mean(values) + critical * standard_error,
        stringsAsFactors = FALSE
      )
    )
  }
  do.call(rbind, rows)
}

cells <- rbind(
  data.frame(
    panel = "dose_response",
    multiplicity = c(1L, 2L, 4L, 8L, 16L),
    trees = 2000L,
    mtry_rule = "fixed_6",
    n_train = 3000L,
    noise_sd = 1,
    method = "ordinary"
  ),
  data.frame(
    panel = "tree_robustness",
    multiplicity = 8L,
    trees = c(800L, 5000L),
    mtry_rule = "fixed_6",
    n_train = 3000L,
    noise_sd = 1,
    method = "ordinary"
  ),
  data.frame(
    panel = "default_mtry",
    multiplicity = 8L,
    trees = 2000L,
    mtry_rule = "grf_default",
    n_train = 3000L,
    noise_sd = 1,
    method = "ordinary"
  ),
  data.frame(
    panel = "larger_sample",
    multiplicity = 8L,
    trees = 2000L,
    mtry_rule = "fixed_6",
    n_train = 10000L,
    noise_sd = 1,
    method = "ordinary"
  ),
  data.frame(
    panel = "noise_robustness",
    multiplicity = 8L,
    trees = 2000L,
    mtry_rule = "fixed_6",
    n_train = 3000L,
    noise_sd = c(0.5, 2),
    method = "ordinary"
  ),
  data.frame(
    panel = "invariant_check",
    multiplicity = 16L,
    trees = 2000L,
    mtry_rule = "fixed_6",
    n_train = 3000L,
    noise_sd = 1,
    method = "class_sampled"
  )
)
tree_scale <- as.numeric(Sys.getenv("RI_ROBUST_TREE_SCALE", "1"))
cells$trees <- pmax(50L, as.integer(round(cells$trees * tree_scale)))
cells$cell_id <- sprintf("cell_%02d", seq_len(nrow(cells)))
cells <- cells[, c(
  "cell_id", "panel", "multiplicity", "trees", "mtry_rule",
  "n_train", "noise_sd", "method"
)]
write.csv(cells, file.path(output_dir, "design_cells.csv"), row.names = FALSE)

max_n_train <- max(cells$n_train)
dimension <- 40L
metric_rows <- list()
subgroup_rows <- list()
metric_cursor <- 0L
subgroup_cursor <- 0L

cat(sprintf(
  "Running %d cells over %d replications with %d threads.\n",
  nrow(cells),
  replications,
  threads
))

for (replication in seq_len(replications)) {
  data_seed <- base_seed + 10000L + replication
  forest_seed <- base_seed + 30000L + replication
  set.seed(data_seed)
  X.train.full <- matrix(
    runif(max_n_train * dimension),
    nrow = max_n_train,
    ncol = dimension
  )
  X.test <- matrix(
    runif(n_test * dimension),
    nrow = n_test,
    ncol = dimension
  )
  colnames(X.train.full) <- colnames(X.test) <- sprintf(
    "x%02d",
    seq_len(dimension)
  )
  W.full <- rbinom(max_n_train, 1, 0.5)
  epsilon.full <- rnorm(max_n_train)

  s1.test <- ifelse(X.test[, "x01"] > 0.5, 1, -1)
  s2.test <- ifelse(X.test[, "x02"] > 0.5, 1, -1)
  tau.test <- 0.20 + 0.65 * s1.test - 0.55 * s2.test +
    0.10 * (X.test[, "x03"] - 0.5)
  subgroup <- paste0(
    "x1_", ifelse(s1.test > 0, "high", "low"),
    "__x2_", ifelse(s2.test > 0, "high", "low")
  )

  for (cell_index in seq_len(nrow(cells))) {
    cell <- cells[cell_index, ]
    index <- seq_len(cell$n_train)
    X.train <- X.train.full[index, , drop = FALSE]
    W <- W.full[index]
    s1.train <- ifelse(X.train[, "x01"] > 0.5, 1, -1)
    s2.train <- ifelse(X.train[, "x02"] > 0.5, 1, -1)
    tau.train <- 0.20 + 0.65 * s1.train - 0.55 * s2.train +
      0.10 * (X.train[, "x03"] - 0.5)
    mu.train <- 1 + X.train[, "x04"] + 0.5 * X.train[, "x05"]^2 -
      0.5 * X.train[, "x06"]
    Y <- mu.train + W * tau.train + cell$noise_sd * epsilon.full[index]
    Y.hat <- mu.train + 0.5 * tau.train
    W.hat <- rep(0.5, cell$n_train)

    x1_design <- make_representation(
      X.train,
      X.test,
      "x01",
      cell$multiplicity
    )
    x2_design <- make_representation(
      X.train,
      X.test,
      "x02",
      cell$multiplicity
    )
    actual_mtry <- if (cell$mtry_rule == "fixed_6") {
      6L
    } else {
      default_mtry(ncol(x1_design$train))
    }
    cat(sprintf(
      "Replication %02d/%02d: %s, panel=%s, m=%d, trees=%d, mtry=%d, n=%d, noise=%.1f, method=%s.\n",
      replication,
      replications,
      cell$cell_id,
      cell$panel,
      cell$multiplicity,
      cell$trees,
      actual_mtry,
      cell$n_train,
      cell$noise_sd,
      cell$method
    ))

    if (cell$method == "ordinary") {
      x1_prediction <- fit_ordinary(
        x1_design,
        Y,
        W,
        Y.hat,
        W.hat,
        cell$trees,
        cell$mtry_rule,
        forest_seed
      )
      x2_prediction <- fit_ordinary(
        x2_design,
        Y,
        W,
        Y.hat,
        W.hat,
        cell$trees,
        cell$mtry_rule,
        forest_seed
      )
    } else {
      x1_prediction <- fit_class_sampled(
        x1_design,
        Y,
        W,
        Y.hat,
        W.hat,
        cell$trees,
        forest_seed
      )
      x2_prediction <- fit_class_sampled(
        x2_design,
        Y,
        W,
        Y.hat,
        W.hat,
        cell$trees,
        forest_seed
      )
    }

    metrics <- comparison_metrics(x2_prediction, x1_prediction, tau.test)
    for (metric in names(metrics)) {
      metric_cursor <- metric_cursor + 1L
      metric_rows[[metric_cursor]] <- data.frame(
        replication = replication,
        data_seed = data_seed,
        forest_seed = forest_seed,
        cell_id = cell$cell_id,
        panel = cell$panel,
        multiplicity = cell$multiplicity,
        trees = cell$trees,
        mtry_rule = cell$mtry_rule,
        actual_mtry = actual_mtry,
        n_train = cell$n_train,
        noise_sd = cell$noise_sd,
        method = cell$method,
        metric = metric,
        value = unname(metrics[[metric]]),
        stringsAsFactors = FALSE
      )
    }

    for (group in sort(unique(subgroup))) {
      keep <- subgroup == group
      x1_mean <- mean(x1_prediction[keep])
      x2_mean <- mean(x2_prediction[keep])
      subgroup_cursor <- subgroup_cursor + 1L
      subgroup_rows[[subgroup_cursor]] <- data.frame(
        replication = replication,
        cell_id = cell$cell_id,
        panel = cell$panel,
        multiplicity = cell$multiplicity,
        trees = cell$trees,
        mtry_rule = cell$mtry_rule,
        actual_mtry = actual_mtry,
        n_train = cell$n_train,
        noise_sd = cell$noise_sd,
        method = cell$method,
        subgroup = group,
        true_mean = mean(tau.test[keep]),
        x1_mean = x1_mean,
        x2_mean = x2_mean,
        sign_reversal = as.integer(x1_mean * x2_mean < 0),
        strong_sign_reversal = as.integer(
          x1_mean * x2_mean < 0 && min(abs(x1_mean), abs(x2_mean)) >= 0.10
        ),
        stringsAsFactors = FALSE
      )
    }
    gc(verbose = FALSE)
  }
}

metrics <- do.call(rbind, metric_rows)
metric_summary <- summarize_long(
  metrics,
  c(
    "cell_id", "panel", "multiplicity", "trees", "mtry_rule",
    "actual_mtry", "n_train", "noise_sd", "method", "metric"
  )
)
subgroups <- do.call(rbind, subgroup_rows)

write.csv(metrics, file.path(output_dir, "robustness_metrics.csv"), row.names = FALSE)
write.csv(
  metric_summary,
  file.path(output_dir, "robustness_metrics_summary.csv"),
  row.names = FALSE
)
write.csv(
  subgroups,
  file.path(output_dir, "robustness_subgroup_signs.csv"),
  row.names = FALSE
)

invariant <- subset(metrics, method == "class_sampled")
if (
  max(abs(invariant$value[invariant$metric == "sign_reversal"])) > 0 ||
    max(abs(invariant$value[invariant$metric == "mean_absolute_prediction_change"])) > 0
) {
  stop("Class-sampled invariance check failed.", call. = FALSE)
}

elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
metadata <- c(
  sprintf("completed_utc: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("grf_version: %s", as.character(packageVersion("grf"))),
  sprintf("allow_unpinned_grf: %s", allow_unpinned_grf),
  sprintf("threads: %d", threads),
  sprintf("replications: %d", replications),
  sprintf("n_test: %d", n_test),
  sprintf("cells: %d", nrow(cells)),
  sprintf("tree_scale: %.4f", tree_scale),
  sprintf("elapsed_seconds: %.0f", elapsed_seconds),
  sprintf(
    "class_sampled_max_prediction_change: %.17g",
    max(abs(invariant$value[invariant$metric == "mean_absolute_prediction_change"]))
  )
)
writeLines(metadata, file.path(output_dir, "run_metadata.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

cat(sprintf("Completed robustness simulation in %.0f seconds.\n", elapsed_seconds))
print(subset(
  metric_summary,
  metric %in% c("sign_reversal", "strong_reversal_010", "reversal_given_true_margin_025")
))
