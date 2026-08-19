#!/usr/bin/env Rscript

started_at <- Sys.time()
script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_flag)) {
  normalizePath(sub("^--file=", "", script_flag[[1L]]))
} else {
  normalizePath("replication/code/r/run_dimension_scaling.R")
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
replications <- as.integer(Sys.getenv("RI_DIM_REPS", "20"))
trees <- as.integer(Sys.getenv("RI_DIM_TREES", "2000"))
n_train <- as.integer(Sys.getenv("RI_DIM_N_TRAIN", "3000"))
n_test <- as.integer(Sys.getenv("RI_DIM_N_TEST", "30000"))
base_seed <- as.integer(Sys.getenv("RI_DIM_SEED", "20260811"))
dimensions <- c(40L, 100L, 250L, 500L)
multiplicity <- 8L
output_dir <- Sys.getenv(
  "RI_OUTPUT_DIR",
  file.path(repository_dir, "rerun_outputs", "dimension_scaling")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

append_feature <- function(design, name, train, target) {
  design$train <- cbind(design$train, setNames(data.frame(train), name))
  design$target <- cbind(design$target, setNames(data.frame(target), name))
  design$train <- as.matrix(design$train)
  design$target <- as.matrix(design$target)
  design
}

monotone_transformations <- function(values) {
  shifted <- values - min(values)
  list(
    standardized = (values - mean(values)) / sd(values),
    complement = min(values) + max(values) - values,
    units = 1000 * values + 250,
    log = log1p(shifted),
    rank = rank(values, ties.method = "average"),
    square = shifted^2,
    exponential = exp(values),
    cube = shifted^3
  )
}

make_representation <- function(X.train, X.test, feature) {
  design <- list(train = X.train, target = X.test)
  n.train <- nrow(X.train)
  values <- c(X.train[, feature], X.test[, feature])
  transformations <- monotone_transformations(values)
  for (name in names(transformations)) {
    transformed <- transformations[[name]]
    design <- append_feature(
      design,
      paste0(feature, "__", name),
      transformed[seq_len(n.train)],
      transformed[-seq_len(n.train)]
    )
  }
  design
}

fit_default <- function(design, Y, W, Y.hat, W.hat, seed) {
  fit <- causal_forest(
    design$train,
    Y,
    W,
    Y.hat = Y.hat,
    W.hat = W.hat,
    num.trees = trees,
    seed = seed,
    num.threads = threads
  )
  as.numeric(predict(fit, design$target)$predictions)
}

fit_invariant <- function(design, Y, W, Y.hat, W.hat, seed, dimension) {
  semantic_mtry <- min(ceiling(sqrt(dimension) + 20), dimension)
  fit <- class_sampled_causal_forest(
    design$train,
    Y,
    W,
    X.target = design$target,
    mtry = semantic_mtry,
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

metrics <- function(x2, x1, truth) {
  reversal <- x1 * x2 < 0
  estimated_010 <- pmin(abs(x1), abs(x2)) >= 0.10
  estimated_025 <- pmin(abs(x1), abs(x2)) >= 0.25
  true_025 <- abs(truth) >= 0.25
  c(
    cate_correlation = cor(x1, x2),
    cate_rank_correlation = cor(x1, x2, method = "spearman"),
    sign_reversal = mean(reversal),
    strong_reversal_010 = mean(reversal & estimated_010),
    strong_reversal_025 = mean(reversal & estimated_025),
    reversal_given_true_margin_025 = conditional_mean(reversal, true_025),
    x1_wrong_sign_true_margin_025 = conditional_mean(x1 * truth < 0, true_025),
    x2_wrong_sign_true_margin_025 = conditional_mean(x2 * truth < 0, true_025),
    mean_absolute_prediction_change = mean(abs(x1 - x2))
  )
}

summarize_results <- function(data) {
  groups <- unique(data[, c("dimension", "method", "metric")])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    keep <- data$dimension == groups$dimension[i] &
      data$method == groups$method[i] & data$metric == groups$metric[i]
    values <- data$value[keep]
    standard_error <- sd(values) / sqrt(length(values))
    critical <- qt(0.975, length(values) - 1L)
    rows[[i]] <- cbind(
      groups[i, ],
      data.frame(
        replications = length(values),
        mean = mean(values),
        sd = sd(values),
        se_mean = standard_error,
        ci95_low = mean(values) - critical * standard_error,
        ci95_high = mean(values) + critical * standard_error
      )
    )
  }
  do.call(rbind, rows)
}

rows <- list()
subgroup_rows <- list()
cursor <- 0L
subgroup_cursor <- 0L

for (replication in seq_len(replications)) {
  data_seed <- base_seed + 10000L + replication
  forest_seed <- base_seed + 30000L + replication
  set.seed(data_seed)
  X.train.full <- matrix(
    runif(n_train * max(dimensions)),
    nrow = n_train,
    ncol = max(dimensions)
  )
  X.test.full <- matrix(
    runif(n_test * max(dimensions)),
    nrow = n_test,
    ncol = max(dimensions)
  )
  full_names <- sprintf("x%03d", seq_len(max(dimensions)))
  colnames(X.train.full) <- colnames(X.test.full) <- full_names
  W <- rbinom(n_train, 1, 0.5)
  epsilon <- rnorm(n_train)

  for (dimension in dimensions) {
    X.train <- X.train.full[, seq_len(dimension), drop = FALSE]
    X.test <- X.test.full[, seq_len(dimension), drop = FALSE]
    s1.train <- ifelse(X.train[, "x001"] > 0.5, 1, -1)
    s2.train <- ifelse(X.train[, "x002"] > 0.5, 1, -1)
    s1.test <- ifelse(X.test[, "x001"] > 0.5, 1, -1)
    s2.test <- ifelse(X.test[, "x002"] > 0.5, 1, -1)
    tau.train <- 0.20 + 0.65 * s1.train - 0.55 * s2.train +
      0.10 * (X.train[, "x003"] - 0.5)
    tau.test <- 0.20 + 0.65 * s1.test - 0.55 * s2.test +
      0.10 * (X.test[, "x003"] - 0.5)
    mu.train <- 1 + X.train[, "x004"] + 0.5 * X.train[, "x005"]^2 -
      0.5 * X.train[, "x006"]
    Y <- mu.train + W * tau.train + epsilon
    Y.hat <- mu.train + 0.5 * tau.train
    W.hat <- rep(0.5, n_train)
    x1_design <- make_representation(X.train, X.test, "x001")
    x2_design <- make_representation(X.train, X.test, "x002")
    actual_mtry <- min(ceiling(sqrt(dimension + multiplicity) + 20), dimension + multiplicity)

    cat(sprintf(
      "Replication %02d/%02d: p=%d, raw_p=%d, default_mtry=%d.\n",
      replication,
      replications,
      dimension,
      dimension + multiplicity,
      actual_mtry
    ))
    x1 <- fit_default(x1_design, Y, W, Y.hat, W.hat, forest_seed)
    x2 <- fit_default(x2_design, Y, W, Y.hat, W.hat, forest_seed)
    estimates <- metrics(x2, x1, tau.test)
    for (metric in names(estimates)) {
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        replication = replication,
        dimension = dimension,
        raw_dimension = dimension + multiplicity,
        actual_mtry = actual_mtry,
        method = "ordinary_default",
        metric = metric,
        value = unname(estimates[[metric]])
      )
    }

    subgroup <- paste0(
      "x1_", ifelse(s1.test > 0, "high", "low"),
      "__x2_", ifelse(s2.test > 0, "high", "low")
    )
    for (group in sort(unique(subgroup))) {
      keep <- subgroup == group
      x1_mean <- mean(x1[keep])
      x2_mean <- mean(x2[keep])
      subgroup_cursor <- subgroup_cursor + 1L
      subgroup_rows[[subgroup_cursor]] <- data.frame(
        replication = replication,
        dimension = dimension,
        subgroup = group,
        true_mean = mean(tau.test[keep]),
        x1_mean = x1_mean,
        x2_mean = x2_mean,
        sign_reversal = as.integer(x1_mean * x2_mean < 0),
        strong_sign_reversal = as.integer(
          x1_mean * x2_mean < 0 && min(abs(x1_mean), abs(x2_mean)) >= 0.10
        )
      )
    }

    if (dimension == max(dimensions)) {
      invariant_x1 <- fit_invariant(
        x1_design, Y, W, Y.hat, W.hat, forest_seed, dimension
      )
      invariant_x2 <- fit_invariant(
        x2_design, Y, W, Y.hat, W.hat, forest_seed, dimension
      )
      invariant_estimates <- metrics(invariant_x2, invariant_x1, tau.test)
      for (metric in names(invariant_estimates)) {
        cursor <- cursor + 1L
        rows[[cursor]] <- data.frame(
          replication = replication,
          dimension = dimension,
          raw_dimension = dimension + multiplicity,
          actual_mtry = min(ceiling(sqrt(dimension) + 20), dimension),
          method = "class_sampled_default",
          metric = metric,
          value = unname(invariant_estimates[[metric]])
        )
      }
    }
    gc(verbose = FALSE)
  }
}

results <- do.call(rbind, rows)
summary <- summarize_results(results)
subgroups <- do.call(rbind, subgroup_rows)
write.csv(results, file.path(output_dir, "dimension_metrics.csv"), row.names = FALSE)
write.csv(summary, file.path(output_dir, "dimension_metrics_summary.csv"), row.names = FALSE)
write.csv(subgroups, file.path(output_dir, "dimension_subgroup_signs.csv"), row.names = FALSE)

invariant <- subset(results, method == "class_sampled_default")
if (max(abs(invariant$value[invariant$metric == "mean_absolute_prediction_change"])) > 0) {
  stop("High-dimensional invariant check failed.", call. = FALSE)
}

elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
writeLines(
  c(
    sprintf("completed_utc: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    sprintf("grf_version: %s", as.character(packageVersion("grf"))),
    sprintf("allow_unpinned_grf: %s", allow_unpinned_grf),
    sprintf("threads: %d", threads),
    sprintf("trees: %d", trees),
    sprintf("replications: %d", replications),
    sprintf("n_train: %d", n_train),
    sprintf("n_test: %d", n_test),
    sprintf("multiplicity: %d", multiplicity),
    sprintf("dimensions: %s", paste(dimensions, collapse = ",")),
    sprintf("elapsed_seconds: %.0f", elapsed_seconds),
    sprintf(
      "class_sampled_max_prediction_change: %.17g",
      max(abs(invariant$value[invariant$metric == "mean_absolute_prediction_change"]))
    )
  ),
  file.path(output_dir, "run_metadata.txt")
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
cat(sprintf("Completed dimension panel in %.0f seconds.\n", elapsed_seconds))
print(subset(
  summary,
  metric %in% c("sign_reversal", "strong_reversal_010", "reversal_given_true_margin_025")
))
