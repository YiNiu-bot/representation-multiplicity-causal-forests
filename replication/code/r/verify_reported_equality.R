#!/usr/bin/env Rscript

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_flag)) {
  normalizePath(sub("^--file=", "", script_flag[[1L]]))
} else {
  normalizePath("replication/code/r/verify_reported_equality.R")
}
repository_dir <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
replication_dir <- file.path(repository_dir, "replication")
rlib <- Sys.getenv("RI_RLIB", "")
if (nzchar(rlib) && dir.exists(rlib)) .libPaths(c(rlib, .libPaths()))

suppressPackageStartupMessages(library(grf))
expected_grf_version <- "2.4.0"
if (as.character(packageVersion("grf")) != expected_grf_version &&
    Sys.getenv("RI_ALLOW_UNPINNED_GRF", "0") != "1") {
  stop(
    sprintf(
      "Expected grf %s, found %s. Set RI_ALLOW_UNPINNED_GRF=1 only for an explicitly non-replication run.",
      expected_grf_version,
      packageVersion("grf")
    ),
    call. = FALSE
  )
}
source(file.path(replication_dir, "code", "r", "class_sampled_forest.R"))

set.seed(20260813L)
n_train <- 400L
n_target <- 200L
n_total <- n_train + n_target
X_all <- cbind(
  income = runif(n_total),
  distance = runif(n_total),
  household_size = sample(1:8, n_total, replace = TRUE)
)
X_train <- X_all[seq_len(n_train), , drop = FALSE]
X_target <- X_all[-seq_len(n_train), , drop = FALSE]
X_augmented_all <- cbind(
  X_all,
  income_alias = X_all[, "income"],
  income_log = log1p(X_all[, "income"]),
  distance_complement = 1 - X_all[, "distance"],
  household_size_rank = rank(
    X_all[, "household_size"], ties.method = "average"
  )
)
X_augmented_train <- X_augmented_all[seq_len(n_train), , drop = FALSE]
X_augmented_target <- X_augmented_all[-seq_len(n_train), , drop = FALSE]

canonical_quotient <- certify_split_classes(X_train, X_target)
augmented_quotient <- certify_split_classes(
  X_augmented_train, X_augmented_target
)
stopifnot(
  identical(canonical_quotient$X.train, augmented_quotient$X.train),
  identical(canonical_quotient$X.target, augmented_quotient$X.target)
)

W <- rbinom(n_train, 1, 0.5)
tau <- as.numeric(X_train[, "income"] > 0.5) +
  0.5 * X_train[, "distance"]
Y <- W * tau + rnorm(n_train)
common <- list(
  Y.hat = rep(mean(Y), n_train),
  W.hat = rep(0.5, n_train),
  num.trees = 400L,
  mtry = 3.25,
  seed = 20260813L,
  num.threads = as.integer(Sys.getenv("RI_THREADS", "1"))
)
canonical <- do.call(
  class_sampled_causal_forest,
  c(list(
    X.train = X_train,
    X.target = X_target,
    Y = Y,
    W = W
  ), common)
)
augmented <- do.call(
  class_sampled_causal_forest,
  c(list(
    X.train = X_augmented_train,
    X.target = X_augmented_target,
    Y = Y,
    W = W
  ), common)
)

canonical_report <- predict(
  canonical,
  estimate.variance = TRUE
)
augmented_report <- predict(
  augmented,
  estimate.variance = TRUE
)
canonical_se <- sqrt(pmax(canonical_report$variance.estimates, 0))
augmented_se <- sqrt(pmax(augmented_report$variance.estimates, 0))

diagnostics <- data.frame(
  grf_version = as.character(packageVersion("grf")),
  seed = 20260813L,
  training_observations = n_train,
  target_observations = n_target,
  trees = 400L,
  raw_dimension_canonical = ncol(X_train),
  raw_dimension_augmented = ncol(X_augmented_train),
  semantic_dimension = canonical$semantic_dimension,
  maximum_absolute_prediction_difference = max(abs(
    canonical_report$predictions - augmented_report$predictions
  )),
  maximum_absolute_variance_difference = max(abs(
    canonical_report$variance.estimates - augmented_report$variance.estimates
  )),
  maximum_absolute_standard_error_difference = max(abs(
    canonical_se - augmented_se
  )),
  stringsAsFactors = FALSE
)
stopifnot(
  diagnostics$maximum_absolute_prediction_difference == 0,
  diagnostics$maximum_absolute_variance_difference == 0,
  diagnostics$maximum_absolute_standard_error_difference == 0
)

output_dir <- Sys.getenv(
  "RI_OUTPUT_DIR",
  file.path(repository_dir, "rerun_outputs", "reported_equality")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  diagnostics,
  file.path(output_dir, "class_sampled_reported_equality.csv"),
  row.names = FALSE
)
print(diagnostics)
