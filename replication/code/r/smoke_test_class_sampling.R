#!/usr/bin/env Rscript

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_flag)) {
  normalizePath(sub("^--file=", "", script_flag[[1L]]))
} else {
  normalizePath("replication/code/r/smoke_test_class_sampling.R")
}
repository_dir <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
replication_dir <- file.path(repository_dir, "replication")
rlib <- Sys.getenv("RI_RLIB", "")
if (nzchar(rlib) && dir.exists(rlib)) .libPaths(c(rlib, .libPaths()))

suppressPackageStartupMessages(library(grf))
source(file.path(replication_dir, "code", "r", "class_sampled_forest.R"))

set.seed(20260808L)
n <- 400L
X <- cbind(
  income = runif(n),
  distance = runif(n),
  household_size = sample(1:8, n, replace = TRUE)
)
X_augmented <- cbind(
  X,
  income_alias = X[, "income"],
  income_log = log1p(X[, "income"]),
  distance_complement = 1 - X[, "distance"],
  household_size_rank = rank(X[, "household_size"], ties.method = "average")
)
canonical_quotient <- certify_split_classes(X)
quotient <- certify_split_classes(X_augmented)
permuted_quotient <- certify_split_classes(
  X_augmented[, rev(seq_len(ncol(X_augmented))), drop = FALSE]
)
replacement_quotient <- certify_split_classes(
  X_augmented[, c(
    "income_alias",
    "distance_complement",
    "household_size_rank"
  ), drop = FALSE]
)
stopifnot(
  quotient$raw_dimension == 7L,
  quotient$semantic_dimension == 3L,
  identical(quotient$X.train, canonical_quotient$X.train),
  identical(permuted_quotient$X.train, canonical_quotient$X.train),
  identical(replacement_quotient$X.train, canonical_quotient$X.train)
)

W <- rbinom(n, 1, 0.5)
tau <- as.numeric(X[, "income"] > 0.5) + 0.5 * X[, "distance"]
Y <- W * tau + rnorm(n)
common <- list(
  Y.hat = rep(mean(Y), n),
  W.hat = rep(0.5, n),
  num.trees = 200L,
  mtry = 3L,
  seed = 20260808L,
  num.threads = as.integer(Sys.getenv("RI_THREADS", "1"))
)
canonical <- do.call(
  class_sampled_causal_forest,
  c(list(X.train = X, Y = Y, W = W), common)
)
corrected <- do.call(
  class_sampled_causal_forest,
  c(list(X.train = X_augmented, Y = Y, W = W), common)
)
stopifnot(max(abs(canonical$predictions - corrected$predictions)) == 0)

cat("Class-sampled GRF smoke test passed.\n")
