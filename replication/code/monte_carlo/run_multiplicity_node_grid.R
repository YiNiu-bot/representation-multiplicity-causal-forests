#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) {
  value <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(value)) return(default)
  stopifnot(length(value) == 1L)
  substring(value, nchar(name) + 4L)
}
mode <- option("mode")
output <- normalizePath(option("output"), mustWork = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
.libPaths(c(option("rlib"), .libPaths()))
suppressPackageStartupMessages(library(grf))
stopifnot(as.character(packageVersion("grf")) == "2.4.0")
source(file.path(dirname(normalizePath(script)), "class_sampled_grf.R"))
n.train <- as.integer(option("n-train", "3000"))
trees <- as.integer(option("trees", "2000"))
node.size <- as.integer(option("node-size", "5"))
base.seed <- as.integer(option("base-seed", "20261007"))

publish <- function(value, path) {
  temporary <- paste0(path, ".", Sys.getpid(), ".tmp")
  saveRDS(value, temporary, compress = "gzip")
  stopifnot(file.rename(temporary, path))
}
basis <- function(X) cbind(intercept = 1,
  psi1 = ifelse(X[, 1L] >= 0.5, 1, -1),
  psi2 = ifelse(X[, 2L] >= 0.5, 1, -1), smooth = X[, 3L] - 0.5)
coefficients <- c(intercept = 0.20, psi1 = 0.65, psi2 = -0.55, smooth = 0.10)
effect <- function(X) as.numeric(basis(X) %*% coefficients)
default.budget <- function(M) min(ceiling(sqrt(M) + 20), M)
budgets <- function(M) c(fixed = 2L, default = default.budget(M),
                         linear = max(floor(M / 3), 1L))
expected.candidates <- function(M, lambda) {
  k <- seq_len(M - 1L)
  exp(-lambda) + sum(k * dpois(k, lambda)) + M * ppois(M - 1L, lambda, lower.tail = FALSE)
}

if (mode == "prepare") {
  set.seed(as.integer(option("eval-seed", "20261008")))
  X <- matrix(runif(as.integer(option("n-test", "10000")) * 6L), ncol = 6L)
  colnames(X) <- paste0("x", 1:6)
  publish(list(X = X, truth = effect(X), basis = basis(X)), file.path(output, "evaluation.rds"))
  writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
  quit(status = 0L)
}
evaluation <- readRDS(file.path(output, "evaluation.rds"))
truth <- evaluation$truth
replication <- as.integer(option("replication", "0"))
reference.path <- file.path(output, "references", sprintf("k%03d_r%03d.rds", node.size, replication))
path.diagnostics <- function(forest, target, signal.columns) {
  points <- unique(as.integer(seq(1L, nrow(target), length.out = min(128L, nrow(target)))))
  X <- target[points, , drop = FALSE]
  indices <- unique(as.integer(seq(1L, forest[["_num_trees"]], length.out = 32L)))
  diagnostics <- lapply(indices, function(index) {
    left <- forest[["_child_nodes"]][[index]][[1L]]
    right <- forest[["_child_nodes"]][[index]][[2L]]
    split.variable <- forest[["_split_vars"]][[index]] + 1L
    split.value <- forest[["_split_values"]][[index]]
    node <- rep(forest[["_root_nodes"]][[index]] + 1L, nrow(X))
    depth <- integer(nrow(X))
    used <- matrix(FALSE, nrow(X), 2L)
    repeat {
      active <- which(left[node] != 0L | right[node] != 0L)
      if (!length(active)) break
      variable <- split.variable[node[active]]
      for (feature in 1:2) {
        seen <- active[variable %in% signal.columns[[feature]]]
        used[seen, feature] <- TRUE
      }
      values <- X[cbind(active, variable)]
      go.left <- values <= split.value[node[active]]
      node[active] <- ifelse(go.left, left[node[active]], right[node[active]]) + 1L
      depth[active] <- depth[active] + 1L
      stopifnot(max(depth) <= length(left))
    }
    c(mean_path_depth = mean(depth), path_used_x1 = mean(used[, 1L]),
      path_used_x2 = mean(used[, 2L]))
  })
  rowMeans(do.call(cbind, diagnostics))
}
fit <- function(design, training, budget, forest.seed, signal.columns) {
  forest <- causal_forest(design$train, training$Y, training$W,
    Y.hat = training$Y.hat, W.hat = rep(0.5, nrow(design$train)),
    num.trees = trees, mtry = budget, seed = forest.seed, num.threads = 1L,
    sample.fraction = 0.5, honesty = TRUE, honesty.fraction = 0.5, min.node.size = node.size)
  prediction <- as.numeric(predict(forest, design$target, num.threads = 1L)$predictions)
  stopifnot(length(prediction) == length(truth), all(is.finite(prediction)))
  list(prediction = prediction, diagnostics = path.diagnostics(forest, design$target, signal.columns))
}
if (mode == "reference") {
  set.seed(base.seed + replication)
  X <- matrix(runif(n.train * 6L), ncol = 6L)
  colnames(X) <- colnames(evaluation$X)
  W <- rbinom(n.train, 1L, 0.5)
  tau <- effect(X)
  mu <- 1 + X[, 4L] + 0.5 * X[, 5L]^2 - 0.5 * X[, 6L]
  training <- list(X = X, Y = mu + W * tau + rnorm(n.train), W = W,
                   Y.hat = mu + 0.5 * tau)
  canonical <- certify_split_classes(X, evaluation$X)
  stopifnot(canonical$semantic_dimension == 6L)
  signal.columns <- lapply(canonical$class_map$class_id[1:2], function(value) value)
  fitted <- lapply(c(2L, 6L), function(budget)
    fit(list(train = canonical$X.train, target = canonical$X.target), training,
        budget, base.seed + 30000L + replication, signal.columns))
  prediction <- do.call(cbind, lapply(fitted, `[[`, "prediction"))
  diagnostics <- do.call(cbind, lapply(fitted, `[[`, "diagnostics"))
  colnames(prediction) <- colnames(diagnostics) <- c("mtry2", "mtry6")
  publish(list(training = training, canonical = canonical, predictions = prediction,
               diagnostics = diagnostics, node_size = node.size,
               replication = replication, n_train = n.train, trees = trees,
               forest_seed = base.seed + 30000L + replication), reference.path)
  cat("Reference complete", replication, "\n")
  quit(status = 0L)
}

if (mode == "worker") {
  multiplicity <- as.integer(option("multiplicity"))
  stopifnot(multiplicity >= 1L)
  reference <- readRDS(reference.path)
  stopifnot(reference$replication == replication, reference$n_train == n.train,
            reference$trees == trees, reference$node_size == node.size)
  canonical <- reference$canonical
  M <- 6L + multiplicity - 1L
  candidate.budgets <- budgets(M)
  predictions <- matrix(NA_real_, length(truth), 6L,
    dimnames = list(NULL, as.vector(outer(c("fixed", "default", "linear"), c("x1", "x2"), paste, sep = "_"))))
  timing.rows <- depth.rows <- list()
  for (feature in 1:2) {
    columns <- c(1:6, rep(canonical$class_map$class_id[feature], multiplicity - 1L))
    design <- list(train = canonical$X.train[, columns, drop = FALSE],
                    target = canonical$X.target[, columns, drop = FALSE])
    colnames(design$train) <- colnames(design$target) <- sprintf("column_%04d", seq_along(columns))
    # Rank coding makes all redundant columns identical as split actions.
    if (identical(option("verify-classes", "0"), "1")) {
      checked <- certify_split_classes(design$train, design$target)
      stopifnot(identical(checked$X.train, canonical$X.train),
                identical(checked$X.target, canonical$X.target))
      rm(checked)
    }
    signal.columns <- lapply(canonical$class_map$class_id[1:2], function(value) which(columns == value))
    for (regime in names(candidate.budgets)) {
      name <- paste(regime, paste0("x", feature), sep = "_")
      started <- proc.time()[["elapsed"]]
      budget <- candidate.budgets[[regime]]
      if (multiplicity == 1L) {
        predictions[, name] <- reference$predictions[, paste0("mtry", budget)]
        diagnostic <- reference$diagnostics[, paste0("mtry", budget)]
      } else {
        fitted <- fit(design, reference$training, budget, reference$forest_seed, signal.columns)
        predictions[, name] <- fitted$prediction
        diagnostic <- fitted$diagnostics
        rm(fitted)
      }
      depth.rows[[name]] <- data.frame(replication, multiplicity, node_size = node.size,
        variant = name, metric = names(diagnostic), value = unname(diagnostic))
      timing.rows[[name]] <- data.frame(replication, multiplicity, node_size = node.size, variant = name,
        supplied_columns = M, mtry = budget, expected_candidates = expected.candidates(M, budget),
        singleton_offer_probability = expected.candidates(M, budget) / M,
        elapsed_seconds = proc.time()[["elapsed"]] - started)
      cat("Complete", replication, multiplicity, name, "\n")
      gc(verbose = FALSE)
    }
    rm(design)
  }
  publish(list(replication = replication, multiplicity = multiplicity, node_size = node.size, predictions = predictions,
    depth = do.call(rbind, depth.rows), timings = do.call(rbind, timing.rows), canonical_verified = identical(option("verify-classes", "0"), "1")),
    file.path(output, "cells", sprintf("k%03d_c%04d_r%03d.rds", node.size, multiplicity, replication)))
  quit(status = 0L)
}

if (mode == "summarize") {
  replications <- as.integer(option("replications"))
  grid <- as.integer(strsplit(option("multiplicities"), ",", fixed = TRUE)[[1L]])
  nodes <- as.integer(strsplit(option("node-sizes"), ",", fixed = TRUE)[[1L]])
  stopifnot(replications >= 3L)
  depth.rows <- list()
  metric.rows <- coefficient.rows <- pair.rows <- timing.rows <- decomposition <- list()
  get.metrics <- function(prediction) c(mse = mean((prediction - truth)^2),
    regret = mean(abs(truth) * ((prediction > 0) != (truth > 0))),
    wrong_sign = mean((prediction > 0) != (truth > 0)))
  qr.basis <- qr(evaluation$basis)
  for (node.size in nodes) for (multiplicity in c(1L, grid)) {
    predictions <- vector("list", replications)
    for (r in seq_len(replications)) {
      ref <- readRDS(file.path(output, "references", sprintf("k%03d_r%03d.rds", node.size, r)))
      if (multiplicity == 1L) {
        reference.columns <- c(1L, 2L, 1L, 1L, 2L, 1L)
        names <- as.vector(outer(c("fixed", "default", "linear"), c("x1", "x2"), paste, sep = "_"))
        base <- ref$predictions[, reference.columns, drop = FALSE]
        colnames(base) <- names
        cell <- list(replication = r, multiplicity = 1L, node_size = node.size, predictions = base,
                     timings = NULL, depth = NULL, canonical_verified = TRUE)
      } else {
        cell <- readRDS(file.path(output, "cells", sprintf("k%03d_c%04d_r%03d.rds", node.size, multiplicity, r)))
      }
      stopifnot(cell$replication == r, cell$multiplicity == multiplicity,
                cell$node_size == node.size, cell$canonical_verified)
      P <- cbind(cell$predictions, class_fixed = ref$predictions[, "mtry2"],
                 class_default = ref$predictions[, "mtry6"], class_linear = ref$predictions[, "mtry2"])
      stopifnot(all(is.finite(P)))
      if (multiplicity == 1L) {
        for (regime in c("fixed", "default", "linear"))
          stopifnot(identical(P[, paste0(regime, "_x1")], P[, paste0(regime, "_x2")]),
                    identical(P[, paste0(regime, "_x1")], P[, paste0("class_", regime)]))
      }
      predictions[[r]] <- P
      key <- paste(node.size, multiplicity, r)
      timing.rows[[key]] <- cell$timings
      depth.rows[[key]] <- cell$depth
      for (variant in colnames(P)) {
        value <- get.metrics(P[, variant])
        metric.rows[[paste(key, variant)]] <- data.frame(replication = r, multiplicity, node_size = node.size,
          variant, metric = names(value), value = unname(value))
        coef <- as.numeric(qr.coef(qr.basis, P[, variant]))
        coefficient.rows[[paste(key, variant)]] <- data.frame(replication = r, multiplicity, node_size = node.size,
          variant, term = names(coefficients), coefficient = coef,
          retained_signal = coef / unname(coefficients))
      }
      for (regime in c("fixed", "default", "linear")) {
        first <- P[, paste0(regime, "_x1")]
        second <- P[, paste0(regime, "_x2")]
        value <- c(sign_disagreement = mean((first > 0) != (second > 0)),
          mean_absolute_prediction_change = mean(abs(second - first)))
        pair.rows[[paste(key, regime)]] <- data.frame(replication = r, multiplicity, node_size = node.size, regime,
          metric = names(value), value = unname(value))
      }
    }
    for (variant in colnames(predictions[[1L]])) {
      P <- vapply(predictions, function(value) value[, variant], numeric(length(truth)))
      mean.prediction <- rowMeans(P)
      plugin <- mean((mean.prediction - truth)^2)
      variance <- mean(rowMeans((P - mean.prediction)^2))
      mse <- mean((P - truth)^2)
      error <- mse - plugin - variance
      stopifnot(abs(error) < 1e-12)
      decomposition[[paste(node.size, multiplicity, variant)]] <- data.frame(multiplicity, node_size = node.size, variant,
        replications, mse, squared_bias_plugin = plugin, variance_population = variance,
        squared_bias_corrected = plugin - variance / (replications - 1L),
        variance_unbiased = variance * replications / (replications - 1L), identity_error = error)
    }
  }
  summarize <- function(data, grouping, value) {
    index <- split(seq_len(nrow(data)), interaction(data[, grouping, drop = FALSE], drop = TRUE))
    result <- lapply(index, function(rows) {
      values <- data[[value]][rows]
      stopifnot(length(values) == replications, all(is.finite(values)))
      cbind(data[rows[1L], grouping, drop = FALSE], data.frame(replications,
        mean = mean(values), mcse = sd(values) / sqrt(replications)))
    })
    do.call(rbind, result)
  }
  metrics <- do.call(rbind, metric.rows)
  coef <- do.call(rbind, coefficient.rows)
  pairs <- do.call(rbind, pair.rows)
  depth <- do.call(rbind, depth.rows)
  results <- list(metrics = metrics, coefficients = coef, comparisons = pairs,
    depth = depth, depth_summary = summarize(depth, c("multiplicity", "node_size", "variant", "metric"), "value"),
    timings = do.call(rbind, timing.rows), bias_variance = do.call(rbind, decomposition),
    metric_summary = summarize(metrics, c("multiplicity", "node_size", "variant", "metric"), "value"),
    signal_summary = summarize(coef, c("multiplicity", "node_size", "variant", "term"), "retained_signal"),
    comparison_summary = summarize(pairs, c("multiplicity", "node_size", "regime", "metric"), "value"))
  for (name in names(results)) write.csv(results[[name]], file.path(output, paste0(name, ".csv")), row.names = FALSE)
  cat("Validated every replication and all decomposition identities\n")
  quit(status = 0L)
}
stop("Unknown mode")
