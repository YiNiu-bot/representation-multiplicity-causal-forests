#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output <- normalizePath(args[1L], mustWork = TRUE)
.libPaths(c(args[2L], .libPaths()))
suppressPackageStartupMessages(library(jsonlite))
config <- fromJSON(file.path(output, "manifest.json"))$configuration
metrics <- read.csv(file.path(output, "metrics.csv"))
comparisons <- read.csv(file.path(output, "comparisons.csv"))
evaluation <- readRDS(file.path(output, "evaluation.rds"))
truth <- evaluation$truth
metric.key <- with(metrics, paste(node_size, multiplicity, replication, variant, metric))
pair.key <- with(comparisons, paste(node_size, multiplicity, replication, regime, metric))
stopifnot(!anyDuplicated(metric.key), !anyDuplicated(pair.key))
metric.lookup <- setNames(metrics$value, metric.key)
pair.lookup <- setNames(comparisons$value, pair.key)
reference <- function(k, r) readRDS(file.path(output, "references", sprintf("k%03d_r%03d.rds", k, r)))
maximum.error <- 0
cell.count <- 0L
for (r in seq_len(config$replications)) {
  first <- reference(config$node_sizes[1L], r)
  for (k in config$node_sizes) {
    ref <- reference(k, r)
    stopifnot(identical(ref$training, first$training),
              identical(ref$canonical, first$canonical), ref$trees == config$trees,
              ref$n_train == config$n_train, ref$forest_seed == first$forest_seed)
    for (c in c(1L, config$multiplicities)) {
      if (c == 1L) {
        P <- ref$predictions[, c(1L, 2L, 1L, 1L, 2L, 1L), drop = FALSE]
        colnames(P) <- as.vector(outer(c("fixed", "default", "linear"), c("x1", "x2"), paste, sep = "_"))
      } else {
        cell <- readRDS(file.path(output, "cells", sprintf("k%03d_c%04d_r%03d.rds", k, c, r)))
        stopifnot(cell$replication == r, cell$node_size == k,
                  cell$multiplicity == c, cell$canonical_verified)
        P <- cell$predictions
        cell.count <- cell.count + 1L
      }
      P <- cbind(P, class_fixed = ref$predictions[, "mtry2"],
                 class_default = ref$predictions[, "mtry6"], class_linear = ref$predictions[, "mtry2"])
      stopifnot(nrow(P) == config$n_test, all(is.finite(P)), all(P != 0))
      for (variant in colnames(P)) {
        prediction <- P[, variant]
        loss <- c(mse = mean((prediction - truth)^2),
          regret = mean(abs(truth) * ((prediction > 0) != (truth > 0))),
          wrong_sign = mean((prediction > 0) != (truth > 0)))
        recorded <- metric.lookup[paste(k, c, r, variant, names(loss))]
        stopifnot(all(is.finite(recorded)))
        maximum.error <- max(maximum.error, abs(recorded - loss))
      }
      for (regime in c("fixed", "default", "linear")) {
        x1 <- P[, paste0(regime, "_x1")]
        x2 <- P[, paste0(regime, "_x2")]
        value <- c(sign_disagreement = mean(x1 * x2 < 0),
                   mean_absolute_prediction_change = mean(abs(x2 - x1)))
        recorded <- pair.lookup[paste(k, c, r, regime, names(value))]
        stopifnot(all(is.finite(recorded)))
        maximum.error <- max(maximum.error, abs(recorded - value))
      }
    }
  }
}
stopifnot(cell.count == config$replications * length(config$node_sizes) * length(config$multiplicities),
          maximum.error < 1e-12)
result <- list(state = "passed", all_prediction_cells_recomputed = cell.count,
               paired_training_data_across_node_sizes = TRUE,
               paired_forest_seeds_across_node_sizes = TRUE,
               zero_predictions = 0L, maximum_loss_disagreement_error = maximum.error)
write_json(result, file.path(output, "PREDICTION_VALIDATION.json"), pretty = TRUE, auto_unbox = TRUE)
print(result)
