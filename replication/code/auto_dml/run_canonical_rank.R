#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)
cache_dir <- normalizePath(argv[[1L]], mustWork = TRUE)
output <- argv[[2L]]
workers <- if (length(argv) > 2L) as.integer(argv[[3L]]) else 3L
script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
package <- dirname(normalizePath(sub("^--file=", "", script_flag[[1L]])))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(randomForest))
stopifnot(getRversion() == "4.0.5", as.character(packageVersion("randomForest")) == "4.6.14")
args <- list(ntree = 1000L)
wanted <- c("make_dictionary", "same_weak_order", "weak_order_classes", "forest_fit", "run_variant")
found <- character()
for (expr in parse(file.path(package, "run_native.R"))) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) && is.symbol(expr[[2L]]) &&
      as.character(expr[[2L]]) %in% wanted) {
    eval(expr, envir = .GlobalEnv)
    found <- c(found, as.character(expr[[2L]]))
  }
}
stopifnot(setequal(found, wanted))
native_forest_fit <- forest_fit

# Transform only forest inputs. run_variant still evaluates the original
# orthogonal score with the original raw dictionary and cached Riesz estimates.
canonicalize <- function(x, maps) {
  x <- as.matrix(x)
  out <- x
  for (name in colnames(x)) {
    map <- maps[[name]]
    index <- match(x[, name], map$values)
    stopifnot(!anyNA(index))
    out[, name] <- map$ranks[index]
  }
  out
}
forest_fit <- function(x, y, ntree, mtry = NULL) {
  if (!is.null(mtry)) stopifnot(mtry == fixed_mtry)
  fit <- native_forest_fit(canonicalize(x, rank_maps), y, ntree, fixed_mtry)
  structure(list(model = fit, maps = rank_maps), class = "canonical_forest")
}
predict.canonical_forest <- function(object, newdata, ...) {
  predict(object$model, newdata = canonicalize(newdata, object$maps), ...)
}

run_cell <- function(cell) {
  dataset <- cell$dataset
  spec <- cell$specification
  stem <- sprintf("%s_spec%d", tolower(dataset), spec)
  target <- file.path(output, stem)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  cache <- readRDS(file.path(cache_dir, paste0(stem, "_theoretical.rds")))
  baseline_cache <- serialize(cache, NULL, version = 2)
  pooled <- rbind(cache$full_dictionary,
                 make_dictionary(rep(0, nrow(cache$x)), cache$x),
                 make_dictionary(rep(1, nrow(cache$x)), cache$x))
  classes <- weak_order_classes(pooled)
  stopifnot(identical(classes$keep, cache$classes$keep), !anyNA(pooled))
  maps <- setNames(vector("list", ncol(pooled)), colnames(pooled))
  detail <- classes$detail
  for (j in seq_len(ncol(pooled))) {
    representative <- detail$representative_index[match(j, detail$member_index)]
    values <- sort(unique(pooled[, j]))
    canonical <- match(pooled[, representative], sort(unique(pooled[, representative])))
    code <- canonical[match(values, pooled[, j])]
    maps[[j]] <- list(values = values, ranks = code)
    stopifnot(identical(as.integer(code[match(pooled[, j], values)]), as.integer(canonical)))
  }
  assign("rank_maps", maps, envir = .GlobalEnv)
  assign("fixed_mtry", max(floor(ncol(pooled) / 3), 1L), envir = .GlobalEnv)
  ranked <- canonicalize(pooled, maps)
  certificate <- detail
  certificate$identical_canonical_values <- vapply(seq_len(nrow(detail)), function(k) {
    identical(ranked[, detail$member_index[k]], ranked[, detail$representative_index[k]])
  }, logical(1))
  stopifnot(all(certificate$identical_canonical_values))
  write.csv(certificate, file.path(target, "certificate.csv"), row.names = FALSE)
  protocol <- data.frame(dataset = dataset, specification = spec,
    full_columns = ncol(pooled), reduced_columns = length(classes$keep),
    mtry = fixed_mtry, trees = args$ntree, forest_seed = 1L, fold_seed = 1L,
    folds = length(cache$folds), c_rr = cache$c_rr,
    rank_support = "observed and both counterfactual arrays pooled; no outcomes",
    rank_ties = "dense integer ranks; equal values tied",
    orientation = "smallest original column index in each weak-order class",
    riesz = "original raw dictionary and cached fold-specific coefficients",
    stringsAsFactors = FALSE)
  write.csv(protocol, file.path(target, "protocol.csv"), row.names = FALSE)
  cat(stem, "started", format(Sys.time()), "\n")
  fits <- list()
  for (variant in c("published_full", "quotient_fixed_mtry")) {
    fit <- run_variant(cache, variant, seed = 1L)
    fit$result$variant <- if (variant == "published_full") "canonical_full" else "canonical_reduced"
    fit$result$mtry <- fixed_mtry
    fits[[fit$result$variant]] <- fit
    saveRDS(fit, file.path(target, paste0(fit$result$variant, ".rds")))
    cat(stem, fit$result$variant, "completed", format(Sys.time()), "\n")
  }
  stopifnot(identical(serialize(cache, NULL, version = 2), baseline_cache))
  full <- fits$canonical_full
  reduced <- fits$canonical_reduced
  score_checks <- lapply(fits, function(fit) {
    independent <- mean(cache$treatment * cache$y - fit$psi) / mean(cache$treatment)
    stopifnot(abs(independent - fit$result$atet) < 1e-8,
              all(is.finite(fit$psi)), all(is.finite(fit$contrast)))
    data.frame(variant = fit$result$variant, score_reconstruction_gap = independent - fit$result$atet)
  })
  write.csv(do.call(rbind, score_checks), file.path(target, "score_checks.csv"), row.names = FALSE)
  diagnostics <- do.call(rbind, lapply(c("all", "treated"), function(subset_name) {
    use <- if (subset_name == "all") rep(TRUE, length(cache$treatment)) else cache$treatment == 1
    data.frame(dataset = dataset, specification = spec, seed = 1L,
      subset = subset_name, observations = sum(use),
      atet_difference = reduced$result$atet - full$result$atet,
      se_difference = reduced$result$se - full$result$se,
      gamma0_rmse = sqrt(mean((reduced$gamma0[use] - full$gamma0[use])^2)),
      gamma1_rmse = sqrt(mean((reduced$gamma1[use] - full$gamma1[use])^2)),
      contrast_rmse = sqrt(mean((reduced$contrast[use] - full$contrast[use])^2)),
      contrast_correlation = cor(reduced$contrast[use], full$contrast[use]),
      contrast_sign_disagreement = mean(sign(reduced$contrast[use]) != sign(full$contrast[use])) )
  }))
  rows <- do.call(rbind, lapply(fits, `[[`, "result"))
  write.csv(rows, file.path(target, "results.csv"), row.names = FALSE)
  write.csv(diagnostics, file.path(target, "diagnostics.csv"), row.names = FALSE)
  list(results = rows, diagnostics = diagnostics, protocol = protocol)
}
cells <- expand.grid(dataset = c("NSW", "PSID", "CPS"), specification = 1:2, stringsAsFactors = FALSE)
answers <- parallel::mclapply(seq_len(nrow(cells)), function(i) run_cell(cells[i, ]),
                             mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(!any(vapply(answers, inherits, logical(1), "try-error")))
for (name in c("results", "diagnostics", "protocol")) {
  write.csv(do.call(rbind, lapply(answers, `[[`, name)), file.path(output, paste0(name, ".csv")), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
cat("All six canonical-rank cells completed. Native results were not changed.\n")
