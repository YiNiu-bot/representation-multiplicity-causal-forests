argv <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(argv[[1L]], mustWork = TRUE)
output <- argv[[2L]]
mode <- if (length(argv) >= 3L) argv[[3L]] else "all"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(randomForest))
stopifnot(getRversion() == "4.0.3", R.version$arch == "x86_64",
          as.character(packageVersion("randomForest")) == "4.6.14")
RNGversion("4.0.3")

load_functions <- function(path, wanted) {
  found <- character()
  for (expr in parse(path)) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% wanted) {
      eval(expr, envir = .GlobalEnv)
      found <- c(found, as.character(expr[[2L]]))
    }
  }
  stopifnot(setequal(found, wanted))
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
code <- dirname(normalizePath(script))
load_functions(file.path(code, "class_sampling.R"),
  c("lex_less", "class_state", "encode_classes"))
load_functions(file.path(code, "run_paired_coordinates.R"), "score_fit")

prepare <- function(dataset) {
  stem <- paste0(tolower(dataset), "_spec3")
  target <- file.path(output, stem)
  dir.create(target, showWarnings = FALSE)
  original_directory <- file.path(root, "paired_coordinates", stem)
  cache <- readRDS(file.path(original_directory, "cache.rds"))
  base <- readRDS(file.path(original_directory, "fit_full_seed1.rds"))
  p_x <- ncol(cache$x)
  polynomial_names <- c("1.0.0.0", "0.1.0.0", "0.0.1.0", "0.0.0.1")
  polynomial <- match(polynomial_names, colnames(cache$x))
  original <- c(1L, 2L, 6L, 7L)
  stopifnot(p_x == 171L, ncol(cache$full) == 344L,
            identical(polynomial, c(47L, 52L, 67L, 102L)),
            !anyNA(polynomial), !anyNA(cache$full),
            identical(base$rng_start, cache$recovered$rng_start))
  differences <- vapply(seq_along(original), function(j) {
    min(max(abs(cache$x[, original[j]] - cache$x[, polynomial[j]])),
        max(abs(cache$x[, original[j]] + cache$x[, polynomial[j]])))
  }, numeric(1))
  stopifnot(all(differences < 1e-9))
  drop <- c(2L + polynomial, 2L + p_x + polynomial)
  counterparts <- c(2L + original, 2L + p_x + original)
  keep <- setdiff(seq_len(ncol(cache$full)), drop)
  pooled <- rbind(cache$full, cache$zero, cache$one)
  exact <- class_state(pooled)
  reduced <- pooled[, keep, drop = FALSE]
  reduced_state <- class_state(reduced)

  # The construction supplies these eight equivalences. No global rounding is used.
  repaired <- pooled
  repaired[, drop] <- pooled[, counterparts]
  repaired_state <- class_state(repaired)
  checks <- c(
    original_exact_classes_344 = ncol(exact$canonical) == 344L,
    exactly_eight_removed = length(drop) == 8L && length(keep) == 336L,
    kept_values_unchanged = identical(reduced, pooled[, keep, drop = FALSE]),
    repaired_classes_336 = ncol(repaired_state$canonical) == 336L,
    reduced_classes_336 = ncol(reduced_state$canonical) == 336L,
    repaired_reduced_canonical = identical(repaired_state$canonical, reduced_state$canonical),
    repaired_permutation = identical(repaired_state$canonical,
      class_state(repaired[, rev(seq_len(ncol(repaired))), drop = FALSE])$canonical),
    repaired_reversal = identical(repaired_state$canonical,
      class_state(-repaired)$canonical))
  duplicated_input <- cbind(repaired, repaired[, ncol(repaired)])
  colnames(duplicated_input)[ncol(duplicated_input)] <- "added_duplicate"
  checks <- c(checks, repaired_duplicate = identical(repaired_state$canonical,
    class_state(duplicated_input)$canonical))
  stopifnot(all(checks))
  cache$manual_keep <- keep
  cache$manual_drop <- drop
  cache$manual_counterparts <- counterparts
  cache$manual_canonical <- repaired_state$canonical
  cache$manual_reduced_canonical <- reduced_state$canonical
  cache$manual_checks <- checks
  saveRDS(cache, file.path(target, "cache.rds"))
  write.csv(data.frame(dataset = dataset,
    covariate = rep(c("age", "education", "re74", "re75"), 2L),
    treatment_interaction = rep(c(FALSE, TRUE), each = 4L),
    removed_dictionary_column = drop, retained_dictionary_column = counterparts,
    polynomial_x_column = rep(polynomial, 2L),
    maximum_affine_difference = rep(differences, 2L)),
    file.path(target, "removed_columns.csv"), row.names = FALSE)
  write.csv(data.frame(column = colnames(pooled),
    exact_class = unname(exact$class_ids),
    repaired_class = unname(repaired_state$class_ids),
    removed = seq_len(ncol(pooled)) %in% drop),
    file.path(target, "classes.csv"), row.names = FALSE)
  message("PREPARED ", stem, " removing ", paste(drop, collapse = ","),
          " columns 344 -> 336")
  cache
}

run_fold <- function(cache, variant, fold) {
  stem <- paste0(tolower(cache$dataset), "_spec3")
  path <- file.path(output, stem, sprintf("fold_%s_%d.rds", variant, fold))
  if (file.exists(path)) return(data.frame(dataset = cache$dataset,
    variant = variant, fold = fold, cached = TRUE))
  n <- length(cache$y)
  if (variant %in% c("class_full", "class_reduced")) {
    encoded <- if (variant == "class_full") cache$manual_canonical else
      cache$manual_reduced_canonical
    observed <- encoded[seq_len(n), , drop = FALSE]
    zero <- encoded[n + seq_len(n), , drop = FALSE]
    one <- encoded[2L * n + seq_len(n), , drop = FALSE]
  } else {
    keep <- if (variant == "control_full") seq_len(ncol(cache$full)) else cache$manual_keep
    observed <- cache$full[, keep, drop = FALSE]
    zero <- cache$zero[, keep, drop = FALSE]
    one <- cache$one[, keep, drop = FALSE]
  }
  mtry <- if (variant %in% c("reduced_fixed", "control_full")) 114L else 112L
  test <- cache$folds[[fold]]
  train <- setdiff(seq_len(n), test)
  assign(".Random.seed", cache$recovered$rng_start[[fold]], .GlobalEnv)
  rng_start <- .Random.seed
  start <- proc.time()
  message("START ", stem, " ", variant, " fold ", fold)
  forest <- randomForest(x = observed[train, , drop = FALSE], y = cache$y[train],
    ntree = 1000L, mtry = mtry, nodesize = 5L, replace = TRUE, na.action = na.omit)
  rng_end <- .Random.seed
  gamma_observed <- unname(predict(forest, newdata = observed[test, , drop = FALSE]))
  gamma0 <- unname(predict(forest, newdata = zero[test, , drop = FALSE]))
  gamma1 <- unname(predict(forest, newdata = one[test, , drop = FALSE]))
  probe <- unique(c(head(test, 4L), tail(test, 4L)))
  single <- function(matrix) vapply(probe, function(i) as.numeric(predict(forest,
    newdata = matrix[i, , drop = FALSE])), numeric(1))
  checks <- c(single_observed = identical(single(observed), gamma_observed[match(probe, test)]),
    single_zero = identical(single(zero), gamma0[match(probe, test)]),
    single_one = identical(single(one), gamma1[match(probe, test)]))
  if (variant == "control_full") {
    base <- readRDS(file.path(root, "paired_coordinates", stem, "fit_full_seed1.rds"))
    checks <- c(checks,
      original_observed_predictions = identical(gamma_observed, unname(base$gamma_observed[test])),
      original_zero_predictions = identical(gamma0, unname(base$gamma0[test])),
      original_one_predictions = identical(gamma1, unname(base$gamma1[test])),
      original_forest_rng = identical(rng_end, base$rng_end[[fold]]))
  }
  stopifnot(all(checks), all(is.finite(c(gamma_observed, gamma0, gamma1))))
  fit <- list(dataset = cache$dataset, variant = variant, fold = fold, test = test,
    gamma_observed = gamma_observed, gamma0 = gamma0, gamma1 = gamma1,
    rng_start = rng_start, rng_end = rng_end, checks = checks,
    columns = ncol(observed), mtry = mtry, elapsed = (proc.time() - start)[["elapsed"]])
  temporary <- paste0(path, ".tmp")
  saveRDS(fit, temporary)
  stopifnot(file.rename(temporary, path))
  message("DONE ", stem, " ", variant, " fold ", fold,
          " seconds ", sprintf("%.1f", fit$elapsed))
  data.frame(dataset = cache$dataset, variant = variant, fold = fold, cached = FALSE)
}

assemble <- function(cache, variant) {
  stem <- paste0(tolower(cache$dataset), "_spec3")
  n <- length(cache$y)
  gamma0 <- gamma1 <- gamma_observed <- alpha <- psi <- rep(NA_real_, n)
  rng_start <- rng_end <- checks <- vector("list", 5L)
  elapsed <- 0
  for (fold in 1:5) {
    fit <- readRDS(file.path(output, stem, sprintf("fold_%s_%d.rds", variant, fold)))
    test <- cache$folds[[fold]]
    stopifnot(identical(fit$test, test))
    gamma0[test] <- fit$gamma0
    gamma1[test] <- fit$gamma1
    gamma_observed[test] <- fit$gamma_observed
    rho <- cache$recovered$rho[[fold]]
    alpha[test] <- vapply(test, function(i) as.numeric(cache$full[i, ] %*% rho), numeric(1))
    psi[test] <- cache$treatment[test] * gamma0[test] +
      alpha[test] * (cache$y[test] - gamma_observed[test])
    rng_start[[fold]] <- fit$rng_start
    rng_end[[fold]] <- fit$rng_end
    checks[[fold]] <- fit$checks
    elapsed <- elapsed + fit$elapsed
  }
  result <- score_fit(cache, psi)
  stopifnot(all(is.finite(psi)), all(unlist(checks)),
            identical(rng_start, cache$recovered$rng_start))
  row <- data.frame(dataset = cache$dataset, specification = 3L,
    variant = variant, seed = 1L, observations = n, original_columns = 344L,
    columns = 336L, mtry = fit$mtry, trees = 1000L, folds = 5L,
    atet = result$atet, se = result$se, elapsed = elapsed)
  combined <- list(row = row, gamma0 = gamma0, gamma1 = gamma1,
    gamma_observed = gamma_observed, alpha = alpha, psi = psi,
    contrast = gamma1 - gamma0, rng_start = rng_start, rng_end = rng_end,
    single_prediction_checks = vapply(checks, all, logical(1)))
  saveRDS(combined, file.path(output, stem, paste0("fit_", variant, "_seed1.rds")))
  row
}

datasets <- c("NSW", "PSID", "CPS")
caches <- lapply(datasets, function(dataset) {
  path <- file.path(output, paste0(tolower(dataset), "_spec3"), "cache.rds")
  if (file.exists(path)) readRDS(path) else prepare(dataset)
})
if (mode == "prepare") quit(status = 0L)
variants <- c("reduced_default", "reduced_fixed", "class_full", "class_reduced")
jobs <- expand.grid(cell = seq_along(caches), variant = variants,
                    fold = 1:5, stringsAsFactors = FALSE)
controls <- data.frame(cell = seq_along(caches), variant = "control_full", fold = 1L)
jobs <- rbind(jobs, controls)
jobs <- jobs[order(-jobs$cell, jobs$fold), ]
workers <- min(nrow(jobs), max(parallel::detectCores() - 1L, 1L))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
message("Running ", nrow(jobs), " fold jobs on ", workers, " workers")
answers <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  run_fold(caches[[job$cell]], job$variant, job$fold)
}, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(!any(vapply(answers, inherits, logical(1), "try-error")),
          !any(vapply(answers, is.null, logical(1))))
write.csv(do.call(rbind, answers), file.path(output, "completed_fold_jobs.csv"), row.names = FALSE)
rows <- list()
for (cache in caches) {
  for (variant in variants) rows[[length(rows) + 1L]] <- assemble(cache, variant)
}
write.csv(do.call(rbind, rows), file.path(output, "all_fits.csv"), row.names = FALSE)
message("All specification 3 fits completed")
