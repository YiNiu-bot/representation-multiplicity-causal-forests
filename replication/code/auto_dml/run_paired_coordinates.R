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

dictionary <- function(treatment, x) {
  out <- cbind(1, treatment, x, x * treatment)
  colnames(out) <- paste0("column_", seq_len(ncol(out)))
  out
}

prepare_cell <- function(dataset, spec) {
  stem <- sprintf("%s_spec%d", tolower(dataset), spec)
  target <- file.path(output, stem)
  dir.create(target, showWarnings = FALSE)
  original <- file.path(root, paste0(stem, "_x86_64"))
  input <- readRDS(file.path(original, "original_inputs.rds"))
  recovered <- readRDS(file.path(original, "literal_result.rds"))
  y <- input[[1L]]
  treatment <- input[[2L]]
  x <- input[[3L]]
  n <- length(y)
  full <- dictionary(treatment, x)
  zero <- dictionary(rep(0, n), x)
  one <- dictionary(rep(1, n), x)
  pooled <- rbind(full, zero, one)
  state <- class_state(pooled)
  keep <- which(!duplicated(unname(state$class_ids)))
  reduced_state <- class_state(pooled[, keep, drop = FALSE])
  permutation <- class_state(pooled[, rev(seq_len(ncol(pooled))), drop = FALSE])
  reversal <- class_state(-pooled[, rev(seq_len(ncol(pooled))), drop = FALSE])
  duplicate <- cbind(pooled, pooled[, ncol(pooled)])
  colnames(duplicate)[ncol(duplicate)] <- "extra_duplicate"
  duplicate_state <- class_state(duplicate)
  checks <- c(
    reduced_canonical = identical(state$canonical, reduced_state$canonical),
    permuted_canonical = identical(state$canonical, permutation$canonical),
    reversed_canonical = identical(state$canonical, reversal$canonical),
    duplicated_canonical = identical(state$canonical, duplicate_state$canonical),
    original_encoding = identical(encode_classes(pooled, state), state$canonical),
    reduced_encoding = identical(encode_classes(pooled[, keep, drop = FALSE],
                                                reduced_state), state$canonical))
  set.seed(1L)
  folds <- split(sample(n, n, replace = FALSE), as.factor(1:5))
  checks <- c(checks, original_fold_rng = identical(.Random.seed,
                                                   recovered$rng_start[[1L]]))
  stopifnot(all(checks), length(recovered$rho) == 5L)
  cache <- list(dataset = dataset, specification = spec, y = y,
    treatment = treatment, x = x, full = full, zero = zero, one = one,
    keep = keep, canonical = state$canonical,
    reduced_canonical = reduced_state$canonical,
    class_ids = state$class_ids, folds = folds,
    recovered = recovered, geometry_checks = checks)
  saveRDS(cache, file.path(target, "cache.rds"))
  write.csv(data.frame(column = colnames(pooled), class_id = state$class_ids),
            file.path(target, "classes.csv"), row.names = FALSE)
  message("PREPARED ", stem, " original columns ", ncol(full),
          " classes ", ncol(state$canonical))
  cache
}

score_fit <- function(cache, psi) {
  y <- cache$y
  w <- cache$treatment
  order <- unlist(cache$folds, use.names = FALSE)
  p <- psi[order]
  est <- mean(p)
  tbar <- mean(w)
  tybar <- mean(y * w)
  centered <- p - est
  v11 <- mean(centered^2)
  v12 <- mean(centered * w[order] * y[order])
  v13 <- mean(centered * w[order])
  v22 <- var(w * y)
  v23 <- cov(w * y, w)
  v33 <- var(w)
  variance <- matrix(c(v11, v12, v13, v12, v22, v23, v13, v23, v33), 3L)
  gradient <- matrix(c(-1 / tbar, 1 / tbar, (est - tybar) / tbar^2), 1L)
  list(atet = (tybar - est) / tbar,
       se = as.numeric(sqrt((gradient %*% variance %*% t(gradient)) / length(w))))
}

run_fit <- function(cache, variant, seed) {
  stem <- sprintf("%s_spec%d", tolower(cache$dataset), cache$specification)
  path <- file.path(output, stem, sprintf("fit_%s_seed%d.rds", variant, seed))
  if (file.exists(path)) return(readRDS(path)$row)
  n <- length(cache$y)
  p <- ncol(cache$full)
  keep <- if (variant %in% c("reduced_default", "reduced_fixed")) cache$keep else seq_len(p)
  canonical <- variant %in% c("class_full", "class_reduced")
  if (canonical) {
    encoded <- if (variant == "class_full") cache$canonical else cache$reduced_canonical
    observed <- encoded[seq_len(n), , drop = FALSE]
    zero <- encoded[n + seq_len(n), , drop = FALSE]
    one <- encoded[2L * n + seq_len(n), , drop = FALSE]
  } else {
    observed <- cache$full[, keep, drop = FALSE]
    zero <- cache$zero[, keep, drop = FALSE]
    one <- cache$one[, keep, drop = FALSE]
  }
  mtry <- max(floor(if (variant == "reduced_fixed") p / 3 else ncol(observed) / 3), 1L)
  if (seed != 1L) {
    set.seed(seed)
    invisible(sample(n, n, replace = FALSE))
  }
  gamma0 <- gamma1 <- gamma_observed <- alpha <- psi <- rep(NA_real_, n)
  rng_start <- rng_end <- vector("list", 5L)
  prediction_checks <- logical(5L)
  start <- proc.time()
  message("START ", stem, " ", variant, " seed ", seed)
  for (fold in seq_len(5L)) {
    test <- cache$folds[[fold]]
    train <- setdiff(seq_len(n), test)
    if (seed == 1L) assign(".Random.seed", cache$recovered$rng_start[[fold]], .GlobalEnv)
    rng_start[[fold]] <- .Random.seed
    forest <- randomForest(x = observed[train, , drop = FALSE], y = cache$y[train],
      ntree = 1000L, mtry = mtry, nodesize = 5L, replace = TRUE, na.action = na.omit)
    rng_end[[fold]] <- .Random.seed
    if (variant == "full" && seed == 1L) {
      stopifnot(identical(rng_end[[fold]], cache$recovered$rng_end[[fold]]))
    }
    gamma_observed[test] <- predict(forest, newdata = observed[test, , drop = FALSE])
    gamma0[test] <- predict(forest, newdata = zero[test, , drop = FALSE])
    gamma1[test] <- predict(forest, newdata = one[test, , drop = FALSE])
    probe <- unique(c(head(test, 4L), tail(test, 4L)))
    single <- vapply(probe, function(i) as.numeric(predict(forest,
      newdata = observed[i, , drop = FALSE])), numeric(1))
    single0 <- vapply(probe, function(i) as.numeric(predict(forest,
      newdata = zero[i, , drop = FALSE])), numeric(1))
    prediction_checks[[fold]] <- identical(unname(single), unname(gamma_observed[probe])) &&
      identical(unname(single0), unname(gamma0[probe]))
    stopifnot(prediction_checks[[fold]])
    rho <- cache$recovered$rho[[fold]]
    alpha[test] <- vapply(test, function(i) as.numeric(cache$full[i, ] %*% rho), numeric(1))
    psi[test] <- cache$treatment[test] * gamma0[test] +
      alpha[test] * (cache$y[test] - gamma_observed[test])
    rm(forest)
  }
  result <- score_fit(cache, psi)
  stopifnot(all(is.finite(psi)), all(is.finite(gamma1 - gamma0)), all(prediction_checks))
  if (variant == "full" && seed == 1L) {
    stopifnot(abs(result$atet - cache$recovered$result[[3L]]) < 1e-8,
              abs(result$se - cache$recovered$result[[4L]]) < 1e-8)
  }
  row <- data.frame(dataset = cache$dataset, specification = cache$specification,
    variant = variant, seed = seed, observations = n, original_columns = p,
    columns = ncol(observed), mtry = mtry, trees = 1000L, folds = 5L,
    atet = result$atet, se = result$se, elapsed = (proc.time() - start)[["elapsed"]])
  fit <- list(row = row, gamma0 = gamma0, gamma1 = gamma1,
    gamma_observed = gamma_observed, alpha = alpha, psi = psi,
    contrast = gamma1 - gamma0, rng_start = rng_start, rng_end = rng_end,
    single_prediction_checks = prediction_checks)
  saveRDS(fit, path)
  message("DONE ", stem, " ", variant, " seed ", seed,
          " ATET ", sprintf("%.2f", result$atet), " SE ", sprintf("%.2f", result$se))
  row
}

cells <- expand.grid(dataset = c("NSW", "PSID", "CPS"), specification = 1:3,
                     stringsAsFactors = FALSE)
if (mode == "pilot") cells <- cells[c(1L, 9L), ]
caches <- lapply(seq_len(nrow(cells)), function(i) {
  stem <- sprintf("%s_spec%d", tolower(cells$dataset[i]), cells$specification[i])
  path <- file.path(output, stem, "cache.rds")
  if (file.exists(path)) readRDS(path) else prepare_cell(cells$dataset[i], cells$specification[i])
})
if (mode == "prepare") quit(status = 0L)
variants <- c("full", "reduced_default", "reduced_fixed", "class_full", "class_reduced",
              rep("full", 4L))
seeds <- c(rep(1L, 5L), 2:5)
jobs <- expand.grid(cell = seq_along(caches), variant_index = seq_along(variants))
if (mode == "pilot") jobs <- jobs[jobs$variant_index == 1L, ]
if (mode == "paired") {
  pilot_cells <- cells$dataset[jobs$cell] == "NSW" & cells$specification[jobs$cell] == 1L |
    cells$dataset[jobs$cell] == "CPS" & cells$specification[jobs$cell] == 3L
  jobs <- jobs[!(pilot_cells & jobs$variant_index == 1L), ]
}
jobs <- jobs[order(-cells$specification[jobs$cell], -match(cells$dataset[jobs$cell],
                                                c("NSW", "PSID", "CPS"))), ]
workers <- min(nrow(jobs), max(parallel::detectCores() - 1L, 1L))
if (mode == "paired") workers <- min(workers, parallel::detectCores() - 2L)
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
message("Running ", nrow(jobs), " fits on ", workers, " workers")
answers <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  run_fit(caches[[job$cell]], variants[job$variant_index], seeds[job$variant_index])
}, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(!any(vapply(answers, inherits, logical(1), "try-error")),
          !any(vapply(answers, is.null, logical(1))))
write.csv(do.call(rbind, answers), file.path(output, paste0(mode, "_fits.csv")), row.names = FALSE)
message("All requested fits completed")
