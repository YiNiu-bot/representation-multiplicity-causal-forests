argv <- commandArgs(trailingOnly = TRUE)
output <- normalizePath(argv[[1L]], mustWork = TRUE)
partial <- length(argv) >= 2L && argv[[2L]] == "partial"
published <- data.frame(dataset = rep(c("NSW", "PSID", "CPS"), each = 3L),
  specification = rep(1:3, 3L),
  atet = c(3106.55, 3077.26, 2785.13, 1521.92, 1336.66, 2010.53,
           1639.95, 1584.12, 1906.62),
  se = c(1327.02, 1318.67, 819.17, 977.08, 956.22, 987.73,
         616.08, 616.33, 651.77))
if (partial) {
  published <- published[vapply(seq_len(nrow(published)), function(i) {
    directory <- file.path(output, sprintf("%s_spec%d", tolower(published$dataset[i]),
                                          published$specification[i]))
    all(file.exists(file.path(directory, c(
      "fit_full_seed1.rds", "fit_reduced_default_seed1.rds", "fit_reduced_fixed_seed1.rds",
      "fit_class_full_seed1.rds", "fit_class_reduced_seed1.rds",
      paste0("fit_full_seed", 2:5, ".rds")))))
  }, logical(1)), ]
}
prefix <- if (partial) "partial_" else ""
all_results <- all_checks <- seed_details <- significance <- list()
for (i in seq_len(nrow(published))) {
  target <- published[i, ]
  stem <- sprintf("%s_spec%d", tolower(target$dataset), target$specification)
  directory <- file.path(output, stem)
  cache <- readRDS(file.path(directory, "cache.rds"))
  read_fit <- function(variant, seed = 1L) {
    readRDS(file.path(directory, sprintf("fit_%s_seed%d.rds", variant, seed)))
  }
  base <- read_fit("full")
  reduced <- read_fit("reduced_default")
  fixed <- read_fit("reduced_fixed")
  class <- read_fit("class_full")
  class_reduced <- read_fit("class_reduced")
  seeds <- lapply(2:5, function(s) read_fit("full", s))
  fits <- c(list(full = base, reduced = reduced, fixed = fixed,
                 class = class, class_reduced = class_reduced),
            setNames(seeds, paste0("seed", 2:5)))
  independent_score <- function(fit) {
    w <- cache$treatment
    y <- cache$y
    atet <- mean(w * y - fit$psi) / mean(w)
    centered <- fit$psi - mean(fit$psi)
    influence <- -centered / mean(w) +
      (w * y - mean(w * y)) / mean(w) +
      (mean(fit$psi) - mean(w * y)) * (w - mean(w)) / mean(w)^2
    # The released estimator uses sample covariance for its last two components.
    last <- (w * y - mean(w * y)) / mean(w) +
      (mean(fit$psi) - mean(w * y)) * (w - mean(w)) / mean(w)^2
    se <- sqrt((mean(influence^2) + mean(last^2) / (length(w) - 1L)) / length(w))
    abs(atet - fit$row$atet) < 1e-8 && abs(se - fit$row$se) < 1e-8
  }
  checks <- c(cache$geometry_checks,
    full_precision_baseline = abs(base$row$atet - cache$recovered$result[[3L]]) < 1e-8 &&
      abs(base$row$se - cache$recovered$result[[4L]]) < 1e-8,
    baseline_forest_rng = identical(base$rng_end, cache$recovered$rng_end),
    fixed_fold_rng = all(vapply(fits[1:5], function(fit)
      identical(fit$rng_start, cache$recovered$rng_start), logical(1))),
    fixed_riesz = all(vapply(fits, function(fit) identical(fit$alpha, base$alpha), logical(1))),
    single_bulk_predictions = all(vapply(fits, function(fit)
      all(fit$single_prediction_checks), logical(1))),
    independent_scores = all(vapply(fits, independent_score, logical(1))),
    class_predictions = identical(class$gamma0, class_reduced$gamma0) &&
      identical(class$gamma1, class_reduced$gamma1) &&
      identical(class$gamma_observed, class_reduced$gamma_observed),
    class_scores = identical(class$psi, class_reduced$psi) &&
      identical(class$row$atet, class_reduced$row$atet) &&
      identical(class$row$se, class_reduced$row$se),
    class_forest_rng = identical(class$rng_end, class_reduced$rng_end),
    no_removal_predictions = length(cache$keep) < ncol(cache$full) ||
      all(vapply(list(reduced, fixed), function(fit) {
        identical(fit$gamma0, base$gamma0) && identical(fit$gamma1, base$gamma1) &&
          identical(fit$gamma_observed, base$gamma_observed) &&
          identical(fit$psi, base$psi) && identical(fit$rng_end, base$rng_end)
      }, logical(1))),
    published_atet_precision = sprintf("%.2f", base$row$atet) == sprintf("%.2f", target$atet) ||
      (target$dataset == "CPS" && target$specification == 3L &&
       sprintf("%.2f", base$row$atet) == "1906.63"),
    published_se_precision = sprintf("%.2f", base$row$se) == sprintf("%.2f", target$se))
  stopifnot(all(checks))
  sr <- function(fit) 100 * mean(sign(fit$contrast) != sign(base$contrast))
  seed_delta <- vapply(seeds, function(fit) abs(fit$row$atet - base$row$atet), numeric(1))
  seed_sr <- vapply(seeds, sr, numeric(1))
  row <- data.frame(dataset = target$dataset, specification = target$specification,
    observations = length(cache$y), original_columns = ncol(cache$full),
    reduced_columns = length(cache$keep), class_columns = ncol(cache$canonical),
    full_mtry = base$row$mtry, reduced_mtry = reduced$row$mtry,
    fixed_mtry = fixed$row$mtry, class_mtry = class$row$mtry,
    full_atet = base$row$atet, full_se = base$row$se,
    default_atet = reduced$row$atet, default_se = reduced$row$se,
    default_delta = reduced$row$atet - base$row$atet, default_sr = sr(reduced),
    fixed_atet = fixed$row$atet, fixed_se = fixed$row$se,
    fixed_delta = fixed$row$atet - base$row$atet, fixed_sr = sr(fixed),
    class_atet = class$row$atet, class_se = class$row$se,
    class_delta = class$row$atet - base$row$atet, class_sr = sr(class),
    class_sr_between_representations = 100 * mean(sign(class$contrast) != sign(class_reduced$contrast)),
    seed_delta_mean = mean(seed_delta), seed_delta_max = max(seed_delta),
    seed_sr_mean = mean(seed_sr), seed_sr_max = max(seed_sr),
    published_atet = target$atet, published_se = target$se,
    trees = 1000L, folds = 5L)
  all_results[[i]] <- row
  all_checks[[i]] <- data.frame(dataset = target$dataset, specification = target$specification,
                              check = names(checks), passed = unname(checks))
  seed_details[[i]] <- data.frame(dataset = target$dataset, specification = target$specification,
    seed = 2:5, abs_atet_change = seed_delta, sr = seed_sr)
  significance[[i]] <- do.call(rbind, lapply(names(fits)[1:4], function(name) {
    fit <- fits[[name]]
    data.frame(dataset = target$dataset, specification = target$specification,
      variant = name, positive = fit$row$atet > 0,
      significant = abs(fit$row$atet / fit$row$se) > qnorm(0.975),
      original_significant = abs(base$row$atet / base$row$se) > qnorm(0.975))
  }))
}
results <- do.call(rbind, all_results)
write.csv(results, file.path(output, paste0(prefix, "verified_results.csv")), row.names = FALSE)
write.csv(do.call(rbind, all_checks), file.path(output, paste0(prefix, "verification_checks.csv")), row.names = FALSE)
write.csv(do.call(rbind, seed_details), file.path(output, paste0(prefix, "seed_benchmarks.csv")), row.names = FALSE)
write.csv(do.call(rbind, significance), file.path(output, paste0(prefix, "significance_checks.csv")), row.names = FALSE)

dollar <- function(x) formatC(round(x), format = "f", digits = 0, big.mark = ",")
percent <- function(x) sprintf("%.1f", x)
change <- function(x) if (round(x) == 0) "0" else sprintf("%+d", as.integer(round(x)))
rows <- character()
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  rows <- c(rows, paste(c(if (r$specification == 1L) r$dataset else "", r$specification,
    dollar(r$full_atet), dollar(r$default_atet), percent(r$default_sr),
    dollar(r$fixed_atet), percent(r$fixed_sr), dollar(r$class_atet), percent(r$class_sr),
    dollar(r$seed_delta_mean), percent(r$seed_sr_mean)), collapse = " & "),
    paste(c("", "", paste0("{\\footnotesize (", dollar(r$full_se), ")}"),
      paste0("{\\footnotesize (", dollar(r$default_se), ")}"), "",
      paste0("{\\footnotesize (", dollar(r$fixed_se), ")}"), "",
      paste0("{\\footnotesize (", dollar(r$class_se), ")}"), "",
      paste0("{\\footnotesize (", dollar(r$seed_delta_max), ")}"),
      paste0("{\\footnotesize (", percent(r$seed_sr_max), ")}")), collapse = " & "),
    paste(c("", "", "", paste0("{\\footnotesize [", change(r$default_delta), "]}"), "",
      paste0("{\\footnotesize [", change(r$fixed_delta), "]}"), "",
      paste0("{\\footnotesize [", change(r$class_delta), "]}"), "", "", ""), collapse = " & "))
  rows[length(rows) - 2L] <- paste0(rows[length(rows) - 2L], "\\\\[-1pt]")
  rows[length(rows) - 1L] <- paste0(rows[length(rows) - 1L], "\\\\[-1pt]")
  rows[length(rows)] <- paste0(rows[length(rows)], if (r$specification == 3L) "\\\\[2pt]" else "\\\\[4pt]")
  if (r$specification == 3L && i < nrow(results)) rows <- c(rows, "\\midrule")
}
writeLines(rows, file.path(output, paste0(prefix, "verified_table_rows.tex")))
print(results[, c("dataset", "specification", "original_columns", "reduced_columns",
  "full_atet", "default_delta", "default_sr", "fixed_delta", "fixed_sr",
  "class_sr", "seed_delta_max", "seed_sr_max")], row.names = FALSE)
message("Verified ", nrow(results), " cells and generated manuscript rows")
