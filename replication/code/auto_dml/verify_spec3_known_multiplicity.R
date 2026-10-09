argv <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(argv[[1L]], mustWork = TRUE)
output <- normalizePath(argv[[2L]], mustWork = TRUE)
old_output <- file.path(root, "paired_coordinates")
combined <- read.csv(file.path(old_output, "verified_results.csv"))
checks_out <- details <- significance <- list()

independent_score <- function(cache, fit) {
  w <- cache$treatment
  y <- cache$y
  atet <- mean(w * y - fit$psi) / mean(w)
  centered <- fit$psi - mean(fit$psi)
  influence <- -centered / mean(w) +
    (w * y - mean(w * y)) / mean(w) +
    (mean(fit$psi) - mean(w * y)) * (w - mean(w)) / mean(w)^2
  last <- (w * y - mean(w * y)) / mean(w) +
    (mean(fit$psi) - mean(w * y)) * (w - mean(w)) / mean(w)^2
  se <- sqrt((mean(influence^2) + mean(last^2) / (length(w) - 1L)) / length(w))
  abs(atet - fit$row$atet) < 1e-8 && abs(se - fit$row$se) < 1e-8
}

for (dataset in c("NSW", "PSID", "CPS")) {
  stem <- paste0(tolower(dataset), "_spec3")
  directory <- file.path(output, stem)
  cache <- readRDS(file.path(directory, "cache.rds"))
  original <- readRDS(file.path(old_output, stem, "cache.rds"))
  input <- readRDS(file.path(root, paste0(stem, "_x86_64"), "original_inputs.rds"))
  base <- readRDS(file.path(old_output, stem, "fit_full_seed1.rds"))
  read_fit <- function(variant) readRDS(file.path(directory,
    paste0("fit_", variant, "_seed1.rds")))
  reduced <- read_fit("reduced_default")
  fixed <- read_fit("reduced_fixed")
  class <- read_fit("class_full")
  class_reduced <- read_fit("class_reduced")
  fits <- list(reduced = reduced, fixed = fixed, class = class, class_reduced = class_reduced)
  control <- readRDS(file.path(directory, "fold_control_full_1.rds"))
  control_test <- cache$folds[[1L]]
  same_predictions <- function(a, b) {
    identical(a$gamma0, b$gamma0) && identical(a$gamma1, b$gamma1) &&
      identical(a$gamma_observed, b$gamma_observed)
  }
  checks <- c(cache$manual_checks,
    original_input_x_unchanged = identical(cache$x, input[[3L]]) && identical(cache$x, original$x),
    original_y_unchanged = identical(cache$y, input[[1L]]) && identical(cache$y, original$y),
    original_treatment_unchanged = identical(cache$treatment, input[[2L]]) &&
      identical(cache$treatment, original$treatment),
    original_dictionary_unchanged = identical(cache$full, original$full) &&
      identical(cache$zero, original$zero) && identical(cache$one, original$one),
    original_folds_unchanged = identical(cache$folds, original$folds),
    baseline_full_precision = abs(base$row$atet - cache$recovered$result[[3L]]) < 1e-8 &&
      abs(base$row$se - cache$recovered$result[[4L]]) < 1e-8,
    parallel_control_predictions = identical(control$gamma0, unname(base$gamma0[control_test])) &&
      identical(control$gamma1, unname(base$gamma1[control_test])) &&
      identical(control$gamma_observed, unname(base$gamma_observed[control_test])),
    parallel_control_rng = identical(control$rng_end, base$rng_end[[1L]]),
    removed_column_indices = identical(cache$manual_drop,
      c(49L, 54L, 69L, 104L, 220L, 225L, 240L, 275L)),
    fixed_start_rng = all(vapply(fits, function(fit)
      identical(fit$rng_start, base$rng_start), logical(1))),
    fixed_riesz_predictions = all(vapply(fits, function(fit)
      identical(fit$alpha, base$alpha), logical(1))),
    single_bulk_predictions = all(vapply(fits, function(fit)
      all(fit$single_prediction_checks), logical(1))),
    independent_score_formula = all(vapply(fits, function(fit)
      independent_score(cache, fit), logical(1))),
    class_full_reduced_predictions = same_predictions(class, class_reduced),
    class_full_reduced_scores = identical(class$psi, class_reduced$psi) &&
      identical(class$row$atet, class_reduced$row$atet) &&
      identical(class$row$se, class_reduced$row$se),
    class_full_reduced_rng = identical(class$rng_end, class_reduced$rng_end),
    all_reduced_columns_336 = all(vapply(fits, function(fit) fit$row$columns == 336L, logical(1))),
    default_mtry_112 = reduced$row$mtry == 112L,
    fixed_mtry_114 = fixed$row$mtry == 114L,
    class_mtry_112 = class$row$mtry == 112L && class_reduced$row$mtry == 112L,
    baseline_rng_unchanged = identical(base$rng_end, cache$recovered$rng_end))
  stopifnot(all(checks))
  sr <- function(fit) 100 * mean(sign(fit$contrast) != sign(base$contrast))
  idx <- which(combined$dataset == dataset & combined$specification == 3L)
  stopifnot(length(idx) == 1L)
  combined$reduced_columns[idx] <- 336L
  combined$class_columns[idx] <- 336L
  combined$reduced_mtry[idx] <- 112L
  combined$class_mtry[idx] <- 112L
  for (pair in list(c("default", "reduced"), c("fixed", "fixed"), c("class", "class"))) {
    prefix <- pair[[1L]]
    fit <- fits[[pair[[2L]]]]
    combined[idx, paste0(prefix, "_atet")] <- fit$row$atet
    combined[idx, paste0(prefix, "_se")] <- fit$row$se
    combined[idx, paste0(prefix, "_delta")] <- fit$row$atet - base$row$atet
    combined[idx, paste0(prefix, "_sr")] <- sr(fit)
  }
  combined$class_sr_between_representations[idx] <-
    100 * mean(sign(class$contrast) != sign(class_reduced$contrast))
  checks_out[[dataset]] <- data.frame(dataset = dataset, specification = 3L,
    check = names(checks), passed = unname(checks))
  details[[dataset]] <- do.call(rbind, lapply(names(fits), function(name) {
    fit <- fits[[name]]
    data.frame(dataset = dataset, variant = name, atet = fit$row$atet, se = fit$row$se,
      delta = fit$row$atet - base$row$atet, sr = sr(fit),
      mean_absolute_contrast_change = mean(abs(fit$contrast - base$contrast)),
      maximum_absolute_contrast_change = max(abs(fit$contrast - base$contrast)))
  }))
  significance[[dataset]] <- do.call(rbind, lapply(names(fits), function(name) {
    fit <- fits[[name]]
    data.frame(dataset = dataset, specification = 3L, variant = name,
      positive = fit$row$atet > 0,
      significant = abs(fit$row$atet / fit$row$se) > qnorm(0.975),
      baseline_significant = abs(base$row$atet / base$row$se) > qnorm(0.975))
  }))
}
write.csv(combined, file.path(output, "verified_results.csv"), row.names = FALSE)
write.csv(do.call(rbind, checks_out), file.path(output, "verification_checks.csv"), row.names = FALSE)
write.csv(do.call(rbind, details), file.path(output, "spec3_diagnostics.csv"), row.names = FALSE)
write.csv(do.call(rbind, significance), file.path(output, "significance_checks.csv"), row.names = FALSE)

dollar <- function(x) formatC(round(x), format = "f", digits = 0, big.mark = ",")
percent <- function(x) sprintf("%.1f", x)
change <- function(x) if (round(x) == 0) "0" else sprintf("%+d", as.integer(round(x)))
rows <- character()
for (i in seq_len(nrow(combined))) {
  r <- combined[i, ]
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
  if (r$specification == 3L && i < nrow(combined)) rows <- c(rows, "\\midrule")
}
writeLines(rows, file.path(output, "verified_table_rows.tex"))
print(combined[combined$specification == 3L, c("dataset", "original_columns", "reduced_columns",
  "full_atet", "default_atet", "default_delta", "default_sr", "fixed_atet", "fixed_delta",
  "fixed_sr", "class_atet", "class_sr", "seed_delta_max", "seed_sr_max")], row.names = FALSE)
message("Verified all three specification 3 comparisons")
