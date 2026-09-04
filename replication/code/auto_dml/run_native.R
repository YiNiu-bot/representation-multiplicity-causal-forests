#!/usr/bin/env Rscript

# Paired replication of the random-forest Auto-DML application in
# Chernozhukov, Newey, and Singh (Econometrica, 2022).

parse_args <- function(x) {
  out <- list(
    package_root = NA_character_,
    output_dir = NA_character_,
    phase = "all",
    tuning = "theoretical",
    datasets = c("NSW", "PSID", "CPS"),
    specs = 1:2,
    seeds = 1L,
    ntree = 1000L,
    rebuild_cache = FALSE
  )
  for (arg in x) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Arguments must have the form --name=value")
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- pieces[[1]]
    value <- paste(pieces[-1], collapse = "=")
    if (!key %in% names(out)) stop("Unknown argument: ", key)
    out[[key]] <- value
  }
  out$datasets <- strsplit(paste(out$datasets, collapse = ","), ",", fixed = TRUE)[[1]]
  out$specs <- as.integer(strsplit(paste(out$specs, collapse = ","), ",", fixed = TRUE)[[1]])
  out$seeds <- as.integer(strsplit(paste(out$seeds, collapse = ","), ",", fixed = TRUE)[[1]])
  out$ntree <- as.integer(out$ntree)
  out$rebuild_cache <- tolower(as.character(out$rebuild_cache)) %in% c("true", "1", "yes")
  out
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
script_dir <- dirname(script_path)
args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.na(args$package_root)) stop("Provide --package_root=/path/to/rrr_lasso_NSW_blackbox")
args$package_root <- normalizePath(args$package_root, mustWork = TRUE)
if (is.na(args$output_dir)) args$output_dir <- file.path(script_dir, "output")
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(args$output_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(args$datasets %in% c("NSW", "PSID", "CPS"))) stop("Unknown data set")
if (!all(args$specs %in% 1:3)) stop("Specifications must be 1, 2, or 3")
if (!args$phase %in% c("baseline", "paired", "seed_noise", "all")) {
  stop("phase must be baseline, paired, seed_noise, or all")
}
if (!args$tuning %in% c("theoretical", "cross_validated")) {
  stop("tuning must be theoretical or cross_validated")
}
if (anyNA(args$seeds) || any(args$seeds < 1L)) stop("Seeds must be positive integers")

suppressPackageStartupMessages({
  library(nnet)
  library(randomForest)
})
stopifnot(getRversion() == "4.0.5", as.character(packageVersion("randomForest")) == "4.6.14")

source(file.path(args$package_root, "primitives.R"), local = TRUE)
source(file.path(args$package_root, "stage0.R"), local = TRUE)
source(file.path(args$package_root, "cv.R"), local = TRUE)
source(file.path(args$package_root, "specifications_intersection.R"), local = TRUE)

names74 <- c(
  "treat", "age", "education", "black", "hispanic", "married",
  "nodegree", "re74", "re75", "re78"
)

read_lalonde <- function(stem) {
  out <- read.table(file.path(args$package_root, paste0(stem, ".txt")))
  names(out) <- names74
  out
}

treated <- read_lalonde("nswre74_treated")
datasets <- list(
  NSW = rbind(treated, read_lalonde("nswre74_control")),
  PSID = rbind(treated, read_lalonde("psid_controls")),
  CPS = rbind(treated, read_lalonde("cps_controls"))
)

base_names <- c(
  "age", "education", "married", "black", "hispanic", "re74", "re75",
  "age_sq", "education_sq", "re74_sq", "re75_sq"
)
spec_names <- list(
  `1` = base_names,
  `2` = c(base_names, "re74_zero", "re75_zero", "nodegree")
)

published <- data.frame(
  dataset = rep(c("NSW", "PSID", "CPS"), each = 3L),
  specification = rep(1:3, times = 3L),
  published_atet = c(
    3106.55, 3077.26, 2785.13,
    1521.92, 1336.66, 2010.53,
    1639.95, 1584.12, 1906.62
  ),
  published_se = c(
    1327.02, 1318.67, 819.17,
    977.08, 956.22, 987.73,
    616.08, 616.33, 651.77
  )
)

make_dictionary <- function(treatment, x) {
  x <- as.matrix(x)
  out <- cbind(1, treatment, x, x * treatment)
  x_names <- colnames(x)
  colnames(out) <- c("intercept", "treatment", x_names, paste0("treatment_x_", x_names))
  out
}

same_weak_order <- function(x, y) {
  rx <- rank(x, ties.method = "min")
  ry <- rank(y, ties.method = "min")
  if (identical(rx, ry)) return("increasing")
  if (identical(rx, rank(-y, ties.method = "min"))) return("decreasing")
  NA_character_
}

weak_order_classes <- function(mat) {
  mat <- as.matrix(mat)
  p <- ncol(mat)
  parent <- seq_len(p)
  find_root <- function(i) {
    while (parent[[i]] != i) i <- parent[[i]]
    i
  }
  unite <- function(i, j) {
    ri <- find_root(i)
    rj <- find_root(j)
    if (ri != rj) parent[[rj]] <<- ri
  }
  nonconstant <- which(apply(mat, 2, function(z) length(unique(z)) > 1L))
  if (length(nonconstant) > 1L) {
    for (ii in seq_len(length(nonconstant) - 1L)) {
      i <- nonconstant[[ii]]
      for (jj in (ii + 1L):length(nonconstant)) {
        j <- nonconstant[[jj]]
        if (!is.na(same_weak_order(mat[, i], mat[, j]))) unite(i, j)
      }
    }
  }
  roots <- vapply(seq_len(p), find_root, integer(1))
  groups <- split(seq_len(p), roots)
  representatives <- vapply(groups, min, integer(1))
  representatives <- sort(representatives)
  detail <- do.call(rbind, lapply(seq_along(groups), function(g) {
    members <- groups[[g]]
    representative <- min(members)
    data.frame(
      class_id = g,
      representative_index = representative,
      representative = colnames(mat)[[representative]],
      member_index = members,
      member = colnames(mat)[members],
      direction = vapply(
        members,
        function(j) if (j == representative) "representative" else same_weak_order(mat[, representative], mat[, j]),
        character(1)
      ),
      class_size = length(members),
      stringsAsFactors = FALSE
    )
  }))
  list(keep = representatives, detail = detail)
}

make_folds <- function(n, seed = 1L) {
  set.seed(seed)
  split(sample(n, n, replace = FALSE), as.factor(1:5))
}

build_cache <- function(dataset_name, spec) {
  message("Building cache for ", dataset_name, ", specification ", spec)
  data <- get_data_intersection(datasets[[dataset_name]], spec)
  y <- data[[1]]
  treatment <- data[[2]]
  x <- as.matrix(data[[3]])
  if (spec <= 2L) {
    colnames(x) <- spec_names[[as.character(spec)]]
  } else {
    colnames(x) <- paste0("z", seq_len(ncol(x)))
  }
  full_dictionary <- make_dictionary(treatment, x)
  classes <- if (spec <= 2L) {
    weak_order_classes(rbind(full_dictionary,
      make_dictionary(rep(0, nrow(x)), x),
      make_dictionary(rep(1, nrow(x)), x)))
  } else {
    list(
      keep = seq_len(ncol(full_dictionary)),
      detail = data.frame(
        class_id = seq_len(ncol(full_dictionary)),
        representative_index = seq_len(ncol(full_dictionary)),
        representative = colnames(full_dictionary),
        member_index = seq_len(ncol(full_dictionary)),
        member = colnames(full_dictionary),
        direction = "representative",
        class_size = 1L,
        stringsAsFactors = FALSE
      )
    )
  }

  dict <- b2
  p <- ncol(full_dictionary)
  p0 <- if (p > 60L) ceiling(p / 40) else ceiling(p / 4)
  c_values <- c(5 / 4, 1, 3 / 4, 1 / 2) * 0.5
  if (args$tuning == "theoretical") {
    # The paper sets c1 = 1. The released code divides by two because its
    # objective writes the penalty as 2 r ||D rho||_1.
    c_rr <- 0.5
  } else {
    set.seed(1)
    c_rr <- get_cv_rr(y, treatment, x, p0, 0, 0.2, 10, dict, 1, 1, c_values)
  }
  folds <- make_folds(nrow(x), 1L)
  set.seed(1)
  invisible(sample(nrow(x), nrow(x), replace = FALSE))
  published_rng_state <- .Random.seed

  rho <- vector("list", length(folds))
  for (fold in seq_along(folds)) {
    train <- setdiff(seq_len(nrow(x)), folds[[fold]])
    rho[[fold]] <- RMD_stable(
      y[train], treatment[train], x[train, , drop = FALSE],
      p0, 0, 0.2, 10, dict, 1, 1, c_rr
    )
  }
  list(
    dataset = dataset_name,
    specification = spec,
    y = y,
    treatment = treatment,
    x = x,
    full_dictionary = full_dictionary,
    classes = classes,
    folds = folds,
    rho = rho,
    c_rr = c_rr,
    p0 = p0,
    published_rng_state = published_rng_state
  )
}

load_cache <- function(dataset_name, spec) {
  path <- file.path(
    cache_dir,
    paste0(tolower(dataset_name), "_spec", spec, "_", args$tuning, ".rds")
  )
  if (args$rebuild_cache || !file.exists(path)) {
    cache <- build_cache(dataset_name, spec)
    saveRDS(cache, path, compress = "xz")
  } else {
    cache <- readRDS(path)
  }
  cache
}

forest_fit <- function(x, y, ntree, mtry = NULL) {
  forest_args <- list(
    x = x,
    y = y,
    clas_nodesize = 1,
    reg_nodesize = 5,
    ntree = ntree,
    na.action = na.omit,
    replace = TRUE
  )
  if (!is.null(mtry)) forest_args$mtry <- mtry
  do.call(randomForest, forest_args)
}

run_variant <- function(cache, variant, seed = NULL, published_rng = FALSE) {
  p_full <- ncol(cache$full_dictionary)
  keep <- switch(
    variant,
    published_full = seq_len(p_full),
    quotient_default = cache$classes$keep,
    quotient_fixed_mtry = cache$classes$keep,
    stop("Unknown variant")
  )
  fixed_mtry <- if (variant == "quotient_fixed_mtry") max(floor(p_full / 3), 1L) else NULL
  mtry_used <- if (is.null(fixed_mtry)) max(floor(length(keep) / 3), 1L) else fixed_mtry
  if (published_rng) {
    assign(".Random.seed", cache$published_rng_state, envir = .GlobalEnv)
  } else {
    set.seed(seed)
  }

  n <- nrow(cache$x)
  psi <- gamma0 <- gamma1 <- gamma_observed <- rep(NA_real_, n)
  for (fold in seq_along(cache$folds)) {
    test <- cache$folds[[fold]]
    train <- setdiff(seq_len(n), test)
    b_train <- cache$full_dictionary[train, keep, drop = FALSE]
    forest <- forest_fit(b_train, cache$y[train], args$ntree, fixed_mtry)

    b_observed <- cache$full_dictionary[test, , drop = FALSE]
    b_zero <- make_dictionary(rep(0, length(test)), cache$x[test, , drop = FALSE])
    b_one <- make_dictionary(rep(1, length(test)), cache$x[test, , drop = FALSE])
    gamma_observed[test] <- predict(forest, newdata = b_observed[, keep, drop = FALSE])
    gamma0[test] <- predict(forest, newdata = b_zero[, keep, drop = FALSE])
    gamma1[test] <- predict(forest, newdata = b_one[, keep, drop = FALSE])
    alpha <- as.vector(b_observed %*% cache$rho[[fold]])
    psi[test] <- cache$treatment[test] * gamma0[test] +
      alpha * (cache$y[test] - gamma_observed[test])
  }

  treatment_mean <- mean(cache$treatment)
  ty_mean <- mean(cache$treatment * cache$y)
  psi_mean <- mean(psi)
  atet <- (ty_mean - psi_mean) / treatment_mean
  centered_psi <- psi - psi_mean
  variance_matrix <- matrix(c(
    mean(centered_psi^2),
    mean(centered_psi * cache$treatment * cache$y),
    mean(centered_psi * cache$treatment),
    mean(centered_psi * cache$treatment * cache$y),
    var(cache$treatment * cache$y),
    cov(cache$treatment * cache$y, cache$treatment),
    mean(centered_psi * cache$treatment),
    cov(cache$treatment * cache$y, cache$treatment),
    var(cache$treatment)
  ), nrow = 3)
  gradient <- matrix(c(
    -1 / treatment_mean,
    1 / treatment_mean,
    (psi_mean - ty_mean) / treatment_mean^2
  ), nrow = 1)
  se <- sqrt(as.numeric(gradient %*% variance_matrix %*% t(gradient)) / n)

  list(
    result = data.frame(
      dataset = cache$dataset,
      specification = cache$specification,
      variant = variant,
      seed = if (published_rng) 1L else seed,
      rng_scheme = if (published_rng) "published" else "fixed_folds_seeded_forests",
      observations = n,
      treated = sum(cache$treatment == 1),
      controls = sum(cache$treatment == 0),
      dictionary_columns = length(keep),
      removed_columns = p_full - length(keep),
      mtry = mtry_used,
      ntree = args$ntree,
      c_rr = cache$c_rr,
      atet = atet,
      se = se,
      ci_low = atet - 1.96 * se,
      ci_high = atet + 1.96 * se,
      stringsAsFactors = FALSE
    ),
    gamma0 = gamma0,
    gamma1 = gamma1,
    contrast = gamma1 - gamma0,
    psi = psi
  )
}

append_rows <- function(path, rows) {
  out <- do.call(rbind, rows)
  write.csv(out, path, row.names = FALSE)
  out
}

selected_cells <- expand.grid(
  dataset = args$datasets,
  specification = args$specs,
  stringsAsFactors = FALSE
)

all_classes <- list()
baseline_rows <- list()
paired_rows <- list()
diagnostic_rows <- list()
seed_noise_rows <- list()
row_baseline <- row_paired <- row_diagnostic <- row_seed_noise <- row_class <- 0L

for (cell in seq_len(nrow(selected_cells))) {
  dataset_name <- selected_cells$dataset[[cell]]
  spec <- selected_cells$specification[[cell]]
  cache <- load_cache(dataset_name, spec)

  class_detail <- cache$classes$detail
  class_detail$dataset <- dataset_name
  class_detail$specification <- spec
  row_class <- row_class + 1L
  all_classes[[row_class]] <- class_detail[, c(
    "dataset", "specification", "class_id", "representative_index",
    "representative", "member_index", "member", "direction", "class_size"
  )]

  if (args$phase %in% c("baseline", "all")) {
    message("Published reproduction: ", dataset_name, ", specification ", spec)
    fit <- run_variant(cache, "published_full", published_rng = TRUE)
    row_baseline <- row_baseline + 1L
    baseline_rows[[row_baseline]] <- fit$result
  }

  if (args$phase %in% c("paired", "seed_noise", "all") && spec <= 2L) {
    full_seed_reference <- NULL
    full_seed_reference_value <- NA_integer_
    for (seed in args$seeds) {
      variants_to_run <- if (args$phase == "seed_noise") {
        "published_full"
      } else {
        c("published_full", "quotient_default", "quotient_fixed_mtry")
      }
      message(
        if (args$phase == "seed_noise") "Seed-noise run: " else "Paired run: ",
        dataset_name, ", specification ", spec, ", seed ", seed
      )
      fits <- lapply(variants_to_run, function(variant) run_variant(cache, variant, seed = seed))
      names(fits) <- variants_to_run
      for (variant in names(fits)) {
        row_paired <- row_paired + 1L
        paired_rows[[row_paired]] <- fits[[variant]]$result
      }
      reference <- fits$published_full
      if (is.null(full_seed_reference)) {
        full_seed_reference <- reference
        full_seed_reference_value <- seed
      } else {
        for (subset_name in c("all", "treated")) {
          use <- if (subset_name == "all") rep(TRUE, length(cache$treatment)) else cache$treatment == 1
          ref_contrast <- full_seed_reference$contrast[use]
          alt_contrast <- reference$contrast[use]
          row_seed_noise <- row_seed_noise + 1L
          seed_noise_rows[[row_seed_noise]] <- data.frame(
            dataset = dataset_name,
            specification = spec,
            reference_seed = full_seed_reference_value,
            comparison_seed = seed,
            subset = subset_name,
            observations = sum(use),
            atet_difference = reference$result$atet - full_seed_reference$result$atet,
            se_difference = reference$result$se - full_seed_reference$result$se,
            gamma0_rmse = sqrt(mean((reference$gamma0[use] - full_seed_reference$gamma0[use])^2)),
            gamma1_rmse = sqrt(mean((reference$gamma1[use] - full_seed_reference$gamma1[use])^2)),
            contrast_rmse = sqrt(mean((alt_contrast - ref_contrast)^2)),
            contrast_correlation = cor(alt_contrast, ref_contrast),
            contrast_sign_disagreement = mean(sign(alt_contrast) != sign(ref_contrast)),
            stringsAsFactors = FALSE
          )
        }
      }
      for (variant in intersect(c("quotient_default", "quotient_fixed_mtry"), names(fits))) {
        comparison <- fits[[variant]]
        for (subset_name in c("all", "treated")) {
          use <- if (subset_name == "all") rep(TRUE, length(cache$treatment)) else cache$treatment == 1
          ref_contrast <- reference$contrast[use]
          alt_contrast <- comparison$contrast[use]
          row_diagnostic <- row_diagnostic + 1L
          diagnostic_rows[[row_diagnostic]] <- data.frame(
            dataset = dataset_name,
            specification = spec,
            seed = seed,
            variant = variant,
            subset = subset_name,
            observations = sum(use),
            atet_difference = comparison$result$atet - reference$result$atet,
            se_difference = comparison$result$se - reference$result$se,
            gamma0_rmse = sqrt(mean((comparison$gamma0[use] - reference$gamma0[use])^2)),
            gamma1_rmse = sqrt(mean((comparison$gamma1[use] - reference$gamma1[use])^2)),
            contrast_rmse = sqrt(mean((alt_contrast - ref_contrast)^2)),
            contrast_correlation = cor(alt_contrast, ref_contrast),
            contrast_sign_disagreement = mean(sign(alt_contrast) != sign(ref_contrast)),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

classes_out <- do.call(rbind, all_classes)
write.csv(classes_out, file.path(args$output_dir, "equivalence_classes.csv"), row.names = FALSE)

if (length(baseline_rows)) {
  baseline_out <- append_rows(
    file.path(args$output_dir, "published_reproduction.csv"), baseline_rows
  )
  baseline_out <- merge(baseline_out, published, by = c("dataset", "specification"), all.x = TRUE)
  baseline_out$atet_difference <- baseline_out$atet - baseline_out$published_atet
  baseline_out$se_difference <- baseline_out$se - baseline_out$published_se
  baseline_out$rounded_atet_match <- round(baseline_out$atet, 2) == baseline_out$published_atet
  baseline_out$rounded_se_match <- round(baseline_out$se, 2) == baseline_out$published_se
  write.csv(
    baseline_out,
    file.path(args$output_dir, "published_reproduction_validation.csv"),
    row.names = FALSE
  )
}

if (length(paired_rows)) {
  paired_out <- append_rows(file.path(args$output_dir, "paired_seed_results.csv"), paired_rows)
  reference <- paired_out[paired_out$variant == "published_full", c(
    "dataset", "specification", "seed", "atet", "se"
  )]
  names(reference)[4:5] <- c("reference_atet", "reference_se")
  paired_with_reference <- merge(
    paired_out, reference,
    by = c("dataset", "specification", "seed"), all.x = TRUE
  )
  paired_with_reference$atet_difference <- paired_with_reference$atet - paired_with_reference$reference_atet
  paired_with_reference$se_difference <- paired_with_reference$se - paired_with_reference$reference_se
  summary_out <- do.call(rbind, lapply(
    split(paired_with_reference, interaction(
      paired_with_reference$dataset,
      paired_with_reference$specification,
      paired_with_reference$variant,
      drop = TRUE
    )),
    function(z) data.frame(
      dataset = z$dataset[[1]],
      specification = z$specification[[1]],
      variant = z$variant[[1]],
      seeds = nrow(z),
      mean_atet = mean(z$atet),
      sd_atet = sd(z$atet),
      min_atet = min(z$atet),
      max_atet = max(z$atet),
      mean_se = mean(z$se),
      mean_atet_difference = mean(z$atet_difference),
      sd_atet_difference = sd(z$atet_difference),
      max_abs_atet_difference = max(abs(z$atet_difference)),
      aggregate_sign_changes = sum(sign(z$atet) != sign(z$reference_atet)),
      inference_changes = sum(
        (z$ci_low > 0 | z$ci_high < 0) !=
          ((z$reference_atet - 1.96 * z$reference_se) > 0 |
             (z$reference_atet + 1.96 * z$reference_se) < 0)
      ),
      stringsAsFactors = FALSE
    )
  ))
  write.csv(summary_out, file.path(args$output_dir, "paired_seed_summary.csv"), row.names = FALSE)

  if (length(diagnostic_rows)) {
    diagnostics_out <- append_rows(
      file.path(args$output_dir, "paired_prediction_diagnostics.csv"), diagnostic_rows
    )
    diagnostic_summary <- do.call(rbind, lapply(
      split(diagnostics_out, interaction(
        diagnostics_out$dataset,
        diagnostics_out$specification,
        diagnostics_out$variant,
        diagnostics_out$subset,
        drop = TRUE
      )),
      function(z) data.frame(
        dataset = z$dataset[[1]],
        specification = z$specification[[1]],
        variant = z$variant[[1]],
        subset = z$subset[[1]],
        seeds = nrow(z),
        mean_contrast_sign_disagreement = mean(z$contrast_sign_disagreement),
        max_contrast_sign_disagreement = max(z$contrast_sign_disagreement),
        mean_contrast_rmse = mean(z$contrast_rmse),
        mean_contrast_correlation = mean(z$contrast_correlation),
        mean_gamma0_rmse = mean(z$gamma0_rmse),
        mean_gamma1_rmse = mean(z$gamma1_rmse),
        stringsAsFactors = FALSE
      )
    ))
    write.csv(
      diagnostic_summary,
      file.path(args$output_dir, "paired_prediction_diagnostic_summary.csv"),
      row.names = FALSE
    )
  }

  if (length(seed_noise_rows)) {
    seed_noise_out <- append_rows(
      file.path(args$output_dir, "seed_noise_diagnostics.csv"), seed_noise_rows
    )
    seed_noise_summary <- do.call(rbind, lapply(
      split(seed_noise_out, interaction(
        seed_noise_out$dataset,
        seed_noise_out$specification,
        seed_noise_out$subset,
        drop = TRUE
      )),
      function(z) data.frame(
        dataset = z$dataset[[1]],
        specification = z$specification[[1]],
        subset = z$subset[[1]],
        seed_comparisons = nrow(z),
        mean_abs_atet_difference = mean(abs(z$atet_difference)),
        max_abs_atet_difference = max(abs(z$atet_difference)),
        mean_contrast_sign_disagreement = mean(z$contrast_sign_disagreement),
        max_contrast_sign_disagreement = max(z$contrast_sign_disagreement),
        mean_contrast_rmse = mean(z$contrast_rmse),
        mean_contrast_correlation = mean(z$contrast_correlation),
        mean_gamma0_rmse = mean(z$gamma0_rmse),
        mean_gamma1_rmse = mean(z$gamma1_rmse),
        stringsAsFactors = FALSE
      )
    ))
    write.csv(
      seed_noise_summary,
      file.path(args$output_dir, "seed_noise_diagnostic_summary.csv"),
      row.names = FALSE
    )
  }
}

writeLines(capture.output(sessionInfo()), file.path(args$output_dir, "session_info.txt"))
writeLines(capture.output(str(args)), file.path(args$output_dir, "run_configuration.txt"))
message("Completed. Outputs: ", normalizePath(args$output_dir))
