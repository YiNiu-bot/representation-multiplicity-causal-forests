#!/usr/bin/env Rscript

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_flag)) {
  normalizePath(sub("^--file=", "", script_flag[[1L]]))
} else {
  normalizePath("replication/code/r/run_fixed_budget_study.R")
}
repository_dir <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
replication_dir <- file.path(repository_dir, "replication")
rlib <- Sys.getenv("RI_RLIB", "")
if (nzchar(rlib) && dir.exists(rlib)) .libPaths(c(rlib, .libPaths()))

suppressPackageStartupMessages({
  library(grf)
  library(haven)
})
expected_grf_version <- "2.4.0"
allow_unpinned_grf <- Sys.getenv("RI_ALLOW_UNPINNED_GRF", "0") == "1"
if (as.character(packageVersion("grf")) != expected_grf_version &&
    !allow_unpinned_grf) {
  stop(
    sprintf(
      "Expected grf %s, found %s. Set RI_ALLOW_UNPINNED_GRF=1 only for a non-replication run.",
      expected_grf_version,
      packageVersion("grf")
    ),
    call. = FALSE
  )
}

source_dir <- Sys.getenv(
  "RI_SOURCE_DIR",
  file.path(replication_dir, "data", "kenya")
)
output_dir <- Sys.getenv(
  "RI_OUTPUT_DIR",
  file.path(repository_dir, "rerun_outputs", "fixed_budget")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_data <- c(
  "GE_HH-Analysis_AllHHs.dta",
  "GE_HH-Survey-BL_Analysis_AllHHs.dta"
)
missing_data <- required_data[!file.exists(file.path(source_dir, required_data))]
if (length(missing_data)) {
  stop(
    sprintf(
      "Missing required Kenya data under %s: %s. See replication/data/README.md.",
      source_dir,
      paste(missing_data, collapse = ", ")
    ),
    call. = FALSE
  )
}
expected_data_md5 <- c(
  "GE_HH-Analysis_AllHHs.dta" = "6841329ee31e3fee6413651bf8a1c5ff",
  "GE_HH-Survey-BL_Analysis_AllHHs.dta" = "35b0b1384ea3cf82e38bd77c73442f3c"
)
observed_data_md5 <- unname(tools::md5sum(file.path(source_dir, required_data)))
names(observed_data_md5) <- required_data
bad_md5 <- names(expected_data_md5)[
  tolower(observed_data_md5[names(expected_data_md5)]) != expected_data_md5
]
allow_data_mismatch <- Sys.getenv("RI_ALLOW_DATA_MISMATCH", "0") == "1"
if (length(bad_md5) && !allow_data_mismatch) {
  stop(
    sprintf(
      "Input checksum mismatch for: %s. Set RI_ALLOW_DATA_MISMATCH=1 only for a non-replication run.",
      paste(bad_md5, collapse = ", ")
    ),
    call. = FALSE
  )
}

threads <- as.integer(Sys.getenv("RI_THREADS", "15"))
real_trees <- as.integer(Sys.getenv("RI_REAL_TREES", "2000"))
simulation_trees <- as.integer(Sys.getenv("RI_SIM_TREES", "1500"))
base_seed <- 20260806L

weighted_mean <- function(x, w) sum(w * x) / sum(w)

cluster_mean <- function(x, w, cluster) {
  estimate <- weighted_mean(x, w)
  influence <- w * (x - estimate) / sum(w)
  cluster_influence <- rowsum(influence, cluster, reorder = FALSE)[, 1]
  g <- length(cluster_influence)
  list(
    estimate = estimate,
    se = sqrt(g / (g - 1) * sum(cluster_influence^2))
  )
}

top_share <- function(score, share, id = seq_along(score)) {
  stopifnot(length(score) == length(id))
  k <- floor(length(score) * share)
  selected <- rep(FALSE, length(score))
  if (k > 0L) selected[order(-score, id)[seq_len(k)]] <- TRUE
  selected
}

comparison_metrics <- function(tau, selected, reference_tau, reference_selected) {
  intersection <- sum(selected & reference_selected)
  union <- sum(selected | reference_selected)
  c(
    cate_correlation = cor(tau, reference_tau),
    cate_rank_correlation = cor(tau, reference_tau, method = "spearman"),
    top_set_overlap = intersection / union,
    assignment_switching = mean(selected != reference_selected)
  )
}

median_impute <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

zero_impute <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x
}

clone_feature <- function(X, feature, multiplicity) {
  stopifnot(feature %in% colnames(X), multiplicity >= 1L)
  if (multiplicity == 1L) return(X)
  copies <- matrix(
    X[, feature],
    nrow = nrow(X),
    ncol = multiplicity - 1L
  )
  colnames(copies) <- sprintf("%s__copy_%03d", feature, seq_len(ncol(copies)))
  cbind(X, copies)
}

balanced_clone <- function(X, multiplicity) {
  stopifnot(multiplicity >= 1L)
  if (multiplicity == 1L) return(X)
  out <- X
  original_names <- colnames(X)
  for (feature in original_names) {
    copies <- matrix(
      X[, feature],
      nrow = nrow(X),
      ncol = multiplicity - 1L
    )
    colnames(copies) <- sprintf("%s__copy_%03d", feature, seq_len(ncol(copies)))
    out <- cbind(out, copies)
  }
  out
}

add_monotone_encodings <- function(X, features) {
  out <- X
  for (feature in features) {
    x <- X[, feature]
    shifted <- x - min(x) + 1
    out <- cbind(
      out,
      setNames(data.frame(log1p(shifted)), paste0(feature, "__log")),
      setNames(data.frame(rank(x, ties.method = "average")), paste0(feature, "__rank"))
    )
  }
  as.matrix(out)
}

deduplicate_exact_columns <- function(X) {
  keep <- rep(TRUE, ncol(X))
  if (ncol(X) <= 1L) return(X)
  for (j in 2:ncol(X)) {
    for (k in seq_len(j - 1L)) {
      if (keep[k] && identical(unname(X[, j]), unname(X[, k]))) {
        keep[j] <- FALSE
        break
      }
    }
  }
  X[, keep, drop = FALSE]
}

split_signature <- function(x) {
  missing <- is.na(x)
  observed <- x[!missing]
  levels <- sort(unique(observed))
  forward <- rep(NA_integer_, length(x))
  forward[!missing] <- match(observed, levels)
  reverse <- forward
  reverse[!missing] <- length(levels) + 1L - forward[!missing]
  missing_code <- length(levels) + 1L
  forward[missing] <- missing_code
  reverse[missing] <- missing_code
  forward_key <- paste(forward, collapse = ",")
  reverse_key <- paste(reverse, collapse = ",")
  paste0(paste(as.integer(missing), collapse = ""), "|", min(forward_key, reverse_key))
}

rank_quotient <- function(X) {
  signatures <- vapply(seq_len(ncol(X)), function(j) split_signature(X[, j]), character(1))
  X[, !duplicated(signatures), drop = FALSE]
}

cat("Preparing randomized cash-transfer data.\n")
analysis <- read_dta(
  file.path(source_dir, "GE_HH-Analysis_AllHHs.dta"),
  col_select = c(
    "hhid_key", "eligible", "baselined", "treat", "village_code",
    "hhweight_EL", "p2_consumption_pc_wins_PPP"
  )
)
baseline <- read_dta(
  file.path(source_dir, "GE_HH-Survey-BL_Analysis_AllHHs.dta"),
  col_select = c(
    "hhid_key", "hhsize1_BL", "female_BL", "haschildhh_BL",
    "has_child_sch_BL", "numchild3_BL", "numchild6_BL", "widowed_BL",
    "age60_BL", "emp_BL", "selfemp_BL", "any_livestock",
    "land_ownland_BL", "own_land_acres_BL", "h1_4_radiotv_BL",
    "num_meals_yest", "num_meals_yest_protein"
  )
)

idx <- match(analysis$hhid_key, baseline$hhid_key)
d <- cbind(analysis, baseline[idx, setdiff(names(baseline), "hhid_key")])
keep <- d$eligible == 1 & d$baselined == 1 & !is.na(d$treat) &
  !is.na(d$village_code) & !is.na(d$hhweight_EL) & d$hhweight_EL > 0 &
  !is.na(d$p2_consumption_pc_wins_PPP)
d <- d[keep, ]

X_real <- cbind(
  hh_size = median_impute(d$hhsize1_BL),
  female_head = zero_impute(d$female_BL),
  has_children = zero_impute(d$haschildhh_BL),
  child_in_school = zero_impute(d$has_child_sch_BL),
  child_under_3 = as.numeric(zero_impute(d$numchild3_BL) > 0),
  child_under_6 = as.numeric(zero_impute(d$numchild6_BL) > 0),
  widow = zero_impute(d$widowed_BL),
  elder = zero_impute(d$age60_BL),
  employed = zero_impute(d$emp_BL),
  self_employed = zero_impute(d$selfemp_BL),
  livestock = zero_impute(d$any_livestock),
  owns_land = zero_impute(d$land_ownland_BL),
  owns_quarter_acre = as.numeric(zero_impute(d$own_land_acres_BL) >= 0.25),
  owns_radio_tv = as.numeric(zero_impute(d$h1_4_radiotv_BL) > 0),
  meals = median_impute(d$num_meals_yest),
  protein_meals = median_impute(d$num_meals_yest_protein)
)
Y_real <- as.numeric(d$p2_consumption_pc_wins_PPP)
W_real <- as.numeric(d$treat)
weights_real <- as.numeric(d$hhweight_EL)
clusters_real <- as.integer(factor(d$village_code))
mtry_mean_real <- ncol(X_real)
budget_share_real <- 0.5

e_hat <- rep(weighted_mean(W_real, weights_real), length(W_real))
y_forest <- regression_forest(
  X_real,
  Y_real,
  sample.weights = weights_real,
  clusters = clusters_real,
  num.trees = real_trees,
  mtry = mtry_mean_real,
  seed = base_seed - 1L,
  num.threads = threads
)
y_hat <- as.numeric(predict(y_forest)$predictions)
gamma <- (W_real - e_hat) * (Y_real - y_hat) / (e_hat * (1 - e_hat))

real_specs <- list(
  list(name = "baseline", kind = "baseline", seed = base_seed),
  list(name = "baseline_seed_2", kind = "baseline", seed = base_seed + 1L),
  list(name = "baseline_seed_3", kind = "baseline", seed = base_seed + 2L),
  list(name = "baseline_seed_4", kind = "baseline", seed = base_seed + 3L)
)
for (multiplicity in c(2L, 4L, 8L, 16L, 32L, 64L, 128L)) {
  real_specs[[length(real_specs) + 1L]] <- list(
    name = sprintf("clone_hh_size_%03d", multiplicity),
    kind = "clone",
    feature = "hh_size",
    multiplicity = multiplicity,
    seed = base_seed
  )
}
for (feature in c("female_head", "employed", "owns_land", "owns_radio_tv", "protein_meals")) {
  real_specs[[length(real_specs) + 1L]] <- list(
    name = sprintf("clone_%s_128", feature),
    kind = "clone",
    feature = feature,
    multiplicity = 128L,
    seed = base_seed
  )
}
real_specs[[length(real_specs) + 1L]] <- list(
  name = "routine_raw_log_rank",
  kind = "monotone",
  seed = base_seed
)
real_specs[[length(real_specs) + 1L]] <- list(
  name = "balanced_clone_008",
  kind = "balanced",
  multiplicity = 8L,
  seed = base_seed
)

real_fits <- vector("list", length(real_specs))
for (i in seq_along(real_specs)) {
  spec <- real_specs[[i]]
  Xm <- switch(
    spec$kind,
    baseline = X_real,
    clone = clone_feature(X_real, spec$feature, spec$multiplicity),
    monotone = add_monotone_encodings(
      X_real,
      c("hh_size", "meals", "protein_meals")
    ),
    balanced = balanced_clone(X_real, spec$multiplicity)
  )
  stopifnot(mtry_mean_real <= ncol(Xm))
  cat(sprintf(
    "Real-data fit %02d/%02d: %s (p=%d, mtry=%d).\n",
    i, length(real_specs), spec$name, ncol(Xm), mtry_mean_real
  ))
  fit <- causal_forest(
    Xm,
    Y_real,
    W_real,
    Y.hat = y_hat,
    W.hat = e_hat,
    sample.weights = weights_real,
    clusters = clusters_real,
    num.trees = real_trees,
    mtry = mtry_mean_real,
    seed = spec$seed,
    num.threads = threads
  )
  tau <- as.numeric(predict(fit)$predictions)
  selected <- top_share(tau, budget_share_real, d$hhid_key)
  value <- cluster_mean(as.numeric(selected) * gamma, weights_real, clusters_real)
  real_fits[[i]] <- list(
    name = spec$name,
    kind = spec$kind,
    seed = spec$seed,
    p = ncol(Xm),
    mtry = mtry_mean_real,
    tau = tau,
    selected = selected,
    value = value$estimate,
    value_se = value$se
  )
  gc(verbose = FALSE)
}

reference_real <- real_fits[[1L]]
real_rows <- vector("list", length(real_fits))
for (i in seq_along(real_fits)) {
  fit <- real_fits[[i]]
  metrics <- comparison_metrics(
    fit$tau,
    fit$selected,
    reference_real$tau,
    reference_real$selected
  )
  paired <- cluster_mean(
    as.numeric(fit$selected - reference_real$selected) * gamma,
    weights_real,
    clusters_real
  )
  real_rows[[i]] <- data.frame(
    scenario = fit$name,
    kind = fit$kind,
    seed = fit$seed,
    n = nrow(X_real),
    clusters = length(unique(clusters_real)),
    p = fit$p,
    mtry = fit$mtry,
    budget_share = budget_share_real,
    cate_correlation = metrics[["cate_correlation"]],
    cate_rank_correlation = metrics[["cate_rank_correlation"]],
    top_set_overlap = metrics[["top_set_overlap"]],
    assignment_switching = metrics[["assignment_switching"]],
    policy_value = fit$value,
    policy_value_se = fit$value_se,
    policy_value_difference = paired$estimate,
    policy_value_difference_se = paired$se,
    stringsAsFactors = FALSE
  )
}
real_results <- do.call(rbind, real_rows)
write.csv(real_results, file.path(output_dir, "real_data_results.csv"), row.names = FALSE)
save_objects <- identical(Sys.getenv("RI_SAVE_INDIVIDUAL_OBJECTS", "0"), "1")
if (save_objects) {
  saveRDS(real_fits, file.path(output_dir, "real_data_fits.rds"))
}

cat("Running the adaptive-forest numerical illustration.\n")
set.seed(base_seed)
n_train <- 5000L
n_test <- 30000L
p0 <- 20L
mtry_mean_sim <- 6L
budget_share_sim <- 0.2
X_train <- matrix(runif(n_train * p0), nrow = n_train, ncol = p0)
X_test <- matrix(runif(n_test * p0), nrow = n_test, ncol = p0)
colnames(X_train) <- colnames(X_test) <- sprintf("x%02d", seq_len(p0))

response_functions <- function(X) {
  group_a <- X[, "x01"] > 0.70 & X[, "x03"] > 0.35
  group_b <- X[, "x02"] > 0.70 & X[, "x03"] <= 0.35
  tau <- 0.05 + 1.40 * group_a + 1.05 * group_b + 0.25 * (X[, "x04"] - 0.5)
  mu <- 1 + X[, "x05"] + 0.5 * X[, "x06"]^2 - 0.5 * X[, "x07"]
  list(mu = mu, tau = tau, group_a = group_a, group_b = group_b)
}

train_signal <- response_functions(X_train)
test_signal <- response_functions(X_test)
W_train <- rbinom(n_train, 1, 0.5)
Y_train <- train_signal$mu + W_train * train_signal$tau + rnorm(n_train)
Y_hat_train <- train_signal$mu + 0.5 * train_signal$tau
W_hat_train <- rep(0.5, n_train)

simulation_specs <- list(
  list(name = "baseline", kind = "baseline"),
  list(name = "clone_effect_a_064", kind = "clone", feature = "x01", multiplicity = 64L),
  list(name = "clone_effect_b_064", kind = "clone", feature = "x02", multiplicity = 64L),
  list(name = "clone_irrelevant_064", kind = "clone", feature = "x20", multiplicity = 64L)
)
simulation_seeds <- base_seed + 100L + 0:4
simulation_rows <- list()
simulation_predictions <- list()

optimal_selected <- top_share(test_signal$tau, budget_share_sim, seq_len(n_test))
optimal_value <- mean(optimal_selected * test_signal$tau)

for (spec in simulation_specs) {
  Xtr <- switch(
    spec$kind,
    baseline = X_train,
    clone = clone_feature(X_train, spec$feature, spec$multiplicity),
    monotone = add_monotone_encodings(X_train, spec$feature)
  )
  Xte <- switch(
    spec$kind,
    baseline = X_test,
    clone = clone_feature(X_test, spec$feature, spec$multiplicity),
    monotone = add_monotone_encodings(X_test, spec$feature)
  )
  for (seed in simulation_seeds) {
    cat(sprintf(
      "Simulation fit: %s, seed=%d (p=%d, mtry=%d).\n",
      spec$name, seed, ncol(Xtr), mtry_mean_sim
    ))
    fit <- causal_forest(
      Xtr,
      Y_train,
      W_train,
      Y.hat = Y_hat_train,
      W.hat = W_hat_train,
      num.trees = simulation_trees,
      mtry = mtry_mean_sim,
      seed = seed,
      num.threads = threads
    )
    tau_hat <- as.numeric(predict(fit, Xte)$predictions)
    selected <- top_share(tau_hat, budget_share_sim, seq_len(n_test))
    value <- mean(selected * test_signal$tau)
    key <- paste(spec$name, seed, sep = "__")
    simulation_predictions[[key]] <- list(tau = tau_hat, selected = selected)
    simulation_rows[[length(simulation_rows) + 1L]] <- data.frame(
      scenario = spec$name,
      seed = seed,
      n_train = n_train,
      n_test = n_test,
      p = ncol(Xtr),
      mtry = mtry_mean_sim,
      budget_share = budget_share_sim,
      cate_rmse = sqrt(mean((tau_hat - test_signal$tau)^2)),
      cate_correlation = cor(tau_hat, test_signal$tau),
      policy_value = value,
      policy_regret = optimal_value - value,
      share_group_a_treated = mean(test_signal$group_a[selected]),
      share_group_b_treated = mean(test_signal$group_b[selected]),
      stringsAsFactors = FALSE
    )
    gc(verbose = FALSE)
  }
}

simulation_results <- do.call(rbind, simulation_rows)
baseline_predictions <- simulation_predictions[[paste("baseline", simulation_seeds[1], sep = "__")]]

# The exact-clone quotient removes repeated columns before calling the unchanged
# learner. Its input matrix is exactly the canonical representation.
X_clone_a <- clone_feature(X_train, "x01", 64L)
X_clone_a_q <- deduplicate_exact_columns(X_clone_a)
stopifnot(identical(unname(X_clone_a_q), unname(X_train)))
X_routine_a <- add_monotone_encodings(X_train, "x01")
X_base_rank_q <- rank_quotient(X_train)
X_clone_rank_q <- rank_quotient(X_clone_a)
X_routine_rank_q <- rank_quotient(X_routine_a)
stopifnot(
  identical(unname(X_clone_rank_q), unname(X_base_rank_q)),
  identical(unname(X_routine_rank_q), unname(X_base_rank_q))
)
quotient_fit <- causal_forest(
  X_clone_rank_q,
  Y_train,
  W_train,
  Y.hat = Y_hat_train,
  W.hat = W_hat_train,
  num.trees = simulation_trees,
  mtry = mtry_mean_sim,
  seed = simulation_seeds[1],
  num.threads = threads
)
quotient_tau <- as.numeric(predict(quotient_fit, X_test)$predictions)
quotient_difference <- max(abs(quotient_tau - baseline_predictions$tau))
routine_quotient_fit <- causal_forest(
  X_routine_rank_q,
  Y_train,
  W_train,
  Y.hat = Y_hat_train,
  W.hat = W_hat_train,
  num.trees = simulation_trees,
  mtry = mtry_mean_sim,
  seed = simulation_seeds[1],
  num.threads = threads
)
routine_quotient_tau <- as.numeric(predict(routine_quotient_fit, X_test)$predictions)
routine_quotient_difference <- max(abs(routine_quotient_tau - baseline_predictions$tau))

X_real_base_q <- rank_quotient(X_real)
X_real_clone_q <- rank_quotient(clone_feature(X_real, "hh_size", 128L))
X_real_routine_q <- rank_quotient(add_monotone_encodings(
  X_real,
  c("hh_size", "meals", "protein_meals")
))
stopifnot(
  identical(unname(X_real_clone_q), unname(X_real_base_q)),
  identical(unname(X_real_routine_q), unname(X_real_base_q))
)

write.csv(
  simulation_results,
  file.path(output_dir, "simulation_results.csv"),
  row.names = FALSE
)
simulation_summary <- aggregate(
  cbind(
    cate_rmse,
    cate_correlation,
    policy_value,
    policy_regret,
    share_group_a_treated,
    share_group_b_treated
  ) ~ scenario,
  data = simulation_results,
  FUN = mean
)
write.csv(
  simulation_summary,
  file.path(output_dir, "simulation_summary.csv"),
  row.names = FALSE
)
write.csv(
  rbind(
    data.frame(
      sample = "simulation",
      correction = "rank_quotient_exact_clones",
      input_columns_before = ncol(X_clone_a),
      input_columns_after = ncol(X_clone_rank_q),
      canonical_columns = ncol(X_base_rank_q),
      input_equals_canonical = identical(unname(X_clone_rank_q), unname(X_base_rank_q)),
      maximum_prediction_difference = quotient_difference
    ),
    data.frame(
      sample = "simulation",
      correction = "rank_quotient_monotone_encodings",
      input_columns_before = ncol(X_routine_a),
      input_columns_after = ncol(X_routine_rank_q),
      canonical_columns = ncol(X_base_rank_q),
      input_equals_canonical = identical(unname(X_routine_rank_q), unname(X_base_rank_q)),
      maximum_prediction_difference = routine_quotient_difference
    ),
    data.frame(
      sample = "real_data",
      correction = "rank_quotient_exact_clones",
      input_columns_before = ncol(clone_feature(X_real, "hh_size", 128L)),
      input_columns_after = ncol(X_real_clone_q),
      canonical_columns = ncol(X_real_base_q),
      input_equals_canonical = identical(unname(X_real_clone_q), unname(X_real_base_q)),
      maximum_prediction_difference = NA_real_
    ),
    data.frame(
      sample = "real_data",
      correction = "rank_quotient_monotone_encodings",
      input_columns_before = ncol(add_monotone_encodings(
        X_real,
        c("hh_size", "meals", "protein_meals")
      )),
      input_columns_after = ncol(X_real_routine_q),
      canonical_columns = ncol(X_real_base_q),
      input_equals_canonical = identical(unname(X_real_routine_q), unname(X_real_base_q)),
      maximum_prediction_difference = NA_real_
    )
  ),
  file.path(output_dir, "quotient_check.csv"),
  row.names = FALSE
)

candidate_probability <- function(P, m, q) {
  if (P - m < q) return(1)
  1 - exp(lchoose(P - m, q) - lchoose(P, q))
}

capped_poisson_pmf <- function(P, lambda) {
  stopifnot(P >= 1L, lambda > 0)
  if (P == 1L) return(1)
  probabilities <- numeric(P)
  probabilities[1L] <- ppois(1, lambda)
  if (P > 2L) {
    interior <- 2L:(P - 1L)
    probabilities[interior] <- dpois(interior, lambda)
  }
  probabilities[P] <- ppois(P - 1L, lambda, lower.tail = FALSE)
  stopifnot(abs(sum(probabilities) - 1) < 1e-12)
  probabilities
}

candidate_probability_random <- function(P, m, lambda) {
  q <- seq_len(P)
  probabilities <- capped_poisson_pmf(P, lambda)
  sum(probabilities * vapply(q, function(qi) {
    candidate_probability(P, m, qi)
  }, numeric(1)))
}

opportunity <- do.call(rbind, lapply(
  c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L),
  function(r) {
    P <- 15L + r
    probabilities <- capped_poisson_pmf(P, mtry_mean_real)
    data.frame(
      multiplicity = r,
      raw_columns = P,
      mtry_poisson_mean = mtry_mean_real,
      expected_candidate_count = sum(seq_len(P) * probabilities),
      cloned_family_opportunity = candidate_probability_random(P, r, mtry_mean_real),
      ordinary_family_opportunity = candidate_probability_random(P, 1L, mtry_mean_real)
    )
  }
))
write.csv(opportunity, file.path(output_dir, "candidate_opportunity.csv"), row.names = FALSE)

metadata <- c(
  sprintf("completed_utc: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("grf_version: %s", packageVersion("grf")),
  sprintf("allow_unpinned_grf: %s", allow_unpinned_grf),
  sprintf("allow_data_mismatch: %s", allow_data_mismatch),
  sprintf("threads: %d", threads),
  sprintf("real_trees: %d", real_trees),
  sprintf("simulation_trees: %d", simulation_trees),
  sprintf("real_n: %d", nrow(X_real)),
  sprintf("real_clusters: %d", length(unique(clusters_real))),
  sprintf("real_mtry_poisson_mean: %d", mtry_mean_real),
  sprintf("simulation_mtry_poisson_mean: %d", mtry_mean_sim),
  sprintf("quotient_clone_max_prediction_difference: %.17g", quotient_difference),
  sprintf("quotient_monotone_max_prediction_difference: %.17g", routine_quotient_difference)
)
writeLines(metadata, file.path(output_dir, "run_metadata.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
if (save_objects) {
  saveRDS(
    simulation_predictions,
    file.path(output_dir, "simulation_predictions.rds")
  )
}

cat("Completed representation-multiplicity study.\n")
print(real_results)
print(simulation_summary)
