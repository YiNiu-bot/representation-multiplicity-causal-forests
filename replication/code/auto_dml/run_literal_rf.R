#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)
package <- normalizePath(argv[[1L]], mustWork = TRUE)
output <- normalizePath(argv[[2L]], mustWork = TRUE)
dataset_selection <- if (length(argv) >= 3L) argv[[3L]] else "all"
spec_selection <- if (length(argv) >= 4L) as.integer(argv[[4L]]) else 1:3
suppressPackageStartupMessages({
  library(nnet)
  library(randomForest)
})
stopifnot(getRversion() == "4.0.3", R.version$arch == "x86_64",
          as.character(packageVersion("randomForest")) == "4.6.14")
for (file in c("primitives.R", "stage0.R", "stage1.R", "stage2.R",
               "specifications_intersection.R")) source(file.path(package, file))
native_stage1 <- get_stage1
read_data <- function(stem) {
  data <- read.table(file.path(package, paste0(stem, ".txt")))
  names(data) <- c("treat", "age", "education", "black", "hispanic",
                   "married", "nodegree", "re74", "re75", "re78")
  data
}
published <- data.frame(dataset = rep(c("NSW", "PSID", "CPS"), each = 3L),
  specification = rep(1:3, 3L),
  atet = c(3106.55, 3077.26, 2785.13, 1521.92, 1336.66, 2010.53,
           1639.95, 1584.12, 1906.62),
  se = c(1327.02, 1318.67, 819.17, 977.08, 956.22, 987.73,
         616.08, 616.33, 651.77))
selected <- if (dataset_selection == "all") c("NSW", "PSID", "CPS") else
  dataset_selection
stopifnot(all(selected %in% c("NSW", "PSID", "CPS")),
          all(spec_selection %in% 1:3))
cells <- expand.grid(dataset = selected, specification = spec_selection,
                     stringsAsFactors = FALSE)
arch <- R.version$arch
writeLines(capture.output(sessionInfo()), file.path(output,
  paste0("literal_rf_session_", arch, ".txt")))

run_cell <- function(i) {
  cell <- cells[i, ]
  stem <- paste0(tolower(cell$dataset), "_spec", cell$specification, "_", arch)
  directory <- file.path(output, stem)
  dir.create(directory, showWarnings = FALSE)
  if (file.exists(file.path(directory, "result.csv"))) {
    return(read.csv(file.path(directory, "result.csv"), stringsAsFactors = FALSE))
  }
  control <- c(NSW = "nswre74_control", PSID = "psid_controls",
               CPS = "cps_controls")[[cell$dataset]]
  raw <- rbind(read_data("nswre74_treated"), read_data(control))
  data <- get_data_intersection(raw, cell$specification)
  saveRDS(data, file.path(directory, "original_inputs.rds"))
  p <- length(b2(data[[2L]][1L], data[[3L]][1L, ]))
  p0 <- if (p > 60L) ceiling(p / 40) else ceiling(p / 4)
  capture <- new.env(parent = emptyenv())
  capture$fold <- 0L
  capture$rho <- vector("list", 5L)
  capture$rng_start <- vector("list", 5L)
  capture$rng_end <- vector("list", 5L)
  recorder <- function(...) {
    capture$fold <- capture$fold + 1L
    fold <- capture$fold
    capture$rng_start[[fold]] <- .Random.seed
    fit <- native_stage1(...)
    capture$rng_end[[fold]] <- .Random.seed
    capture$rho[[fold]] <- get("rho_hat", envir = environment(fit[[1L]]))
    message("STAGE1 ", stem, " fold ", fold)
    saveRDS(list(rho = capture$rho[[fold]],
                 rng_start = capture$rng_start[[fold]],
                 rng_end = capture$rng_end[[fold]]),
            file.path(directory, paste0("fold", fold, "_stage1.rds")))
    fit
  }
  assign("get_stage1", recorder, envir = .GlobalEnv)
  RNGversion("4.0.3")
  set.seed(1L)
  start <- proc.time()
  message("START ", stem, " original RF pipeline")
  result <- rrr(data[[1L]], data[[2L]], data[[3L]], p0,
                0, 0.2, 10, b2, 1, 2, 0, 0.5, 0.5)
  target <- published[published$dataset == cell$dataset &
    published$specification == cell$specification, ]
  row <- data.frame(dataset = cell$dataset, specification = cell$specification,
    r_version = as.character(getRversion()), arch = arch,
    forest_version = as.character(packageVersion("randomForest")),
    nnet_version = as.character(packageVersion("nnet")),
    treated = result[[1L]], controls = result[[2L]], columns = p,
    mtry = floor(p / 3), trees = 1000L, penalty = 0.5,
    atet = result[[3L]], se = result[[4L]],
    published_atet = target$atet, published_se = target$se,
    elapsed = (proc.time() - start)[["elapsed"]],
    atet_match = sprintf("%.2f", result[[3L]]) == sprintf("%.2f", target$atet),
    se_match = sprintf("%.2f", result[[4L]]) == sprintf("%.2f", target$se))
  saveRDS(list(row = row, result = result, rho = capture$rho,
    rng_start = capture$rng_start, rng_end = capture$rng_end,
    final_rng = .Random.seed), file.path(directory, "literal_result.rds"))
  write.csv(row, file.path(directory, "result.csv"), row.names = FALSE)
  message("DONE ", stem, " ATET ", sprintf("%.2f", row$atet),
          " SE ", sprintf("%.2f", row$se), " exact ", row$atet_match && row$se_match)
  row
}
rows <- parallel::mclapply(seq_len(nrow(cells)), run_cell,
  mc.cores = min(9L, parallel::detectCores() - 1L),
  mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(!any(vapply(rows, inherits, logical(1), "try-error")))
result <- do.call(rbind, rows)
write.csv(result, file.path(output, paste0("literal_rf_", dataset_selection,
  "_", paste(spec_selection, collapse = ""), "_", arch, "_results.csv")),
  row.names = FALSE)
print(result, row.names = FALSE)
