canonical_split_column <- function(x) {
  x <- as.numeric(x)
  missing <- is.na(x)
  observed <- x[!missing]
  levels <- sort(unique(observed))

  forward <- rep(NA_integer_, length(x))
  reverse <- rep(NA_integer_, length(x))
  if (length(observed)) {
    forward[!missing] <- match(observed, levels)
    reverse[!missing] <- length(levels) + 1L - forward[!missing]
  }

  missing_code <- length(levels) + 1L
  missing_key <- paste0(as.integer(missing), collapse = "")
  forward_key <- paste(ifelse(missing, missing_code, forward), collapse = ",")
  reverse_key <- paste(ifelse(missing, missing_code, reverse), collapse = ",")
  use_forward <- forward_key <= reverse_key
  list(
    certificate = paste0(
      missing_key,
      "|",
      if (use_forward) forward_key else reverse_key
    ),
    values = as.numeric(if (use_forward) forward else reverse)
  )
}

split_certificate <- function(x) {
  canonical_split_column(x)$certificate
}

certify_split_classes <- function(X.train, X.target = NULL) {
  X.train <- as.matrix(X.train)
  storage.mode(X.train) <- "double"
  if (is.null(colnames(X.train))) {
    colnames(X.train) <- sprintf("x%03d", seq_len(ncol(X.train)))
  }

  if (!is.null(X.target)) {
    X.target <- as.matrix(X.target)
    storage.mode(X.target) <- "double"
    stopifnot(
      ncol(X.target) == ncol(X.train),
      identical(colnames(X.target), colnames(X.train))
    )
    declared <- rbind(X.train, X.target)
  } else {
    declared <- X.train
  }

  encoded <- lapply(
    seq_len(ncol(declared)),
    function(j) canonical_split_column(declared[, j])
  )
  certificates <- vapply(encoded, `[[`, character(1), "certificate")
  class_certificates <- sort(unique(certificates))
  class_id <- match(certificates, class_certificates)
  representatives <- match(class_certificates, certificates)
  canonical_declared <- vapply(
    representatives,
    function(j) encoded[[j]]$values,
    numeric(nrow(declared))
  )
  if (length(representatives) == 1L) {
    canonical_declared <- matrix(canonical_declared, ncol = 1L)
  }
  colnames(canonical_declared) <- sprintf(
    "semantic_%04d",
    seq_along(representatives)
  )
  n.train <- nrow(X.train)

  class_map <- data.frame(
    raw_column = seq_len(ncol(X.train)),
    raw_name = colnames(X.train),
    class_id = class_id,
    representative = seq_len(ncol(X.train)) %in% representatives,
    representative_name = colnames(X.train)[representatives[class_id]],
    stringsAsFactors = FALSE
  )

  list(
    X.train = canonical_declared[seq_len(n.train), , drop = FALSE],
    X.target = if (is.null(X.target)) {
      NULL
    } else {
      canonical_declared[-seq_len(n.train), , drop = FALSE]
    },
    class_map = class_map,
    representatives = representatives,
    raw_dimension = ncol(X.train),
    semantic_dimension = length(representatives)
  )
}

class_sampled_causal_forest <- function(
    X.train,
    Y,
    W,
    X.target = NULL,
    mtry,
    ...
) {
  quotient <- certify_split_classes(X.train, X.target)
  semantic_mtry <- min(as.integer(mtry), quotient$semantic_dimension)
  fit <- grf::causal_forest(
    quotient$X.train,
    Y,
    W,
    mtry = semantic_mtry,
    ...
  )
  predictions <- if (is.null(quotient$X.target)) {
    as.numeric(predict(fit)$predictions)
  } else {
    as.numeric(predict(fit, quotient$X.target)$predictions)
  }

  list(
    fit = fit,
    predictions = predictions,
    class_map = quotient$class_map,
    raw_dimension = quotient$raw_dimension,
    fitted_dimension = quotient$semantic_dimension,
    semantic_dimension = quotient$semantic_dimension,
    mtry = semantic_mtry
  )
}

drop_exact_columns <- function(X.train, X.target = NULL) {
  X.train <- as.matrix(X.train)
  if (!is.null(X.target)) {
    X.target <- as.matrix(X.target)
    stopifnot(
      ncol(X.target) == ncol(X.train),
      identical(colnames(X.target), colnames(X.train))
    )
  }

  keep <- rep(TRUE, ncol(X.train))
  if (ncol(X.train) > 1L) {
    for (j in 2:ncol(X.train)) {
      for (k in seq_len(j - 1L)) {
        same_train <- keep[k] &&
          identical(unname(X.train[, j]), unname(X.train[, k]))
        same_target <- is.null(X.target) ||
          identical(unname(X.target[, j]), unname(X.target[, k]))
        if (same_train && same_target) {
          keep[j] <- FALSE
          break
        }
      }
    }
  }

  list(
    X.train = X.train[, keep, drop = FALSE],
    X.target = if (is.null(X.target)) NULL else X.target[, keep, drop = FALSE],
    keep = keep
  )
}
