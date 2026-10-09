canonical_split_column <- function(x) {
  x <- as.numeric(x)
  if (anyNA(x)) {
    stop(
      "Class-sampled certification requires a complete declared array.",
      call. = FALSE
    )
  }
  levels <- sort(unique(x))
  forward <- match(x, levels)
  reverse <- length(levels) + 1L - forward
  forward_key <- paste(forward, collapse = ",")
  reverse_key <- paste(reverse, collapse = ",")
  use_forward <- forward_key <= reverse_key
  list(
    certificate = if (use_forward) forward_key else reverse_key,
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
    if (!is.null(X.target) && !is.null(colnames(X.target))) {
      colnames(X.train) <- colnames(X.target)
    } else {
      colnames(X.train) <- sprintf("x%03d", seq_len(ncol(X.train)))
    }
  }

  if (!is.null(X.target)) {
    X.target <- as.matrix(X.target)
    storage.mode(X.target) <- "double"
    if (is.null(colnames(X.target))) {
      colnames(X.target) <- colnames(X.train)
    }
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
    certificate = certificates,
    class_certificate = class_certificates[class_id],
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
  candidate_poisson_mean <- as.numeric(mtry)
  if (length(candidate_poisson_mean) != 1L ||
      !is.finite(candidate_poisson_mean) || candidate_poisson_mean <= 0) {
    stop("mtry must be one positive finite number.", call. = FALSE)
  }
  fit <- grf::causal_forest(
    quotient$X.train,
    Y,
    W,
    mtry = candidate_poisson_mean,
    ...
  )
  predictions <- if (is.null(quotient$X.target)) {
    as.numeric(predict(fit)$predictions)
  } else {
    as.numeric(predict(fit, quotient$X.target)$predictions)
  }

  result <- list(
    fit = fit,
    predictions = predictions,
    canonical_target = quotient$X.target,
    class_map = quotient$class_map,
    raw_dimension = quotient$raw_dimension,
    fitted_dimension = quotient$semantic_dimension,
    semantic_dimension = quotient$semantic_dimension,
    mtry = candidate_poisson_mean
  )
  class(result) <- c("class_sampled_causal_forest", "list")
  result
}

predict.class_sampled_causal_forest <- function(
    object,
    newdata = NULL,
    estimate.variance = FALSE,
    ...
) {
  if (!is.null(newdata)) {
    stop(
      paste(
        "Predictions are certified only on the array declared at fitting.",
        "Refit with the new points supplied as X.target."
      ),
      call. = FALSE
    )
  }
  if (is.null(object$canonical_target)) {
    predict(object$fit, estimate.variance = estimate.variance, ...)
  } else {
    predict(
      object$fit,
      object$canonical_target,
      estimate.variance = estimate.variance,
      ...
    )
  }
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
