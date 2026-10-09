lex_less <- function(x, y) {
  difference <- which(x != y)
  length(difference) > 0L && x[difference[[1L]]] < y[difference[[1L]]]
}

class_state <- function(pooled) {
  pooled <- as.matrix(pooled)
  stopifnot(!anyNA(pooled), !anyDuplicated(colnames(pooled)))
  codes <- matrix(0L, nrow(pooled), ncol(pooled), dimnames = dimnames(pooled))
  maps <- setNames(vector("list", ncol(pooled)), colnames(pooled))
  for (j in seq_len(ncol(pooled))) {
    values <- sort(unique(pooled[, j]))
    ranks <- match(pooled[, j], values)
    reversed <- length(values) + 1L - ranks
    code <- if (lex_less(reversed, ranks)) reversed else ranks
    codes[, j] <- code
    maps[[j]] <- list(values = values, ranks = code[match(values, pooled[, j])])
  }
  unique_columns <- which(!duplicated(codes, MARGIN = 2L))
  order <- integer()
  for (j in unique_columns) {
    position <- length(order) + 1L
    while (position > 1L && lex_less(codes[, j], codes[, order[position - 1L]])) {
      position <- position - 1L
    }
    order <- append(order, j, after = position - 1L)
  }
  canonical <- codes[, order, drop = FALSE]
  colnames(canonical) <- sprintf("class_%03d", seq_len(ncol(canonical)))
  class_ids <- vapply(seq_len(ncol(codes)), function(j) {
    matches <- which(vapply(seq_len(ncol(canonical)), function(k) {
      identical(codes[, j], canonical[, k])
    }, logical(1)))
    stopifnot(length(matches) == 1L)
    matches[[1L]]
  }, integer(1))
  names(class_ids) <- colnames(pooled)
  list(canonical = canonical, maps = maps, class_ids = class_ids,
       mtry = max(floor(ncol(canonical) / 3L), 1L))
}

encode_classes <- function(x, state) {
  x <- as.matrix(x)
  ids <- state$class_ids[colnames(x)]
  stopifnot(!anyNA(ids), setequal(unique(ids), seq_len(ncol(state$canonical))))
  out <- matrix(0L, nrow(x), ncol(state$canonical))
  for (id in seq_len(ncol(out))) {
    j <- which(ids == id)[[1L]]
    map <- state$maps[[colnames(x)[[j]]]]
    index <- match(x[, j], map$values)
    stopifnot(!anyNA(index))
    out[, id] <- map$ranks[index]
  }
  colnames(out) <- colnames(state$canonical)
  out
}
