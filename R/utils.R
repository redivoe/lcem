lca_patterns <- function(X, counts) {
  X <- as.data.frame(X)
  stopifnot(nrow(X) > 0L, ncol(X) > 0L, !anyNA(X))
  if (is.null(counts)) counts <- rep(1, nrow(X))
  stopifnot(length(counts) == nrow(X), all(is.finite(counts)), all(counts > 0))
  X[] <- lapply(X, function(x) if (is.factor(x)) x else factor(x))
  levels <- lapply(X, levels)
  x <- matrix(unlist(lapply(X, as.integer)), nrow(X), ncol(X))
  key <- do.call(paste, c(as.data.frame(x), sep = ":"))
  keep <- !duplicated(key)
  index <- match(key, key[keep])
  out <- list(x = x[keep, , drop = FALSE], patterns = X[keep, , drop = FALSE],
              counts = as.numeric(rowsum(counts, index, reorder = FALSE)),
              index = index, levels = levels)
  return(out)
}

row_probs <- function(x) {
  # An empty class with a flat prior has no preferred item probabilities.
  x[rowSums(x) == 0, ] <- 1
  return(x / rowSums(x))
}

lca_estep <- function(x, pi, theta, counts) {
  point_loglik <- matrix(log(pi), nrow(x), length(pi), byrow = TRUE)
  for (j in seq_len(ncol(x))){
    point_loglik <- point_loglik + t(log(theta[[j]][, x[, j], drop = FALSE]))
  }
  m <- apply(point_loglik, 1L, max)
  tau <- exp(point_loglik - m)
  s <- rowSums(tau)
  loglik <- sum(counts * (m + log(s)))
  tau <- tau / s

  return(list("tau" = tau,
              "loglik" = loglik))
}

lca_logprior <- function(pi, theta, delta, alpha) {
  term <- function(p, a) {
    keep <- a > 1
    return(sum((a[keep] - 1) * log(p[keep])))
  }
  return(term(pi, delta) + sum(vapply(seq_along(theta), function(j)
    term(theta[[j]], alpha[[j]]), numeric(1))))
}
