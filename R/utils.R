lca_patterns_simple <- function(X, counts) {
  X <- as.data.frame(X)
  stopifnot(nrow(X) > 0L, ncol(X) > 0L, !anyNA(X))
  if (is.null(counts)) counts <- rep(1, nrow(X))
  stopifnot(length(counts) == nrow(X), all(is.finite(counts)), all(counts > 0))
  X[] <- lapply(X, function(x) if (is.factor(x)) x else factor(x))
  levels <- lapply(X, levels)
  
  x <- matrix(unlist(lapply(X, as.integer), use.names = FALSE), nrow(X), ncol(X))
  
  key <- do.call(paste, c(as.data.frame(x), sep = ":"))
  keep <- !duplicated(key)
  index <- match(key, key[keep])
  out <- list(x = x[keep, , drop = FALSE], patterns = X[keep, , drop = FALSE],
              counts = as.numeric(rowsum(counts, index, reorder = FALSE)),
              index = index, levels = levels)
  return(out)
}


lca_patterns <- function(X, counts = NULL) {
  X <- as.data.frame(X)
  stopifnot(nrow(X) > 0L, ncol(X) > 0L, !anyNA(X))
  if (!is.null(counts)){
    stopifnot(length(counts) == nrow(X), all(is.finite(counts)), all(counts > 0))
  }
  levels <- vector("list", ncol(X))
  names(levels) <- names(X)
  key <- numeric(nrow(X))
  stride <- 1
  for (j in seq_len(ncol(X))) {
    f <- if (is.factor(X[[j]])) X[[j]] else factor(X[[j]])
    levels[[j]] <- levels(f)
    L <- nlevels(f)
    # Keep numeric keys exact when category products become large.
    if (stride * L > 2^52) {
      key <- match(key, unique(key)) - 1L
      stride <- max(key) + 1
    }
    if (stride * L > 2^52) {
      key <- paste(key, as.integer(f), sep = ":")
      key <- match(key, unique(key)) - 1L
      stride <- max(key) + 1
    } else {
      key <- key + (as.integer(f) - 1L) * stride
      stride <- stride * L
    }
  }
  keep <- !duplicated(key)
  index <- match(key, key[keep])
  # Build factors and the integer matrix only for distinct patterns.
  patterns <- X[keep, , drop = FALSE]
  patterns[] <- Map(function(x, lev) if (is.factor(x)) x else factor(x, levels = lev),
                    patterns, levels)
  x <- matrix(unlist(lapply(patterns, as.integer), use.names = FALSE),
              nrow(patterns), ncol(patterns))
  
  if (is.null(counts)){
    counts <- as.numeric(tabulate(index, nrow(patterns)))
  }else{
    counts <- as.numeric(rowsum(counts, index, reorder = FALSE))
    
  }
  
  return(list(x = x, patterns = patterns, counts = counts, index = index, levels = levels))
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
