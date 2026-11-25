#' EM estimation of a latent class model
#'
#' Fit a latent class analysis model: a finite mixture model for categorical
#' variables assumed to be conditionally independent given the latent class.
#' Estimation uses the EM algorithm on distinct response patterns, weighted
#' by their frequencies.
#'
#' @param X Data frame or matrix of categorical variables, without missing
#'   values. Non-factor columns are converted to factors; existing factor
#'   levels, including unused levels, are retained.
#' @param G Positive integer giving the number of latent classes (groups).
#' @param alpha Dirichlet concentration parameters for the class-specific
#'   category probabilities. Supply either:
#'   \itemize{
#'     \item A scalar, used for every category, variable and class.
#'     \item A list with one entry per variable, in the same order as
#'       \code{X}. Each entry can be a scalar, a vector with one concentration
#'       per category (shared across classes), or a matrix with \code{G} rows
#'       and one column per category.
#'   }
#'   Category order follows the factor levels of each variable.
#' @param delta Dirichlet concentration parameters for the class proportions:
#'   a scalar or a vector of length \code{G}.
#' @param counts Optional vector of positive row frequencies. Defaults to one
#'   per row.
#' @param nstarts Number of random starts. The fit with the largest final
#'   log posterior is returned.
#' @param max_iter Maximum number of EM iterations per start.
#' @param eps Relative convergence tolerance for the log posterior.
#'
#' @details All concentration parameters must be at least one. The defaults
#'   \code{alpha = delta = 1} give uniform priors and maximum-likelihood EM.
#'   Other choices give maximum a posteriori (MAP) updates proportional to
#'   expected counts plus concentration minus one.
#'
#' @return A list of class \code{lca_em}, including:
#' \describe{
#'   \item{pi}{Estimated class proportions.}
#'   \item{theta}{A list of class-by-category probability matrices, one per
#'     variable. Rows correspond to classes and columns to categories.}
#'   \item{tau}{Posterior class probabilities for each distinct pattern.
#'     For input rows, use \code{tau[index, , drop = FALSE]}.}
#'   \item{patterns, counts, index, levels}{Distinct rows, their frequencies,
#'     the input-row-to-pattern mapping, and category labels.}
#'   \item{loglik, logpost}{Data log likelihood and log posterior, with fixed
#'     Dirichlet normalising constants omitted from the latter.}
#'   \item{log_post_seq, lp_starts}{Log posterior values for the selected run,
#'     including its initial value, and the final value for each random start.}
#'   \item{iter, converged}{Number of EM iterations and convergence indicator
#'     for the selected fit.}
#' }
#' @examples
#' set.seed(1)
#' X <- data.frame(a = sample(letters[1:3], 200, TRUE),
#'                 b = sample(0:1, 200, TRUE),
#'                 c = sample(1:4, 200, TRUE))
#' fit <- lca_em(X, G = 2)
#' fit$pi
#' fit$theta
#' head(predict(fit, type = "posterior"))
#' @seealso \code{\link{predict.lca_em}}, \code{\link{logLik.lca_em}}
#' @importFrom stats logLik predict
#' @export
lca_em <- function(X, G, alpha = 1, delta = 1, counts = NULL, nstarts = 5,
                   max_iter = 1000, eps = 1e-8) {
  call <- match.call()
  
  stopifnot(length(G) == 1L,
            G >= 1,
            G == as.integer(G),
            nstarts >= 1,
            max_iter >= 1,
            eps > 0)
  G <- as.integer(G)
  
  dat <- lca_patterns(X, counts)
  x <- dat$x
  counts <- dat$counts
  L <- lengths(dat$levels)
  p <- length(L)

  # Dirichlet priors
  if (!is.list(alpha)){
    alpha <- rep(list(alpha), p)
  }
  stopifnot(length(alpha) == p)
  for(j in seq_len(p)){
    # alpha[[j]] must be of size G x L[j]
    a <- alpha[[j]]
    if (is.null(dim(a))) {
      stopifnot(length(a) %in% c(1L, L[j]))
      a <- matrix(rep(a, each = G), nrow = G, ncol = L[j])
    }
    stopifnot(identical(dim(a), c(G, unname(L[j]))),
              all(is.finite(a)),
              all(a >= 1))
    alpha[[j]] <- a
  }
  
  stopifnot(length(delta) %in% c(1L, G))
  if(length(delta) == 1L){
    delta <- rep(delta, G)
  }
  stopifnot(all(is.finite(delta)), all(delta >= 1))
  
  
  lp_starts <- numeric(nstarts)
  best_rep <- NULL

  for (r in seq_len(nstarts)) {
    # initial parameters and E-step
    pi <- rep(1 / G, G)
    theta <- lapply(L, function(Lj) row_probs(matrix(stats::rexp(G * Lj), nrow = G, ncol = Lj)))
    
    estep <- lca_estep(x, pi, theta, counts)
    
    log_post_seq <- numeric(max_iter + 1L)
    log_post_seq[1L] <- estep$loglik + lca_logprior(pi, theta, delta, alpha)
    converged <- FALSE
    
    for (iter in seq_len(max_iter)) {
      # M-step: weighted expected counts plus prior minus one
      w <- estep$tau * counts
      pi <- colSums(w) + delta - 1
      pi <- pi / sum(pi)
      
      for(j in seq_len(p)){
        weighted_freqs <- matrix(0, nrow = G, ncol = L[j])
        for (l in seq_len(L[j])){
          weighted_freqs[, l] <- colSums(w[x[, j] == l, , drop = FALSE])
        }
        theta[[j]] <- row_probs(weighted_freqs + alpha[[j]] - 1)
      }

      # E-step and convergence
      estep <- lca_estep(x, pi, theta, counts)
      log_post_seq[iter + 1L] <- estep$loglik + lca_logprior(pi, theta, delta, alpha)
      converged <- abs(log_post_seq[iter + 1L] - log_post_seq[iter]) <= eps * (1 + abs(log_post_seq[iter]))
      if (converged) break
    }
    
    log_post_seq <- log_post_seq[seq_len(iter + 1L)]
    lp_starts[r] <- log_post_seq[iter + 1L]
    
    if (is.null(best_rep) || lp_starts[r] > best_rep$logpost){
      best_rep <- list(pi = pi, theta = theta, tau = estep$tau, loglik = estep$loglik,
                       logpost = lp_starts[r], log_post_seq = log_post_seq, iter = iter,
                       converged = converged)
    }
  }

  classes <- paste0("class", seq_len(G))
  names(best_rep$pi) <- classes
  names(best_rep$theta) <- names(alpha) <- names(dat$levels)
  for (j in seq_len(p)){
    dimnames(best_rep$theta[[j]]) <- dimnames(alpha[[j]]) <- list(classes, dat$levels[[j]])
  }
  colnames(best_rep$tau) <- classes
  out <- c(best_rep, dat[c("patterns", "counts", "index", "levels")],
           list(G = G, prior = list(alpha = alpha, delta = delta),
                lp_starts = lp_starts, npar = G - 1 + G * sum(L - 1), call = call))
  class(out) <- "lca_em"
  return(out)
}
