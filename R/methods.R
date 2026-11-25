#' Methods for fitted latent class models
#'
#' Print a fit, predict class membership, or extract the data log likelihood.
#'
#' @param x,object A fitted model returned by \code{\link{lca_em}}.
#' @param newdata Optional data frame or matrix with the training column names
#'   and known category levels. If \code{NULL}, use the original input rows.
#' @param type Either \code{"class"} for assignments or \code{"posterior"}
#'   for class membership probabilities.
#' @param ... Further arguments (unused).
#' @return \code{predict} returns class labels in \code{1:G} or a matrix with
#'   one row per input row and one column per class. \code{logLik} returns a
#'   \code{logLik} object with degrees of freedom and sample size;
#'   \code{print} invisibly returns the fitted model.
#' @seealso \code{\link{lca_em}}
#' @examples
#' set.seed(1)
#' fit <- lca_em(data.frame(a = rep(0:1, 10), b = rep(letters[1:4], 5)), G = 2)
#' print(fit)
#' head(predict(fit))
#' head(predict(fit, type = "posterior"))
#' logLik(fit)
#' @name lca_methods
NULL

#' @rdname lca_methods
#' @export
print.lca_em <- function(x, ...) {
  cat("Latent class analysis:", x$G, "classes,", sum(x$counts), "observations,",
      nrow(x$patterns), "patterns\n")
  cat("Log likelihood:", format(x$loglik), " | Iterations:", x$iter,
      " | Converged:", x$converged, "\n")
  print(x$pi)
  invisible(x)
}

#' @rdname lca_methods
#' @export
predict.lca_em <- function(object, newdata = NULL,
                           type = c("class", "posterior"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    tau <- object$tau[object$index, , drop = FALSE]
  } else {
    newdata <- as.data.frame(newdata)[names(object$levels)]
    x <- matrix(unlist(Map(match, newdata, object$levels)), nrow(newdata))
    if (anyNA(x)) stop("Unknown levels or missing values in newdata.")
    tau <- lca_estep(x, object$pi, object$theta, rep(1, nrow(x)))$tau
    if (any(!is.finite(tau))) stop("A new pattern has zero probability under this fit.")
    colnames(tau) <- names(object$pi)
  }
  if (type == "class") return(max.col(tau, ties.method = "first"))
  return(tau)
}

#' @rdname lca_methods
#' @export
logLik.lca_em <- function(object, ...) {
  return(structure(object$loglik, df = object$npar, nobs = sum(object$counts),
                   class = "logLik"))
}
