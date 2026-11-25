# lcem

A small, dependency-free R package for latent class analysis of categorical
data. Choose any positive number of classes `G`; each variable can have a
different number of levels. Estimation uses soft EM on distinct response
patterns with frequency weights, rather than repeating calculations for
identical observations.

## Install and fit

```r
install.packages("lcem_0.1.1.tar.gz", repos = NULL, type = "source")
library(lcem)

set.seed(7)
X <- data.frame(
  a = sample(letters[1:3], 500, replace = TRUE),
  b = sample(0:1, 500, replace = TRUE),
  c = sample(1:5, 500, replace = TRUE)
)
fit <- lca_em(X, G = 3)
fit
fit$pi                         # class proportions
fit$theta                      # one class-by-level matrix per variable
head(predict(fit))              # class assignments for the input rows
head(predict(fit, type = "posterior"))
predict(fit, newdata = X[1:5, ], type = "posterior")
```

Non-factor columns are converted to factors. Factor levels are preserved,
including unused levels. Use factors when you want a specific category order.
Missing values are not supported.

## Patterns

Patterns are ordered by their first appearance in the input. Only observed
patterns are stored; the full Cartesian product of levels is never created.

```r
fit$patterns                   # distinct rows, with original category labels
fit$counts                     # frequency of each pattern
fit$tau                          # posterior membership for each pattern
fit$index                      # maps input rows to patterns
fit$tau[fit$index, , drop = FALSE]  # membership for each input row

# Refit from an already tabulated data set:
fit2 <- lca_em(fit$patterns, G = 3, counts = fit$counts)
```

`counts` must contain positive frequencies. Duplicate rows in a tabulated
input are combined too. When fitting tabulated input, `predict(fit)` returns
one result per supplied row, without expanding its frequency.

## Dirichlet priors

Class proportions have a `Dirichlet(delta)` prior. Each class-specific item
probability vector has its own `Dirichlet(alpha)` prior, independent of the
others. The defaults `alpha = 1, delta = 1` are uniform densities on the
simplexes, so MAP estimation reduces exactly to maximum likelihood.

```r
# Symmetric shrinkage toward equal probabilities:
fit <- lca_em(X, G = 3, alpha = 2, delta = 2)

# Separate priors by variable, category and class:
alpha <- list(
  a = c(2, 3, 2),               # category concentrations, shared by classes
  b = matrix(c(4, 2, 2, 2, 4, 2), nrow = 3),
  c = 1                        # uniform for every class and category
)
fit <- lca_em(X, G = 3, alpha = alpha, delta = c(2, 2, 3))
fit$levels                     # order used for categories and prior columns
```

`alpha` is a scalar or a list in input-column order. A list element is a
scalar, a vector with one entry per category, or a `G` by number-of-categories
matrix. Names are labels, not a mechanism for matching or reordering priors.
`delta` is a scalar or a length-`G` vector.

Concentrations must be at least 1. Values below 1 can give an unbounded
posterior at zero probabilities, so they are excluded from this MAP estimator.

For pattern `s`, let `n_s` be its frequency and `r_sg` its posterior class
probability. With `N_g = sum_s n_s r_sg` and
`N_jgk = sum_s n_s r_sg I(x_sj = k)`, the M-step is

```text
pi_g       = (N_g + delta_g - 1) / (N + sum(delta) - G)
theta_jgk  = (N_jgk + alpha_jgk - 1) /
             (N_g + sum_k(alpha_jgk) - K_j)
```

An empty class with uniform item priors receives uniform item probabilities;
its item probabilities do not affect the objective. The E-step uses log
probabilities and log-sum-exp, retaining exact zero estimates without adding
arbitrary smoothing constants.

## Fitting controls and results

```r
set.seed(7)
fit <- lca_em(X, G = 3, nstarts = 10, max_iter = 2000, eps = 1e-9)
fit$converged
fit$log_post_seq                      # log posterior for the selected start
fit$lp_starts                     # final log posterior for every start
logLik(fit)                    # data log likelihood, with df and nobs
```

The best final posterior across random starts is returned. Stopping uses
`abs(new - old) <= eps * (1 + abs(old))` on the log posterior. `log_post_seq` includes
the initial value; `iter` counts completed EM updates. `logpost` omits fixed
Dirichlet normalising constants and equals `loglik` under uniform priors.
The returned responsibilities correspond to the final returned parameters.

EM can reach local optima, labels are arbitrary, and a requested `G` need not
be identifiable. More starts can help with local optima. No standard errors,
automatic class selection, or missing-data model are included. `AIC()` and
`BIC()` work through `logLik()`; with informative priors they evaluate the MAP
fit rather than the maximum-likelihood fit.

## Code and references

The implementation uses small procedural functions, direct updates, and named
lists, following the style of
[hcoref](https://github.com/redivoe/partition-estimation-record-linkage/tree/main/hcoref)
and [linkbrl](https://github.com/redivoe/linkbrl). The observed-pattern frequency
interface follows the approach described for `blca.em` in the
[BayesLCA manual](https://cran.r-project.org/web/packages/BayesLCA/BayesLCA.pdf).
The implementation here was written independently.

The `DESCRIPTION` maintainer entry is a placeholder to replace before
publication. To rebuild documentation, run `roxygen2::roxygenise()` in the
package directory. The package itself needs only base R and `stats`.
