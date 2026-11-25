# lcem

Latent class analysis for categorical data with any number of classes and
categories per variable. Estimation uses EM on distinct response patterns,
weighted by their frequencies, with Dirichlet priors.

## Installation

```r
# install.packages("remotes")
remotes::install_github("redivoe/lcem")
```

## Example

```r
library(lcem)

set.seed(1)
X <- data.frame(a = sample(letters[1:3], 200, replace = TRUE),
                b = sample(0:1, 200, replace = TRUE),
                c = sample(1:4, 200, replace = TRUE))

fit <- lca_em(X, G = 2)
fit$pi                                      # class proportions
fit$theta                                   # class-specific category probabilities
head(predict(fit, type = "posterior"))      # posterior class probabilities
```
