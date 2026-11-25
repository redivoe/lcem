X <- expand.grid(a = letters[1:2], b = LETTERS[1:3], c = 0:1)
counts <- c(16, 3, 5, 11, 2, 8, 4, 12, 6, 1, 14, 7)

test_that("one class gives the closed-form MAP", {
  fit <- lca_em(X, G = 1, counts = counts, alpha = 2, nstarts = 1)
  expect_equal(unname(fit$pi), 1)
  for (j in seq_along(X)) {
    n <- as.numeric(tapply(counts, X[[j]], sum))
    expect_equal(as.numeric(fit$theta[[j]]), (n + 1) / sum(n + 1))
  }
})

test_that("weighted patterns and expanded rows give the same fit", {
  set.seed(1)
  weighted <- lca_em(X, G = 2, counts = counts, alpha = 2, nstarts = 1)
  set.seed(1)
  expanded <- lca_em(X[rep(seq_len(nrow(X)), counts), ], G = 2,
                     alpha = 2, nstarts = 1)
  expect_equal(weighted$pi, expanded$pi)
  expect_equal(weighted$theta, expanded$theta)
  expect_equal(weighted$loglik, expanded$loglik)
})

test_that("multiple starts return the best fit", {
  set.seed(1)
  multi <- lca_em(X, G = 2, counts = counts, alpha = 2,
                  nstarts = 4, max_iter = 1)
  set.seed(1)
  single <- lapply(1:4, function(i)
    lca_em(X, G = 2, counts = counts, alpha = 2, nstarts = 1, max_iter = 1))
  scores <- vapply(single, function(fit) fit$logpost, numeric(1))
  expect_true(which.max(scores) != length(scores))
  best <- single[[which.max(scores)]]
  expect_equal(multi$lp_starts, scores)
  expect_equal(multi$logpost, best$logpost)
  expect_equal(multi$pi, best$pi)
  expect_equal(multi$theta, best$theta)
})
