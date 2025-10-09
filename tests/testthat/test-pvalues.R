test_that("pvalues_one return p_values in (0, 1)", {
  pvals <- pvalues_(object, eba = 2, peba = 2, unbiased = 2, trad = "sb", pols = 2)
  expect_true(all(pvals < 1))
  expect_true(all(pvals > 0))
})
