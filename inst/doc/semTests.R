## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## ----fit----------------------------------------------------------------------
library("semTests")
library("lavaan")

model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"
fit <- cfa(model, HolzingerSwineford1939)
pvalues(fit)

## ----battery------------------------------------------------------------------
pvalues(fit, c("STD_ML", "SB_ML", "SS_ML", "PEBA4_ML", "SB", "PEBA4_UG_RLS"))
fit_mlm <- cfa(model, HolzingerSwineford1939, estimator = "MLM")
as.numeric(fitMeasures(fit_mlm, "pvalue.scaled"))

## ----nested-------------------------------------------------------------------
restricted <- "
  visual  =~ x1 + a*x2 + a*x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"
fit0 <- cfa(restricted, HolzingerSwineford1939)
pvalues_nested(fit0, fit, tests = c("SB_ML", "PALL_UG_ML"))

## ----units--------------------------------------------------------------------
rescaled <- HolzingerSwineford1939
rescaled$x1 <- 10 * rescaled$x1
fit_rescaled <- cfa(model, rescaled)
fit0_rescaled <- cfa(restricted, rescaled)
sapply(c(weighted = "weighted", delta = "delta"), function(A) {
  c(
    original = pvalues_nested(fit0, fit, tests = "SB_ML", A.method = A)[[1]],
    rescaled = pvalues_nested(fit0_rescaled, fit_rescaled,
      tests = "SB_ML", A.method = A
    )[[1]]
  )
})

## ----estimators---------------------------------------------------------------
pvalues(cfa(model, HolzingerSwineford1939, estimator = "GLS"))
pvalues(cfa(model, HolzingerSwineford1939, estimator = "ULS"))

