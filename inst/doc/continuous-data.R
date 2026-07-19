## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## ----setup--------------------------------------------------------------------
library("semTests")
library("lavaan")

model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"
data <- HolzingerSwineford1939

## ----ml-fit-------------------------------------------------------------------
fit_ml <- cfa(model, data, estimator = "MLM")
pvalues(fit_ml)

## ----ml-battery---------------------------------------------------------------
pvalues(
  fit_ml,
  c("STD_ML", "SB_ML", "SS_ML", "ALL_ML", "PEBA4_ML")
)

## ----classical-options--------------------------------------------------------
pvalues(
  fit_ml,
  c("SB_UG_ML", "PEBA4_UG_ML", "PEBA4_RLS", "PEBA4_UG_RLS")
)

## ----provenance---------------------------------------------------------------
result <- pvalues(fit_ml, "PEBA4_ML")
attr(result, "semtests")

## ----nested-fit---------------------------------------------------------------
constrained <- "
  visual  =~ x1 + a*x2 + a*x3
  textual =~ x4 + b*x5 + b*x6
  speed   =~ x7 + x8 + x9
"

m1 <- cfa(model, data, estimator = "MLM")
m0 <- cfa(constrained, data, estimator = "MLM")

pvalues_nested(m0, m1)

## ----nested-methods-----------------------------------------------------------
pvalues_nested(m0, m1, method = "2000", tests = c("SB_ML", "PALL_ML"))

## ----least-squares------------------------------------------------------------
fit_gls <- cfa(model, data, estimator = "GLS")
fit_uls <- cfa(
  model, data,
  estimator = "ULS", test = "satorra.bentler"
)

pvalues(fit_gls, "PEBA4")
pvalues(fit_uls, "PEBA4")

