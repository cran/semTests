## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)
options(warn = -1)

## ----setup--------------------------------------------------------------------
library("semTests")
library("lavaan")

## ----basic--------------------------------------------------------------------
model <- "visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9"
fit <- cfa(model, HolzingerSwineford1939, estimator = "MLM")
pvalues(fit)

## ----grammar------------------------------------------------------------------
pvalues(fit, c("SB_ML", "pEBA4_ML", "pOLS2_ML"))

## ----nested-------------------------------------------------------------------
constrained <- "visual  =~ x1 + x2 + x3
                textual =~ a*x4 + a*x5 + a*x6
                speed   =~ x7 + x8 + x9"
m1 <- cfa(model, HolzingerSwineford1939, estimator = "MLM")
m0 <- cfa(constrained, HolzingerSwineford1939, estimator = "MLM")
pvalues_nested(m0, m1)

## ----estimators---------------------------------------------------------------
gls <- cfa(model, HolzingerSwineford1939, estimator = "GLS")
pvalues(gls, "pEBA4")

uls <- cfa(model, HolzingerSwineford1939, estimator = "ULS",
           test = "satorra.bentler")
pvalues(uls, "pEBA4")

## ----fiml---------------------------------------------------------------------
HS <- HolzingerSwineford1939
set.seed(1)
HS$x1[sample(nrow(HS), 60)] <- NA
fit_fiml <- cfa(model, HS, missing = "fiml", estimator = "MLR")
pvalues(fit_fiml)

