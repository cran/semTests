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

ordered_names <- paste0("x", 1:9)
ordinal_data <- HolzingerSwineford1939
ordinal_data[ordered_names] <- lapply(
  ordinal_data[ordered_names],
  function(x) {
    ordered(cut(
      x,
      breaks = c(-Inf, quantile(x, c(.33, .67)), Inf),
      include.lowest = TRUE
    ))
  }
)

## ----ordinal-fit--------------------------------------------------------------
fit_ordinal <- cfa(
  model, ordinal_data,
  ordered = ordered_names,
  estimator = "WLSMV"
)

pvalues(fit_ordinal, c("STD", "SB", "SS", "ALL", "PEBA4"))

## ----ordinal-provenance-------------------------------------------------------
attr(pvalues(fit_ordinal, "PEBA4"), "semtests")

## ----mixed-fit----------------------------------------------------------------
mixed_names <- paste0("x", 1:6)
mixed_data <- HolzingerSwineford1939
mixed_data[mixed_names] <- lapply(
  mixed_data[mixed_names],
  function(x) {
    ordered(cut(
      x,
      breaks = c(-Inf, quantile(x, c(.33, .67)), Inf),
      include.lowest = TRUE
    ))
  }
)

fit_mixed <- cfa(
  model, mixed_data,
  ordered = mixed_names,
  estimator = "WLSMV"
)
pvalues(fit_mixed, "PEBA4")

## ----pairwise-fit-------------------------------------------------------------
pairwise_data <- ordinal_data
set.seed(314)
for (variable in ordered_names) {
  pairwise_data[sample(nrow(pairwise_data), 35), variable] <- NA
}

fit_pairwise <- cfa(
  model, pairwise_data,
  ordered = ordered_names,
  estimator = "WLSMV",
  missing = "pairwise"
)
pvalues(fit_pairwise, c("SB", "SS", "PEBA4"))

## ----categorical-nested-------------------------------------------------------
constrained <- "
  visual  =~ x1 + a*x2 + a*x3
  textual =~ x4 + b*x5 + b*x6
  speed   =~ x7 + x8 + x9
"

m1 <- cfa(
  model, ordinal_data,
  ordered = ordered_names, estimator = "WLSMV"
)
m0 <- cfa(
  constrained, ordinal_data,
  ordered = ordered_names, estimator = "WLSMV"
)

pvalues_nested(
  m0, m1,
  method = "2000", A.method = "delta",
  tests = c("SB", "SS", "PALL")
)

