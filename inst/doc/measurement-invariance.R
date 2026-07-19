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
mi_data <- HolzingerSwineford1939
mi_data[ordered_names] <- lapply(
  mi_data[ordered_names],
  function(x) {
    ordered(cut(
      x,
      breaks = c(-Inf, quantile(x, c(.25, .50, .75)), Inf),
      include.lowest = TRUE
    ))
  }
)

## ----configural---------------------------------------------------------------
configural <- cfa(
  model, mi_data,
  ordered = ordered_names,
  group = "school", estimator = "WLSMV"
)

pvalues(configural, c("SB", "SS", "ALL", "PEBA4"))

## ----thresholds---------------------------------------------------------------
equal_thresholds <- cfa(
  model, mi_data,
  ordered = ordered_names,
  group = "school", group.equal = "thresholds", estimator = "WLSMV"
)

pvalues_nested(
  equal_thresholds, configural,
  tests = c("SB", "SS", "ALL", "PEBA4")
)

## ----loadings-----------------------------------------------------------------
equal_loadings <- cfa(
  model, mi_data,
  ordered = ordered_names,
  group = "school",
  group.equal = c("thresholds", "loadings"),
  estimator = "WLSMV"
)

pvalues_nested(
  equal_loadings, equal_thresholds,
  tests = c("SB", "SS", "ALL", "PEBA4")
)

## ----provenance---------------------------------------------------------------
attr(
  pvalues_nested(
    equal_loadings, equal_thresholds,
    tests = "PEBA4"
  ),
  "semtests"
)

