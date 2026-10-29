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
quartiles <- function(x) {
  ordered(cut(x, c(-Inf, quantile(x, c(.25, .50, .75)), Inf)))
}

items <- paste0("x", 1:9)
ordinal_data <- HolzingerSwineford1939
ordinal_data[items] <- lapply(ordinal_data[items], quartiles)

## ----fit----------------------------------------------------------------------
fit <- cfa(model, ordinal_data, ordered = items, estimator = "WLSMV")
pvalues(fit, c("STD", "SB", "SS", "ALL", "PEBA4"))

## ----lavaan-------------------------------------------------------------------
as.numeric(fitMeasures(fit, "pvalue.scaled"))

## ----mixed--------------------------------------------------------------------
mixed_data <- HolzingerSwineford1939
mixed_data[items[1:6]] <- lapply(mixed_data[items[1:6]], quartiles)
fit_mixed <- cfa(model, mixed_data, ordered = items[1:6], estimator = "WLSMV")
pvalues(fit_mixed)

## ----pairwise-----------------------------------------------------------------
pairwise_data <- ordinal_data
set.seed(314)
for (item in items) {
  pairwise_data[sample(nrow(pairwise_data), 35), item] <- NA
}
fit_pairwise <- cfa(model, pairwise_data,
  ordered = items, estimator = "WLSMV", missing = "pairwise"
)
pvalues(fit_pairwise, c("SB", "SS", "PEBA4"))

## ----invariance-fits----------------------------------------------------------
configural <- cfa(model, ordinal_data,
  ordered = items, group = "school", estimator = "WLSMV"
)
equal_thresholds <- cfa(model, ordinal_data,
  ordered = items, group = "school", estimator = "WLSMV",
  group.equal = "thresholds"
)
equal_loadings <- cfa(model, ordinal_data,
  ordered = items, group = "school", estimator = "WLSMV",
  group.equal = c("thresholds", "loadings")
)

## ----invariance-tests---------------------------------------------------------
thresholds_test <- pvalues_nested(equal_thresholds, configural,
  tests = c("SB", "SS", "PALL")
)
loadings_test <- pvalues_nested(equal_loadings, equal_thresholds,
  tests = c("SB", "SS", "PALL")
)
thresholds_test
loadings_test

## ----configural---------------------------------------------------------------
pvalues(configural)

