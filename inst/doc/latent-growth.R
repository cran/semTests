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

growth_data <- Demo.growth

## ----missing------------------------------------------------------------------
set.seed(20260718)
driver <- as.numeric(scale(growth_data$t1))
growth_data$t3[runif(nrow(growth_data)) < plogis(qlogis(.20) + .8 * driver)] <- NA
growth_data$t4[runif(nrow(growth_data)) < plogis(qlogis(.20) + .8 * driver)] <- NA
colMeans(is.na(growth_data[paste0("t", 1:4)]))

## ----linear-------------------------------------------------------------------
linear <- growth(
  "i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
   s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4",
  growth_data,
  missing = "fiml", estimator = "MLR"
)

pvalues(linear, c("SB", "SS", "ALL", "PEBA4"))

## ----basis--------------------------------------------------------------------
basis <- growth(
  "i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
   s =~ 0*t1 + 1*t2 + t3 + t4",
  growth_data,
  missing = "fiml", estimator = "MLR"
)

pvalues_nested(linear, basis, tests = c("SB", "SS", "ALL", "PEBA2"))

