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

## ----missing-data-------------------------------------------------------------
fiml_data <- HolzingerSwineford1939
set.seed(20260717)
driver <- as.numeric(scale(fiml_data$x9))

missing_x1 <- runif(nrow(fiml_data)) <
  plogis(qlogis(.20) + .8 * driver)
missing_x5 <- runif(nrow(fiml_data)) <
  plogis(qlogis(.20) - .8 * driver)

fiml_data$x1[missing_x1] <- NA
fiml_data$x5[missing_x5] <- NA
colMeans(is.na(fiml_data[paste0("x", 1:9)]))

## ----fiml-fit-----------------------------------------------------------------
fit_observed <- cfa(
  model, fiml_data,
  missing = "fiml",
  estimator = "MLR"
)

pvalues(fit_observed, c("SB", "SS", "PEBA4"))

## ----random-x-fit-------------------------------------------------------------
model_x <- "
  visual =~ x1 + x2 + x3
  visual ~ ageyr
"

random_x_data <- fiml_data
random_x_data$ageyr[seq(7, nrow(random_x_data), by = 23)] <- NA

fit_random_x <- cfa(
  model_x, random_x_data,
  missing = "fiml", estimator = "MLR",
  fixed.x = FALSE, conditional.x = FALSE
)

pvalues(fit_random_x, c("SB", "SS", "PEBA2"))

## ----information-study, echo = FALSE------------------------------------------
information_coverage <- data.frame(
  Method = c(
    "Expected, model-based",
    "Expected + sandwich",
    "Observed H1, model-based",
    "Observed H1 + sandwich",
    "Full observed Hessian, model-based",
    "Full observed Hessian + sandwich"
  ),
  `Correct, normal` = c(94.48, 94.13, 94.67, 94.39, 94.73, 94.72),
  `Correct, t5` = c(79.41, 92.45, 79.98, 92.52, 80.14, 93.64),
  `Misspecified, normal` = c(92.10, 91.53, 92.72, 92.46, 93.46, 94.53),
  `Misspecified, t5` = c(75.21, 89.28, 76.38, 89.98, 77.51, 93.00),
  check.names = FALSE
)
knitr::kable(
  information_coverage,
  digits = 2,
  caption = paste(
    "Mean 95% coverage across sample sizes, missingness mechanisms,",
    "and focal parameters. Misspecified cells use centered coverage."
  )
)

## ----expected-fit-------------------------------------------------------------
fit_expected <- cfa(
  model, fiml_data,
  missing = "fiml",
  estimator = "MLR",
  information = "expected"
)

max(abs(coef(fit_observed) - coef(fit_expected)))

## ----convention-table---------------------------------------------------------
tests <- c("SB", "SS", "PEBA4")

comparison <- cbind(
  observed_fit__observed_convention = as.numeric(
    pvalues(fit_observed, tests, fiml.convention = "observed")
  ),
  observed_fit__lavaan_convention = as.numeric(
    pvalues(fit_observed, tests, fiml.convention = "lavaan")
  ),
  expected_fit__observed_convention = as.numeric(
    pvalues(fit_expected, tests, fiml.convention = "observed")
  ),
  expected_fit__lavaan_convention = as.numeric(
    pvalues(fit_expected, tests, fiml.convention = "lavaan")
  )
)
rownames(comparison) <- tests
format(comparison, scientific = TRUE, digits = 4)

## ----nested-fits--------------------------------------------------------------
constrained <- "
  visual  =~ x1 + a*x2 + a*x3
  textual =~ x4 + b*x5 + b*x6
  speed   =~ x7 + x8 + x9
"

m1 <- cfa(
  model, fiml_data,
  missing = "fiml", estimator = "MLR"
)
m0 <- cfa(
  constrained, fiml_data,
  missing = "fiml", estimator = "MLR"
)

## ----nested-comparison--------------------------------------------------------
nested_tests <- c("SB", "SS", "PEBA2")
nested <- cbind(
  observed_delta = as.numeric(pvalues_nested(
    m0, m1,
    tests = nested_tests,
    A.method = "delta", fiml.convention = "observed"
  )),
  observed_exact = as.numeric(pvalues_nested(
    m0, m1,
    tests = nested_tests,
    A.method = "exact", fiml.convention = "observed"
  )),
  lavaan_delta = as.numeric(pvalues_nested(
    m0, m1,
    tests = nested_tests,
    A.method = "delta", fiml.convention = "lavaan"
  ))
)
rownames(nested) <- nested_tests
format(nested, scientific = TRUE, digits = 4)

## ----multigroup-fits----------------------------------------------------------
mg1 <- cfa(
  model, fiml_data,
  group = "school",
  missing = "fiml", estimator = "MLR"
)
mg0 <- cfa(
  model, fiml_data,
  group = "school", group.equal = "loadings",
  missing = "fiml", estimator = "MLR"
)

pvalues(mg1, c("SB", "SS", "PEBA4"))
pvalues_nested(mg0, mg1, tests = c("SB", "SS", "PEBA4"))

