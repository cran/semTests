## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(semTests)

## ----models, echo=FALSE, out.width="100%", fig.alt="Two one-factor models with four indicators. The unrestricted model has loadings a and b on x2 and x3; the restricted model uses the same loading on both paths.", fig.cap="The same factor model before and after imposing equality. Blue paths identify the restriction. Residuals are uncorrelated; intercepts and variances are omitted from the drawing."----
knitr::include_graphics(if (knitr::is_latex_output()) {
  "figures/restriction-models.pdf"
} else {
  "figures/restriction-models.svg"
})

## ----fits---------------------------------------------------------------------
d <- lavaan::HolzingerSwineford1939
model <- "f =~ x1 + a*x2 + b*x3 + x4"
fit <- lavaan::cfa(model, d, meanstructure = TRUE)
restricted <- lavaan::cfa(paste(model, "a == b", sep = "\n"),
                          d, meanstructure = TRUE)
lavaan::coef(fit)[c("a", "b")]

## ----score--------------------------------------------------------------------
score <- semtest(restricted, fit)
score

## ----template-----------------------------------------------------------------
template <- lavaan::cfa(model, d, meanstructure = TRUE, do.fit = FALSE)
semtest(restricted, template)

## ----criterion----------------------------------------------------------------
criterion <- semtest(restricted, fit, statistic = "criterion")
criterion

## ----comparison, echo=FALSE---------------------------------------------------
comparison <- data.frame(
  Test = c("Score", "Criterion difference"),
  Statistic = c(score$statistic, criterion$statistic),
  `p-value` = c(score$p_value, criterion$p_value),
  check.names = FALSE
)
knitr::kable(comparison, digits = 4)

## ----weights------------------------------------------------------------------
model3 <- "f =~ x1 + a*x2 + b*x3 + c*x4"
fit3 <- lavaan::cfa(model3, d, meanstructure = TRUE)
equal3 <- lavaan::cfa(paste(model3, "a == b", "b == c", sep = "\n"),
                      d, meanstructure = TRUE)
weights <- c("curvature", "expected_curvature", "score_covariance")
do.call(rbind, lapply(weights, function(w) {
  test <- semtest(equal3, fit3, weight = w)
  data.frame(weight = w, statistic = test$statistic, df = test$df,
             reference = test$reference, p_value = test$p_value)
}))

## ----control------------------------------------------------------------------
inference_control()

## ----control-override---------------------------------------------------------
semtest(equal3, fit3, curvature = "expected")
semtest(equal3, fit3, control = inference_control(curvature = "expected"))

## ----calibrations-------------------------------------------------------------
test3 <- semtest(equal3, fit3)
test3$eigenvalues
pvalues(test3, c("ALL", "SB", "SS", "PEBA2"))

## ----global-------------------------------------------------------------------
global <- semtest(fit)
global
pvalues(global, c("SB", "ALL"))

## ----fiml---------------------------------------------------------------------
d_missing <- d
set.seed(1)
d_missing$x2[sample(nrow(d), 60)] <- NA
fit_fiml <- lavaan::cfa(model, d_missing, meanstructure = TRUE,
                        missing = "fiml")
pvalues(semtest(fit_fiml), c("STD", "SB", "ALL"))

## ----weights-dwls-------------------------------------------------------------
fit_dwls <- lavaan::cfa(model3, d, meanstructure = TRUE, estimator = "DWLS",
                        ordered = FALSE)
equal_dwls <- lavaan::cfa(paste(model3, "a == b", "b == c", sep = "\n"), d,
                          meanstructure = TRUE, estimator = "DWLS",
                          ordered = FALSE)
semtest(equal_dwls, fit_dwls, weight = "score_covariance")
semtest(equal_dwls, fit_dwls, weight = "score_covariance",
        weight_effect = "fixed")

## ----inspect------------------------------------------------------------------
score$diagnostics$nuisance_stationarity_norm

