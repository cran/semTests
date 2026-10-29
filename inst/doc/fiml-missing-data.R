## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## ----data---------------------------------------------------------------------
library("semTests")
library("lavaan")

model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"

fiml_data <- HolzingerSwineford1939
set.seed(20260717)
driver <- as.numeric(scale(fiml_data$x9))
missing_x1 <- runif(nrow(fiml_data)) < plogis(qlogis(.20) + .8 * driver)
missing_x5 <- runif(nrow(fiml_data)) < plogis(qlogis(.20) - .8 * driver)
fiml_data$x1[missing_x1] <- NA
fiml_data$x5[missing_x5] <- NA
colMeans(is.na(fiml_data[c("x1", "x5")]))

## ----fit----------------------------------------------------------------------
fit <- cfa(model, fiml_data, missing = "fiml")
pvalues(fit, c("SB", "SS", "PEBA4"))

## ----information--------------------------------------------------------------
fit_expected <- cfa(model, fiml_data, missing = "fiml", information = "expected")
lavaan_sb <- function(fit) {
  lavTest(fit, test = "satorra.bentler")$satorra.bentler$pvalue
}
rbind(
  semTests = c(
    observed = pvalues(fit, "SB")[[1]],
    expected = pvalues(fit_expected, "SB")[[1]]
  ),
  lavaan = c(lavaan_sb(fit), lavaan_sb(fit_expected))
)

## ----predictor----------------------------------------------------------------
predictor_data <- fiml_data
predictor_data$ageyr[seq(7, nrow(predictor_data), by = 23)] <- NA
fit_predictor <- cfa("visual =~ x1 + x2 + x3
                      visual ~ ageyr",
  predictor_data,
  missing = "fiml", fixed.x = FALSE, conditional.x = FALSE
)
pvalues(fit_predictor, c("SB", "PALL"))

## ----growth-data--------------------------------------------------------------
growth_data <- Demo.growth
set.seed(20260718)
driver <- as.numeric(scale(growth_data$t1))
for (wave in c("t3", "t4")) {
  dropped <- runif(nrow(growth_data)) < plogis(qlogis(.20) + .8 * driver)
  growth_data[dropped, wave] <- NA
}
colMeans(is.na(growth_data[paste0("t", 1:4)]))

## ----growth-------------------------------------------------------------------
linear <- growth("i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
                  s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4",
  growth_data,
  missing = "fiml"
)
basis <- growth("i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
                 s =~ 0*t1 + 1*t2 + t3 + t4",
  growth_data,
  missing = "fiml"
)
growth_test <- pvalues_nested(linear, basis, tests = c("SB", "SS", "PALL"))
growth_test
coef(basis)[c("s=~t3", "s=~t4")]

## ----exact, error = TRUE------------------------------------------------------
try({
pvalues_nested(linear, basis, A.method = "exact")
})

## ----groups-------------------------------------------------------------------
configural <- cfa(model, fiml_data, group = "school", missing = "fiml")
equal_loadings <- cfa(model, fiml_data,
  group = "school", group.equal = "loadings", missing = "fiml"
)
pvalues_nested(equal_loadings, configural)

