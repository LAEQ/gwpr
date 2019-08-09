context("Geographical weight panel regression")

library(testthat)

source("../../R/bandwith_option.R")

test_that("Bandwith options: 1", {
  load(file = "../../data/us_data/Data.Rdata")
  load(file = "../../data/us_data/us_dmat.rda")

  Formula <- gsp ~ pcap + pc + emp + unemp
  effect  <- "individual"
  model   <- "within"
  index   <- c("state", "year")

  result <- gpr(Formula, Produc, effect, model, index)

  load(file = "../../data/us_data/gpr/gpr_expected.rda")
  coef_expected <- c(pcap = -0.2555487, pc = 0.1633191, emp = 34.9584965,unemp = -209.3940985)

  diff_coef <- setdiff(result$coefficients, expected_gpr$coefficients)
  diff_vcov <- setdiff(result$vcov, expected_gpr$vcov)

  expect_equal(length(diff_coef), 0)
  expect_equal(length(diff_vcov), 0)
})




