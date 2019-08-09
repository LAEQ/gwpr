context("Compute dmat")

library(testthat)

source("../../R/bandwith_option.R")

test_that("Bandwith options: 1", {
  load(file = "../../data/us_data/Data.Rdata")
  load(file = "../../data/us_data/qx_qy_expected.rda")

  equation <- gsp ~ pcap + pc + emp + unemp
  index <- c("state", "year")
  effect <- "individual"
  model <- "within"

  result_QX_QY <- compute_QX_QY(Produc, equation, index, model, effect)

  expect_equal(result_QX_QY@k, 4)
  expect_equal(result_QX_QY@n, 48)
  expect_equal(result_QX_QY@t, 17)
  expect_equal(result_QX_QY@QX, QX)
  expect_equal(result_QX_QY@QY, Qy)
})




