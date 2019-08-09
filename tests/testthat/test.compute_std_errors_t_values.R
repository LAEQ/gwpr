context("Compute standard errors and T values")

library(testthat)

source("../../R/bandwith_option.R")

test_that("compute_std_errors_T_values", {
  load(file = "../../data/us_data/qx_qy_expected.rda")
  load(file = "../../data/us_data/gwpr_expected.rda")
  load(file = "../../data/us_data/SEsMat_TVsMat_expected.rda")

  qx_qy <- new("QXQY", QX = QX, QY = Qy, n = n, t = t, k = k, x = x, y = y)
  gwpr  <- new('GWPR', ListC = ListC, CoefsMat = CoefsMat, wmat = wmat, HatMat = HatMat, yHat = yHat, Resid = Resid)
  result <- compute_std_errors_T_values(qx_qy, gwpr)

  expect_equal(result$SEsMat, SEsMat)
  expect_equal(result$TVsMat, TVsMat)
})
