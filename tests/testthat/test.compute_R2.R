context("Compute R2")

library(testthat)

source("../../R/bandwith_option.R")

test_that("Computation of LocalR2Matrice", {
  load(file = "../../data/us_data/qx_qy_expected.rda")
  load(file = "../../data/us_data/gwpr_expected.rda")
  load(file = "../../data/us_data/R2_expected.rda")

  qx_qy  <- new("QXQY", QX = QX, QY = Qy, n = n, t = t, k = k, x = x, y = y)
  gwpr   <- new('GWPR', ListC = ListC, CoefsMat = CoefsMat, wmat = wmat, HatMat = HatMat, yHat = yHat, Resid = Resid)
  result <- compute_R2(qx_qy, gwpr)

  expect_equal(result, R2)
})
