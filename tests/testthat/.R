context("Compute gwpr")

library(testthat)

source("../../R/bandwith_option.R")

test_that("Bandwith options: 1", {
  load(file = "../../data/us_data/Data.Rdata")
  load(file = "../../data/us_data/qx_qy_expected.rda")
  load(file = "../../data/us_data/dmat_expected.rda")
  load(file = "../../data/us_data/bandwith_expected.rda")
  load(file = "../../data/us_data/gwpr_expected.rda")

  equation <- gsp ~ pcap + pc + emp + unemp
  index <- c("state", "year")
  kernel <- "bisquare"
  adaptive <- TRUE

  result_QX_QY <- new("QXQY", QX = QX, QY = Qy, n = n, t = t, k = k, x = x, y = y)

  result_gwpr <- compute_gwpr(result_QX_QY, dmat_expected, bandwidth_expected, kernel, adaptive)


  expect_equal(result_gwpr@CoefsMat, CoefsMat)
  expect_equal(result_gwpr@HatMat, HatMat)
  expect_equal(result_gwpr@ListC, ListC)
  expect_equal(result_gwpr@wmat, wmat)
  expect_equal(result_gwpr@Resid, Resid)
  expect_equal(result_gwpr@yHat, yHat)
})
