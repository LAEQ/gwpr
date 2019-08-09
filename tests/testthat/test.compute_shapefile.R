context("Compute shape file")

library(testthat)

source("../../R/bandwith_option.R")

test_that("compute_shapefile", {
  load(file = "../../data/us_data/Data.Rdata")
  load(file = "../../data/us_data/qx_qy_expected.rda")
  load(file = "../../data/us_data/gwpr_expected.rda")
  load(file = "../../data/us_data/R2_expected.rda")
  load(file = "../../data/us_data/LocalR2Mat_expected.rda")
  load(file = "../../data/us_data/SEsMat_TVsMat_expected.rda")

  USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
  result_shapefile <- new('ShapeFile')
  result_shapefile@LocalR2Mat <- LocalR2Mat
  result_shapefile@R2 <- R2
  result_shapefile@SEsMat <- SEsMat
  result_shapefile@TVsMat <- TVsMat

  qx_qy  <- new("QXQY", QX = QX, QY = Qy, n = n, t = t, k = k, x = x, y = y)
  gwpr   <- new('GWPR', ListC = ListC, CoefsMat = CoefsMat, wmat = wmat, HatMat = HatMat, yHat = yHat, Resid = Resid)
  result <- compute_shapefile(USStates, result_shapefile, gwpr, qx_qy)

  colnames(CoefsMat) <- paste(colnames(x), "_Coef", sep = "")
  colnames(SEsMat) <- paste(colnames(x), "_SE", sep = "")
  colnames(TVsMat) <- paste(colnames(x), "_TV", sep = "")
  colnames(LocalR2Mat) <- c("Local_RSquared")

  expect_equal(length(result[[1]]), 4)
  expect_equal(result[[1]][[1]], LocalR2Mat)
  expect_equal(result[[1]][[2]], CoefsMat)
  expect_equal(result[[1]][[3]], SEsMat)
  expect_equal(result[[1]][[4]], TVsMat)
  expect_equal(result[[3]], R2)
})
