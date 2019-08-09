context("Geographical weight layer panel regression")

library(testthat)

source("../../R/bandwith_option.R")

test_that("gwlpr with us data", {
  load(file = "../../data/us_data/Data.Rdata")
  load(file = "../../data/us_data/gwlpr_expected.rda")
  formula <- gsp ~ pcap + pc + emp + unemp
  sequence <- seq(1, 48, 1)
  id <- "state"
  index <- c("state", "year")

  result <- gwlpr(SpatialDataFrame = USStates, dataset = Produc, formula = formula, id = id, index = index, sequence = sequence)

  expect_equal(result[[1]][[1]], LocalR2Mat)
  expect_equal(result[[1]][[2]], CoefsMat)
  expect_equal(result[[1]][[3]], SEsMat)
  expect_equal(result[[1]][[4]], TVsMat)
  expect_equal(all.equal(names, names(result[[2]]), order = FALSE), TRUE)
  expect_equal(result[[3]], R2)
})
