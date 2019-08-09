context("Compute dmat")

library(testthat)

source("../../R/bandwith_option.R")

test_that("Bandwith options: 1", {
  load(file = "../../data/us_data/Data.Rdata")
  load(file = "../../data/us_data/dmat_expected.rda")

  USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))

  result <- compute_dmat(coordinates(USStates), p = 2, longlat = F)

  expect_equal(result, dmat_expected)
})




