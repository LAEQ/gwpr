context("Bandwith optimisation")

library(testthat)

source("../../R/bandwith_option.R")

# test_that("Bandwith optimisation: 1", {
#   Equation <- Y_FR ~ Chomag + FaMono + FaibSc + ImgRec + P65 + Menag1
#   load(file = "../../data/data_preparation/data_prepared.rda")
#   load(file = "../../data/bandwith_option/dmat.rda")
#
#   sequence_bw <- seq(50, 200, 10)
#
#   result <- bandwidth_optimisation(Equation, data_prepared, dmat, sequence_bw)
#   expect_equal(as.numeric(result), 90)
# })
#
#
# test_that("Bandwith optimisation: 1", {
#   Equation <- Y_FR ~ Chomag + FaMono + FaibSc + ImgRec + P65 + Menag1
#   load(file = "../../data/data_preparation/data_prepared.rda")
#   load(file = "../../data/bandwith_option/dmat.rda")
#
#   sequence_bw <- seq(85, 95, 1)
#
#   result <- bandwidth_optimisation(Equation, data_prepared, dmat, sequence_bw)
#   expect_equal(as.numeric(result), 90)
# })


test_that("Bandwith optimisation: Product USA", {
  load(file = "../../data/us_data/Data.Rdata")
  load(file = "../../data/us_data/data_preparation/produc_result.rda")
  load(file = "../../data/us_data/us_dmat.rda")

  Equation <- gsp ~ pcap + pc + emp + unemp
  sequence_bw <- seq(1, 48, 1)

  result <- bandwidth_optimisation(Equation, expected, us_dmat, sequence_bw)
  result_expected = c(Bandwidth=6)
  expect_equal(result, result_expected)
})
