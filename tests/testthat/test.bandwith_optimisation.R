context("Bandwith optimisation")

library(testthat)

source("../../R/bandwith_option.R")

# test_that("Bandwith optimisation: 1", {
#   Equation <- Y_FR ~ Chomag + FaMono + FaibSc + ImgRec + P65 + Menag1
#   load(file = "../../data/mtl_data/data_prepared.rda")
#   load(file = "../../data/mtl_data/dmat.rda")
#
#   sequence_bw <- seq(50, 200, 10)
#
#   result <- bandwidth_optimisation(Equation, data_prepared, dmat, sequence_bw)
#   expect_equal(as.numeric(result), 90)
# })

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
#
#
test_that("Bandwith optimisation: Product USA", {
  load(file = "../../data/us_data/dataAVG_expected.rda")
  load(file = "../../data/us_data/us_dmat.rda")
  load(file = "../../data/us_data/bandwidth_expected.rda")

  Equation <- gsp ~ pcap + pc + emp + unemp
  sequence_bw <- seq(1, 48, 1)

  result <- bandwidth_optimisation(Equation, dataAVG_expected, us_dmat, sequence_bw)

  expect_equal(result[[1]], BwOpt)
  expect_equal(result[[2]], Opt)
})
