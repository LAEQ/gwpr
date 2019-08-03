# context("Bandwith options")

library(testthat)


test_that("Bandwith options: 1", {
  Equation <- Y_FR ~ Chomag + FaMono + FaibSc + ImgRec + P65 + Menag1
  load(file = "../../data/data_preparation/data_prepared.rda")
  load(file = "../../data/bandwith_option/dmat.rda")

  result <- bandwidth_optimisation(Equation, data_prepared, dmat)
  expect_equal(result, 5)
})
