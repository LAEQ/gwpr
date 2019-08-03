context("Data preparation")

library(testthat)

test_that("Datapreparation: fails if formula is invalid", {
  expect_that(data_preparation(1, 1, 1), throws_error())
})

test_that("Datapreparation: fails if dataframe and formula are incompatible", {
  formula <-  y ~ x + b
  load(file = "../../data/PanelMtl.rda")
  expect_that(data_preparation(PanelMtl, formula, "CTNAME86"), throws_error())
})

test_that("Data preparation: valid test", {
  Equation <- Y_FR ~ Chomag + FaMono + FaibSc + ImgRec + P65 + Menag1
  load(file = "../../data/PanelMtl.rda")

  result <- data_preparation(PanelMtl, Equation, "CTNAME86")
  load(file = "../../data/data_preparation/expected.rda")

  expect_equal(result, expected)
})
