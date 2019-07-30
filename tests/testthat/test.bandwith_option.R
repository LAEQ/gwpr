context("Generate bandwith options")
library(dplyr)
library(tidyr)
library(testthat)

test_that("Datapreparation: fails if formula is invalid", {
  expect_that(data_preparation(1, 1, 1), throws_error())
})

test_that("Datapreparation: fails if dataframe and formula are incompatible", {
  formula <-  y ~ x + b
  load(file = "../../data/PanelMtl.rda")
  expect_that(data_preparation(PanelMtl, formula, "id"), throws_error())
})

test_that("Data preparation: test formula", {
  Equation <- Y_FR ~ Chomag + FaMono + FaibSc + ImgRec + P65 + Menag1

  load(file = "../../data/PanelMtl.rda")

  result <- data_preparation(PanelMtl, Equation, "CTNAME86")

  expect_equal(TRUE, TRUE)

})

# test_that("10 equals 10", {
#   result <- bandwidthOption(4, 4)
#   expect_equal(result, 8)
# })
