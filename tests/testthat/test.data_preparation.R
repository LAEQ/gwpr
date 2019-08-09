context("Data preparation")

library(testthat)

source("../../R/bandwith_option.R")

test_that("Datapreparation: fails if formula is invalid", {
  expect_that(data_preparation(1, 1, 1), throws_error())
})

test_that("Datapreparation: fails if dataframe and formula are incompatible", {
  formula <-  y ~ x + b
  load(file = "../../data/us_data/Data.Rdata")
  expect_that(data_preparation(PanelMtl, formula, "CTNAME86"), throws_error())
})

test_that("Data preparation: valid test", {
  Equation <- gsp ~ pcap + pc + emp + unemp
  load(file = "../../data/us_data/Data.Rdata")

  result <- data_preparation(Produc, Equation, "state")

  load(file="../../data/us_data/dataAVG_expected.rda")

  # https://github.com/tidyverse/dplyr/issues/2751
  expect_equal(as.data.frame(result[,-1]), as.data.frame(dataAVG_expected))
})
