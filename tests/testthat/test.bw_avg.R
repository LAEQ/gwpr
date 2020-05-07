# Context("Bandwith options")

library(testthat)

test_that("Bandwith option: Version avec donnees moyennes / adaptive bw ", {
  # path <- file.path(system.file("data", package = "gwpr"), "Data.rda")

  test <- readRDS(file.path("..", "..", "data", "test"))

  expect_equal(test, 1)


})

