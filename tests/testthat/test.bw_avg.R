# Context("Bandwith options")

library(testthat)

test_that("Bandwith option: Version avec donnees moyennes / adaptive bw ", {
  # path <- file.path(system.file("data", package = "gwpr"), "Data.rda")

  test <- readRDS(system.file("test", package = "gwpr"))

  expect_equal(test, 1)


})

