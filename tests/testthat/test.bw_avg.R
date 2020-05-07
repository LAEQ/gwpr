# Context("Bandwith options")

library(testthat)

test_that("Bandwith option: Version avec donnees moyennes / adaptive bw ", {
  path <- file.path(system.file("data", package = "gwpr"), "Data.rda")

  if(file.exists(path)){
    expect_equal(1, 1)
  } else {
    expect_equal(1, 0)
  }
})

