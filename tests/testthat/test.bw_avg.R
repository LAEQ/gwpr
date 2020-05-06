# Context("Bandwith options")

library(testthat)

test_that("Bandwith option: Version avec donnees moyennes / adaptive bw ", {
  wd <- getwd()
  path <- file.path(getwd(), "..","..", "data", "Data.rda")
  if(file.exists(path)){
    expect_equal(1, 1)
  } else {
    expect_equal(1, 0)
  }
})

