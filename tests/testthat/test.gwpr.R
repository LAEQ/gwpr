context("Geographical weight panel regression")

library(testthat)

source("../../R/bandwith_option.R")

test_that("Bandwith options: 1", {
  load(file = "../../data/us_data/Data.Rdata")
  load(file = "../../data/us_data/us_dmat.rda")
  Formula <- gsp ~ pcap + pc + emp + unemp
  bwselection <- c(Bandwidth=90)

  result <- gwpr(Produc, Formula, c("state", "year"), bwselection, us_dmat)


  expect_equal(result, TRUE)
})




