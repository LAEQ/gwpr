context("Geographical weight panel regression")

library(testthat)

source("../../R/bandwith_option.R")

# https://github.com/r-lib/testthat/issues/473

test_that("Bandwith options: 1", {
  load(file = "../../data/us_data/Data.Rdata")
  load(file = "../../data/us_data/us_dmat.rda")
  USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
  Formula <- gsp ~ pcap + pc + emp + unemp
  bwselection <- c(Bandwidth=90)

  result <- gwpr(USStates, Produc, Formula, c("state", "year"),
                 bwselection, us_dmat, kernel = "bisquare", effect = "individual", model = "within", adaptive = T)
  load(file = "../../data/us_data/gwpr/r2_expected.rda")

  expect_equal(result[3], r2_expected)
})




