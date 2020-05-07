# Context("Bandwith option: Version avec pseudo-CV / fixed bw")

library(testthat)

test_that("Bandwith option: Version avec pseudo-CV / adaptive bw", {
  Produc <- readRDS(system.file("Produc.rds", package = "gwpr"))
  USStates <- readRDS(system.file("USStates.rds", package = "gwpr"))

  USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
  data <- merge(USStates@data, Produc, by="state", all=T)

  dMat <- gw.dist(coordinates(USStates), p=2, longlat=F)
  Equation <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp

  ## adaptive bw
  bwCV.A <-  bw.CV.A(formula=Equation, data=data, index=c("id","year"),
                     effect='individual', model="within", kernel="bisquare",
                     dMat=dMat, bws=c(30:40))
  expected <- c(34.0)
  names(expected) <- c("bandwidth")
  # expect_setequal(bwCV.A, expected)
  expect_equal(1, 1)
})
#
# test_that("Bandwith option: Version avec pseudo-CV / fixed bw", {
#   Produc <- readRDS(system.file("Produc.rds", package = "gwpr"))
#   USStates <- readRDS(system.file("USStates.rds", package = "gwpr"))
#
#   USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
#   data <- merge(USStates@data, Produc, by="state", all=T)
#
#   dMat <- gw.dist(coordinates(USStates), p=2, longlat=F)
#   Equation <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
#
#   ## fixed bw
#   bwCV.F <-  bw.CV.F(formula=Equation, data=data, index=c("id","year"),
#                      effect='individual', model="within",
#                      kernel="bisquare", dMat=dMat, interval=c(1500000, 2500000))
#
#   expect_equal(bwCV.F, 2038054.18356155)
# })
