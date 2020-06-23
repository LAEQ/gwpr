# Context("GWPLM: adaptive bw")

library(testthat)


test_that("Bandwith option:Version avec donnees moyennes / adaptive bw ", {
  Produc <- readRDS(system.file("Produc.rds", package = "gwpr"))
  USStates <- readRDS(system.file("USStates.rds", package = "gwpr"))

  USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
  data <- merge(USStates@data, Produc, by="state", all=T)

  dMat <- GWmodel::gw.dist(sp::coordinates(USStates), p=2, longlat=F)
  Equation <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp

  bwCV.A <- readRDS(system.file("bwCV.A.rds", package = "gwpr"))

  result <- gwplm(SpDF=USStates, data=data, index=c("id", "year"),
                     formula=Equation, bw=bwCV.A, kernel="bisquare",
                     adaptive=T, effect="individual", model="within",
                     dMat=dMat)


  USgwplm_A <- readRDS(system.file("USgwplm.A.rds", package = "gwpr"))

  expect_equal(result,  USgwplm_A)
})


test_that("Bandwith option:Version avec donnees moyennes / adaptive bw ", {
  Produc <- readRDS(system.file("Produc.rds", package = "gwpr"))
  USStates <- readRDS(system.file("USStates.rds", package = "gwpr"))

  USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
  data <- merge(USStates@data, Produc, by="state", all=T)

  dMat <- GWmodel::gw.dist(sp::coordinates(USStates), p=2, longlat=F)
  Equation <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp

  bwCV.F <- readRDS(system.file("bwCV.F.rds", package = "gwpr"))

  result <- gwplm(SpDF=USStates, data=data, index=c("id", "year"),
                     formula=Equation, bw=bwCV.F, kernel="bisquare",
                     adaptive=F, effect="individual", model="within",
                     dMat=dMat)


  USgwplm_B <- readRDS(system.file("USgwplm.B.rds", package = "gwpr"))

  expect_equal(result, USgwplm_B)
})
