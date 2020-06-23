# Context("Bandwidth options")

library(testthat)

test_that("Bandwith option: averages with adaptive bandwith", {
  Produc <- readRDS(system.file("Produc.rds", package = "gwpr"))
  USStates <- readRDS(system.file("USStates.rds", package = "gwpr"))
  USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
  data <- merge(USStates@data, Produc, by="state", all=T)

  dMat <- GWmodel::gw.dist(sp::coordinates(USStates), p=2, longlat=F)
  Equation <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp

  bwAVG.A <- bw.avg(formula=Equation, data=data, SDF=USStates, index=c("id","year"), approach="AICc",
                    kernel="bisquare", adaptive=T, p=2, longlat=F, dMat=dMat)
  expect_equal(bwAVG.A, 37)
})


