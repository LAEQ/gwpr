# Context("Bandwith options")

library(testthat)

test_that("Bandwith option: Version avec donnees moyennes / adaptive bw ", {
  load(file = '../../inst/Data.Rdata')
  USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
  data <- merge(USStates@data, Produc, by="state", all=T)

  dMat <- GWmodel::gw.dist(sp::coordinates(USStates),  p=2, longlat=F)
  Equation <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp

  bwAVG.A <- bw.avg(formula=Equation, data=data, SDF=USStates, index=c("id","year"), approach="AICc",  kernel="bisquare", adaptive=T, p=2, longlat=F, dMat=dMat)

  expect_equal(bwAVG.A, 37)
})

# test_that("Bandwith option: Version avec donnees moyennes / fixed bw", {
#   load(file = '../../data/Data.Rdata')
#   USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
#   data <- merge(USStates@data, Produc, by="state", all=T)
#
#   dMat <- GWmodel::gw.dist(sp::coordinates(USStates), p=2, longlat=F)
#   Equation <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
#
#   bwAVG.F <- bw.avg(formula=Equation, data=data, SDF=USStates, index=c("id","year"), approach="AICc",  kernel="bisquare", adaptive=F, p=2, longlat=F, dMat=dMat)
#
#   expect_equal(bwAVG.F,2077193.92924246)
# })
