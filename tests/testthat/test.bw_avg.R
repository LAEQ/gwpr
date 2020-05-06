# Context("Bandwith options")

library(testthat)

test_that("Bandwith option: Version avec donnees moyennes / adaptive bw ", {
  path <- file.path(getwd(), "..","..", "data", "Data.rda")
  print(path)

  # base::load(file = path)
  # USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
  # data <- merge(USStates@data, Produc, by="state", all=T)
  #
  # dMat <- GWmodel::gw.dist(sp::coordinates(USStates),  p=2, longlat=F)
  # Equation <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
  #
  # bwAVG.A <- bw.avg(formula=Equation, data=data, SDF=USStates, index=c("id","year"), approach="AICc",  kernel="bisquare", adaptive=T, p=2, longlat=F, dMat=dMat)
  #
  # expect_equal(bwAVG.A, 37)

  expect_equal(path, "test/path")

})

