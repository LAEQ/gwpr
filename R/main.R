source("R/bandwith_option.R")



load(file = "data/us_data/Data.Rdata")
load(file = "data/us_data/us_dmat.rda")
USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
Formula <- gsp ~ pcap + pc + emp + unemp
bwselection <- c(Bandwidth=90)



result <- gwpr(USStates, Produc, Formula, c("state", "year"), bwselection, us_dmat, kernel = "bisquare", adaptive = T)
