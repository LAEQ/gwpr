# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
library(plyr)
library(dplyr)
library(tidyr)
library(foreach)
library(iterators)
library(doParallel)
registerDoParallel(cores=2)

data_preparation <- function(dataset, formula, id){
  if( !is.formula(formula)) stop("You must provide a valid formula")

  formula_vars <- all.vars(formula)
  vars <- c(formula_vars, id)

  if( ! all(vars %in% colnames(dataset))) stop("Your formula cannot match the dataset")

  dataAVG <- dataset %>%
    select(vars) %>%
    group_by(.dots = id) %>%
    summarise_at(formula_vars, mean)
}

bandwidth_optimisation <- function(formula, data, dmdmat) {
  lowerBw <- 50
  upperBW <- 200
  step    <- 10

  bandwith_list = seq(lowerBw, upperBW, step)
  CVsMat <- matrix(NA, nrow = length(bandwith_list), ncol = 2)
  colnames(CVsMat) <- c("Bandwidth", "CV")
  ResidCV2 <- c()
  adaptive <- TRUE


  foreach(bandwidth = bandwith_list, i = icount()) %do% {
    for(j in 1:672){
      W.j <- gw.weight(as.numeric(DMAT[j,]), bw = bandwidth, kernel = "bisquare", adaptive = adaptive)
      W.j[j] <- 0
      dataAVG$wgt <- W.j
      lm.j <- lm(formula = Equation, data= dataAVG, weights = wgt)
      ResidCV2[j] <- (lm.j$residuals[j])^2
    }

    dataAVG$wgt <- NULL
    CVsMat[i,1] <- bandwidth
    CVsMat[i,2] <- sum(ResidCV2)
  }

  BwOpt <- subset(CVsMat, CVsMat[,2]==min(CVsMat[,2], na.rm=TRUE))[1,1]
  Opt <- subset(CVsMat, CVsMat[,2]==min(CVsMat[,2], na.rm=TRUE))[1,2]


  return(2+2)
}
