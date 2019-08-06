# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#  http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#  Install Package:      'Ctrl + Shift + B'
#  Check Package:       'Ctrl + Shift + E'
#  Test Package:       'Ctrl + Shift + T'
library(plyr)
library(dplyr)
library(tidyr)
library(foreach)
library(iterators)
library(doParallel)
library(GWmodel)
library(plm)


registerDoParallel(cores=2)

#' Prepare the data ...
#'
#' @param dataset
#' @param formula
#' @param id
#'
#' @export
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

#' Bandwith optimization
#'
#' Desc to come
#'
#' @param formula (desc to come)
#' @param data (desc to come)
#' @param dmat (desc to come)
#' @param sequence (desc to come)
#'
#' @export
#'
bandwidth_optimisation <- function(formula, data, dmat, sequence) {
  CVsMat <- matrix(NA, nrow = length(sequence), ncol = 2)
  colnames(CVsMat) <- c("Bandwidth", "CV")
  ResidCV2 <- c()
  adaptive <- TRUE


  foreach(bandwidth = sequence, i = icount()) %do% {
    for(j in 1:length(data)){
      W.j <- gw.weight(as.numeric(dmat[j,]), bw = bandwidth, kernel = "bisquare", adaptive = adaptive)
      W.j[j] <- 0
      data$wgt <- W.j
      lm.j <- lm(formula = formula, data = data, weights = wgt)
      ResidCV2[j] <- (lm.j$residuals[j])^2
    }

    data$wgt <- NULL
    CVsMat[i,1] <- bandwidth
    CVsMat[i,2] <- sum(ResidCV2)
  }

  BwOpt <- subset(CVsMat, CVsMat[,2]==min(CVsMat[,2], na.rm=TRUE))[1,1]

  return(BwOpt)
}

# gwpr(SpDF=USStates, data=Produc, index=c("state", "year"), formula=Equation,
#      bw = bwselection, kernel="bisquare", adaptive=T, effect="individual", model="within", dmat=dmat)

#' Geopgraphical weight panel regression
#'
#' This function implements basic geographically weighted panel regression (GWPR)
#'
#' @param data      A data frame for the panel data
#' @param formula   A regression model formula6
#' @param index     A vector of the indexes
#' @param bandwith  A named vector e.g c(Bandwith=90)
#' @param dmat
#'
#' @export
gwpr <- function(dataset, formula, indexes, bandwidth, dmat){
  if( ! all(indexes %in% colnames(dataset))) stop("Indexes are missing in the dataframe")

  formula_vars <- all.vars(formula)
  if( ! all(formula_vars %in% colnames(dataset))) stop("Your formula does not match the dataset")


  #Panel data preparation
  Pdata <- pdata.frame(dataset, index = indexes, drop.index = FALSE, row.names = FALSE, stringsAsFactors = default.stringsAsFactors())


  # Namin y, x ???
  y <- model.response(model.frame(formula = formula, data = Pdata))
  x <- model.matrix(formula, data = Pdata)



  #To parametarise
  kernel    <- "bisquare"
  adaptive  <- TRUE
  effect    <- "individual"
  model     <- "within"

  if(model == "within"){
    x <- x[, -1]
  }

  # What if we passe more than 2 index
  # naming n, t, k
  n <- length(unique(Pdata[, indexes[1]]))
  t <- length(unique(Pdata[, indexes[2]]))
  K <- ncol(x)

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Following section adapted from spfeml
  # (Giovanni Millo, Sep 2017, https://github.com/cran/splm/blob/master/R/spfeml.R)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  #Check whether the panel is balanced
  balanced <- n * t == dim(Pdata)[1]
  if(balanced == FALSE){
    stop("Estimation method unavailable for unbalanced panels")
  }

  # First t observations
  indic   <- seq(1, n)
  inde    <- as.numeric(rep(indic, each = t))

  # Observations: 1, t + 1, 2t + 1 ...
  indic1  <- seq(1, t)
  inde1   <- rep(indic1, n) ####takes observations 1,  t+1, 2t+1...

  # Within transformation: Calculating the mean of Y and X for each time period
  if (effect %in% c("time", "twoways")){
    ytms<-tapply(y,inde,mean)
    tpms<-function(q) tapply(q,inde,mean)
    xtms<-apply(x,2,tpms)
    ytm<-rep(ytms,each=t)
    xtm<-matrix(NA,(n*t),K)
    for (k in 1:K) xtm[,k]<-rep(xtms[,k],each=t)
  }

  # Calculating the mean of Y and X for each individual
  if (effect %in% c("individual", "twoways")){
    ysms <- tapply(y,inde1,mean)
    spms <- function(q) tapply(q,inde1,mean)
    xsms <- apply(x,2,spms)
    ysm  <- rep(ysms,n)
    xsm  <- matrix(NA,(n*t),K)

    for (k in 1:K){
      xsm[,k]<-rep(xsms[,k],n)
    }
  }


  # Calculating the within-transformed Y and X (Qy and QX) according to the selected model
  # Issue: if effect not time, individual, two_ways ???
  if (effect == "time"){
    Qy <- y - ytm
    QX <- x - xtm
  }

  if(effect == "individual"){
    Qy <- y - ysm
    QX <- x - xsm
  }

  if (effect=="twoways"){
    Qy <- y - ysm - ytm + rep(mean(y),n*t)
    xmm <- matrix(NA,n*t,K)
    for (k in 1:K){
      xmm[,k] <- rep(mean(x[,k]), n * t)
    }
    QX <- x - xsm - xtm + xmm
  }


  # Computation of the GWPR model
  ListC     <- list()
  CoefsMat  <- matrix(NA, nrow = n, ncol = k)
  wmat      <- matrix(NA, nrow = n, ncol = n*t)
  HatMat    <- matrix(NA, nrow = n*t, ncol = n*t)
  yHat      <- c()
  Resid     <- c()

  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n){
    dist.vi <- dmat[i,]
    W.i <- gw.weight(dist.vi, bw = bandwidth, kernel = kernel, adaptive = adaptive)
    W.i <- diag(W.i)
    W.i <- kronecker(W.i,diag(t))
    wmat[i,] <- diag(W.i)
    QXw.i <- QX*wmat[i,]
    ListC[[i]] <- (solve(t(QXw.i)%*%QX))%*%t(QXw.i)
    CoefsMat[i,] <- ListC[[i]]%*%Qy
    for(j in 1:t){
      HatMat[(i-1)*t+j,] <- QX[(i-1)*t+j,]%*%ListC[[i]]
      yHat.temp <- CoefsMat[i,]%*%t(QX)
      yHat[(i-1)*t+j] <- yHat.temp[(i-1)*t+j]
      Resid[(i-1)*t+j] <- Qy[(i-1)*t+j] - yHat[(i-1)*t+j]
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  return(FALSE)
}


#' Return Linear Models for Panel Data
#'
#' @param formula   a formula
#' @param data      a data.frame
#' @param effect    effect introduced in the model c("individual", "time", "twoways", "nested")
#' @param model     one of c("pooling", "within", "between", "random", "fd", "ht")
#' @param index     list of indexes
#'
#' @export
gpr <- function(formula, data, effect, model, index){
  PanelModel <- plm(formula = formula, data = data, effect = effect, model = model, index = indexes)
  return(PanelModel)
}
