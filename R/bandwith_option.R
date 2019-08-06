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
#' @param SpDF      A spacial polygone DataFrame
#' @param data      A data frame for the panel data
#' @param formula   A regression model formula6
#' @param index     A vector of the indexes
#' @param bandwith  A named vector e.g c(Bandwith=90)
#' @param dmat
#' @param kernel    function chosen as follow: gaussian, exponential, bisquare, tricube, boxcar
#'
#' @export
gwpr <- function(SpDF, dataset, formula, indexes, bandwidth, dmat, kernel, effect, model, adaptive = F){
  if( ! all(indexes %in% colnames(dataset))) stop("Indexes are missing in the dataframe")
  if( ! all(all.vars(formula) %in% colnames(dataset))) stop("Your formula does not match the dataset")
  if( ! effect %in% c("individual", "twoways", "time")) stop(paste(effect, " is not supported. (invidual, twoways, time)"))

  #########################################################################################
  #Panel data preparation
  Pdata <- pdata.frame(dataset, index = indexes, drop.index = FALSE, row.names = FALSE, stringsAsFactors = default.stringsAsFactors())

  # Namin y, x ???
  y <- model.response(model.frame(formula = formula, data = Pdata))
  x <- model.matrix(formula, data = Pdata)


  # To parametarise
  # kernel    <- "bisquare"
  # adaptive  <- TRUE
  # effect    <- "individual"
  # model     <- "within"

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
  #########################################################################################

  #########################################################################################
  # Validation ? if effect not in c(time, twoways, individual)
  #########################################################################################

  # Within transformation: Calculating the mean of Y and X for each time period (Y, X ???)
  if (effect %in% c("time", "twoways")){
    ytms  <- tapply(y, inde, mean)
    tpms  <- function(q) tapply(q, inde, mean)
    xtms  <- apply(x, 2, tpms)
    ytm   <- rep(ytms, each = t)
    xtm   <- matrix(NA, (n * t), K)

    for (k in 1:K){
      xtm[, k] <- rep(xtms[, k], each = t)
    }
  }

  # Calculating the mean of Y and X for each individual
  if (effect %in% c("individual", "twoways")){
    ysms <- tapply(y, inde1, mean)
    spms <- function(q) tapply(q, inde1, mean)
    xsms <- apply(x, 2, spms)
    ysm  <- rep(ysms, n)
    xsm  <- matrix(NA, (n * t), K)

    for (k in 1:K){
      xsm[,k]<-rep(xsms[, k], n)
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
    Qy <- y - ysm - ytm + rep(mean(y), n * t)
    xmm <- matrix(NA, n*t, K)
    for (k in 1:K){
      xmm[,k] <- rep(mean(x[, k]), n * t)
    }
    QX <- x - xsm - xtm + xmm
  }



  #########################################################################################

  ######################################################
  # Computation of the GWPR model
  # @param dmat
  # @param bandwith
  # @param kernel
  # @param adaptive
  # @param QX
  # @param QY
  #
  # @return wmat, yHat, Resid, HatMat
  ListC     <- list()
  CoefsMat  <- matrix(NA, nrow = n, ncol = k)
  wmat      <- matrix(NA, nrow = n, ncol = n * t)
  HatMat    <- matrix(NA, nrow = n*t, ncol = n * t)
  yHat      <- c()
  Resid     <- c()

  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n){
    dist.vi       <- dmat[i,]
    W.i           <- gw.weight(dist.vi, bw = bandwidth, kernel = kernel, adaptive = adaptive)
    W.i           <- diag(W.i)
    W.i           <- kronecker(W.i, diag(t))
    wmat[i,]      <- diag(W.i)
    QXw.i         <- QX * wmat[i, ]
    ListC[[i]]    <- (solve(t(QXw.i) %*% QX)) %*% t(QXw.i)
    CoefsMat[i,]  <- ListC[[i]] %*% Qy

    for(j in 1:t){
      HatMat[(i - 1) * t + j, ] <- QX[(i - 1) * t + j, ]%*% ListC[[i]]
      yHat.temp <- CoefsMat[i,] %*% t(QX)
      yHat[(i - 1) * t + j] <- yHat.temp[(i-1) * t + j]
      Resid[(i - 1) * t + j] <- Qy[(i - 1) * t + j] - yHat[(i - 1) * t + j]
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  ######################################################


  ######################################################
  # Local R-Squared
  # @param wmat
  # @param Qy
  # @param yHat
  #
  # @return LocalR2Mat
  LocalR2Mat <- matrix(NA, nrow = n, ncol = 1)
  for (i in 1:n){
    TSSw.i <- t((Qy) * wmat[i,]) %*% (Qy)
    RSSw.i <- t((Qy - yHat) * wmat[i,]) %*% (Qy - yHat)
    LocalR2Mat[i,] <- (TSSw.i - RSSw.i) / TSSw.i
  }
  ######################################################

  ######################################################
  # Global R-Squared
  # @param yHat
  # @param Qy
  #
  # @return R2
  R2 <- cor(yHat, Qy) ^ 2
  ######################################################

  ######################################################
  # Standard errors and T values (matrices)
  # @param Resid
  # @param HtMat
  #
  # @return SesMat, TVsMat
  RSS <- t(Resid) %*% (Resid)
  v1 <- sum(diag(HatMat))
  v2 <- sum(HatMat ^ 2) # A ete verifie que c'est equivalent a trace(StS)

  # Selon l'equation adaptee au format panel de la GWR (IMPORTANT - is it fixed ?)
  sigma2 <- as.numeric(RSS / (n * t - 2 * v1 + v2))

  SEsMat <- matrix(nrow = n, ncol = k)
  TVsMat <- matrix(nrow = n, ncol = k)

  for (i in 1:n){
    C.i <- ListC[[i]]
    SEs.i <- c()
    SEs.i <- sqrt(diag(C.i %*% (t(C.i))) * sigma2)
    for (j in 1:k){
      SEsMat[i,j] <- SEs.i[j]
      TVsMat[i,j] <- CoefsMat[i, j] / SEsMat[i, j]
    }
  }
  ######################################################

  #Rename objects
  colnames(CoefsMat) <- paste(colnames(x), "_Coef", sep = "")
  colnames(SEsMat) <- paste(colnames(x), "_SE", sep = "")
  colnames(TVsMat) <- paste(colnames(x), "_TV", sep = "")
  colnames(LocalR2Mat) <- c("Local_RSquared")

  env <- new.env(parent=globalenv())
  assign("Coefficients matrix", CoefsMat, envir = env)
  assign("Std.Errors matrix", SEsMat, envir = env)
  assign("T-values matrix", TVsMat, envir = env)
  assign("Local R-squared", LocalR2Mat, envir = env)
  assign("Global R-squared", R2, envir = env)
  assign("Spatial panel data frame", SpDF, envir = env)

  # Creating the shapefile with the GWPR results(coefficients, standard errors and T values)
  listMat <- list(LocalR2Mat, CoefsMat, SEsMat, TVsMat)
  for(i in 1:length(listMat)){
    newVars <- colnames(listMat[[i]])

    for(j in 1:length(newVars)){
      SpDF$x <- listMat[[i]][, j]
      last <- length(names(SpDF))
      names(SpDF)[last] <- newVars[j]
    }
  }

  res <- list(listMat, SpDF, R2)

  return(res)
}

#' Geopgraphical weight panel regression
#'
#' This function implements basic geographically weighted panel regression (GWPR)
#'
#' @param SpDF      A spacial polygone DataFrame
#' @param data      A data frame for the panel data
#' @param formula   A regression model formula6
#' @param index     A vector of the indexes
#' @param bandwith  A named vector e.g c(Bandwith=90)
#' @param dmat
#' @param kernel    function chosen as follow: gaussian, exponential, bisquare, tricube, boxcar
#'
#' @export
gwpr <- function(SpDF, dataset, formula, indexes, bandwidth, dmat, kernel, effect, model, adaptive = F){
  if( ! all(indexes %in% colnames(dataset))) stop("Indexes are missing in the dataframe")
  if( ! all(all.vars(formula) %in% colnames(dataset))) stop("Your formula does not match the dataset")
  if( ! effect %in% c("individual", "twoways", "time")) stop(paste(effect, " is not supported. (invidual, twoways, time)"))

  #########################################################################################
  #Panel data preparation
  Pdata <- pdata.frame(dataset, index = indexes, drop.index = FALSE, row.names = FALSE, stringsAsFactors = default.stringsAsFactors())

  # Namin y, x ???
  y <- model.response(model.frame(formula = formula, data = Pdata))
  x <- model.matrix(formula, data = Pdata)


  # To parametarise
  # kernel    <- "bisquare"
  # adaptive  <- TRUE
  # effect    <- "individual"
  # model     <- "within"

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
  #########################################################################################

  #########################################################################################
  # Validation ? if effect not in c(time, twoways, individual)
  #########################################################################################

  # Within transformation: Calculating the mean of Y and X for each time period (Y, X ???)
  if (effect %in% c("time", "twoways")){
    ytms  <- tapply(y, inde, mean)
    tpms  <- function(q) tapply(q, inde, mean)
    xtms  <- apply(x, 2, tpms)
    ytm   <- rep(ytms, each = t)
    xtm   <- matrix(NA, (n * t), K)

    for (k in 1:K){
      xtm[, k] <- rep(xtms[, k], each = t)
    }
  }

  # Calculating the mean of Y and X for each individual
  if (effect %in% c("individual", "twoways")){
    ysms <- tapply(y, inde1, mean)
    spms <- function(q) tapply(q, inde1, mean)
    xsms <- apply(x, 2, spms)
    ysm  <- rep(ysms, n)
    xsm  <- matrix(NA, (n * t), K)

    for (k in 1:K){
      xsm[,k]<-rep(xsms[, k], n)
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
    Qy <- y - ysm - ytm + rep(mean(y), n * t)
    xmm <- matrix(NA, n*t, K)
    for (k in 1:K){
      xmm[,k] <- rep(mean(x[, k]), n * t)
    }
    QX <- x - xsm - xtm + xmm
  }

  #########################################################################################

  ######################################################
  # Computation of the GWPR model
  # @param dmat
  # @param bandwith
  # @param kernel
  # @param adaptive
  # @param QX
  # @param QY
  #
  # @return wmat, yHat, Resid, HatMat
  ListC     <- list()
  CoefsMat  <- matrix(NA, nrow = n, ncol = k)
  wmat      <- matrix(NA, nrow = n, ncol = n * t)
  HatMat    <- matrix(NA, nrow = n*t, ncol = n * t)
  yHat      <- c()
  Resid     <- c()

  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n){
    dist.vi       <- dmat[i,]
    W.i           <- gw.weight(dist.vi, bw = bandwidth, kernel = kernel, adaptive = adaptive)
    W.i           <- diag(W.i)
    W.i           <- kronecker(W.i, diag(t))
    wmat[i,]      <- diag(W.i)
    QXw.i         <- QX * wmat[i, ]
    ListC[[i]]    <- (solve(t(QXw.i) %*% QX)) %*% t(QXw.i)
    CoefsMat[i,]  <- ListC[[i]] %*% Qy

    for(j in 1:t){
      HatMat[(i - 1) * t + j, ] <- QX[(i - 1) * t + j, ]%*% ListC[[i]]
      yHat.temp <- CoefsMat[i,] %*% t(QX)
      yHat[(i - 1) * t + j] <- yHat.temp[(i-1) * t + j]
      Resid[(i - 1) * t + j] <- Qy[(i - 1) * t + j] - yHat[(i - 1) * t + j]
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  ######################################################


  ######################################################
  # Local R-Squared
  # @param wmat
  # @param Qy
  # @param yHat
  #
  # @return LocalR2Mat
  LocalR2Mat <- matrix(NA, nrow = n, ncol = 1)
  for (i in 1:n){
    TSSw.i <- t((Qy) * wmat[i,]) %*% (Qy)
    RSSw.i <- t((Qy - yHat) * wmat[i,]) %*% (Qy - yHat)
    LocalR2Mat[i,] <- (TSSw.i - RSSw.i) / TSSw.i
  }
  ######################################################

  ######################################################
  # Global R-Squared
  # @param yHat
  # @param Qy
  #
  # @return R2
  R2 <- cor(yHat, Qy) ^ 2
  ######################################################

  ######################################################
  # Standard errors and T values (matrices)
  # @param Resid
  # @param HatMat
  #
  # @return SesMat, TVsMat
  RSS <- t(Resid) %*% (Resid)
  v1 <- sum(diag(HatMat))
  v2 <- sum(HatMat ^ 2) # A ete verifie que c'est equivalent a trace(StS)

  # Selon l'equation adaptee au format panel de la GWR (IMPORTANT - is it fixed ?)
  sigma2 <- as.numeric(RSS / (n * t - 2 * v1 + v2))

  SEsMat <- matrix(nrow = n, ncol = k)
  TVsMat <- matrix(nrow = n, ncol = k)

  for (i in 1:n){
    C.i <- ListC[[i]]
    SEs.i <- c()
    SEs.i <- sqrt(diag(C.i %*% (t(C.i))) * sigma2)
    for (j in 1:k){
      SEsMat[i,j] <- SEs.i[j]
      TVsMat[i,j] <- CoefsMat[i, j] / SEsMat[i, j]
    }
  }
  ######################################################

  #Rename objects
  colnames(CoefsMat) <- paste(colnames(x), "_Coef", sep = "")
  colnames(SEsMat) <- paste(colnames(x), "_SE", sep = "")
  colnames(TVsMat) <- paste(colnames(x), "_TV", sep = "")
  colnames(LocalR2Mat) <- c("Local_RSquared")

  env <- new.env(parent=globalenv())
  assign("Coefficients matrix", CoefsMat, envir = env)
  assign("Std.Errors matrix", SEsMat, envir = env)
  assign("T-values matrix", TVsMat, envir = env)
  assign("Local R-squared", LocalR2Mat, envir = env)
  assign("Global R-squared", R2, envir = env)
  assign("Spatial panel data frame", SpDF, envir = env)

  # Creating the shapefile with the GWPR results(coefficients, standard errors and T values)
  listMat <- list(LocalR2Mat, CoefsMat, SEsMat, TVsMat)
  for(i in 1:length(listMat)){
    newVars <- colnames(listMat[[i]])

    for(j in 1:length(newVars)){
      SpDF$x <- listMat[[i]][, j]
      last <- length(names(SpDF))
      names(SpDF)[last] <- newVars[j]
    }
  }

  res <- list(listMat, SpDF, R2)

  return(res)
}

calcul_QX_QY <- function(dataset, formula, index, model, effect){
  #########################################################################################
  #Panel data preparation
  Pdata <- pdata.frame(dataset, index = index, drop.index = FALSE, row.names = FALSE, stringsAsFactors = default.stringsAsFactors())

  # Namin y, x ???
  y <- model.response(model.frame(formula = formula, data = Pdata))
  x <- model.matrix(formula, data = Pdata)


  # To parametarise
  # kernel    <- "bisquare"
  # adaptive  <- TRUE
  # effect    <- "individual"
  # model     <- "within"

  if(model == "within"){
    x <- x[, -1]
  }

  # What if we passe more than 2 index
  # naming n, t, k
  n <- length(unique(Pdata[, index[1]]))
  t <- length(unique(Pdata[, index[2]]))
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
  #########################################################################################

  #########################################################################################
  # Validation ? if effect not in c(time, twoways, individual)
  #########################################################################################

  # Within transformation: Calculating the mean of Y and X for each time period (Y, X ???)
  if (effect %in% c("time", "twoways")){
    ytms  <- tapply(y, inde, mean)
    tpms  <- function(q) tapply(q, inde, mean)
    xtms  <- apply(x, 2, tpms)
    ytm   <- rep(ytms, each = t)
    xtm   <- matrix(NA, (n * t), K)

    for (k in 1:K){
      xtm[, k] <- rep(xtms[, k], each = t)
    }
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
    xmm <- matrix(NA, n*t, K)
    for (k in 1:K){
      xmm[,k] <- rep(mean(x[, k]), n * t)
    }
    QX <- x - xsm - xtm + xmm
  }

  return(data.frame(QX = QX, Qy = Qy))
}
my_gwpr <- function (dmat, bandwith, kernel, adaptive, data_QXY){
  ListC     <- list()
  CoefsMat  <- matrix(NA, nrow = n, ncol = k)
  wmat      <- matrix(NA, nrow = n, ncol = n * t)
  HatMat    <- matrix(NA, nrow = n*t, ncol = n * t)
  yHat      <- c()
  Resid     <- c()

  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n){
    dist.vi       <- dmat[i,]
    W.i           <- gw.weight(dist.vi, bw = bandwidth, kernel = kernel, adaptive = adaptive)
    W.i           <- diag(W.i)
    W.i           <- kronecker(W.i, diag(t))
    wmat[i,]      <- diag(W.i)
    QXw.i         <- QX * wmat[i, ]
    ListC[[i]]    <- (solve(t(QXw.i) %*% QX)) %*% t(QXw.i)
    CoefsMat[i,]  <- ListC[[i]] %*% Qy

    for(j in 1:t){
      HatMat[(i - 1) * t + j, ] <- QX[(i - 1) * t + j, ]%*% ListC[[i]]
      yHat.temp <- CoefsMat[i,] %*% t(QX)
      yHat[(i - 1) * t + j] <- yHat.temp[(i-1) * t + j]
      Resid[(i - 1) * t + j] <- Qy[(i - 1) * t + j] - yHat[(i - 1) * t + j]
    }
    setTxtProgressBar(pb, i)
  }

  close(pb)

  # return wmat, Resid,
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
