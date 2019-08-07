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
bandwidth_optimisation <- function(formula, data, dmat, sequence, kernel = "bisquare", adaptive = TRUE, verbose = FALSE) {
  CVsMat <- matrix(NA, nrow = length(sequence), ncol = 2)
  colnames(CVsMat) <- c("Bandwidth", "CV")
  ResidCV2 <- c()
  adaptive <- TRUE


  foreach(bandwidth = sequence, i = icount(), .verbose = verbose) %do% {
    for(j in 1:length(data)){
      W.j <- gw.weight(as.numeric(dmat[j,]), bw = bandwidth, kernel = kernel, adaptive = adaptive)
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
  # A discuter avec Jeremy
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

compute_QX_QY <- function(dataset, formula, index, model, effect){
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

  result <- new("QXQY", QX = QX, QY = Qy, n = n, t = t, k = k, x = x, y = y)

  return(result)
}

######################################################
# Computation of the GWPR model
# @param dmat
# @param bandwith
# @param kernel
# @param adaptive
# @param QX
# @param QY
# @param n (total index[1])
# @param t (total index[2])
#
# @return wmat, yHat, Resid, HatMat
compute_gwpr <- function(qx_qy, dmat, bandwidth, kernel, adaptive){
  n <- qx_qy@n
  t <- qx_qy@t
  k <- qx_qy@k
  QX <- qx_qy@QX
  Qy <- qx_qy@QY

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

  result <- new('GWPR', ListC = ListC, CoefsMat = CoefsMat, wmat = wmat, HatMat = HatMat, yHat = yHat, Resid = Resid)

  return(result)
}

#' Compute a numeric distance matrix or vector (from GW tools)
#' @param dataset A numeric matrix of two columns giving the coordinates of the data points
#' @param p       the power of the Minkowski distance, default is 2, i.e. the Euclidean distance
#' @param longlat TRUE, great circle distances will be calculated
#'
#' @return Returns a numeric distance matrix or vector
#'
compute_dmat <- function(dataset, p, longlat){
  dmat <- gw.dist(dataset, p = p, longlat = longlat)

  return(dmat)
}


######################################################
# Local R-Squared
# @param wmat
# @param Qy
# @param yHat
#
# @return matrix of ....
compute_localR2Mat <- function(QXQY, GWPR){
  n <- QXQY@n
  Qy <- QXQY@QY
  wmat <- GWPR@wmat
  yHat <- GWPR@yHat

  LocalR2Mat <- matrix(NA, nrow = n, ncol = 1)

  for (i in 1:n){
    TSSw.i <- t((Qy) * wmat[i,]) %*% (Qy)
    RSSw.i <- t((Qy - yHat) * wmat[i,]) %*% (Qy - yHat)
    LocalR2Mat[i,] <- (TSSw.i - RSSw.i) / TSSw.i
  }

  return(LocalR2Mat)
}


# Global R-Squared
# @param yHat
# @param Qy
#
# @return R2
compute_R2 <- function(QX_QY, GWPR){
  yHat <- GWPR@yHat
  Qy <- QX_QY@QY

  return(cor(yHat, Qy) ^ 2)
}

#' Standard errors and T values (matrices)
#' @param Resid
#' @param HatMat
#'
#' @return SesMat, TVsMat
#'
compute_std_errors_T_values <- function(QX_QY, GWPR){
  n         <- QX_QY@n
  t         <- QX_QY@t
  k         <- QX_QY@k
  Resid     <- GWPR@Resid
  HatMat    <- GWPR@HatMat
  CoefsMat  <- GWPR@CoefsMat
  ListC     <- GWPR@ListC

  RSS <- t(Resid) %*% (Resid)
  v1  <- sum(diag(HatMat))
  v2  <- sum(HatMat ^ 2) # A ete verifie que c'est equivalent a trace(StS)

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

  return(list(SEsMat = SEsMat, TVsMat = TVsMat))
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

#' Compute Spacial Polygon Dataframes
#' @param SpDF        a spacial polygon data frame
#' @param Shapefile   a customized shape file
#' @param GWPR        class result from gwpr method
#' @param QX_QY       class result from QX_QY method
compute_shapefile <- function(SpDF, result_shapefile, result_gwpr, result_QX_QY){
  CoefsMat   <- result_gwpr@CoefsMat
  SEsMat     <- result_shapefile@SEsMat
  TVsMat     <- result_shapefile@TVsMat
  LocalR2Mat <- result_shapefile@LocalR2Mat
  x          <- result_QX_QY@x

  colnames(CoefsMat) <- paste(colnames(x), "_Coef", sep = "")
  colnames(SEsMat) <- paste(colnames(x), "_SE", sep = "")
  colnames(TVsMat) <- paste(colnames(x), "_TV", sep = "")
  colnames(LocalR2Mat) <- c("Local_RSquared")


  listMat <- list(LocalR2Mat, CoefsMat, SEsMat, TVsMat)
  for(i in 1:length(listMat)){
    newVars <- colnames(listMat[[i]])

    for(j in 1:length(newVars)){
      SpDF$x <- listMat[[i]][, j]
      last <- length(names(SpDF))
      names(SpDF)[last] <- newVars[j]
    }
  }

  res <- list(listMat, SpDF, result_shapefile@R2)

  return(res)
}

compute_gpwr_plot <- function(Dataset, P=c("0.05","0.01","0.001")){
  # T values
  P <- match.arg(P)
  n <- nrow(Dataset)
  if(P=="0.05") {TValue = qt(1-(0.05/2), n-2)}
  if(P=="0.01") {TValue = qt(1-(0.05/2), n-2)}
  if(P=="0.001"){TValue = qt(1-(0.05/2), n-2)}
  TVLimits <- c(qt(1 - (0.05 / 2), n - 2), qt( 1 - (0.01 / 2), n - 2),qt( 1 - (0.001 / 2), n - 2))

  #Data preparation
  Dataset$OID <- 1:nrow(Dataset)
  Dataset$LocalR2 <- Dataset$Local_RSquared

  # Note: Fortify is deprecated (https://www.rdocumentation.org/packages/ggplot2/versions/3.2.0/topics/fortify)
  FDataset <- fortify(Dataset, region='OID')
  FDataset <- merge(FDataset, Dataset, by.x = "id", by.y = "OID", all.x = TRUE)

  Vars      <- names(Dataset)
  TNames    <- Vars[str_detect(Vars, regex("TV"))]
  CoefNames <- Vars[str_detect(Vars, regex("Coef"))]

  R2Plot <- ggplot(data = FDataset, mapping = aes(x = long, y = lat, group = group)) +
            geom_polygon(aes(fill = LocalR2), colour = "black") +
            scale_fill_gradient(high = "black", low = "white") +
            labs(x = NULL, y = NULL, title = "Geographically Weighted Regression", subtitle = "R squared")+
            theme_void() +
            coord_equal()


  #Mapping coefficients
  CoeffPlots <- list()
  for (Name in CoefNames){
    NameT <- str_replace(Name, "_Coef", "_TV")
    FDataset$Binary <- ifelse(abs(FDataset[[NameT]])>=TValue, "1", "0")
    FDatasetB <- subset(FDataset,FDataset$Binary=="0")

    Plot <- ggplot(data = FDataset, mapping= aes(x=long,y=lat,group=group)) +
      geom_polygon(aes_string(fill=Name), colour="black")+
      scale_fill_gradient(high="#99000d",low="#fff5f0")+
      geom_polygon(data = FDatasetB, mapping= aes(x=long,y=lat,group=group), fill="gray",colour="black")+
      xlab(Name)+ylab("")+
      labs(caption = paste("Gray: not significant at P=",P, sep=""))+
      theme_void()+
      coord_equal()
    CoeffPlots[[length(CoeffPlots)+1]] <- Plot
  }

  #Mapping Significance
  TPlots <- list()
  cols <- c( "Negative: P<0.001" = "#2166AC",
             "Negative: P=[0.001-0.01[" = "#67A9CF",
             "Negative: P=[0.01-0.05[" = "#D1E5F0",
             "Not significant" = "#E6E6E6",
             "Positive: P=[0.01-0.05[" = "#FDDBC7",
             "Positive: P=[0.001-0.01[" = "#EF8A62",
             "Positive: P<0.001" = "#B2182B")
  for (Name in TNames){
    FDataset[[paste(Name,"SignLevel",sep="_")]]  <-
      cut(FDataset[[paste(Name,sep="_")]],
          breaks = c(-Inf, -TVLimits[3], -TVLimits[2], -TVLimits[1], TVLimits[1], TVLimits[2], TVLimits[3], Inf),
          labels = c("Negative: P<0.001", "Negative: P=[0.001-0.01[","Negative: P=[0.01-0.05[", "Not significant","Positive: P=[0.01-0.05[","Positive: P=[0.001-0.01[", "Positive: P<0.001"))
    Plot <- ggplot(data = FDataset, mapping= aes(x=long,y=lat,group=group)) +
      geom_polygon(aes_string(fill=paste(Name,"SignLevel",sep="_")), colour="black")+
      scale_fill_manual(values = cols)+
      xlab("TITI")+ylab("")+
      theme(axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank()
      )+
      theme_void()+
      coord_equal()
    TPlots[[length(TPlots)+1]] <- Plot
  }

  #Mapping the most local significant variable
  AbsTNames  <- c()
  for(Name in TNames){
    FDataset[[paste("Abs",Name,sep="_")]] <- abs(FDataset[[Name]])
    AbsTNames[[length(AbsTNames)+1]] <- paste("Abs",Name,sep="_")
  }
  FDataset$MostSignif <- as.factor(AbsTNames[apply(FDataset[AbsTNames],1,which.max)])
  PlotMostSignif <- ggplot(data = FDataset, mapping= aes(x=long,y=lat,group=group)) +
    geom_polygon(aes(fill=MostSignif), colour="black")+
    xlab("")+ylab("")+
    labs(x = NULL, y = NULL,
         title = "Geographically Weighted Regression",
         subtitle = paste("The most significant variable at the level of P =",P, sep=" "))+
    theme_void()+
    coord_equal()

  #Etape 6 : cartographier le nombre de variable significative
  FDataset$NSignif <- rowSums((FDataset[,TNames] > TValue | FDataset[,TNames] < -TValue ))
  FDataset$NSignif <- as.character(FDataset$NSignif)
  PlotNSignif <- ggplot(data = FDataset, mapping= aes(x=long,y=lat,group=group)) +
    geom_polygon(aes(fill=NSignif), colour="black")+
    #scale_fill_manual(values = Colors)+
    scale_fill_brewer(type = "seq", palette = "OrRd", direction = 1)+
    labs(x = NULL, y = NULL,
         title = "Geographically Weighted Regression",
         subtitle = paste("Number of significant variables at the level of P =",P, sep=" "))+
    theme_void()+
    coord_equal()

  return (list(R2Plot,CoeffPlots,TPlots,PlotMostSignif,PlotNSignif))
}
