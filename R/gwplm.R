#' This function implements basic geographically weighted panel regression (GWPR)
#' The function requires the plm packages
#' Je pense que ce package est utile seulement si l'on imprime les effets individuels/temporels (indisponible pour le moment)
#' @param SpDF:         large SpatialPolygonsdataFrame
#' @param data:         dataFrame for the Panel data
#' @param index:        List for the indexes : (c(" ID, Time"))
#' @param formula:      Regression model formula : Y ~ X1 + ... + Xk
#' @param bw:           bandwidth to be used (see GWmodel::gwr.basic)
#' @param kernel:       gaussian,exponential, bisquare, tricube, boxcar (see GWmodel::gw.weight)
#' @param adaptive:     TRUE or FALSE (see GWmodel::gw.weight)
#' @param dMat:         a distance matrix or vector (Optional parameter, see GWmodel::gw.weight)
#' @param effect:       the effects introduced in the model, one of "individual", "time", or "twoways" (see plm::plm)
#' @param model:        one of "pooling", "within", "between", "random", "fd", or "ht" (see plm::plm)
#'
#' @return
#'
#' @export
gwplm <- function(SpDF, data, index, formula, bw, kernel, adaptive=F, dMat,
                  effect=c("individual", "time", "twoways", "nested"),
                  model = c("within", "random", "ht", "between", "pooling", "fd")){

  # Description:
  # This function implements basic geographically weighted panel regression (GWPR)
  # The function requires the plm and GWmodel packages

  # Arguments of the function
  #SpDF:         large SpatialPolygonsdataFrame
  #data:         dataFrame for the Panel data
  #index:        List for the indexes : (c(" ID, Time"))
  #formula:      Regression model formula : Y ~ X1 + ... + Xk
  #bw:           bandwidth to be used (see GWmodel::gwr.basic)
  #kernel:       gaussian,exponential, bisquare, tricube, boxcar (see GWmodel::gw.weight)
  #adaptive:     TRUE or FALSE (see GWmodel::gw.weight)
  #dMat:         a distance matrix or vector (Optional parameter, see GWmodel::gw.weight)
  #effect:       the effects introduced in the model, one of "individual", "time", or "twoways" (see plm::plm)
  #model:        one of "pooling", "within", "between", "random", "fd", or "ht" (see plm::plm)

  #Panel data preparation
  Pdata <- plm::pdata.frame(data, index = index, drop.index = FALSE, row.names = FALSE,
                            stringsAsFactors = default.stringsAsFactors())
  PanelModel <- plm::plm(formula = formula, effect=effect, model=model, data=Pdata, index=index)
  y <- plm::pmodel.response(PanelModel)
  x <- model.matrix(PanelModel)
  n <- length(unique(Pdata[,index[1]]))
  t <- length(unique(Pdata[,index[2]]))
  K <- ncol(x)

  #Global panel model
  cat("***********************************************************************\n")
  cat("*                  Results of Global Panel Regression                 *\n")
  cat("***********************************************************************\n")

  print(summary(PanelModel))
  cat("R-square value: ", round(plm::r.squared(PanelModel, type="rss"), 5))
  cat("\nAdjusted R-square value: ", round(plm::r.squared(PanelModel, dfcor=TRUE), 5))

  #Computation of the GWPR model
  ###############################

  #Computation of local weighted panel regressions
  wMat <- matrix(NA, nrow = n, ncol = n*t)
  listC <- list()
  coefsMat <- matrix(NA, nrow = n, ncol = K)
  hatMat <- matrix(NA, nrow = n*t, ncol = n*t)
  yHat <- c()
  resid <- c()

  cat("\nComputing ", n, " local weighted panel regressions\n", sep="")
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n){
    W.i <- GWmodel::gw.weight(as.numeric(dMat[i,]),bw=bw,kernel=kernel,adaptive=adaptive)
    W.i <- diag(W.i)
    W.i <- kronecker(W.i,diag(t))
    wMat[i,] <- diag(W.i)
    C.i <- (solve(t(x)%*%W.i%*%x))%*%t(x)%*%W.i
    listC[[i]] <- C.i
    Pdata[['wgt']] <- wMat[i,]
    plm.i <- plm::plm(formula = formula, effect=effect, model=model, data=Pdata, index=index, weights = wgt)
    coefsMat[i,] <- plm.i[['coefficients']]

    for(j in 1:t){
      hatMat[(i-1)*t+j,] <- x[(i-1)*t+j,]%*%C.i
      yHat[(i-1)*t+j] <- predict(plm.i)[(i-1)*t+j]
      resid[(i-1)*t+j] <- plm.i[['residuals']][(i-1)*t+j]
    }
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  #Local R-Squared
  localR2 <- numeric()

  for (i in 1:n){
    TSSw.i <- t((y-mean(y))*wMat[i,])%*%(y-mean(y))
    RSSw.i <- t((y-yHat)*wMat[i,])%*%(y-yHat)
    localR2[i] <- (TSSw.i - RSSw.i)/TSSw.i
  }

  #Global R-Squared
  TSS <- t(y-mean(y))%*%(y-mean(y))
  RSS <- t(resid)%*%(resid)
  v1 <- sum(diag(hatMat))
  v2 <- sum(hatMat^2)
  edf <- n*t -2*v1 + v2 # effective degrees of freedom
  R2 <- (TSS - RSS)/TSS
  adj.R2 <- 1 - (1 - R2)*(n*t - 1)/(edf - 1)

  #Standard errors and T values (matrices)
  sigma2 <- as.numeric(RSS/edf)
  SEsMat <- matrix(nrow = n, ncol = K)
  TVsMat <- matrix(nrow = n, ncol = K)

  for (i in 1:n){
    C.i <- listC[[i]]
    SEs.i <- c()
    SEs.i <- sqrt(diag(C.i%*%(t(C.i)))*sigma2)
    for (k in 1:K){
      SEsMat[i,k] <- SEs.i[k]
      TVsMat[i,k] <- coefsMat[i,k] / SEsMat[i,k]
    }
  }

  #Rename objects
  colnames(coefsMat) <- paste(colnames(x), "_coef", sep = "")
  colnames(SEsMat) <- paste(colnames(x), "_SE", sep = "")
  colnames(TVsMat) <- paste(colnames(x), "_TV", sep = "")

  env <- new.env(parent=globalenv())
  assign("Coefficients matrix", coefsMat, envir=env)
  assign("Std.Errors matrix", SEsMat, envir=env)
  assign("T-values matrix", TVsMat, envir=env)
  assign("Local R-squared", localR2, envir=env)
  assign("Global R-squared", R2, envir=env)
  assign("Adjusted global R-squared", adj.R2, envir=env)
  assign("Spatial panel data frame", SpDF, envir=env)

  #Creating the shapefile with the GWPR results(coefficients, standard errors and T values)
  listMat <- list(coefsMat, SEsMat, TVsMat)
  for(i in 1:length(listMat)){
    newVars <- colnames(listMat[[i]])
    for(j in 1:length(newVars)){
      SpDF[['x']] <- listMat[[i]][,j]
      last <- length(names(SpDF))
      names(SpDF)[last] <- newVars[j]
    }
  }
  SpDF[['localR2']] <- localR2

  cat("\n***********************************************************************\n")
  cat("*         Results of Geographically Weighted Panel Regression         *\n")
  cat("***********************************************************************\n")
  cat("*********************Model Calibration information*********************\n")
  cat("Kernel Function:", kernel, "\n")
  if(adaptive==TRUE) cat("Adaptive bandwidth:", bw, "(number of nearest neighbours)\n") else cat("Fixed bandwidth:", round(bw,3), "(distance)\n")
  cat("Panel model used:", model, "\n")
  cat("Type of effects:", effect, "\n")

  cat("\n***************Summary of GWPR coefficient estimates:*****************\n")
  print(summary(coefsMat))
  cat("\n*********************Summary of GWPR T values:************************\n")
  print(summary(TVsMat))
  cat("\n**************Summary of GWPR local R squared values:*****************\n")
  print(summary(localR2))

  cat("\n************************Diagnostic information************************\n")
  cat("Number of observations (n):", n,"\n")
  cat("Number of time periods (t):", t,"\n")
  cat("Number of data points (n*t):", n*t,"\n")
  cat("Effective degrees of freedom (n*t-2trace(S) + trace(S'S)):", round(edf, 3),"\n")
  cat("Residual sum of squares:", round(RSS, 5),"\n")
  cat("R-squared value:", round(R2, 5),"\n")
  cat("Adjusted R-squared value:", round(adj.R2, 5),"\n")

  res <- list(listMat, SpDF, R2, adj.R2)
  return(res)
}
