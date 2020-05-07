#' A function for bandwith selection to calibrate a GWPR model, based on the mean over time of the data.
#' The formula for calculating the CV value for Panel is based on the study of:
#' YU, Danlin. 2010. Exploring spatiotemporally varying regressed relationships: the geographically weighted panel regression analysis.
#' In : Proceedings of the Joint International Conference on Theory, data Handling and modelling in GeoSpatial Information Science. p. 134-
#'
#' @param formula Regression model formula : Y ~ X1 + ... + Xk
#' @param data         dataFrame for the Panel data
#' @param SDF          large SpatialPolygonsdataFrame on which is based the data
#' @param index        List for the indexes : (c(" ID, Time"))
#' @param approach     score used to optimize the bandwidth (see GWmodel::bw.gwr)
#' @param kernel       gaussian,exponential, bisquare, tricube, boxcar (see GWmodel::gw.weight)
#' @param adaptive     TRUE or FALSE (see GWmodel::gw.weight)
#' @param p	           the power of the Minkowski distance, default is 2, i.e. the Euclidean distance (see GWmodel::bw.gwr)
#' @param longlat      if TRUE, great circle distances will be calculated (see GWmodel::bw.gwr)
#' @param dMat         a distance matrix or vector (Optional parameter, see GWmodel::gw.weight)
#'
#' @return double
#'
#' @export
bw.avg <- function(formula, data, SDF, index, approach=c("CV","AICc"), kernel="bisquare", adaptive=FALSE, p=2, longlat=FALSE, dMat){
  #Extraction of parameters to be used
  pdata       <- plm::pdata.frame(data, index = index, drop.index = FALSE, row.names = FALSE, stringsAsFactors = default.stringsAsFactors())
  panelModel  <- plm::plm(formula = formula, model="pooling", data=pdata, index=index)
  n           <- length(unique(pdata[,index[1]]))
  t           <- length(unique(pdata[,index[2]]))
  K           <- length(all.vars(formula)[-1])
  y           <- plm::pmodel.response(panelModel) #allows for keeping variable transformations (log, ^x, ...)
  x           <- stats::model.matrix(panelModel)[,-1]

  #Verification that the panel dataset is balanced
  if(dim(pdata)[1]!=n*t) stop("Averaging process available only for balanced panels")

  #Averaging the dataset over time
  idx       <- as.numeric(rep(seq(1,n),each=t))
  yAVG      <- tapply(y,idx,mean)
  tpms      <- function(q) tapply(q,idx,mean)
  xAVG      <- apply(x,2,tpms)
  dataAVG   <- as.data.frame(cbind(yAVG,xAVG))
  colnames(dataAVG) <- c(formula[[2]], colnames(x))
  SDF@data  <- dataAVG

  #Optimization process with transformed data (left without time dimension).
  formulaAVG  <- formula(dataAVG) #Creation of a new formula to avoid transforming variables a second time (ex. log(log(X)), ...)
  bw          <- GWmodel::bw.gwr(formula=formulaAVG, data=SDF, approach=approach, kernel=kernel, adaptive=adaptive, p=p, longlat=longlat, dMat=dMat)
  return(bw)
}

#' A function for bandwith selection to calibrate a GWPR model, based on the mean over time of the data.

#' Arguments of the function
#' @param formula Regression model formula : Y ~ X1 + ... + Xk
#' @param data dataFrame for the Panel data
#' @param index List for the indexes : (c(" ID, Time"))
#' @param effect the effects introduced in the model, one of "individual", "time", or "twoways" (see plm::plm)
#' @param model one of "pooling", "within", "between", "random", "fd", or "ht" (see plm::plm)
#' @param kernel gaussian,exponential, bisquare, tricube, boxcar (see GWmodel::gw.weight)
#' @param dMat a distance matrix or vector (Optional parameter, see GWmodel::gw.weight)
#' @param bws bandwidths to be used for calculations of CV-score
#'
#' @return list(bandwidth = double)
#'
#' @export
bw.CV.A <- function(formula, data, index, effect=c("individual", "time", "twoways", "nested"),
                    model=c("within", "random", "ht", "between", "pooling", "fd"),
                    kernel="bisquare", dMat, bws){


  #Data preparation
  Pdata <- plm::pdata.frame(data, index = index, drop.index = FALSE, row.names = FALSE, stringsAsFactors = default.stringsAsFactors())
  n <- length(unique(Pdata[,index[1]]))
  t <- length(unique(Pdata[,index[2]]))

  CVsMat <- matrix(NA, nrow = length(bws), ncol = 2)
  colnames(CVsMat) <- c("Bandwidth", "CV-score")
  resid <- c()

  # Faire ici le critere de mesure du data pour, si n*t > C, impression d'un message tel quel
  # "Long data, calculations may take a long time, optimisation on time-averaged data suggested
  # by choosing 'short = TRUE'." et auquel cas moyenner les donnees et les passer dans bw.avg

  #Optimization Process
  cat("\nOptimization Process:", length(bws), "iterations \n")
  pb <- txtProgressBar(min = 0, max = length(bws), style = 3)

  for(i in 1:length(bws)){
    bw <- bws[i]

    for(j in 1:n){
      W.j <- GWmodel::gw.weight(as.numeric(dMat[j,]), bw=bw, kernel=kernel, adaptive=T)
      W.j[j] <- 0
      W.j <- diag(W.j)
      W.j <- kronecker(W.j,diag(t))
      W.j <- diag(W.j)
      Pdata$wgt <- W.j
      plm.j <- plm::plm(formula=formula, model=model, data=Pdata, index=index, weights=wgt)

      for(l in 1:t){
        resid[(j-1)*t+l] <- plm::pmodel.response(plm.j)[(j-1)*t+l] - predict(plm.j)[(j-1)*t+l]
      }
    }

    CVsMat[i,1] <- bw
    CVsMat[i,2] <- t(resid)%*%(resid)
    Pdata$wgt <- NULL
    resid <- NULL
    setTxtProgressBar(pb, i)
  }
  close(pb)

  bw <- subset(CVsMat, CVsMat[,2]==min(CVsMat[,2], na.rm=TRUE))
  print(as.data.frame(CVsMat))
  cat("Optimal value for the bandwidth (number of neighbors): ", bw[,'Bandwidth'], sep="")
  cat("\nwith a CV score value: ", bw[,'CV-score'], sep="")

  return(bw[,'Bandwidth'])
}
