library(plm)
# library(GWmodel)

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
  pdata <- pdata.frame(data, index = index, drop.index = FALSE, row.names = FALSE, stringsAsFactors = default.stringsAsFactors())
  panelModel <- plm(formula = formula, model="pooling", data=pdata, index=index)
  n <- length(unique(pdata[,index[1]]))
  t <- length(unique(pdata[,index[2]]))
  K <- length(all.vars(formula)[-1])
  y <- pmodel.response(panelModel) #allows for keeping variable transformations (log, ^x, ...)
  x <- model.matrix(panelModel)[,-1]

  #Verification that the panel dataset is balanced
  if(dim(pdata)[1]!=n*t) stop("Averaging process available only for balanced panels")

  #Averaging the dataset over time
  idx <- as.numeric(rep(seq(1,n),each=t))
  yAVG <- tapply(y,idx,mean)
  tpms <- function(q) tapply(q,idx,mean)
  xAVG <- apply(x,2,tpms)
  dataAVG <- as.data.frame(cbind(yAVG,xAVG))
  colnames(dataAVG) <- c(formula[[2]], colnames(x))
  SDF@data <- dataAVG

  #Optimization process with transformed data (left without time dimension).
  formulaAVG <- formula(dataAVG) #Creation of a new formula to avoid transforming variables a second time (ex. log(log(X)), ...)
  bw <- GWmodel::bw.gwr(formula=formulaAVG, data=SDF, approach=approach, kernel=kernel,
               adaptive=adaptive, p=p, longlat=longlat, dMat=dMat)
  return(bw)

  return(1)
}
