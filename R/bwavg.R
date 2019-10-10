

bwAVG <- function(formula, data, SDF, index, approach=c("CV","AICc"), kernel="bisquare", 
                  adaptive=FALSE, p=2, longlat=FALSE, dMat){
  
  #packages used in the function
  require(plm)
  require(GWmodel)
  
  #Extraction of parameters to be used
  pdata <- pdata.frame(data, index = index, drop.index = FALSE, row.names = FALSE, 
                       stringsAsFactors = default.stringsAsFactors())
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
  bw <- bw.gwr(formula=formulaAVG, data=SDF, approach=approach, kernel=kernel,
               adaptive=adaptive, p=p, longlat=longlat, dMat=dMat)
  return(bw)
}
