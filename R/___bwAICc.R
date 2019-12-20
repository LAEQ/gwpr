###############################################################################################
# Author: Maxime Gaboriault-Boudreau (functions: bwopt.avg, bw.CV.A, bw.CV.F)
# December 2019
# Functions for calculating and mapping a geographically weighted panel regression (GWPR)
# This reseach project was funded by the LAEQ (http://laeq.ucs.inrs.ca/) 
# and supervised by Philippe Apparicio
###############################################################################################

#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#
# Les fonctions qui suivent ne fonctionnent pas. Les resultats convergent,
# mais systematiquement vers les bornes de l'intervalle de bandwidths propose.
#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#

bw.AICc.A <- function(formula, data, effect=c("individual", "time", "twoways", "nested"), index, 
                      model=c("within", "random", "ht", "between", "pooling", "fd"), 
                      kernel="bisquare", dMat, Lw.bw, Up.bw, step){
  
  Pdata <- pdata.frame(data, index, drop.index=F, row.names=F, stringsAsFactors=default.stringsAsFactors())
  PanelModel <- plm(formula=formula, effect=effect, model=model, data=Pdata, index=index)
  y <- pmodel.response(PanelModel)
  x <- model.matrix(PanelModel)
  n <- length(unique(Pdata[,index[1]]))
  t <- length(unique(Pdata[,index[2]]))
  K <- ncol(x)
  
  seq <- seq(from=Lw.bw, to=Up.bw, by=step)
  AICcMat <- matrix(NA, nrow=length(seq), ncol=2)
  colnames(AICcMat) <- c("Bandwidth", "AICc-score")
  hatMat <- matrix(NA, nrow = n*t, ncol = n*t)
  resid <- c()
  cat("\nOptimization Process:", length(seq), "iterations \n")
  pb <- txtProgressBar(min = 0, max = length(seq), style = 3)
  
  for(z in 1:length(seq)){
    bw <- seq[z]
    for (i in 1:n){
      W.i <- gw.weight(as.numeric(dMat[i,]),bw=bw,kernel=kernel,adaptive=T)
      W.i <- diag(W.i)
      W.i <- kronecker(W.i,diag(t))
      C.i <- (solve(t(x)%*%W.i%*%x))%*%t(x)%*%W.i
      Pdata$wgt <- as.numeric(diag(W.i))
      plm.i <- plm(formula = formula, effect=effect, model=model, data=Pdata, index=index, weights=wgt)
      
      for(j in 1:t){
        hatMat[(i-1)*t+j,] <- x[(i-1)*t+j,]%*%C.i
        resid[(i-1)*t+j] <- pmodel.response(plm.i)[(i-1)*t+j] - predict(plm.i)[(i-1)*t+j]
      }
    }
    AICcMat[z,1] <- bw
    sigmaHat <- sd(resid)
    # En utilisant n*t plutot que n dans l'equation qui suit, ca converge vers la bw minimale 
    # (voir fonction bw.AICc.F) alors qu'avec seulement n ca converge vers la valeur maximale
    AICcMat[z,2] <- 2*n*log(sigmaHat)+n*log(2*pi)+n*((n+sum(diag(hatMat)))/(n-2-sum(diag(hatMat)))) 
    Pdata$wgt <- NULL                                                                               
    resid <- NULL
    setTxtProgressBar(pb, z)
  }
  close(pb)
  
  bw <- subset(AICcMat, AICcMat[,2]==min(AICcMat[,2], na.rm=TRUE))
  print(as.data.frame(AICcMat))
  cat("Optimal value for the bandwidth (number of neighbors): ", bw[,'Bandwidth'], sep="")
  cat("\nwith an AICc score value: ", bw[,'AICc-score'], sep="")
  return(bw[,'Bandwidth'])
}
# FIN DE LA FONCTION -----------------------------------------------------------------------------

bw.AICc.F <- function(formula, data, index, effect=c("individual", "time", "twoways", "nested"), 
                      model=c("within", "random", "ht", "between", "pooling", "fd"),
                      kernel="bisquare", dMat, interval){
  
  #Data preparation
  Pdata <- pdata.frame(data, index, drop.index=F, row.names=F, stringsAsFactors=default.stringsAsFactors())
  PanelModel <- plm(formula=formula, effect=effect, model=model, data=Pdata, index=index)
  y <- pmodel.response(PanelModel)
  x <- model.matrix(PanelModel)
  n <- length(unique(Pdata[,index[1]]))
  t <- length(unique(Pdata[,index[2]]))
  hatMat <- matrix(NA, nrow = n*t, ncol = n*t)
  resid <- c()

  #Optimization Process
  AICc <- function(bw, n, kernel, adaptive){
    
    for (i in 1:n){
      W.i <- gw.weight(as.numeric(dMat[i,]),bw=bw,kernel=kernel,adaptive=T)
      W.i <- diag(W.i)
      W.i <- kronecker(W.i,diag(t))
      C.i <- (solve(t(x)%*%W.i%*%x))%*%t(x)%*%W.i
      Pdata$wgt <- as.numeric(diag(W.i))
      plm.i <- plm(formula = formula, effect=effect, model=model, data=Pdata, index=index, weights=wgt)
      
      for(j in 1:t){
        hatMat[(i-1)*t+j,] <- x[(i-1)*t+j,]%*%C.i
        resid[(i-1)*t+j] <- pmodel.response(plm.i)[(i-1)*t+j] - predict(plm.i)[(i-1)*t+j]
      }
    }
    
    sigmaHat <- sd(resid)
    AICc <- 2*n*t*log(sigmaHat)+n*t*log(2*pi)+n*t*((n*t+sum(diag(hatMat)))/(n*t-2-sum(diag(hatMat))))
    return(AICc)
  }
  
  bw.AICc <- optimize(AICc, interval=interval, n=n, kernel=kernel, adaptive=F, maximum=F)
  bw <- as.numeric(bw.AICc)
  names(bw) <- c('Bandwidth', 'AICc')
  cat("Optimal value for the bandwidth (distance): ", bw[['Bandwidth']], sep="")
  cat("\nwith an AICc of: ", bw[['AICc']], sep="")
  
  return(bw[['Bandwidth']])
}