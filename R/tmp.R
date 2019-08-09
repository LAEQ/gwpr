# # Libraries
# library(Rcpp)
# library(foreign)
# library(rgdal)
# library(GWmodel)
# library(plm)
# library(lmtest)
# library(ggplot2)
# library(ggpubr)
# library(stringr)
# library(spdep)
# library(stargazer)
# library(plyr)
#
# # source("R/bandwith_option.R")
#
#
# # US
# load("data/us_data/Data.Rdata")
# index <- c("state", "year")
# model <- "within"
# kernel <- "bisquare"
# adaptive <- TRUE
# effect <- "individual"
#
#
# formula <- gsp ~ pcap + pc + emp + unemp
# dmat    <- gw.dist(coordinates(USStates), p=2, longlat=F)
#
# Pdata <- pdata.frame(Produc, index = index, drop.index = FALSE, row.names = FALSE, stringsAsFactors = default.stringsAsFactors())
# y     <- model.response(model.frame(formula ,data = Produc))
# x     <- model.matrix(formula, data = Produc)
#
# if (model == "within") x <- x[,-1]
# n <- length(unique(Pdata[,index[1]]))
# t <- length(unique(Pdata[,index[2]]))
# K <- ncol(x)
#
# balanced<-n*t==dim(Produc)[1]
# if(!balanced) stop("Estimation method unavailable for unbalanced panels")
# indic<-seq(1,n)
# inde<-as.numeric(rep(indic,each=t)) ####takes the first t observations
# indic1<-seq(1,t)
# inde1<-rep(indic1,n) ####takes observations 1,  t+1, 2t+1...
#
#
# #Within transformation
# ##Calculating the mean of Y and X for each time period
# if (effect %in% c("time", "twoways")){
#   ytms<-tapply(y,inde,mean)
#   tpms<-function(q) tapply(q,inde,mean)
#   xtms<-apply(x,2,tpms)
#   ytm<-rep(ytms,each=t)
#   xtm<-matrix(NA,(n*t),K)
#   for (k in 1:K) xtm[,k]<-rep(xtms[,k],each=t)
# }
# ##Calculating the mean of Y and X for each individual
# if (effect %in% c("individual", "twoways")){
#   ysms<-tapply(y,inde1,mean)
#   spms<-function(q) tapply(q,inde1,mean)
#   xsms<-apply(x,2,spms)
#   ysm<-rep(ysms,n)
#   xsm<-matrix(NA,(n*t),K)
#   for (k in 1:K) xsm[,k]<-rep(xsms[,k],n)
# }
# ##Calculating the within-transformed Y and X (Qy and QX) according to the selected model
# if (effect=="time"){
#   Qy<-y-ytm
#   QX<-x-xtm
# }
# if(effect=="individual"){
#   Qy<-y-ysm
#   QX<-x-xsm
# }
# if (effect=="twoways"){
#   Qy<-y - ysm - ytm + rep(mean(y),n*t)
#   xmm<-matrix(NA,n*t,K)
#   for (k in 1:K) xmm[,k]<-rep(mean(x[,k]),n*t)
#   QX<-x - xsm - xtm + xmm
# }
#
# ListC <- list()
# CoefsMat <- matrix(NA, nrow = n, ncol = k)
# wmat <- matrix(NA, nrow = n, ncol = n*t)
# HatMat <- matrix(NA, nrow = n*t, ncol = n*t)
# yHat <- c()
# Resid <- c()
#
# load(file="data/us_data/bandwidth_expected.rda")
#
# bw <- BwOpt
#
# cat("\nComputing ", n, " local weighted panel regressions\n", sep="")
# pb <- txtProgressBar(min = 0, max = n, style = 3)
# for (i in 1:n){
#   dist.vi <- dmat[i,]
#   W.i <- gw.weight(dist.vi,bw=bw,kernel=kernel,adaptive=adaptive)
#   W.i <- diag(W.i)
#   W.i <- kronecker(W.i,diag(t))
#   wmat[i,] <- diag(W.i)
#   QXw.i <- QX*wmat[i,]
#   ListC[[i]] <- (solve(t(QXw.i)%*%QX))%*%t(QXw.i)
#   CoefsMat[i,] <- ListC[[i]]%*%Qy
#   for(j in 1:t){
#     HatMat[(i-1)*t+j,] <- QX[(i-1)*t+j,]%*%ListC[[i]]
#     tmp <- CoefsMat[i,]%*%t(QX)
#     yHat[(i-1)*t+j] <- tmp[(i-1)*t+j]
#     Resid[(i-1)*t+j] <- Qy[(i-1)*t+j] - yHat[(i-1)*t+j]
#   }
#   setTxtProgressBar(pb, i)
# }
# close(pb)
#
#
#
# LocalR2Mat <- matrix(NA, nrow = n, ncol = 1)
# for (i in 1:n){
#   TSSw.i <- t((Qy)*wmat[i,])%*%(Qy)
#   RSSw.i <- t((Qy-yHat)*wmat[i,])%*%(Qy-yHat)
#   LocalR2Mat[i,] <- (TSSw.i - RSSw.i)/TSSw.i
# }
#
# R2 <- cor(yHat, Qy)^2
#
# #Standard errors and T values (matrices)
# RSS <- t(Resid)%*%(Resid)
# v1 <- sum(diag(HatMat))
# v2 <- sum(HatMat^2) # A ?t? v?rifi? que c'est ?quivalent ? trace(StS)
# sigma2 <- as.numeric(RSS/(n*t - 2*v1 + v2)) # Selon l'?quation adapt?e au format panel de la GWR
# SEsMat <- matrix(nrow = n, ncol = k)
# TVsMat <- matrix(nrow = n, ncol = k)
#
# for (i in 1:n){
#   C.i <- ListC[[i]]
#   SEs.i <- c()
#   SEs.i <- sqrt(diag(C.i%*%(t(C.i)))*sigma2)
#   for (j in 1:k){
#     SEsMat[i,j] <- SEs.i[j]
#     TVsMat[i,j] <- CoefsMat[i,j] / SEsMat[i,j]
#   }
# }
#
# save(ListC, CoefsMat, wmat, HatMat, yHat, Resid, file="data/us_data/gwpr_expected_2.rda")
#
# #Rename objects
# colnames(CoefsMat) <- paste(colnames(x), "_Coef", sep = "")
# colnames(SEsMat) <- paste(colnames(x), "_SE", sep = "")
# colnames(TVsMat) <- paste(colnames(x), "_TV", sep = "")
# colnames(LocalR2Mat) <- c("Local_RSquared")
#
# #Creating the shapefile with the GWPR results(coefficients, standard errors and T values)
# listMat <- list(LocalR2Mat,CoefsMat,SEsMat,TVsMat)
# for(i in 1:length(listMat)){
#   newVars <- colnames(listMat[[i]])
#   for(j in 1:length(newVars)){
#     USStates$x <- listMat[[i]][,j]
#     last <- length(names(USStates))
#     names(USStates)[last] <- newVars[j]
#   }
# }
#
#
#
# # save(BwOpt, Opt, file="data/us_data/bandwidth_expected.rda")
# # save(dataAVG_expected, file = "data/us_data/dataAVG_expected.rda")
# # my_dataAVG <- data_preparation(dataset = Produc, formula = formula, id = "state")
#
#
# # bwselection <- bwopt_cv(data=Produc, index=c("state", "year"), formula=Equation, kernel="bisquare",
# #                         adaptive=T, dmat=dmat, lowerBW=1, upperBW=48, step=1)
# #
# #
# # my_avg = data_preparation(Produc, Equation, "state")
# # my_bw <- bandwidth_optimisation(formula, my_dataAVG, dmat, seq(1, 48, 1))/
# #
# #
# # save(my_avg, file="my_avg.rda")
# #
# #
# # bwselection == my_bw
#
#
#
#
# # Mtl
# # load("data/mtl_data/PanelMtl.rda")
# # load("data/mtl_data/dmat.rda")
# # Equation <- Y_FR ~ Chomag + FaMono + FaibSc + ImgRec + P65 + Menag1
# #
# # bwselection <- bwopt_cv(data=PanelMtl, index=c("CTNAME86", "AN"), formula=Equation, kernel="bisquare",
# #                         adaptive=T, dmat=dmat,  lowerBW=50, upperBW=200, step=10)
# #
# # my_avg = data_preparation(PanelMtl, Equation, "CTNAME86")
# # my_bw <- bandwidth_optimisation(Equation, my_avg, dmat, seq(50, 200, 10))
#
