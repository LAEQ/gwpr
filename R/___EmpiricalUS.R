# Libraries
library(Rcpp)
library(foreign)
library(rgdal)
library(GWmodel)
library(plm)
library(ggplot2)
library(ggpubr)
library(stringr)
library(spdep)
library(stargazer)
library(plyr)

# Working directory
rm(list=ls())
setwd("Z:/___articles/shp")

# Importing functions for the GWPR
source("Z:/___articles/article1_GWPR/___gwplm.R")
source("Z:/___articles/article1_GWPR/___bwOpt.R")
source("Z:/___articles/article1_GWPR/___bwAICc.R")
source("Z:/___articles/article1_GWPR/___PlotGWPR.R")

##############################################
# Data preparation
#############################################

load("Data.Rdata")
# L'ordre des Etats dans le dataset ("Produc") n'est pas le même que celui du shapefile.
# On commence donc par reordonner et fusionner les fichiers pour que les weights soient
# attribues aux bons Etats danms les prochaines etapes.
USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
data <- merge(USStates@data, Produc, by="state", all=T)

dMat <- gw.dist(coordinates(USStates), p=2, longlat=F)
Equation <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp


#############################################
# Bandwidth selection (CV-score)
#############################################

# Version avec donnees moyennes (sur le temps)
## adaptive bw
bwAVG.A <- bw.avg(formula=Equation, data=data, SDF=USStates, index=c("id","year"), approach="AICc", 
                  kernel="bisquare", adaptive=T, p=2, longlat=F, dMat=dMat)
## fixed bw
bwAVG.F <- bw.avg(formula=Equation, data=data, SDF=USStates, index=c("id","year"), approach="AICc", 
                  kernel="bisquare", adaptive=F, p=2, longlat=F, dMat=dMat)

# Version avec pseudo-CV
## adaptive bw
bwCV.A <-  bw.CV.A(formula=Equation, data=data, index=c("id","year"), effect='individual', model="within", 
                   kernel="bisquare", dMat=dMat, bws=c(30:40))
## fixed bw
bwCV.F <-  bw.CV.F(formula=Equation, data=data, index=c("id","year"), effect='individual', model="within", 
                   kernel="bisquare", dMat=dMat, interval=c(1500000, 2500000))

# Version avec pseudo-AICc
# !!! CONVERGENT VERS LES BORNES DES INTERVALLES, DONC LE CALCUL DES AICc DOIT ETRE ERRONE !!!
## adaptive bw
bwAICc.A <- bw.AICc.A(formula=Equation, data=data, index=c("id","year"), effect='individual', model="within", 
                      kernel="bisquare", dMat=dMat, Lw.bw=40, Up.bw=48, step=1)
## fixed bw
bwAICc.F <- bw.AICc.F(formula=Equation, data=data, index=c("id","year"), effect='individual', model="within", 
                      kernel="bisquare", dMat=dMat, interval=c(750000, 2500000))

#############################################
# GWPR model
#############################################

# L'application empirique avec laquelle nous comparons (panel classique) est celle de Baltagi (2005, p.25)
# qui reprend l'exemple de Munnell (1990)

# individual FE
## adaptive bw
USgwplm.A <- gwplm(SpDF=USStates, data=data, index=c("id", "year"), formula=Equation, bw=bwCV.A, kernel="bisquare", 
                   adaptive=T, effect="individual", model="within", dMat=dMat)
## fixed bw
USgwplm.B <- gwplm(SpDF=USStates, data=data, index=c("id", "year"), formula=Equation, bw=bwCV.F, kernel="bisquare", 
                   adaptive=F, effect="individual", model="within", dMat=dMat)


#############################################
# Cartography of the models 
#############################################

###Plots <- PlotGWPR(USgwplm[[2]], P="0.05")
#Mapping local coefficients
###ggarrange(plotlist = list(Plots[[2]][[1]], Plots[[2]][[2]], Plots[[2]][[3]], Plots[[2]][[4]]), ncol = 2, nrow = 2)
#Mapping local TVs 
###ggarrange(plotlist = list(Plots[[3]][[1]], Plots[[3]][[2]], Plots[[3]][[3]], Plots[[3]][[4]]), ncol = 2, nrow = 2)
#Mapping the most significant variable
###ggarrange(plotlist = list(Plots[[5]]), ncol = 1, nrow = 1)
#Mapping the number of significant variables
###ggarrange(plotlist = list(Plots[[4]]), ncol = 1, nrow = 1)
#Mapping R squared
###Plots[[1]]
