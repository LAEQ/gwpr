# Geographically Weighted Panel Regression (GWPR)

[![Build status](https://ci.appveyor.com/api/projects/status/axu6kxqja7p96r54/branch/master?svg=true)](https://ci.appveyor.com/project/davidmaignan/gwpr/branch/master)


Library aims to offer explicit equations for geographically weighted panel regression (GWPR) and R functions to optimize a bandwidth selection in the panel data case and to produce GWPR results and illustrations.

## Requirements

R (>= 4.0)

## Installation

```bash
devtools.install_github("LAEQ/gwpr")

??gwpr
```

## Examples

Some examples of the main features


```bash
library(gwpr)
data(USStates)

# L'ordre des Etats dans le dataset ("Produc") n'est pas le m?me que celui du shapefile.
# On commence donc par reordonner et fusionner les fichiers pour que les weights soient
# attribues aux bons Etats danms les prochaines etapes.
USStates@data$id <- c(1:length(unique(USStates@data[,"state"])))
data <- merge(USStates@data, Produc, by="state", all=T)

dMat <- gw.dist(coordinates(USStates), p=2, longlat=F)
Equation <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp


# Bandwidth selection (CV-score)

# Version avec donnees moyennes (sur le temps)
# adaptive bw
bwAVG.A <- bw.avg(formula=Equation, data=data, SDF=USStates, index=c("id","year"), approach="AICc", 
                  kernel="bisquare", adaptive=T, p=2, longlat=F, dMat=dMat)
# fixed bw
bwAVG.F <- bw.avg(formula=Equation, data=data, SDF=USStates, index=c("id","year"), approach="AICc", 
                  kernel="bisquare", adaptive=F, p=2, longlat=F, dMat=dMat)

# Version avec pseudo-CV
# adaptive bw
bwCV.A <-  bw.CV.A(formula=Equation, data=data, index=c("id","year"), effect='individual', model="within", 
                   kernel="bisquare", dMat=dMat, bws=c(30:40))
# fixed bw
bwCV.F <-  bw.CV.F(formula=Equation, data=data, index=c("id","year"), effect='individual', model="within", 
                   kernel="bisquare", dMat=dMat, interval=c(1500000, 2500000))

# GWPR model

# L'application empirique avec laquelle nous comparons (panel classique) est celle de
# Baltagi, B. H. (2005) Econometric Analysis of Panel Data (third ed.) qui reprend l'exemple de Munnell (1990)

# individual FE
# adaptive bw
USgwplm.A <- gwplm(SpDF=USStates, data=data, index=c("id", "year"), formula=Equation, bw=bwCV.A, kernel="bisquare", 
                   adaptive=T, effect="individual", model="within", dMat=dMat)
# fixed bw
USgwplm.F <- gwplm(SpDF=USStates, data=data, index=c("id", "year"), formula=Equation, bw=bwCV.F, kernel="bisquare", 
                   adaptive=F, effect="individual", model="within", dMat=dMat)

```

## Authors
- Maxime Gaboriault - Creator
- Philippe Apparicio - Author
- David Maignan - Compiler



