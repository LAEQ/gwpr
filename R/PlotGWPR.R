PlotGWPR<-function(Dataset, P=c("0.05","0.01","0.001")){
  require(ggplot2)
  require(ggpubr)
  require(stringr)
  require(rgeos)
  require(maptools)

  #T values
  P <- match.arg(P)
  n <- nrow(Dataset)
  if(P=="0.05") {TValue = qt(1-(0.05/2), n-2)}
  if(P=="0.01") {TValue = qt(1-(0.05/2), n-2)}
  if(P=="0.001"){TValue = qt(1-(0.05/2), n-2)}
  TVLimits <- c(qt(1-(0.05/2), n-2), qt(1-(0.01/2), n-2),qt(1-(0.001/2), n-2))

  #Data preparation
  Dataset$OID <- 1:nrow(Dataset)
  Dataset$LocalR2 <- Dataset$Local_RSquared
  FDataset <- fortify(Dataset, region='OID')
  FDataset <- merge(FDataset,Dataset, by.x ="id",by.y="OID",all.x=TRUE)

  Vars <- names(Dataset)
  TNames <-c()
  CoefNames <-c()


  for (NameVar in Vars){
    nchars <- nchar(NameVar)
    if(substr(NameVar, nchars-4, nchars) == "_Coef"){
      CoefNames[[length(CoefNames)+1]] = NameVar
    }

    if(substr(NameVar, nchars-2, nchars) == "_TV"){
      TNames[[length(TNames)+1]] = NameVar
    }
  }

  R2Plot <- ggplot(data = FDataset, mapping= aes(x=long,y=lat,group=group)) +
    geom_polygon(aes(fill=LocalR2), colour="black")+
    scale_fill_gradient(high="black",low="white")+
    labs(x = NULL, y = NULL,
         title = "Geographically Weighted Regression",
         subtitle = "R squared")+
    theme_void()+
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
    CoeffPlots[[length(CoeffPlots)+1]]<-Plot
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
    FDataset[[paste(Name,"SignLevel",sep="_")]] <-
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
    TPlots[[length(TPlots)+1]]<-Plot
  }

  #Mapping the most local significant variable
  AbsTNames <-c()
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

  #?tape 6 : cartographier le nombre de variable significative
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
