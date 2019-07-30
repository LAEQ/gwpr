# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

library(plyr)

hello <- function() {
  print("Hello, world!")
}



data_preparation <- function(dataset, formula, id){
  if(! is.formula(formula)) stop("You must provide a valid formula")

  vars <- c(all.vars(formula), id)

  if(! all(vars %in%  colnames(dataset))) stop("Your formula cannot match the dataset")


  #todo: validate Equation and id fields are in the data set
  vars <- c("CTNAME86", all.vars(formula))

  dataAVG <- dataset %>%
    select(vars) %>%
    group_by(CTNAME86) %>%
    summarise(Y_FR = mean(Y_FR), Chomag = mean(Chomag),
              FaMono = mean(FaMono), FaibSc = mean(FaibSc),
              ImgRec = mean(ImgRec), P65= mean(P65),
              Menag1 = mean(Menag1))


}

bandwidthOption <- function(a, b) {
  return (a + b)
}
