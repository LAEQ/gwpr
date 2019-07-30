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
library(dplyr)
library(tidyr)


data_preparation <- function(dataset, formula, id){
  if( !is.formula(formula)) stop("You must provide a valid formula")

  formula_vars <- all.vars(formula)
  vars <- c(formula_vars, id)

  if( ! all(vars %in% colnames(dataset))) stop("Your formula cannot match the dataset")

  dataAVG <- dataset %>%
    select(vars) %>%
    group_by(.dots = id) %>%
    summarise_at(formula_vars, mean)
}

bandwidthOption <- function(a, b) {
  return (a + b)
}
