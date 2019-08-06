\name{data_preparation}
\alias{data_preparation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
title data preparation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data_preparation(dataset, formula, id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dataset, formula, id)
{
    if (!is.formula(formula))
        stop("You must provide a valid formula")
    formula_vars <- all.vars(formula)
    vars <- c(formula_vars, id)
    if (!all(vars \%in\% colnames(dataset)))
        stop("Your formula cannot match the dataset")
    dataAVG <- dataset \%>\% select(vars) \%>\% group_by(.dots = id) \%>\%
        summarise_at(formula_vars, mean)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
