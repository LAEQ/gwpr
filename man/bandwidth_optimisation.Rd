\name{bandwidth_optimisation}
\alias{bandwidth_optimisation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
title bandwith optimisation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bandwidth_optimisation(formula, data, dmat, sequence = sequence_bw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{dmat}{
%%     ~~Describe \code{dmat} here~~
}
  \item{sequence}{
%%     ~~Describe \code{sequence} here~~
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
function (formula, data, dmat, sequence = sequence_bw)
{
    CVsMat <- matrix(NA, nrow = length(sequence), ncol = 2)
    colnames(CVsMat) <- c("Bandwidth", "CV")
    ResidCV2 <- c()
    adaptive <- TRUE
    foreach(bandwidth = sequence, i = icount()) \%do\% {
        for (j in 1:672) {
            W.j <- gw.weight(as.numeric(dmat[j, ]), bw = bandwidth,
                kernel = "bisquare", adaptive = adaptive)
            W.j[j] <- 0
            data$wgt <- W.j
            lm.j <- lm(formula = formula, data = data, weights = wgt)
            ResidCV2[j] <- (lm.j$residuals[j])^2
        }
        data$wgt <- NULL
        CVsMat[i, 1] <- bandwidth
        CVsMat[i, 2] <- sum(ResidCV2)
    }
    BwOpt <- subset(CVsMat, CVsMat[, 2] == min(CVsMat[, 2], na.rm = TRUE))[1,
        1]
    return(BwOpt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
