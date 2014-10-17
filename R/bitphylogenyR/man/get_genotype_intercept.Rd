\name{get_genotype_intercept}
\alias{get_genotype_intercept}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get common genotypes between two populations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_genotype_intercept(g1, g2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g1}{
%%     ~~Describe \code{g1} here~~
}
  \item{g2}{
%%     ~~Describe \code{g2} here~~
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
function (g1, g2) 
{
    tt = vector("list", dim(g1)[1])
    for (i in 1:dim(g1)[1]) {
        for (j in 1:dim(g2)[1]) {
            if (sum(g1[i, ] == g2[j, ]) == dim(g1)[2]) {
                tt[[i]] = c(tt[[i]], j)
            }
        }
    }
    ind1 <- c()
    for (i in 1:length(tt)) {
        if (!is.null(tt[[i]])) {
            ind1 <- c(ind1, i)
        }
    }
    ind2 <- c()
    for (i in 1:length(tt)) {
        if (!is.null(tt[[i]])) {
            ind2 <- c(ind2, tt[[i]][1])
        }
    }
    return(list(ind1 = ind1, ind2 = ind2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
