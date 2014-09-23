\name{baseline}
\alias{baseline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform baseline analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
baseline(x, K, method = "hc", true_label = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{true_label}{
%%     ~~Describe \code{true_label} here~~
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
function (x, K, method = "hc", true_label = NA) 
{
    if (method == "hc") {
        res <- get_label_hc(x, K)
    }
    if (method == "kc") {
        res <- get_label_kc(x, K)
    }
    if (length(true_label) > 1) {
        vmeasure <- vmeasureR(res$label, true_label)
        res$vmeasure <- vmeasure
    }
    genotype <- res$genotype
    mst <- get_mst(genotype)
    res$mst <- mst
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
