\name{plot_sankey_mft}
\alias{plot_sankey_mft}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_sankey_mft%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_sankey_mft(fh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fh}{
%%     ~~Describe \code{fh} here~~
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
function (fh) 
{
    treefreq <- read.csv(fh)
    mft <- treefreq[which.max(treefreq[, "freq"]), 1]
    fn <- paste("nodes-", mft, ".gdl", sep = "")
    df <- gdl2df(fn)
    plot_sankey(df$nodemat, df$edgemat, df$genomat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
