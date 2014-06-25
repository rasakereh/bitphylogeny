\name{bitphyloR}
\alias{bitphyloR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function to run BitPhylogeny analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bitphyloR(fin, fout, n = 100, b = 10, t = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fin}{
    data file
}
  \item{fout}{
    result directory 
}
\item{n}{
  number of mcmc samples
}
\item{b}{
  number of burnin samples
}
\item{t}{
  thinning window
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The function doesn't return anything. Instead, it makes a result
  directory and saves MCMC traces. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ke Yuan
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
function (fin, fout, n = 100, b = 10, t = 1) 
{
    python.exec("from bitphylogeny.run import run_analysis")
    python.call("run_analysis", fin, fout, n, b, t)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
