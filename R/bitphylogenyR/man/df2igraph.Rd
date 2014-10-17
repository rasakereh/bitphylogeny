\name{df2igraph}
\alias{df2igraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert sankey data frame to igraph object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
df2igraph(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
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
function (df) 
{
    nodemat <- df$nodemat
    edgemat <- df$edgemat
    nn <- nodemat[, 1]
    g <- graph.empty() + vertices(nn)
    ee <- unlist(t(edgemat[, 1:2]))
    g <- g + edges(sapply(ee, function(i) which(nn == i)))
    return(g)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
