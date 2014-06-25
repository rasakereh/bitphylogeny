\name{get_mst}
\alias{get_mst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get_mst
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_mst(genotype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype}{
%%     ~~Describe \code{genotype} here~~
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
function (genotype) 
{
    n <- dim(genotype)[1]
    adj <- hamming.distance(as.matrix(genotype))
    gr <- graph.adjacency(adjmatrix = adj, weighted = T, mode = "undirected")
    mst <- minimum.spanning.tree(gr, algorithm = "prim")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
