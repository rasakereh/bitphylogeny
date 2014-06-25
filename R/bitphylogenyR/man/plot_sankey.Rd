\name{plot_sankey}
\alias{plot_sankey}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_sankey%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_sankey(nodes, edges_mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodes}{
%%     ~~Describe \code{nodes} here~~
}
  \item{edges_mat}{
%%     ~~Describe \code{edges_mat} here~~
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
function (nodes, edges_mat) 
{
    nn <- nodes[, 1]
    g <- graph.empty() + vertices(nn)
    ee <- unlist(t(edges_mat[, 1:2]))
    ee <- edges(sapply(ee, function(i) which(nn == i)))
    g <- g + edges(ee)
    ll <- layout.reingold.tilford(g)
    nodes$y = l[, 1] * 100
    river <- makeRiver(nodes, edges)
    style <- list(edgecol = "col")
    riverplot(river, srt = 90, lty = 1, default_style = style)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
