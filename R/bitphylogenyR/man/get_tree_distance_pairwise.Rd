\name{get_tree_distance_pairwise}
\alias{get_tree_distance_pairwise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute consensus node based shortest path distance between two trees.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_tree_distance_pairwise(tree1, tree2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree1}{
%%     ~~Describe \code{tree1} here~~
}
  \item{tree2}{
%%     ~~Describe \code{tree2} here~~
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
function (tree1, tree2) 
{
    tmp <- get_genotype(tree1)
    t1 <- tmp$t
    g1 <- tmp$g
    tmp <- get_genotype(tree2)
    t2 <- tmp$t
    g2 <- tmp$g
    ind <- get_genotype_intercept(t1, t2)
    ind1 <- ind$ind1
    ind2 <- ind$ind2
    mat1 <- shortest.paths(g1)
    mat2 <- shortest.paths(g2)
    submat1 <- mat1[ind1, ind1]
    submat2 <- mat2[ind2, ind2]
    m = sum(abs(submat1[lower.tri(submat1)] - submat2[lower.tri(submat2)]))
    return(list(m = m, n = length(ind1)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
