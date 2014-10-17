\name{get_tree_distance_common}
\alias{get_tree_distance_common}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute consensus node based shortest path distance between four trees
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_tree_distance_common(tree1, tree2, tree3, tree4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree1}{
%%     ~~Describe \code{tree1} here~~
}
  \item{tree2}{
%%     ~~Describe \code{tree2} here~~
}
  \item{tree3}{
%%     ~~Describe \code{tree3} here~~
}
  \item{tree4}{
%%     ~~Describe \code{tree4} here~~
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
function (tree1, tree2, tree3, tree4) 
{
    tmp <- get_genotype(tree1)
    t1 <- tmp$t
    g1 <- tmp$g
    tmp <- get_genotype(tree2)
    t2 <- tmp$t
    g2 <- tmp$g
    tmp <- get_genotype(tree3)
    t3 <- tmp$t
    g3 <- tmp$g
    tmp <- get_genotype(tree4)
    t4 <- tmp$t
    g4 <- tmp$g
    ind1 <- try(get_genotype_intercept(t1, t2))
    ind2 <- try(get_genotype_intercept(t1, t3))
    ind3 <- try(get_genotype_intercept(t1, t4))
    ind <- Reduce(intersect, list(ind1$ind1, ind2$ind1, ind3$ind1))
    tt <- t1[ind, ]
    ind1 <- get_genotype_intercept(tt, t2)
    ind2 <- get_genotype_intercept(tt, t3)
    ind3 <- get_genotype_intercept(tt, t4)
    mat1 <- shortest.paths(g1)
    mat2 <- shortest.paths(g2)
    mat3 <- shortest.paths(g3)
    mat4 <- shortest.paths(g4)
    submat1 <- mat1[ind, ind]
    submat2 <- mat2[ind1$ind2, ind1$ind2]
    submat3 <- mat3[ind2$ind2, ind2$ind2]
    submat4 <- mat4[ind3$ind2, ind3$ind2]
    m1 <- sum(abs(submat1[lower.tri(submat1)] - submat2[lower.tri(submat2)]))
    m2 <- sum(abs(submat1[lower.tri(submat1)] - submat3[lower.tri(submat3)]))
    m3 <- sum(abs(submat1[lower.tri(submat1)] - submat4[lower.tri(submat4)]))
    return(list(m = c(m1, m2, m3), n = length(ind)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
