\name{get_genotype}
\alias{get_genotype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get genotype from sankey data frame
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_genotype(tree2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (tree2) 
{
    if (is.null(tree2$nodemat)) {
        g2 <- get_mst(tree2)
        t2 <- as.matrix(tree2)
    }
    else {
        g2 <- df2igraph(tree2)
        genomat <- tree2$genomat
        t2 <- unlist(strsplit(genomat[1, 3], " "))
        for (i in 2:dim(genomat)[1]) {
            t2 <- rbind(t2, unlist(strsplit(genomat[i, 3], " ")))
        }
        t2 = t2[, -1]
        rownames(t2) <- NULL
        class(t2) <- "numeric"
        t2 <- t2 > 0.5
    }
    class(t2) <- "numeric"
    return(list(t = t2, g = g2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
