\name{get_label_hc}
\alias{get_label_hc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get_label_hc
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_label_hc(x, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
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
function (x, K) 
{
    dis <- dist(x, "binary")
    hc_cand <- lapply(K, function(ii) cutree(hclust(dis), ii))
    hc_silhouette_res <- sapply(1:length(K), function(ii) summary(silhouette(hc_cand[[ii]], 
        dis))$avg.width)
    idx <- which.max(hc_silhouette_res)
    hc_label <- hc_cand[[idx]]
    clone <- sapply(unique(hc_label), function(i) which(hc_label == 
        i))
    n <- length(clone)
    genotype <- matrix(0, n, dim(x)[2])
    for (i in 1:n) {
        idx <- clone[[i]]
        if (length(idx) == 1) {
            genotype[i, ] = as.matrix(x[idx, ])
        }
        else {
            genotype[i, ] = as.numeric(colMeans(as.matrix(x[idx, 
                ])) > 0.5)
        }
    }
    return(list(hc_label = hc_label, hc_genotype = genotype))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
