\name{get_common_tree_distance_df}
\alias{get_common_tree_distance_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get common consensus node based shortest path distance from synthetic data 
results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_common_tree_distance_df(f1, f2, f3, treeType, f4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f1}{
%%     ~~Describe \code{f1} here~~
}
  \item{f2}{
%%     ~~Describe \code{f2} here~~
}
  \item{f3}{
%%     ~~Describe \code{f3} here~~
}
  \item{treeType}{
%%     ~~Describe \code{treeType} here~~
}
  \item{f4}{
%%     ~~Describe \code{f4} here~~
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
function (f1, f2, f3, treeType, f4) 
{
    folders <- dir(f1, pattern = "mutmat", full.names = TRUE, 
        recursive = TRUE, include.dirs = TRUE)
    if (treeType == "star-clone") {
        folders <- folders[c(2, 4, 6, 8)]
    }
    folders = folders[c(4, 1, 2, 3)]
    pathToSynth <- "extdata/synthetic/"
    treeType <- treeType
    pathToFiles <- system.file(paste(pathToSynth, treeType, sep = ""), 
        package = "bitphylogenyR")
    files <- dir(pathToFiles, pattern = "mutmat", full.names = T)
    true <- gdl2df(file_gdl = f4)
    hc_genotype_files <- dir(path = f2, pattern = "hc_genotype", 
        recursive = T, full.names = T)
    kc_genotype_files <- dir(path = f3, pattern = "kc_genotype", 
        recursive = T, full.names = T)
    errors = c("0", "0.01", "0.02", "0.05")
    bitphylo_trees = vector("list", length = 4)
    for (i in 1:4) {
        fh <- paste(folders[i], "/treescripts/", sep = "")
        treefreq <- read.csv(dir(path = fh, pattern = "tree-freq", 
            recursive = T, full.names = T))
        mft <- treefreq[which.max(treefreq[, "freq"]), 1]
        fn <- paste(fh, "/nodes-", mft, ".gdl", sep = "")
        bitphylo_trees[[i]] <- gdl2df(fn)
    }
    distance = vector("numeric", length = 12)
    matched_clone = distance
    for (i in 1:4) {
        res <- get_tree_distance_common(true, bitphylo_trees[[i]], 
            read.csv(kc_genotype_files[i]), read.csv(hc_genotype_files[i]))
        distance[((i - 1) * 3 + 1):(i * 3)] = res$m
        matched_clone[((i - 1) * 3 + 1):(i * 3)] = res$n
    }
    metrics <- data.frame(distance = distance, method = c("BitPhylogeny", 
        "k-centroids", "hierarchical clustering"), tree_type = treeType, 
        matched_clone = matched_clone, error = c(rep("0", 3), 
            rep("0.01", 3), rep("0.02", 3), rep("0.05", 3)))
    return(metrics)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
