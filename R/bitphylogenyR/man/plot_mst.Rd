\name{plot_mst}
\alias{plot_mst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot_mst}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_mst(genotype, label, flag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype}{
%%     ~~Describe \code{genotype} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{flag}{
%%     ~~Describe \code{flag} here~~
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
function (genotype, label, flag = FALSE) 
{
    reads <- sapply(unique(unlist(label)), function(ii) length(which(label == 
        ii)))
    mst <- get_mst(genotype)
    l = layout.reingold.tilford(graph = mst, root = which.min(rowSums(genotype)))
    nodes = matrix(0, dim(genotype)[1], 3)
    for (ii in 1:dim(nodes)[1]) {
        geno = c()
        for (jj in 1:dim(genotype)[2]) {
            geno <- paste(geno, genotype[ii, jj])
        }
        nodes[ii, 1] <- ii
        nodes[ii, 2] <- reads[ii]
        nodes[ii, 3] <- geno
    }
    colnames(nodes) = c("Nodes", "Read Counts", "Genotype")
    if (flag) {
        pdf(paste(p2, "/", fns[i], ".pdf", sep = ""), width = 15, 
            height = 9)
    }
    par(mfrow = c(1, 2), oma = c(3, 3, 0, 0) + 0.1, mar = c(0, 
        0, 1, 0.5))
    plot(mst, layout = l)
    textplot(nodes, show.rownames = F, cex = 1.25)
    if (flag) {
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
