\name{run_baseline}
\alias{run_baseline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run hc and kc pipline for installed synthetic data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
run_baseline(output, tree_type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
%%     ~~Describe \code{output} here~~
}
  \item{tree_type}{
%%     ~~Describe \code{tree_type} here~~
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
function (output, tree_type) 
{
    filepath <- system.file(paste("extdata/synthetic/", tree_type, 
        "/", sep = ""), package = "bitphylogenyR")
    files <- dir(filepath, pattern = c("mutmat"))
    fp1 <- paste(output, "/cluster/hc/", tree_type, sep = "")
    fp2 <- paste(output, "/cluster/kc/", tree_type, sep = "")
    K <- seq(2, 20, 1)
    for (file in files) {
        cat(sprintf("Processing \%s ...\n", file))
        data <- read.csv(paste(filepath, file, sep = ""))
        tmpname <- substr(file, nchar(file) - 14, nchar(file) - 
            11)
        true_label <- data[, 9]
        data <- data[, -9]
        hcres <- get_label_hc(data, K)
        kcres <- get_label_kc(data, K)
        hc_vmeasure <- vmeasureR(hcres$hc_label, true_label)
        kc_vmeasure <- vmeasureR(kcres$kc_label, true_label)
        write.csv(hcres$hc_label, file = paste(fp1, "/", "hc_labels_", 
            tmpname, ".csv", sep = ""), row.names = F)
        write.csv(hc_vmeasure, file = paste(fp1, "/", "hc_vmeasure_", 
            tmpname, ".csv", sep = ""), row.names = F)
        write.csv(hcres$hc_genotype, file = paste(fp1, "/", "hc_genotype_", 
            tmpname, ".csv", sep = ""), row.names = F)
        write.csv(kcres$kc_label, file = paste(fp2, "/", "kc_labels_", 
            tmpname, ".csv", sep = ""), row.names = F)
        write.csv(kc_vmeasure, file = paste(fp2, "/", "kc_vmeasure_", 
            tmpname, ".csv", sep = ""), row.names = F)
        write.csv(kcres$kc_genotype, file = paste(fp2, "/", "kc_genotype_", 
            tmpname, ".csv", sep = ""), row.names = F)
    }
    cat(sprintf("Plotting trees ...\n"))
    p1 = paste(output, "/cluster/kc/", tree_type, "/", sep = "")
    p2 = paste(output, "/tree/kc/", tree_type, "/", sep = "")
    plot_mst(p1, p2)
    p1 = paste(output, "/cluster/hc/", tree_type, "/", sep = "")
    p2 = paste(output, "/tree/hc/", tree_type, "/", sep = "")
    plot_mst(p1, p2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
