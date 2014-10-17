\name{write_vmeasure_traces}
\alias{write_vmeasure_traces}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
write_vmeasure_traces
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write_vmeasure_traces(filepath, ft)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{
%%     ~~Describe \code{filepath} here~~
}
  \item{ft}{
%%     ~~Describe \code{ft} here~~
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
function (filepath, ft) 
{
    labels <- load_label_traces(filepath)
    true_label <- load_true_label(ft)
    vmeasures <- t(sapply(1:dim(labels)[1], function(i) vmeasureR(labels[i, 
        ], true_label)))
    colnames(vmeasures) <- c("homogeneity", "completeness", "v_measure")
    write.csv(vmeasures, file = paste(filepath, "vmeasure_traces.csv", 
        sep = "/"), row.names = F)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
