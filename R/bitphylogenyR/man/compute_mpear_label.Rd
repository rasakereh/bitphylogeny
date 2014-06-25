\name{compute_mpear_label}
\alias{compute_mpear_label}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{compute_mpear_label
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_mpear_label(label_traces)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{label_traces}{
%%     ~~Describe \code{label_traces} here~~
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
function (label_traces) 
{
    ltmat <- as.matrix(label_traces)
    class(ltmat) <- "integer"
    ltmat <- ltmat + 1
    psm <- comp.psm(ltmat)
    mpear <- maxpear(psm)
    mpear_label <- mpear$cl
    return(mpear_label)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
