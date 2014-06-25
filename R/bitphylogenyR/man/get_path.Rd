\name{get_path}
\alias{get_path}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get_path
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_path(folder, key)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folder}{
%%     ~~Describe \code{folder} here~~
}
  \item{key}{
%%     ~~Describe \code{key} here~~
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
function (folder, key) 
{
    trace <- dir(folder, recursive = T, include.dirs = T, pattern = key)
    filepath <- paste(folder, trace, sep = "/")
    return(filepath)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
