\name{bitphyloR}
\alias{bitphyloR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function to run BitPhylogeny analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bitphyloR(fin, fout, contains_true_label = FALSE, n = 100, b = 10, t = 1, 
mode = "methylation", rand_seed = 1234)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fin}{
  path to data file
}
  \item{fout}{
  path to result directory 
}
  \item{contains_true_label}{
  flag indicating whether the data file contains the true labels  
}
  \item{n}{
number of mcmc samples
}
  \item{b}{
number of burnin samples
}
  \item{t}{
thinning window
}
  \item{mode}{
type of data. Choose from "methylation" and "mutation". The default value is 
"methylation"
}
  \item{rand_seed}{
default 1234
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function doesn't return anything. Instead, it makes a result
  directory and saves MCMC traces. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ke Yuan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## example_file <- system.file('sample_data.csv', package='bitphylogenyR')
## bitphyloR(example_file, './output', T, 200, 50, 5)
## bitphyloR(example_file, './output', T, 200, 50, 5, mode = "mutation")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
