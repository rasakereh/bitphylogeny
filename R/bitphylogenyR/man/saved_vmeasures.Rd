\name{saved_vmeasures}
\alias{saved_vmeasures}
\docType{data}
\title{
A list of V-measures reported in the manuscript
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(saved_vmeasures)}
\format{
  The format is:
List of 4
 $ mcmc_vmeasures :List of 2
  ..$ big_clone  : num [1:10000, 1:4] 0.909 0.909 0.909 0.907 0.906 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:4] "0" "0.01" "0.02" "0.05"
  ..$ small_clone: num [1:10000, 1:4] 0.828 0.825 0.824 0.819 0.824 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:4] "0" "0.01" "0.02" "0.05"
 $ hc_vmeasures   :List of 2
  ..$ big_clone  : Named num [1:4] 1 0.856 0.794 0.595
  .. ..- attr(*, "names")= chr [1:4] "0" "0.01" "0.02" "0.05"
  ..$ small_clone: Named num [1:4] 1 0.816 0.755 0.581
  .. ..- attr(*, "names")= chr [1:4] "0" "0.01" "0.02" "0.05"
 $ kc_vmeasures   :List of 2
  ..$ big_clone  : Named num [1:4] 1 0.9 0.839 0.701
  .. ..- attr(*, "names")= chr [1:4] "0" "0.01" "0.02" "0.05"
  ..$ small_clone: Named num [1:4] 1 0.908 0.828 0.679
  .. ..- attr(*, "names")= chr [1:4] "0" "0.01" "0.02" "0.05"
 $ mpear_vmeasures:List of 2
  ..$ big_clone  : Named num [1:4] 0.909 0.928 0.879 0.825
  .. ..- attr(*, "names")= chr [1:4] "0" "0.01" "0.02" "0.05"
  ..$ small_clone: Named num [1:4] 0.89 0.914 0.879 0.699
  .. ..- attr(*, "names")= chr [1:4] "0" "0.01" "0.02" "0.05"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(saved_vmeasures)
## maybe str(saved_vmeasures) ; plot(saved_vmeasures) ...
}
\keyword{datasets}
