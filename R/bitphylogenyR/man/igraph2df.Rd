\name{igraph2df}
\alias{igraph2df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
igraph2df
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
igraph2df(g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
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
function (g) 
{
    dataFrame <- get.data.frame(g, what = "both")
    nodes <- dataFrame$vertices
    edges <- dataFrame$edges
    d1 <- strsplit(nodes$name, "-")
    depth <- sapply(1:length(d1), function(i) length(d1[[i]])) - 
        1
    nodes$layer <- depth
    nodes$y <- layout.reingold.tilford(g)[, 1]
    nodes <- rename(nodes, c(name = "ID", branch = "x"))
    for (k in 1:nrow(nodes)) {
        code = unlist(strsplit(nodes[k, "ID"], "-"))
        if (length(code) >= 3) {
            index = which(nodes[, "layer"] == (length(code) - 
                1))
            for (id in index) {
                code1 = unlist(strsplit(nodes[id, "ID"], "-"))
                tt = sapply(2:length(code1), function(i) code1[i] == 
                  code[i])
                if (sum(tt) == (length(code1) - 1)) {
                  nodes[k, "x"] = as.numeric(nodes[k, "x"]) + 
                    as.numeric(nodes[id, "x"])
                  break
                }
            }
        }
    }
    palette = gg_color_hue(length(unique(depth)))
    nodes$col <- palette[factor(nodes$layer)]
    nn <- subset(nodes, select = c(ID, x, y, col))
    edges <- rename(edges, c(from = "N1", to = "N2"))
    edges <- subset(edges, select = c(N1, N2, Value))
    edges$direction <- "A"
    edges$col <- "gray90"
    genomat <- subset(nodes, select = c(ID, size, params))
    return(list(nodemat = nn, edgemat = edges, genomat = genomat))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
