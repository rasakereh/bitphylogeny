\name{gdl2df}
\alias{gdl2df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gdl2df%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gdl2df(file_gdl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_gdl}{
%%     ~~Describe \code{file_gdl} here~~
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
function (file_gdl) 
{
    filetext <- read.csv(file_gdl, as.is = T, quote = "")
    filetext <- filetext[-nrow(filetext), ]
    ftsplit = unlist(strsplit(filetext[1], " "))
    ID = unlist(strsplit(ftsplit[23], ":"))[2]
    layer = length(ID)
    x = as.numeric(ftsplit[18])
    y = 1000
    col = 0
    nodes = c(ID, x, y, col, layer)
    edges = c()
    for (ft in 2:length(filetext)) {
        node_sp = unlist(strsplit(filetext[ft], " "))
        if (node_sp[1] == "node:") {
            ID = unlist(strsplit(node_sp[24], ":"))[2]
            layer = length(unlist(strsplit(ID, "-")))
            x = as.numeric(node_sp[19])
            y = 0
            col = 0
            nodes = rbind(nodes, c(ID, x, y, col, layer))
            edge_sp = unlist(strsplit(filetext[ft + 1], " "))
            if (edge_sp[1] == "edge:") {
                ID1 = unlist(strsplit(edge_sp[3], ":"))[2]
                ID2 = unlist(strsplit(edge_sp[4], ":"))[2]
                node_label = unlist(strsplit(filetext[ft], " "))
                value = unlist(strsplit(node_label[4], ":"))[2]
                direction = "A"
                edges = rbind(edges, c(ID1, ID2, value, direction))
            }
        }
    }
    colnames(nodes) = c("ID", "x", "y", "col", "layer")
    rownames(nodes) = NULL
    colnames(edges) = c("ID1", "ID2", "value", "direction")
    for (k in 1:nrow(nodes)) {
        code = unlist(strsplit(nodes[k, "ID"], "-"))
        if (length(code) == 3) {
            index = which(nodes[, "layer"] == (length(code) - 
                1))
            for (id in index) {
                if (unlist(strsplit(nodes[id, "ID"], "-"))[length(code) - 
                  1] == code[length(code) - 1]) {
                  nodes[k, "x"] = as.numeric(nodes[k, "x"]) + 
                    as.numeric(nodes[id, "x"])
                  break
                }
            }
        }
        if (length(code) == 4) {
            index = which(nodes[, "layer"] == (length(code) - 
                1))
            for (id in index) {
                if ((unlist(strsplit(nodes[id, "ID"], "-"))[2] == 
                  code[2]) & (unlist(strsplit(nodes[id, "ID"], 
                  "-"))[3] == code[3])) {
                  nodes[k, "x"] = as.numeric(nodes[k, "x"]) + 
                    as.numeric(nodes[id, "x"])
                  break
                }
            }
        }
    }
    col_pa = c("coral2", "deepskyblue", "lightgreen", "pink2", 
        "black")
    for (i in unique(nodes[, "layer"])) {
        i = as.numeric(i)
        nodes[nodes[, "layer"] == as.character(i), "col"] = col_pa[i]
        if (i != 1) {
            no_i = sum(nodes[, "layer"] == i)
            lower = as.numeric(nodes[1, "y"]) - i * 0.1 * as.numeric(nodes[1, 
                "y"])
            upper = as.numeric(nodes[1, "y"]) + i * 0.1 * as.numeric(nodes[1, 
                "y"])
            nodes[nodes[, "layer"] == i, "y"] = seq(lower, upper, 
                length = no_i)
        }
    }
    nodes = nodes[, -ncol(nodes)]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
