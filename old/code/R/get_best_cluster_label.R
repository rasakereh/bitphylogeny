library(mcclust)

# filename <- ''
# label_traces <- read.csv(filename, header = F)
ltmat <- as.matrix(label_traces)
class(ltmat) <- 'integer'
ltmat <- ltmat + 1
psm <- comp.psm(ltmat)
mpear <- maxpear(psm)
best_label <- mpear$cl

# outfile = 'test.csv'
# write.csv(best_label, file=outfile, row.names = F)
