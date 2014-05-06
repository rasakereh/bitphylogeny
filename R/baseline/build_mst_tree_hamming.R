rm(list=ls())
library(igraph)
library(e1071)

get_mst <- function(genotype){
  n <- dim(genotype)[1]
  adj <- hamming.distance(as.matrix(genotype))
  gr <- graph.adjacency( adjmatrix= adj, weighted= T, mode='undirected')
  mst <- minimum.spanning.tree(gr, algorithm= 'prim')
}

plot_mst <- function(p1,p2){
  fns <- dir(p1, pattern='genotype')
  fns1 <- dir(p1, pattern='label') 
  n <- length(fns1)
  for ( i in 1:n ){
    genotype <- read.csv(paste(p1,fns[i],sep='/'))
    label <- read.csv(paste(p1,fns1[i],sep='/'))
    reads <- sapply(unique(unlist( label ) ), 
                    function(ii) length( which( label==ii )))
    mst <- get_mst(genotype)
    l = layout.reingold.tilford(graph=mst, 
                                root=which.min(rowSums(genotype)))
    nodes = matrix(0, dim(genotype)[1], 3 )
    for (ii in 1:dim(nodes)[1]){
      geno = c()
      for (jj in 1:dim(genotype)[2]){
        geno <- paste( geno, genotype[ii,jj] )
      }
      nodes[ii,1] <- ii
      nodes[ii,2] <- reads[ii]
      nodes[ii,3] <- geno
    }
    colnames(nodes) = c('Nodes','Read Counts','Genotype')
    
    pdf( paste(p2, '/', fns[i], '.pdf',sep=''),  width=15, height=9)
    par(mfrow = c(1, 2), oma = c(3,3,0,0) + 0.1,
        mar = c(0,0,1,0.5) )
    plot( mst, layout=l )
    textplot(nodes, show.rownames = F, cex=1.25)
    dev.off()
  }
}

p1 = '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster/kc/big-clone/'
p2 = '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster-pdf/kc/big-clone/'
plot_mst(p1,p2)

p1 = '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster/kc/small-clone/'
p2 = '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster-pdf/kc/small-clone/'
plot_mst(p1,p2)

p1 = '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster/hc/big-clone/'
p2 = '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster-pdf/hc/big-clone/'
plot_mst(p1,p2)

p1 = '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster/hc/small-clone/'
p2 = '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster-pdf/hc/small-clone/'
plot_mst(p1,p2)
