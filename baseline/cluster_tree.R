rm(list=ls())
library(cluster)
library(e1071)
library(igraph)
library(gplots)

get_label_kc <- function( data, K ){
  
  dis <- dist(data, 'binary')
  
  kc_cand <- lapply(K, function(ii) pam( dis, ii) )
  
  kc_silhouette_res <- sapply(1:length(K), 
                              function(ii) 
                                summary( silhouette(kc_cand[[ii]]$clustering,dis) )$avg.width )
  idx <- which.max( kc_silhouette_res )
  
  kc_label <- kc_cand[[idx]]$clustering
  kc_genotype <- data[kc_cand[[idx]]$medoids,]
  
  return(list(kc_label = kc_label, kc_genotype = kc_genotype))
}

get_label_hc <- function( data, K ){
  
  dis <- dist(data, 'binary')
  
  hc_cand <- lapply(K, function(ii) cutree( hclust( dis ), ii ) )
  
  hc_silhouette_res <- sapply(1:length(K), 
                              function(ii) 
                                summary( silhouette(hc_cand[[ii]] ,dis) )$avg.width )
  idx <- which.max( hc_silhouette_res )
  
  hc_label <- hc_cand[[idx]]
  
  clone <- sapply(unique(hc_label), function(i) which(hc_label==i) )
  
  n <- length(clone)
  
  genotype <- matrix(0, n, dim(data)[2])
  
  for (i in 1:n){
    idx <- clone[[i]]
    if ( length(idx)==1 ){
      genotype[i,] = as.matrix(data[idx,])
    }else{
      genotype[i,] = as.numeric( colMeans(as.matrix(data[idx,])) > 0.5 )
    }
  }
  
  return(list(hc_label = hc_label, hc_genotype = genotype))
  
}

## clustering ##################################################################

filepath <- '../data/synthetic/big-clone/' 
files <- dir(filepath, pattern = c('mutmat'))
fp1 <- './cluster//hc//big-clone'
fp2 <- './cluster//kc//big-clone'
K <- seq(2,20,1)

for (file in files){
  
  cat(sprintf('Processing %s ...\n', file))
  
  data <- read.csv(paste(filepath,file,sep=''))
  tmpname <- substr(file, nchar(file)-14, nchar(file)-11)
  true_label <- data[,9]
  data <- data[,-9]
  
  hcres = get_label_hc(data, K)
  
  kcres = get_label_kc(data,K)
  
  write.csv(hcres$hc_label, file = paste(fp1, '/' ,'hc_labels_', tmpname, '.csv', 
                                         sep=''), row.names=F )
  write.csv(hcres$hc_genotype, file = paste(fp1, '/' ,'hc_genotype_', tmpname, '.csv', 
                                            sep=''), row.names=F )
  write.csv(kcres$kc_label, file = paste(fp2, '/' ,'kc_labels_', tmpname, '.csv', 
                                   sep=''), row.names=F )
  write.csv(kcres$kc_genotype, file = paste(fp2, '/' ,'kc_genotype_', tmpname, '.csv', 
                                         sep=''), row.names=F )
}


filepath <- '../data/synthetic/small-clone/' 
files <- dir(filepath, pattern = c('mutmat'))
fp1 <- './cluster//hc//small-clone'
fp2 <- './cluster//kc//small-clone'
K <- seq(2,20,1)

for (file in files){
  
  cat(sprintf('Processing %s ...\n', file))
  
  data <- read.csv(paste(filepath,file,sep=''))
  tmpname <- substr(file, nchar(file)-14, nchar(file)-11)
  true_label <- data[,9]
  data <- data[,-9]
  
  hcres = get_label_hc(data, K)
  
  kcres = get_label_kc(data,K)
  
  write.csv(hcres$hc_label, file = paste(fp1, '/' ,'hc_labels_', tmpname, '.csv', 
                                   sep=''), row.names=F )
  write.csv(hcres$hc_genotype, file = paste(fp1, '/' ,'hc_genotype_', tmpname, '.csv', 
                                            sep=''), row.names=F )
  write.csv(kcres$kc_label, file = paste(fp2, '/' ,'kc_labels_', tmpname, '.csv', 
                                         sep=''), row.names=F )
  write.csv(kcres$kc_genotype, file = paste(fp2, '/' ,'kc_genotype_', tmpname, '.csv', 
                                            sep=''), row.names=F )
}

## tree building ##################################################################

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

cat(sprintf('Plotting trees ...\n'))
## kc big-clone trees
p1 = './cluster//kc//big-clone/'
p2 = './tree//kc//big-clone/'
plot_mst(p1,p2)


## kc small-clone trees
p1 = './cluster//kc//small-clone/'
p2 = './tree//kc//small-clone/'
plot_mst(p1,p2)

## hc big-clone trees
p1 = './cluster//hc//big-clone/'
p2 = './tree//hc//big-clone/'
plot_mst(p1,p2)

## hc small-clone trees
p1 = './cluster//hc//small-clone/'
p2 = './tree//hc//small-clone/'
plot_mst(p1,p2)
