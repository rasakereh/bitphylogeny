library(cluster)
library(e1071)

get_label_kc <- function( data, K ){
  
  dis <- dist(data, 'binary')
  #is <- hamming.distance(as.matrix(data))
  
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
  #dis <- hamming.distance(as.matrix(data))
  
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


filepath <- '~/Dropbox/cancer-evolution/tssb/data/full_methy/' 
files <- dir(filepath, pattern = c('mutmat'))
fp1 <- '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster/hc/big-clone'
fp2 <- '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster/kc/big-clone'
K <- seq(2,20,1)

for (file in files){
  
  cat(sprintf('Processing %s ...\n', file))
  
  data <- read.csv(paste(filepath,file,sep=''))
  tmpname <- substr(file, nchar(file)-14, nchar(file)-11)
  true_label <- data[,9]
  data <- data[,-9]
  
  hcres = get_label_hc(data, K)
  
#   kcres = get_label_kc(data,K)
  
  write.csv(hcres$hc_label, file = paste(fp1, '/' ,'hc_labels_', tmpname, '.csv', 
                                         sep=''), row.names=F )
  write.csv(hcres$hc_genotype, file = paste(fp1, '/' ,'hc_genotype_', tmpname, '.csv', 
                                            sep=''), row.names=F )
#   write.csv(kcres$kc_label, file = paste(fp2, '/' ,'kc_labels_', tmpname, '.csv', 
#                                    sep=''), row.names=F )
#   write.csv(kcres$kc_genotype, file = paste(fp2, '/' ,'kc_genotype_', tmpname, '.csv', 
#                                          sep=''), row.names=F )
}


filepath <- '~/Dropbox/cancer-evolution/tssb/data/full_methy/small-clone/' 
files <- dir(filepath, pattern = c('mutmat'))
fp1 <- '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster/hc/small-clone'
fp2 <- '~/Dropbox/cancer-evolution/tssb/post-process/dump-cluster/kc/small-clone'
K <- seq(2,20,1)

for (file in files){
  
  cat(sprintf('Processing %s ...\n', file))
  
  data <- read.csv(paste(filepath,file,sep=''))
  tmpname <- substr(file, nchar(file)-14, nchar(file)-11)
  true_label <- data[,9]
  data <- data[,-9]
  
  hcres = get_label_hc(data, K)
  
#   kcres = get_label_kc(data,K)
  
  write.csv(hcres$hc_label, file = paste(fp1, '/' ,'hc_labels_', tmpname, '.csv', 
                                   sep=''), row.names=F )
  write.csv(hcres$hc_genotype, file = paste(fp1, '/' ,'hc_genotype_', tmpname, '.csv', 
                                            sep=''), row.names=F )
#   write.csv(kcres$kc_label, file = paste(fp2, '/' ,'kc_labels_', tmpname, '.csv', 
#                                          sep=''), row.names=F )
#   write.csv(kcres$kc_genotype, file = paste(fp2, '/' ,'kc_genotype_', tmpname, '.csv', 
#                                             sep=''), row.names=F )
}

