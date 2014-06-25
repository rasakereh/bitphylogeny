
run_baseline <- function(output){ 
## clustering ##################################################################

filepath <- system.file('extdata/synthetic/big-clone/',
                        package='bitphylogenyR')
files <- dir(filepath, pattern = c('mutmat'))
fp1 <- paste(output, '/cluster/hc/big-clone', sep='')
fp2 <- paste(output, '/cluster/kc/big-clone', sep='')
K <- seq(2,20,1)

for (file in files){
  
  cat(sprintf('Processing %s ...\n', file))
  
  data <- read.csv(paste(filepath,file,sep=''))
  tmpname <- substr(file, nchar(file)-14, nchar(file)-11)
  true_label <- data[,9]
  data <- data[,-9]
  
  hcres = get_label_hc(data, K)
  
  kcres = get_label_kc(data,K)
  
  write.csv(hcres$hc_label, file = paste(fp1, '/' ,'hc_labels_',
                                         tmpname, '.csv', 
                                         sep=''), row.names=F )
  write.csv(hcres$hc_genotype, file = paste(fp1, '/' ,'hc_genotype_',
                                            tmpname, '.csv', 
                                            sep=''), row.names=F )
  write.csv(kcres$kc_label, file = paste(fp2, '/' ,'kc_labels_',
                                         tmpname, '.csv', 
                                         sep=''), row.names=F )
  write.csv(kcres$kc_genotype, file = paste(fp2, '/' ,'kc_genotype_',
                                            tmpname, '.csv', 
                                            sep=''), row.names=F )
}


filepath <- system.file('extdata/synthetic/small-clone/',
                        package='bitphylogenyR')
files <- dir(filepath, pattern = c('mutmat'))
fp1 <- paste(output, '/cluster/hc/small-clone', sep='')
fp2 <- paste(output, '/cluster/kc/small-clone', sep='')
K <- seq(2,20,1)

for (file in files){
  
  cat(sprintf('Processing %s ...\n', file))
  
  data <- read.csv(paste(filepath,file,sep=''))
  tmpname <- substr(file, nchar(file)-14, nchar(file)-11)
  true_label <- data[,9]
  data <- data[,-9]
  
  hcres = get_label_hc(data, K)
  
  kcres = get_label_kc(data,K)
  
  write.csv(hcres$hc_label, file = paste(fp1, '/' ,'hc_labels_',
                                tmpname, '.csv', 
                                   sep=''), row.names=F )
  write.csv(hcres$hc_genotype, file = paste(fp1, '/' ,'hc_genotype_',
                                   tmpname, '.csv', 
                                            sep=''), row.names=F )
  write.csv(kcres$kc_label, file = paste(fp2, '/' ,'kc_labels_',
                                tmpname, '.csv', 
                                         sep=''), row.names=F )
  write.csv(kcres$kc_genotype, file = paste(fp2, '/' ,'kc_genotype_',
                                   tmpname, '.csv', 
                                            sep=''), row.names=F )
}

## tree building ###############################################################

cat(sprintf('Plotting trees ...\n'))
## kc big-clone trees
p1 = paste(output, '/cluster/kc/big-clone/', sep='')
p2 = paste(output, '/tree/kc/big-clone/', sep='')
plot_mst(p1,p2)


## kc small-clone trees
p1 = paste(output, '/cluster/kc/small-clone/', sep='')
p2 = paste(output, '/tree/kc/small-clone/', sep='')
plot_mst(p1,p2)

## hc big-clone trees
p1 = paste(output, '/cluster/hc/big-clone/', sep='')
p2 = paste(output, '/tree/hc/big-clone/', sep='')
plot_mst(p1,p2)

## hc small-clone trees
p1 = paste(output, '/cluster/hc/small-clone/', sep='')
p2 = paste(output, '/tree/hc/small-clone/', sep='')
plot_mst(p1,p2)
}
