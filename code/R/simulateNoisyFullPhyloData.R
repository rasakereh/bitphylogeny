# TODO: Add comment
# 
# Author: warriour
###############################################################################

construct_genotype <- function(num_snvs,mix){
  # Mannually construct the genotypes 7 clones phylogeny 
  # tree:
  #        clone1
  #          |
  #        clone2
  #          |   
  # clone3 clone4 clone5
  #   |             |
  # clone6        clone7
  # 
  # Inputs: 
  #   num_snvs: number of SNVs
  #   mix:      true clonal composition 
  
  # clone 1
  clone1<-runif(num_snvs)<mix[1]
  
  # clone 2
  clone2 <- clone1 
  clone2[ which(clone1==0) ] <- runif(length(which(clone1==0))) < 
    mix[2]/(1-mix[1])
  
  # clone 3
  clone3 <- clone2 
  clone3[ which(clone2==0) ] <- runif(length(which(clone2==0))) < 
    mix[3]/(1-sum(mix[1:2]))
  
  # clone 4
  clone4 <- clone2 
  clone4[ which(clone3==0) ] <- runif(length(which(clone3 == 0))) < 
    mix[4]/(1-sum(mix[1:3]))
  
  # clone 5
  clone5 <- clone2 
  index <- intersect(which(clone3 == 0), which(clone4 == 0) ) 
  clone5[ index ] <- runif(length(index)) < mix[5]/(1-sum(mix[1:4]))
  
  # clone 6
  clone6 <- clone3 
  index <- intersect(index, which(clone5 == 0)) 
  clone6[ index ] <- runif(length(index)) < mix[6]/(1-sum(mix[1:5]))
  
  # clone 7
  clone7 <- clone5 
  index <- intersect(index, which(clone6 == 0))
  clone7[index] <- 1
  
  clone <- cbind(clone1,clone2,clone3,clone4,clone5,clone6,clone7)
  return(clone)
}

sample_full_reads <- function(n, genotype, epsilon){
  # Sampling of full reads data
  # Inputs:
  #   n:        number of reads
  #   genotype: true binary genotype
  #   epsilon:  noise level
  # Outputs
  #   reads: noisy reads data
  
  p = length(genotype)
  reads <- matrix(0,n,p)
  for (ii in 1:n){
    flip = runif(p) < epsilon
    idx = which(flip == 1)
    tmp = genotype
    tmp[idx] = !tmp[idx]
    reads[ii,] = tmp
  }
  return(reads)
}

num_snvs=8

num_reads=2000

mix<-c(0.3,0.1,0.3,0.07,0.14,0.05,0.04)

snvprop<-c(0.2,0.2,0.2,0.1,0.14,0.1,0.06)

epslist <- c(0,1e-2,2e-2,5e-2)

clone = construct_genotype(num_snvs,snvprop)
while((ncol(t(unique(t(clone))))<7) || (sum(clone[,1])<1) ) {
	clone = construct_genotype(num_snvs,snvprop)
}

reads_in_clone <- rmultinom(1, num_reads, mix)   

for(epsilon in epslist) {
	# sample reads for each clone
	reads <- sapply(1:length(mix), 
	                function(ii) 
	                  sample_full_reads(reads_in_clone[ii], 
	                                    clone[,ii], epsilon) )
	
	mutmat <- vector()
	for ( ii in 1:length(mix) ){
	  mutmat <- rbind(mutmat, reads[[ii]])
	}
	
	
	
	write.csv(cbind(mutmat,rep(1:7,reads_in_clone)), paste('noisy_fullsyn_', num_snvs, '_', 
	                        num_reads, '_', epsilon, '_mutmat.csv', sep=""), 
	          row.names = FALSE)

}

write.csv(clone, paste('noisy_fullsyn_', num_snvs, '_',
                       num_reads, '_genotype.csv', sep = ''),
          row.names = FALSE)
