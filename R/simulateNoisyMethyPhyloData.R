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
  clone1<-rep(0,num_snvs)
  
  # clone 2
   
  clone2 <- runif(num_snvs)<mix 
  clone2 <- (clone1+clone2)%%2	  

  # clone 3
  clone3 <- runif(num_snvs)<mix 
  clone3 <- (clone2+clone3)%%2	  
  
  # clone 4
  clone4 <- runif(num_snvs)<mix 
  clone4 <- (clone2+clone4)%%2	  
  
  # clone 5
  clone5 <- runif(num_snvs)<mix 
  clone5 <- (clone2+clone5)%%2	  
  
  # clone 6
  clone6 <- runif(num_snvs)<mix 
  clone6 <- (clone3+clone6)%%2	  
  
  # clone 7
  clone7 <- runif(num_snvs)<mix 
  clone7 <- (clone5+clone7)%%2	  
  
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

snvprop<-1/5

epslist <- c(0,1e-2,2e-2,5e-2)

clone = construct_genotype(num_snvs,snvprop)
while(ncol(t(unique(t(clone))))<7) {
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
	
	
	
	write.csv(cbind(mutmat,rep(1:7,reads_in_clone)), paste('noisy_full_methy_', num_snvs, '_', 
	                        num_reads, '_', epsilon, '_mutmat.csv', sep=""), 
	          row.names = FALSE)

}

write.csv(clone, paste('noisy_full_methy_', num_snvs, '_',
                       num_reads, '_genotype.csv', sep = ''),
          row.names = FALSE)
