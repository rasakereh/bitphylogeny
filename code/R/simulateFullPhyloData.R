# TODO: Add comment
# 
# Author: warriour
###############################################################################


P=8

N=2000

mix<-c(0.3,0.1,0.3,0.07,0.14,0.05,0.04)

#print(sum(mix))

mutmat<-matrix(rep(0,P*N),nrow=N,byrow=T)

snvprop<-c(0.2,0.2,0.2,0.1,0.14,0.1,0.06)

#snvprop <- rdirichlet(1, rep(1, 7))

#print(sum(snvprop))

construc_p <- function(P,mix){
p1<-runif(P)<mix[1]

p2 <- p1

p2[ which(p1==0) ] <- runif(length(which(p1==0))) < mix[2]/(1-mix[1])

p3 <- p2

p3[ which(p2==0) ] <- runif(length(which(p2==0))) < mix[3]/(1-sum(mix[1:2]))

p4 <- p2


p4[ which(p3==0) ] <- runif(length(which(p3 == 0))) < mix[4]/(1-sum(mix[1:3]))

p5 <- p2

index <- intersect(which(p3 == 0), which(p4 == 0) )

p5[ index ] <- runif(length(index)) < mix[5]/(1-sum(mix[1:4]))

#print(index)

p6 <- p3

index <- intersect(index, which(p5 == 0)) 

p6[ index ] <- runif(length(index)) < mix[6]/(1-sum(mix[1:5]))

print(index)

p7 <- p5

index <- intersect(index, which(p6 == 0))

p7[index] <- 1

#print(index)

p <- cbind(p1,p2,p3,p4,p5,p6,p7)

return (p)
}

p = construc_p(P,snvprop)

print(p)

for(i in 1:N) {

       if( i/N < mix[1] ) {
             mutmat[i,] = p[,1]
 
	} else {  if  (i/N< sum(mix[1:2]) ){
          mutmat[i,] <- p[,2]
        } else {
                

        if (i/N<sum(mix[1:3])){
          mutmat[i,] <- p[,3]
          
        } else {  

        if (i/N<sum(mix[1:4])){
          mutmat[i,] <- p[,4]
          
        } else{ 

        if (i/N<sum(mix[1:5])){
          mutmat[i,] <- p[,5]
          
        } else {

        if (i/N<sum(mix[1:6])){
          mutmat[i,] <- p[,6]
          
        } else{
       
        mutmat[i,] <- p[,7] }}}}}}
         
        
}


#m <- rbind(mutmat[,c(1,3)],mutmat[,c(2,4)])

#freq <- matrix(rep(0,2*P),nrow=P)

#reads <- matrix(rep(0,3*P),nrow=P)

#for (i in 1:P){
  
#  freq[i,] <- c(i, sum(m[which(m[,1]==i),2])/length(which(m[,1]==i)))
  
#  reads[i,] <- c(i, sum(m[which(m[,1]==i),2]), length(which(m[,1]==i)))
#}


#hist(freq[,2],100)

#write.csv(reads, 'groundtsyn_reads.csv', row.names = FALSE)

write.csv(mutmat, 'fullsyn8_2000_mutmat.csv', row.names = FALSE)

#write.csv(freq, 'groundtsyn_freq.csv', row.names = FALSE)



#write.csv(p, 'groundtsyn_genotypes.csv', row.names = FALSE)

#write.csv(mutmat[1:floor(N*mix[1]),], 'groundtsyn_mutmat1.csv', row.names = FALSE)
#write.csv(mutmat[(floor(N*mix[1])+1):floor(N*sum(mix[1:2])),], 'groundtsyn_mutmat2.csv', row.names = FALSE)
#write.csv(mutmat[(floor(N*sum(mix[1:2]))+1):floor(N*sum(mix[1:3])),], 'groundtsyn_mutmat3.csv', row.names = FALSE)
#write.csv(mutmat[(floor(N*sum(mix[1:3]))+1):floor(N*sum(mix[1:4])),], 'groundtsyn_mutmat4.csv', row.names = FALSE)
#write.csv(mutmat[(floor(N*sum(mix[1:4]))+1):floor(N*sum(mix[1:5])),], 'groundtsyn_mutmat5.csv', row.names = FALSE)
#write.csv(mutmat[(floor(N*sum(mix[1:5]))+1):floor(N*sum(mix[1:6])),], 'groundtsyn_mutmat6.csv', row.names = FALSE)
#write.csv(mutmat[(floor(N*sum(mix[1:6]))+1):N,], 'groundtsyn_mutmat7.csv', row.names = FALSE)