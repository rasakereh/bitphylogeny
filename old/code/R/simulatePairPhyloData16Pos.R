# TODO: Add comment
# 
# Author: warriour
###############################################################################


P=50

N=5000

mix<-c(0.3,0.1,0.3,0.07,0.14,0.05,0.04)

print(sum(mix))

datadim = 16

mutmat<-matrix(rep(0,2*datadim*N),nrow=N,byrow=T)

snvprop<-c(0.2,0.2,0.2,0.1,0.14,0.1,0.06)

print(sum(snvprop))



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

print(index)

p6 <- p3

index <- intersect(index, which(p5 == 0)) 

p6[ index ] <- runif(length(index)) < mix[6]/(1-sum(mix[1:5]))

print(index)

p7 <- p5

index <- intersect(index, which(p6 == 0))

p7[index] <- 1

print(index)

p <- cbind(p1,p2,p3,p4,p5,p6,p7)

return (p)
}

p = construc_p(P,snvprop)

summary(p)



for(i in 1:N) {

       if( i/N < mix[1] ) {
	     ind <- ceiling(runif(datadim)*P)
         mutmat[i,1:datadim] = ind
         mutmat[i,(datadim+1):(2*datadim)] = p[ind,1]
 
	} else {  if  (i/N< sum(mix[1:2]) ){
			ind <- ceiling(runif(datadim)*P)
			mutmat[i,1:datadim] = ind
			mutmat[i,(datadim+1):(2*datadim)] = p[ind,2]
        } else {
                

        if (i/N<sum(mix[1:3])){
			ind <- ceiling(runif(datadim)*P)
			mutmat[i,1:datadim] = ind
			mutmat[i,(datadim+1):(2*datadim)] = p[ind,3]
          
        } else {  

        if (i/N<sum(mix[1:4])){
			ind <- ceiling(runif(datadim)*P)
			mutmat[i,1:datadim] = ind
			mutmat[i,(datadim+1):(2*datadim)] = p[ind,4]
          
        } else{ 

        if (i/N<sum(mix[1:5])){
			ind <- ceiling(runif(datadim)*P)
			mutmat[i,1:datadim] = ind
			mutmat[i,(datadim+1):(2*datadim)] = p[ind,5]
          
        } else {

        if (i/N<sum(mix[1:6])){
			ind <- ceiling(runif(datadim)*P)
			mutmat[i,1:datadim] = ind
			mutmat[i,(datadim+1):(2*datadim)] = p[ind,6]
          
        } else{
       
			ind <- ceiling(runif(datadim)*P)
			mutmat[i,1:datadim] = ind
			mutmat[i,(datadim+1):(2*datadim)] = p[ind,7] }}}}}}
         
        
}

m<-mutmat[,c(1,datadim+1)]
for ( i in 2:datadim) {
	m <- rbind(m, mutmat[,c(i,datadim+i)])
}
freq <- matrix(rep(0,2*P),nrow=P)

for (i in 1:P){
  
  freq[i,] <- c(i, sum(m[which(m[,1]==i),2])/length(which(m[,1]==i)))
  
}



write.csv(mutmat, '../data/syn_mutmat_5000_50_16.csv', row.names = FALSE)

write.csv(freq, '../data/syn_freq_5000_50_16.csv', row.names = FALSE)
