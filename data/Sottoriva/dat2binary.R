file_name = dir('~/Dropbox/dp/tssb/data/Sottoriva/',pattern = '.dat')

repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}

for (nfile in file_name){
  
  x<- read.table(nfile, 
                 quote="\"", sep = ',')
  x = as.matrix(x)
  
  xx = sapply(1:length(x), function(ii) strsplit(x[ii,1], " " ))
  
  idx = vector()
  
  for (ii in 1:length(xx)){
    
    if ( length(xx[[ii]]) > 2 ){
      idx = rbind(idx, ii)
    }
    
  }
  
  data = vector('list', length(idx))
  
  for (ii in 2: length(idx)){
    selected = (idx[ii-1] + 1) : (idx[ii] - 1)
    tmp = vector()
    for (jj in selected){
      tmp = rbind(tmp, xx[[jj]])
    }
    data[[ii-1]] = tmp
  }
  
  tmp = vector()
  for (ii in ( idx[ length(idx) ] + 1 ) : length(xx) ){
    tmp = rbind(tmp, xx[[ii]])
  }
  
  data[[length(idx)]] = tmp
  
  tmp_name = substr(nfile, 1, nchar(nfile)-4)
  
  for (ii in 1:length(data)){
    t = data[[ii]]
    t1 = sapply(1:dim(t)[1], function(jj) unlist( strsplit( t[jj,1], '' )  )   )
    t1 = t(t1)
    t2 = cbind(t1,t[,2])
    class(t2) = 'numeric'
    
    tmp = vector()
    
    for (jj in 1:dim(t2)[1]){
      tmp1 = t(as.matrix(t2[jj,-dim(t2)[2]]))
      tmp = rbind(tmp, repmat(tmp1,t2[jj,dim(t2)[2]],1))
    }
    
    
    data[[ii]] = tmp
    
    fn = paste( tmp_name, '_', substr(xx[idx[ii]], 9,9),
      substr(xx[idx[ii]],4,4), '.csv', sep='')
    write.csv(tmp, fn, row.names=F)
  } 

}


