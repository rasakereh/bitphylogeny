## EM for Bernoulli and logistic mixtures
## Author: Ke Yuan 

rdirichlet <- function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}

logsumexp <- function(x, y, dim=1){
  ## todo: Bugs when dim = 2
  x <- as.matrix(x)
  nargin <- length(as.list(match.call())) - 1
  if (nargin == 1){
      rdim <- dim(x)[2:1] 
      dim <- which(rdim!=1)[1] 
  }  
  y <- apply(x, dim, max)
  x <- x - y
  s <- y + log(apply(exp(x), dim, sum))
  return(s)
}

repmat <- function(x, n, m){
  x <- as.matrix(x)
  return(kronecker(matrix(1,n,m),x))
}

sigmoid <- function(x){
  return(1/(1+exp(-x)))
}

sigmoidln <- function(x){
  return(-log(1+exp(-x)))
}

constriant <- function(x){
  if (length(which(x > 0.9999) ) >0 ){
    x[which(x>0.9999)] = 0.9999
  }
  if (length(which(x < 1e-20)) > 0){
    x[which(x<1e-20)] = 1e-20
  }
  return(x)
}

logit <- function(x){
  x <- constriant(x)
  return(-log(1/x-1))
}

logit_llh <- function(x, y){
  x = t(as.matrix(x))
  x = repmat(x, dim(y)[1], 1) 
  llh = y*sigmoidln(x) + (1-y)*sigmoidln(-x)
  return(rowSums(llh, na.rm=T))
}

bernoulli_llh <- function(x,y){
  x = t(as.matrix(x))
  
  x = repmat(x, dim(y)[1], 1) 
  llh = y*log(x) + (1-y)*log(1-x)
  return(rowSums(llh))
}


compute_true_llh <- function(data, func, Epi, Eparam){
  
  N <- dim(data)[1]
  K <- length(Epi)
  
  
  vecEpi <- repmat(Epi, N, 1)
  veclogEpi <- log(vecEpi)
  
  fxln <- matrix(0, N, K)
  
  for (k in 1:K){
    fxln[,k] <- do.call(func, list(Eparam[k,], data))
  }  
  
  return( sum(logsumexp(veclogEpi + fxln, 1)) )
  
}


## EM

expectation_maximization <- function(data, K, func, eps, totem){
  
  ## Initial 
  P <- dim(data)[2]
  N <- dim(data)[1]
  Eparam <- matrix(runif(K*P), K, P)
  a <- rep(0.5, K)
  Epi <- rdirichlet(1, a)
  vecEpi <- repmat(Epi, N, 1)
  veclogEpi <- log(vecEpi)
  fxln <- matrix(0, N, K)
  llh <- rep(0,totem)
  
  for (nem in 1:totem){
    
    ## E-step: compute responsabilies
    for (k in 1:K){
      fxln[,k] <- do.call(func, list(Eparam[k,], data))
    }  
    llh_tmp <- logsumexp(veclogEpi + fxln, 1)
    llh_tmp <- repmat(llh_tmp, 1, K)
    resln <- veclogEpi + fxln - llh_tmp
    
    ## E-step: compute log-marginal
    llh[nem] <- sum(llh_tmp[,1]) 
    res <- exp(resln)
    
    
    ## M-step: compute node params
    Nk <- colSums(res)
    vecNk <- repmat(Nk, 1, P)
    data <- as.matrix(data)
    xkbar <- t(res)%*%data/vecNk
    
    if (func == 'logit_llh'){
      Eparam <- logit(xkbar)
    }
    else {
      Eparam <- constriant(xkbar)
    }
    
    ## M-step: compute weights
    Epi <- Nk/N
    vecEpi <- t(repmat(Epi, 1, N))
    veclogEpi <- log(vecEpi)  
    
    ## stopping condition
    if (nem>2){
      gain = llh[nem] - llh[nem-1]
      cat(sprintf("\rEM steps: %i, log-marginal likelihood: %f, relative gain: %f\r", 
                    nem, llh[nem], abs(gain/llh[nem-1])))
      if (abs( gain/llh[nem-1] ) < eps && gain > 0) {
        cat(sprintf('\nconverged.\n'))
        break
      }
    }
  }
  
  pred_label <- apply(res, 1, which.max)
  
  return( list(pred_label=pred_label, llh=llh[1:nem], Eparam=Eparam, Epi=Epi) )
  
}

bic <- function(l, k, n){
  return( -2*l + k*( log(n)+log(2*pi) ) )
}

aic <- function(l, k, n){
  return( -2*l + 2*k*(k+1)/(n-k-1) )
}


filepath <- '~/Dropbox/cancer-evolution/tssb/data/full_methy/' 
files <- dir(filepath, pattern = c('mutmat'))

bic_bernoulli <- vector('list', length(files) )
aic_bernoulli <- vector('list', length(files) )
bic_logit <- vector('list', length(files) )
aic_logit <- vector('list', length(files) )
em_bernoulli <- vector('list', length(files) )
em_logit <- vector('list', length(files) )
true_llh_bernoulli <- vector('numeric', length(files))
true_llh_logit <- vector('numeric', length(files))

true_Epi <- cbind(0.3,0.1,0.3,0.07,0.14,0.05,0.04)
true_genotype <- read.csv(paste(filepath,
                                'noisy_full_methy_8_2000_genotype.csv',sep=''))
true_Eparam <- t(constriant(as.matrix(true_genotype)))

i = 1

for (file in files){

  cat(sprintf('Processing %s ...\n', file))
  data <- read.csv(paste(filepath,file,sep=''))
  
  tmpname <- substr(file, nchar(file)-15, nchar(file)-11)
  
  true_label <- data[,9]
  data <- data[,-9]
  K <- seq(2,14,1)
  eps <- 1e-5
  totem <- 1000
  
  
  true_llh_logit[i] <- compute_true_llh(data, 
                                     'logit_llh', 
                                     true_Epi, 
                                     logit(true_Eparam))
  true_llh_bernoulli[i] <- compute_true_llh(data, 
                                         'bernoulli_llh', 
                                         true_Epi, 
                                         true_Eparam)
  

  cat(sprintf('Fitting Bernoulli mixture:\n'))
  
  em_bernoulli[[i]] <- lapply(K, function(ii)
    expectation_maximization(data, ii, 'bernoulli_llh', eps, totem) )

  bic_bernoulli[[i]] <- sapply(1:length(K), 
                        function(ii) 
                          bic(max(em_bernoulli[[i]][[ii]]$llh), 
                              K[ii]*dim(data)[2]+K[ii], 
                              dim(data)[1]))

  aic_bernoulli[[i]] <- sapply(1:length(K), 
                        function(ii) 
                          aic(max(em_bernoulli[[i]][[ii]]$llh), 
                              K[ii]*dim(data)[2]+K[ii], 
                              dim(data)[1]))
                        
  cat(sprintf('Fitting logistic mixture:\n'))

  em_logit[[i]] <- lapply(K, function(ii)
    expectation_maximization(data, ii, 'logit_llh', eps, totem) )

  bic_logit[[i]] <- sapply(1:length(K), 
                          function(ii) 
                            bic(max(em_logit[[i]][[ii]]$llh), 
                                K[ii]*dim(data)[2]+K[ii], 
                                dim(data)[1]))
  
  aic_logit[[i]] <- sapply(1:length(K), 
                          function(ii) 
                            aic(max(em_logit[[i]][[ii]]$llh), 
                                K[ii]*dim(data)[2]+K[ii], 
                                dim(data)[1]))
  
  
  
  i = i + 1  
}


for (i in 1:4){
  
  tmpname <- substr(files[i], nchar(files[i])-14, nchar(files[i])-11)
  
  #BIC
  
  idx <- which.min( bic_logit[[i]] )
  
  pred_label_bernoulli <- em_bernoulli[[i]][[idx]]$pred_label
  pred_label_logit     <- em_logit[[i]][[idx]]$pred_label
  
  labels <- data.frame( ber_label = pred_label_bernoulli, 
                        logit_label = pred_label_logit,
                        true_label = true_label )
  
  write.csv( labels, file = paste('simple_mixtures_labels_syn_methy_bic_', tmpname, '.csv', 
                                  sep=''), row.names=F )
  #AIC
  
  idx <- which.min( aic_logit[[i]] )
  
  pred_label_bernoulli <- em_bernoulli[[i]][[idx]]$pred_label
  pred_label_logit     <- em_logit[[i]][[idx]]$pred_label
  
  labels <- data.frame( ber_label = pred_label_bernoulli, 
                        logit_label = pred_label_logit,
                        true_label = true_label )
  
  write.csv( labels, file = paste('simple_mixtures_labels_syn_methy_aic_', tmpname, '.csv', 
                                  sep=''), row.names=F )
}





