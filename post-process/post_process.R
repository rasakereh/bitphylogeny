library(mcclust)

get_best_cluster_label <- function(label_traces){
  ltmat <- as.matrix(label_traces)
  class(ltmat) <- 'integer'
  ltmat <- ltmat + 1
  psm <- comp.psm(ltmat)
  mpear <- maxpear(psm)
  best_label <- mpear$cl
} 

get_path <- function(codename_folder){
  filepath <- paste(codename_folder, dir(codename_folder), sep='/')
  return(filepath)
}
  
load_label_traces <- function(filepath){
  filename <- dir(filepath, pattern = 'label_traces')
  label_traces <- read.csv(paste(filepath,filename,sep='/'), 
                                header = F)
  return(label_traces)
}

write_best_label <- function(filepath){
  label_traces <- load_label_traces(filepath)
  best_label <- get_best_cluster_label(label_traces)
  outfile = 'best_label.csv'
  write.csv(best_label, 
            file=paste(filepath,outfile,sep='/'), 
            row.names = F)
}

find_synth_folders <- function(pathname){
  codename_folders <- dir(pathname)
  filepath <- vector(mode='character',0)
  for (cd in codename_folders){
    samplefolder <- dir(paste(pathname,cd,sep='/'))
    if ( length(grep('noisy', samplefolder) )> 0) {
      filepath <- cbind(filepath, 
                        paste(pathname, cd, samplefolder, sep='/'))
    }
  }
  return(filepath)
}

load_vmeasure_traces <- function( filepaths ){
  fn <- 'v_measure_traces.csv'
  x <- lapply(filepaths, function(i) read.csv( paste(i, fn, sep='/'), header=T, as.is = T ))
  xn <- sapply(filepaths, function(ii) substr(ii, nchar(ii)-14, nchar(ii)-11 ) )
  names(x) <- xn
  return(x)
}

plot_vmeasure <- function(x){
  n <- length(x)
  xn <- names(x)
  vmeasure <- as.numeric(x[[1]]$v_measure)
  for(i in 2:n){vmeasure <- cbind(vmeasure, as.numeric(x[[i]]$v_measure ))}
  colnames(vmeasure) <- xn

  homogeneity <- as.numeric(x[[1]]$homogeneity)
  for(i in 2:n){homogeneity <- cbind(homogeneity, as.numeric(x[[i]]$homogeneity))}
  colnames(homogeneity) <- xn
  
  completeness <- as.numeric(x[[1]]$completeness)
  for(i in 2:n){completeness <- cbind(completeness, as.numeric(x[[i]]$completeness))}
  colnames(completeness) <- xn
  
  pdf(file="./v_measures.pdf")
  par(mfrow=c(2,2)) 
  boxplot(homogeneity[, order(colnames(homogeneity))], xlab= 'Error Rate', ylab='Homogeneity' )
  boxplot(completeness[, order(colnames(completeness))], xlab= 'Error Rate', ylab='Completeness' )
  boxplot(vmeasure[, order(colnames(vmeasure))], xlab= 'Error Rate', ylab='V-measure' )
  dev.off()
}


##Draw v-measure plots
pathname = '/home/yuan03/Downloads/traces-eth/Sottoriva2/pairwise/'
synth_files = find_synth_folders(pathname)
x = load_vmeasure_traces(synth_files)
plot_vmeasure(x)

## Compute best_labels
compute_best_labesl = F
if (compute_best_labesl){
  codename_folders <- dir(pathname)
  for (f1 in codename_folders) {
    write_best_label(get_path(paste(pathname,f1, sep='/')))
  }  
}



