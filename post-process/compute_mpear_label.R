library(mcclust)

get_best_cluster_label <- function(label_traces){
  ltmat <- as.matrix(label_traces)
  class(ltmat) <- 'integer'
  ltmat <- ltmat + 1
  psm <- comp.psm(ltmat)
  mpear <- maxpear(psm)
  best_label <- mpear$cl
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

get_path <- function(folder){
  trace <- dir(folder, pattern = 'mcmc-traces')
  filepath <- paste(folder, trace, sep='/')
  return(filepath)
}

##Draw v-measure plots
args <- commandArgs(trailingOnly = TRUE)
## Compute best_labels
write_best_label(get_path(args[1]))





