rm(list=ls())
load_vmeasures<- function( fp, key ){
  fn <- dir( fp, pattern = key)
  x <- read.csv( paste(fp, fn, sep='/'), header=T, as.is = T)
  return(x)
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

args <- commandArgs(trailingOnly = TRUE)
fp <- paste(args, dir(args[1], pattern = 'mcmc-traces'), sep='/')
vmeasure_traces <- load_vmeasures(fp, 'vmeasure_traces.csv')
mpear_vmeasure <- load_vmeasures(fp, 'best_label_vmeasure.csv')

x<- as.matrix(vmeasure_traces)
class(x) <- 'numeric'
y<- as.matrix(mpear_vmeasure)
class(y) <- 'numeric'


pdf(file=paste (args,'vmeasures.pdf', sep='/'))  
par(cex.lab=1.5, cex.axis=1.5)
boxplot(x, outline=F, ylim=c(0.5,1), 
        cex.main=1.3, border=c('gray60'), col='gray90')
points(c(1,2,3),y, pch=22,cex = 1.5, bg= 'black')

colors1 <- c("gray90",'black')
colors2 <- c("gray",'black')
add_legend("bottomleft", legend=c("traces", 'BitPhylogeny'),
           pch=c(22,22), inset = c(0.1,0.13), col=colors1, 
           pt.bg=colors2,
           horiz=F, bty='n', cex=1.5)
dev.off()