library("bitphylogenyR")


compute_vmeasures('/mnt/home/mcmcrun_new_syn/Sottoriva/final/noyaded-7082/noisy_full_methy_mono_8_2000_0_mutmat.csv','/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/mono-clone/noisy_full_methy_mono_8_2000_0_mutmat.csv')

compute_vmeasures('/mnt/home/mcmcrun_new_syn/Sottoriva/final/overfastidiousness-8958/noisy_full_methy_mono_8_2000_0.01_mutmat.csv/','/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/mono-clone/noisy_full_methy_mono_8_2000_0.01_mutmat.csv')

compute_vmeasures('/mnt/home/mcmcrun_new_syn/Sottoriva/final/wuzzled-3691/noisy_full_methy_mono_8_2000_0.02_mutmat.csv/','/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/mono-clone/noisy_full_methy_mono_8_2000_0.02_mutmat.csv')

compute_vmeasures('/mnt/home/mcmcrun_new_syn/Sottoriva/final/forevers-1577/noisy_full_methy_mono_8_2000_0.05_mutmat.csv/','/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/mono-clone/noisy_full_methy_mono_8_2000_0.05_mutmat.csv')

compute_vmeasures('/mnt/home/mcmcrun_new_syn/Sottoriva/final/zoogeog-3538/noisy_hyper_clones_hyper_8_2000_0_mutmat.csv/','/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/hyper-clone/noisy_small_clones_hyper_8_2000_0_mutmat.csv')

compute_vmeasures('/mnt/home/mcmcrun_new_syn/Sottoriva/final/trilite-4965/noisy_hyper_clones_hyper_8_2000_0.01_mutmat.csv/','/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/hyper-clone/noisy_small_clones_hyper_8_2000_0.01_mutmat.csv')

compute_vmeasures('/mnt/home/mcmcrun_new_syn/Sottoriva/final/testibrachial-9012/noisy_hyper_clones_hyper_8_2000_0.02_mutmat.csv/','/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/hyper-clone/noisy_small_clones_hyper_8_2000_0.02_mutmat.csv')

compute_vmeasures('/mnt/home/mcmcrun_new_syn/Sottoriva/final/tymp-8396/noisy_hyper_clones_hyper_8_2000_0.05_mutmat.csv/','/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/hyper-clone/noisy_small_clones_hyper_8_2000_0.05_mutmat.csv')

####mono clone:

vmeasure_traces_mat = matrix(0,nrow=30000,ncol=4)
mpearvec = vector()

fp <- get_path('/mnt/home/mcmcrun_new_syn/Sottoriva/final/noyaded-7082/noisy_full_methy_mono_8_2000_0_mutmat.csv/', 'mcmc-traces')

vmeasure_traces <- as.matrix(load_vmeasures(fp, 'vmeasure_traces.csv'))
mpear_vmeasure <- as.matrix(load_vmeasures(fp, 'mpear_label_vmeasure.csv'))
class(vmeasure_traces) <- 'numeric'
class(mpear_vmeasure) <- 'numeric'

vmeasure_traces_mat[,1] <- vmeasure_traces
mpearvec[1] <- mpear_vmeasure[3]

fp <- get_path('/mnt/home/mcmcrun_new_syn/Sottoriva/final/overfastidiousness-8958/noisy_full_methy_mono_8_2000_0.01_mutmat.csv/', 'mcmc-traces')

vmeasure_traces <- as.matrix(load_vmeasures(fp, 'vmeasure_traces.csv'))
mpear_vmeasure <- as.matrix(load_vmeasures(fp, 'mpear_label_vmeasure.csv'))
class(vmeasure_traces) <- 'numeric'
class(mpear_vmeasure) <- 'numeric'

vmeasure_traces_mat[,2] <- vmeasure_traces
mpearvec[2] <- mpear_vmeasure[3]

fp <- get_path('/mnt/home/mcmcrun_new_syn/Sottoriva/final/wuzzled-3691/noisy_full_methy_mono_8_2000_0.02_mutmat.csv/', 'mcmc-traces')

vmeasure_traces <- as.matrix(load_vmeasures(fp, 'vmeasure_traces.csv'))
mpear_vmeasure <- as.matrix(load_vmeasures(fp, 'mpear_label_vmeasure.csv'))
class(vmeasure_traces) <- 'numeric'
class(mpear_vmeasure) <- 'numeric'

vmeasure_traces_mat[,3] <- vmeasure_traces
mpearvec[3] <- mpear_vmeasure[3]

fp <- get_path('/mnt/home/mcmcrun_new_syn/Sottoriva/final/forevers-1577/noisy_full_methy_mono_8_2000_0.05_mutmat.csv/', 'mcmc-traces')

vmeasure_traces <- as.matrix(load_vmeasures(fp, 'vmeasure_traces.csv'))
mpear_vmeasure <- as.matrix(load_vmeasures(fp, 'mpear_label_vmeasure.csv'))
class(vmeasure_traces) <- 'numeric'
class(mpear_vmeasure) <- 'numeric'

vmeasure_traces_mat[,4] <- vmeasure_traces
mpearvec[4] <- mpear_vmeasure[3]

mcmc_vmeasures <- list()
colnames(vmeasure_traces_mat) <- c("0", "0.01", "0.02", "0.05")
mcmc_vmeasures$mono_clone <- vmeasure_traces_mat

mpear_vmeasures <- list()
names(mpearvec) <- c("0", "0.01", "0.02", "0.05")
mpear_vmeasures$mono_clone <- mpearvec

####hyper clone:
vmeasure_traces_mat = matrix(0,nrow=30000,ncol=4)
mpearvec = vector()

fp <- get_path('/mnt/home/mcmcrun_new_syn/Sottoriva/final/zoogeog-3538/noisy_small_clones_hyper_8_2000_0_mutmat.csv/', 'mcmc-traces')

vmeasure_traces <- as.matrix(load_vmeasures(fp, 'vmeasure_traces.csv'))
mpear_vmeasure <- as.matrix(load_vmeasures(fp, 'mpear_label_vmeasure.csv'))
class(vmeasure_traces) <- 'numeric'
class(mpear_vmeasure) <- 'numeric'

vmeasure_traces_mat[,1] <- vmeasure_traces
mpearvec[1] <- mpear_vmeasure[3]

fp <- get_path('/mnt/home/mcmcrun_new_syn/Sottoriva/final/trilite-4965/noisy_small_clones_hyper_8_2000_0.01_mutmat.csv/', 'mcmc-traces')

vmeasure_traces <- as.matrix(load_vmeasures(fp, 'vmeasure_traces.csv'))
mpear_vmeasure <- as.matrix(load_vmeasures(fp, 'mpear_label_vmeasure.csv'))
class(vmeasure_traces) <- 'numeric'
class(mpear_vmeasure) <- 'numeric'

vmeasure_traces_mat[,2] <- vmeasure_traces
mpearvec[2] <- mpear_vmeasure[3]

fp <- get_path('/mnt/home/mcmcrun_new_syn/Sottoriva/final/testibrachial-9012/noisy_small_clones_hyper_8_2000_0.02_mutmat.csv/', 'mcmc-traces')

vmeasure_traces <- as.matrix(load_vmeasures(fp, 'vmeasure_traces.csv'))
mpear_vmeasure <- as.matrix(load_vmeasures(fp, 'mpear_label_vmeasure.csv'))
class(vmeasure_traces) <- 'numeric'
class(mpear_vmeasure) <- 'numeric'

vmeasure_traces_mat[,3] <- vmeasure_traces
mpearvec[3] <- mpear_vmeasure[3]

fp <- get_path('/mnt/home/mcmcrun_new_syn/Sottoriva/final/tymp-8396/noisy_small_clones_hyper_8_2000_0.05_mutmat.csv/', 'mcmc-traces')

vmeasure_traces <- as.matrix(load_vmeasures(fp, 'vmeasure_traces.csv'))
mpear_vmeasure <- as.matrix(load_vmeasures(fp, 'mpear_label_vmeasure.csv'))
class(vmeasure_traces) <- 'numeric'
class(mpear_vmeasure) <- 'numeric'

vmeasure_traces_mat[,4] <- vmeasure_traces
mpearvec[4] <- mpear_vmeasure[3]

colnames(vmeasure_traces_mat) <- c("0", "0.01", "0.02", "0.05")
mcmc_vmeasures$hyper_clone <- vmeasure_traces_mat

names(mpearvec) <- c("0", "0.01", "0.02", "0.05")
mpear_vmeasures$hyper_clone <- mpearvec


###baseline methods:

###mono clone:

library(cluster)
hc_vmeasures <- list()
hcvmlist <- vector()

kc_vmeasures <- list()
kcvmlist <- vector()

mono0<-read.csv("/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/mono-clone/noisy_full_methy_mono_8_2000_0_mutmat.csv")
true_label_mono0 <- mono0[,dim(mono0)[2]]
mono0 <- mono0[,-dim(mono0)[2]]
K <- seq(2,20,1)
hc <- get_label_hc(mono0,K)
hcvmlist[1] <- (vmeasureR(hc$hc_label, true_label_mono0))[3]

kc <- get_label_kc(mono0,K)
kcvmlist[1] <- (vmeasureR(kc$kc_label, true_label_mono0))[3]

mst<-get_mst(hc$hc_genotype)
#plot_mst(hc$hc_genotype,hc$hc_label,flag=F)
mst<-get_mst(kc$kc_genotype)
#plot_mst(kc$kc_genotype,kc$kc_label,flag=F)


mono001<-read.csv("/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/mono-clone/noisy_full_methy_mono_8_2000_0.01_mutmat.csv")
true_label_mono001 <- mono001[,dim(mono001)[2]]
mono001 <- mono001[,-dim(mono001)[2]]
K <- seq(2,20,1)
hc <- get_label_hc(mono001,K)
hcvmlist[2] <- (vmeasureR(hc$hc_label, true_label_mono001))[3]

kc <- get_label_kc(mono001,K)
kcvmlist[2] <- (vmeasureR(kc$kc_label, true_label_mono001))[3]

mst<-get_mst(hc$hc_genotype)
#plot_mst(hc$hc_genotype,hc$hc_label,flag=F)
mst<-get_mst(kc$kc_genotype)
#plot_mst(kc$kc_genotype,kc$kc_label,flag=F)


mono002<-read.csv("/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/mono-clone/noisy_full_methy_mono_8_2000_0.02_mutmat.csv")
true_label_mono002 <- mono002[,dim(mono002)[2]]
mono002 <- mono002[,-dim(mono002)[2]]
K <- seq(2,20,1)
hc <- get_label_hc(mono002,K)
hcvmlist[3] <- (vmeasureR(hc$hc_label, true_label_mono002))[3]

kc <- get_label_kc(mono002,K)
kcvmlist[3] <- (vmeasureR(kc$kc_label, true_label_mono002))[3]

mst<-get_mst(hc$hc_genotype)
#plot_mst(hc$hc_genotype,hc$hc_label,flag=F)
mst<-get_mst(kc$kc_genotype)
#plot_mst(kc$kc_genotype,kc$kc_label,flag=F)

mono005<-read.csv("/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/mono-clone/noisy_full_methy_mono_8_2000_0.05_mutmat.csv")
true_label_mono005 <- mono005[,dim(mono005)[2]]
mono005 <- mono005[,-dim(mono005)[2]]
K <- seq(2,20,1)
hc <- get_label_hc(mono005,K)
hcvmlist[4] <- (vmeasureR(hc$hc_label, true_label_mono005))[3]

kc <- get_label_kc(mono005,K)
kcvmlist[4] <- (vmeasureR(kc$kc_label, true_label_mono005))[3]

mst<-get_mst(hc$hc_genotype)
#plot_mst(hc$hc_genotype,hc$hc_label,flag=F)
mst<-get_mst(kc$kc_genotype)
#plot_mst(kc$kc_genotype,kc$kc_label,flag=F)

names(hcvmlist) <- c("0", "0.01", "0.02", "0.05")
hc_vmeasures$mono_clone <- hcvmlist

names(kcvmlist) <- c("0", "0.01", "0.02", "0.05")
kc_vmeasures$mono_clone <- kcvmlist


###hyper clone: (use higher max number of clusters)
hcvmlist <- vector()

kcvmlist <- vector()

hyper0<-read.csv("/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/hyper-clone/noisy_small_clones_hyper_8_2000_0_mutmat.csv")
true_label_hyper0 <- hyper0[,dim(hyper0)[2]]
hyper0 <- hyper0[,-dim(hyper0)[2]]
K <- seq(2,30,1)
hc <- get_label_hc(hyper0,K)
hcvmlist[1] <- (vmeasureR(hc$hc_label, true_label_hyper0))[3]

kc <- get_label_kc(hyper0,K)
kcvmlist[1] <- (vmeasureR(kc$kc_label, true_label_hyper0))[3]

mst<-get_mst(hc$hc_genotype)
#plot_mst(hc$hc_genotype,hc$hc_label,flag=F)
mst<-get_mst(kc$kc_genotype)
#plot_mst(kc$kc_genotype,kc$kc_label,flag=F)


hyper001<-read.csv("/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/hyper-clone/noisy_small_clones_hyper_8_2000_0.01_mutmat.csv")
true_label_hyper001 <- hyper001[,dim(hyper001)[2]]
hyper001 <- hyper001[,-dim(hyper001)[2]]
K <- seq(2,30,1)
hc <- get_label_hc(hyper001,K)
hcvmlist[2] <- (vmeasureR(hc$hc_label, true_label_hyper001))[3]

kc <- get_label_kc(hyper001,K)
kcvmlist[2] <- (vmeasureR(kc$kc_label, true_label_hyper001))[3]

mst<-get_mst(hc$hc_genotype)
#plot_mst(hc$hc_genotype,hc$hc_label,flag=F)
mst<-get_mst(kc$kc_genotype)
#plot_mst(kc$kc_genotype,kc$kc_label,flag=F)


hyper002<-read.csv("/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/hyper-clone/noisy_small_clones_hyper_8_2000_0.02_mutmat.csv")
true_label_hyper002 <- hyper002[,dim(hyper002)[2]]
hyper002 <- hyper002[,-dim(hyper002)[2]]
K <- seq(2,30,1)
hc <- get_label_hc(hyper002,K)
hcvmlist[3] <- (vmeasureR(hc$hc_label, true_label_hyper002))[3]

kc <- get_label_kc(hyper002,K)
kcvmlist[3] <- (vmeasureR(kc$kc_label, true_label_hyper002))[3]

mst<-get_mst(hc$hc_genotype)
#plot_mst(hc$hc_genotype,hc$hc_label,flag=F)
mst<-get_mst(kc$kc_genotype)
#plot_mst(kc$kc_genotype,kc$kc_label,flag=F)


hyper005<-read.csv("/home/warriour/git/phylo-tree/R/bitphylogenyR/inst/extdata/synthetic/hyper-clone/noisy_small_clones_hyper_8_2000_0.05_mutmat.csv")
true_label_hyper005 <- hyper005[,dim(hyper005)[2]]
hyper005 <- hyper005[,-dim(hyper005)[2]]
K <- seq(2,30,1)
hc <- get_label_hc(hyper005,K)
hcvmlist[4] <- (vmeasureR(hc$hc_label, true_label_hyper005))[3]

kc <- get_label_kc(hyper005,K)
kcvmlist[4] <- (vmeasureR(kc$kc_label, true_label_hyper005))[3]

mst<-get_mst(hc$hc_genotype)
#plot_mst(hc$hc_genotype,hc$hc_label,flag=F)
mst<-get_mst(kc$kc_genotype)
#plot_mst(kc$kc_genotype,kc$kc_label,flag=F)


names(hcvmlist) <- c("0", "0.01", "0.02", "0.05")
hc_vmeasures$hyper_clone <- hcvmlist

names(kcvmlist) <- c("0", "0.01", "0.02", "0.05")
kc_vmeasures$hyper_clone <- kcvmlist


####plot:

mon<-gdl2df(file="/home/warriour/git/phylo-tree/old/tree-truth-full-mono-methy2.gdl")
hyp<-gdl2df(file="/home/warriour/git/phylo-tree/old/tree-truth-full-hyper-methy.gdl")

pdf(file="Figure-sim2a.pdf",width=4,height=4)
plot_sankey(mon)
dev.off()
pdf(file="Figure-sim2b.pdf",width=10,height=10)
plot_sankey(hyp)
dev.off()


pdf(file="Figure-sim2c.pdf",width=7,height=4)

par(mfrow=c(1,2), oma = c(3,3,0,0) + 0.1,mar = c(0,0,1,0.5) + 0.1, cex.lab=1.5, cex.axis=1.5)
boxplot(mcmc_vmeasures$mono_clone, outline=F,ylim=c(0,1) ,cex.main=1.3,border=c('gray60'), col='gray90')
points( mpear_vmeasures$mono_clone, pch=22,cex = 1.5, bg= 'black')
points( hc_vmeasures$mono_clone, pch=24,cex = 1.5, bg= 'black')
points( kc_vmeasures$mono_clone, pch=25,cex = 1.5, bg= 'black')


boxplot(mcmc_vmeasures$hyper_clone, outline=F,ylim=c(0.5,1) ,cex.main=1.3,border=c('gray60'), col='gray90')
points( mpear_vmeasures$hyper_clone, pch=22,cex = 1.5, bg= 'black')
points( hc_vmeasures$hyper_clone, pch=24,cex = 1.5, bg= 'black')
points( kc_vmeasures$hyper_clone, pch=25,cex = 1.5, bg= 'black')

colors1 <- c("gray90",'black', 'black', "black")
colors2 <- c("gray",'black', 'black', "black")
add_legend("bottomleft", legend=c("traces", 'BitPhylogeny','k-centroids','hierarchical clustering'),pch=c(22,22,25,24), inset = c(0.1,0.13),col=colors1,pt.bg=colors2,horiz=F, bty='n', cex=1.5)
title(xlab = "error",ylab = "v-measure",outer = TRUE, line = 2.2)

dev.off()

pdf(file="Figure-sim2e.pdf",width=7,height=4)

mono_clone_t <- c(1, 2)
mono_clone_bit <- c(1, 2)
mono_clone_hc <- c(1, 2)
mono_clone_kc <- c(1, 2)
mono_clone_bit <- rbind(mono_clone_bit,c(2, 4), c(2,3), c(1,5) )
mono_clone_hc <- rbind(mono_clone_hc,c(10, 20), c(13,20), c(13,19))
mono_clone_kc <- rbind(mono_clone_kc,c(10.95, 20), c(11.05,20), c(9,20))
hyper_clone_t <- c(5, 18)
hyper_clone_bit <- c(2, 12)
hyper_clone_hc <- c(7, 18)
hyper_clone_kc <- c(7, 18)
hyper_clone_bit <- rbind(hyper_clone_bit,c(3, 27),c(3,21), c(3,26) )
hyper_clone_hc <- rbind(hyper_clone_hc,c(5, 22),c(12,30), c(9,30))
hyper_clone_kc <- rbind(hyper_clone_kc,c(9, 30),c(10, 30), c(12,30))
par(mfrow=c(1,2), oma = c(3,3,0,0) + 0.1,mar = c(0,0,1,0.5) + 0.1, cex.lab=1.5, cex.axis=1.5)
color <- c('blue', 'red', 'red', 'red')
plot(mono_clone_t[2], mono_clone_t[1],pch=3, ylim=c(0,15) ,xlim= c(2,20) ,cex=1.5)
points(mono_clone_bit[,2], mono_clone_bit[,1],pch=0,cex=1.5, col=color)
points(mono_clone_hc[,2], mono_clone_hc[,1],pch=2,cex=1.5, col=color)
points(mono_clone_kc[,2], mono_clone_kc[,1],pch=6,cex=1.5, col=color)
plot(hyper_clone_t[2], hyper_clone_t[1],pch=3, ylim=c(0,15) ,xlim= c(5,32),cex=1.5,yaxt='n')
points(hyper_clone_bit[,2], hyper_clone_bit[,1],pch=0,cex=1.5, col=color)
points(hyper_clone_hc[,2], hyper_clone_hc[,1],pch=2,cex=1.5, col=color)
points(hyper_clone_kc[,2], hyper_clone_kc[,1],pch=6,cex=1.5, col=color)
colors1 <- c("black",'black', 'black', "black")
colors2 <- c("black",'black', 'black', "black")
add_legend("topleft", legend=c("truth", 'BitPhylogeny', 'k-centroids', 'hierarchical clustering'),pch=c(3,0,6,2), inset = c(0.08,0.02), col=colors1,pt.bg=colors2,horiz=F, bty='n', cex=1.5)
add_legend("topleft", legend=c("noiseless", 'noise levels: \n0.01,0.02,0.05'),
inset = c(0.50,0.02), text.col=c('blue','red'),horiz=F, bty='n', cex=1.5)
title(xlab = "number of clones",ylab = "maximum tree depth",outer = TRUE, line = 2.2)

dev.off()

