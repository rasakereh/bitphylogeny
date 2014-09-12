CX_R1<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/demirelief-3049/CX_IRX2P_R1.csv/traces-CX_IRX2P_R1.csv",sep=",",header=T)
CX_R2<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/spatulose-3680/CX_IRX2P_R2.csv/traces-CX_IRX2P_R2.csv",sep=",",header=T)
CX_R6<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/alevin-293/CX_IRX2P_R6.csv/traces-CX_IRX2P_R6.csv",sep=",",header=T)

CX_L3<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/uplinked-8718/CX_IRX2P_L3.csv/traces-CX_IRX2P_L3.csv",sep=",",header=T)
CX_L4<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/subtransverse-8545/CX_IRX2P_L4.csv/traces-CX_IRX2P_L4.csv",sep=",",header=T)
CX_L5<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/marbleness-8646/CX_IRX2P_L5.csv/traces-CX_IRX2P_L5.csv",sep=",",header=T)

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CX_nodenum.pdf")

#par(mfrow=c(1,2))
feature<-2
band<-0.7
g<-density(CX_R1[,feature],bw=band)

plot(g, col=rgb(0,0.7,0.7,1/2), lwd=2, main="CX node numbers", xlim=c(min(CX_R1[,feature],CX_R2[,feature],CX_L3[,feature],CX_L4[,feature],CX_L5[,feature]),25))
lines(density(CX_R2[,feature],bw=band), col=rgb(0,0.8,0.8,1/2),lwd=2)
lines(density(CX_R6[,feature],bw=band), col=rgb(0,0.9,0.9,1/2),lwd=2)

lines(density(CX_L3[,feature],bw=band), col=rgb(0.7,0,0.7,1/2),lwd=2)
lines(density(CX_L4[,feature],bw=band), col=rgb(0.8,0,0.8,1/2),lwd=2)
lines(density(CX_L5[,feature],bw=band), col=rgb(0.9,0,0.9,1/2),lwd=2)

dev.off()

t.test(c(mean(CX_R1[,feature]),mean(CX_R2[,feature]),mean(CX_R6[,feature])),c(mean(CX_L3[,feature]),mean(CX_L4[,feature]),mean(CX_L5[,feature])))


pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CX_bignodenum.pdf")

#pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CX_bignodenum_depth.pdf",width=10,height=5)


band<-0.7
feature<-3
g<-density(CX_R1[,feature],bw=band)

plot(g, col=rgb(0,0.7,0.7,1/2), main="CX big node numbers", xlim=c(min(CX_R1[,feature],CX_R2[,feature],CX_L3[,feature],CX_L4[,feature],CX_L5[,feature]),max(CX_R1[,feature],CX_R2[,feature],CX_L3[,feature],CX_L4[,feature],CX_L5[,feature])),
		ylim=c(0,0.5),lwd=2)
lines(density(CX_R2[,feature],bw=band), col=rgb(0,0.8,0.8,1/2),lwd=2)
lines(density(CX_R6[,feature],bw=band), col=rgb(0,0.9,0.9,1/2),lwd=2)

lines(density(CX_L3[,feature],bw=band), col=rgb(0.7,0,0.7,1/2),lwd=2)
lines(density(CX_L4[,feature],bw=band), col=rgb(0.8,0,0.8,1/2),lwd=2)
lines(density(CX_L5[,feature],bw=band), col=rgb(0.9,0,0.9,1/2),lwd=2)

dev.off()


pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CX_depth.pdf")

feature<-5
band<-0.2
g<-density(CX_R1[,feature],bw=band)

plot(g, col=rgb(0,0.7,0.7,1/2), main="CX max depth", xlim=c(min(CX_R1[,feature],CX_R2[,feature],CX_L3[,feature],CX_L4[,feature],CX_L5[,feature]),max(CX_R1[,feature],CX_R2[,feature],CX_L3[,feature],CX_L4[,feature],CX_L5[,feature])),
		ylim=c(0,2),lwd=2, xaxt='n')
axis(1, at=c(1,2,3,4,5), labels=c(1,2,3,4,5))
lines(density(CX_R2[,feature],bw=band), col=rgb(0,0.8,0.8,1/2),lwd=2)
lines(density(CX_R6[,feature],bw=band), col=rgb(0,0.9,0.9,1/2),lwd=2)

lines(density(CX_L3[,feature],bw=band), col=rgb(0.7,0,0.7,1/2),lwd=2)
lines(density(CX_L4[,feature],bw=band), col=rgb(0.8,0,0.8,1/2),lwd=2)
lines(density(CX_L5[,feature],bw=band), col=rgb(0.9,0,0.9,1/2),lwd=2)

t.test(c(mean(CX_R1[,feature]),mean(CX_R2[,feature]),mean(CX_R6[,feature])),c(mean(CX_L3[,feature]),mean(CX_L4[,feature]),mean(CX_L5[,feature])))

dev.off()

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CX_base_std_rootbias.pdf",height=5,width=15)

par(mfrow=c(1,2))
feature<-6
band<-0.2
g<-density(CX_R1[,feature],bw=band)

plot(g, col=rgb(0,0.7,0.7,1/2), main="CX base value", xlim=c(min(CX_R1[,feature],CX_R2[,feature],CX_L3[,feature],CX_L4[,feature],CX_L5[,feature]),max(CX_R1[,feature],CX_R2[,feature],CX_L3[,feature],CX_L4[,feature],CX_L5[,feature])),
		ylim=c(0,1.6),lwd=2)
lines(density(CX_R2[,feature],bw=band), col=rgb(0,0.8,0.8,1/2),lwd=2)
lines(density(CX_R6[,feature],bw=band), col=rgb(0,0.9,0.9,1/2),lwd=2)

lines(density(CX_L3[,feature],bw=band), col=rgb(0.7,0,0.7,1/2),lwd=2)
lines(density(CX_L4[,feature],bw=band), col=rgb(0.8,0,0.8,1/2),lwd=2)
lines(density(CX_L5[,feature],bw=band), col=rgb(0.9,0,0.9,1/2),lwd=2)

feature<-7
band<-0.2
g<-density(CX_R1[,feature],bw=band)

plot(g, col=rgb(0.7,0,0.7,1/2), main="CX standard deviation", xlim=c(min(CX_R1[,feature],CX_R2[,feature],CX_L3[,feature],CX_L4[,feature],CX_L5[,feature]),max(CX_R1[,feature],CX_R2[,feature],CX_L3[,feature],CX_L4[,feature],CX_L5[,feature])),
		ylim=c(0,1.05),lwd=2)
lines(density(CX_R2[,feature],bw=band), col=rgb(0.8,0,0.8,1/2),lwd=2)
lines(density(CX_R6[,feature],bw=band), col=rgb(0.9,0,0.9,1/2),lwd=2)

lines(density(CX_L3[,feature],bw=band), col=rgb(0.7,0,0.7,1/2),lwd=2)
lines(density(CX_L4[,feature],bw=band), col=rgb(0.8,0,0.8,1/2),lwd=2)
lines(density(CX_L5[,feature],bw=band), col=rgb(0.9,0,0.9,1/2),lwd=2)


dev.off()

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CX_mass_barplot.pdf",width=9,height=5)

feature<-23

massesR<-matrix(c(mean(CX_R1[,feature]),mean(CX_R1[,feature+1]),mean(CX_R1[,feature+2]),mean(CX_R1[,feature+3]),mean(CX_R1[,feature+4]),mean(CX_R2[,feature]),mean(CX_R2[,feature+1]),mean(CX_R2[,feature+2]),mean(CX_R2[,feature+3]),mean(CX_R2[,feature+4]),
		mean(CX_R6[,feature]),mean(CX_R6[,feature+1]),mean(CX_R6[,feature+2]),mean(CX_R6[,feature+3]),mean(CX_R6[,feature+4])),nrow=5)
massesL<-matrix(c(mean(CX_L3[,feature]),mean(CX_L3[,feature+1]),mean(CX_L3[,feature+2]),mean(CX_L3[,feature+3]),mean(CX_L3[,feature+4]),mean(CX_L4[,feature]),mean(CX_L4[,feature+1]),mean(CX_L4[,feature+2]),mean(CX_L4[,feature+3]),mean(CX_L4[,feature+4]),
				mean(CX_L5[,feature]),mean(CX_L5[,feature+1]),mean(CX_L5[,feature+2]),mean(CX_L5[,feature+3]),mean(CX_L5[,feature+4])),nrow=5)
masses<-cbind(massesL,massesR)
masses<-apply(masses,2,rev)
par(mar=c(5.1,4.1,4.1,6.1))
barplot(masses,col=rev(c("red","blue","green","pink","orange")),names.arg=c("CX_L3", "CX_L4", "CX_L5", "CX_R1", "CX_R2", "CX_R6"),main="mean posterior layer-wise tumor masses",
		args.legend = list(x=11.5, y=1), bty = "n")


dev.off()


pdf("CX_total_branch_length.pdf")

band<-0.3

CX_R1branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/demirelief-3049/CX_IRX2P_R1.csv/traces-CX_IRX2P_R1.csv_branch_traces.csv",sep=",",header=F)
tb1<-0
tb1list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb1list[i]<-sum(unique(unlist(CX_R1branch[i,])))
}
tb1<-sum(tb1list)
plot(density(tb1list,bw=band),ylim=c(0,0.5), xlim=c(0,12), col=rgb(0,0.7,0.7,1/2),lwd=2)

CX_R2branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/spatulose-3680/CX_IRX2P_R2.csv/traces-CX_IRX2P_R2.csv_branch_traces.csv",sep=",",header=F)
tb2<-0
tb2list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb2list[i]<-sum(unique(unlist(CX_R2branch[i,])))
}
tb2<-sum(tb2list)
lines(density(tb2list,bw=band), col=rgb(0,0.8,0.8,1/2),lwd=2)

CX_R6branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/alevin-293/CX_IRX2P_R6.csv/traces-CX_IRX2P_R6.csv_branch_traces.csv",sep=",",header=F)
tb6<-0
tb6list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb6list[i]<-sum(unique(unlist(CX_R6branch[i,])))
}
tb6<-sum(tb6list)
lines(density(tb6list,bw=band), col=rgb(0,0.9,0.9,1/2),lwd=2)


CX_L3branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/uplinked-8718/CX_IRX2P_L3.csv/traces-CX_IRX2P_L3.csv_branch_traces.csv",sep=",",header=F)
tb3<-0
tb3list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb3list[i]<-sum(unique(unlist(CX_L3branch[i,])))
}
tb3<-sum(tb3list)
lines(density(tb3list,bw=band), col=rgb(0.7,0,0.7,1/2),lwd=2)

CX_L4branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/subtransverse-8545/CX_IRX2P_L4.csv/traces-CX_IRX2P_L4.csv_branch_traces.csv",sep=",",header=F)
tb4<-0
tb4list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb4list[i]<-sum(unique(unlist(CX_L4branch[i,])))
}
tb4<-sum(tb4list)
lines(density(tb4list,bw=band), col=rgb(0.8,0,0.8,1/2),lwd=2)

CX_L5branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/marbleness-8646/CX_IRX2P_L5.csv/traces-CX_IRX2P_L5.csv_branch_traces.csv",sep=",",header=F)
tb5<-0
tb5list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb5list[i]<-sum(unique(unlist(CX_L5branch[i,])))
}
tb5<-sum(tb5list)
lines(density(tb5list,bw=band), col=rgb(0.9,0,0.9,1/2),lwd=2)

t.test(c(mean(tb1list),mean(tb2list),mean(tb6list)),c(mean(tb3list),mean(tb4list),mean(tb5list)))


dev.off()

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CX_nodenum_depth.pdf")


feature1<-2
feature2<-5
plot(mean(CX_R1[,feature1]),mean(CX_R1[,feature2]),xlim=c(4,12),ylim=c(2.5,3.3),col=rgb(0,0.7,0.7,1/2),lwd=4, cex=3, ylab=colnames(CX_R1)[feature2],xlab=colnames(CX_R1)[feature1])
points(mean(CX_R2[,feature1]),mean(CX_R2[,feature2]),col=rgb(0,0.7,0.7,1/2),lwd=4,cex=3)
points(mean(CX_R6[,feature1]),mean(CX_R6[,feature2]),col=rgb(0,0.8,0.8,1/2),lwd=4,cex=3)

points(mean(CX_L3[,feature1]),mean(CX_L3[,feature2]),col=rgb(0.7,0,0.7,1/2),lwd=4,pch=2,cex=3)
points(mean(CX_L4[,feature1]),mean(CX_L4[,feature2]),col=rgb(0.8,0,0.8,1/2),lwd=4,pch=2,cex=3)
points(mean(CX_L5[,feature1]),mean(CX_L5[,feature2]),col=rgb(0.9,0,0.9,1/2),lwd=4,pch=2,cex=3)


dev.off()