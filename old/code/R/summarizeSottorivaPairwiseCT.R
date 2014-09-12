CT_R1<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/incorrespondent-6142/CT_IRX2P_R1.csv/traces-CT_IRX2P_R1.csv",sep=",",header=T)
CT_R4<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/scolia-5868/CT_IRX2P_R4.csv/traces-CT_IRX2P_R4.csv",sep=",",header=T)
CT_R5<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/neuromerous-9376/CT_IRX2P_R5.csv/traces-CT_IRX2P_R5.csv",sep=",",header=T)
CT_R6<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/proaccelerin-3347/CT_IRX2P_R6.csv/traces-CT_IRX2P_R6.csv",sep=",",header=T)

CT_L2<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/oceanity-4938/CT_IRX2P_L2.csv/traces-CT_IRX2P_L2.csv",sep=",",header=T)
CT_L3<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/nonprolific-2900/CT_IRX2P_L3.csv/traces-CT_IRX2P_L3.csv",sep=",",header=T)
CT_L7<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/buoy-3083/CT_IRX2P_L7.csv/traces-CT_IRX2P_L7.csv",sep=",",header=T)
CT_L8<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/epistyle-3189/CT_IRX2P_L8.csv/traces-CT_IRX2P_L8.csv",sep=",",header=T)

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CT_nodenum.pdf")

feature<-2
g<-density(CT_R1[,feature])

plot(g, col=rgb(0,0,0.25,1/4), main="CT node numbers", xlim=c(min(CT_R1[,feature],CT_R4[,feature],CT_R5[,feature],CT_R6[,feature],CT_L2[,feature],CT_L3[,feature],CT_L7[,feature],CT_L8[,feature]),max(CT_R1[,feature],CT_R4[,feature],CT_R5[,feature],CT_R6[,feature],CT_L2[,feature],CT_L3[,feature],CT_L7[,feature],CT_L8[,feature])))
lines(density(CT_R4[,feature]), col=rgb(0,0,0.5,1/4))
lines(density(CT_R5[,feature]), col=rgb(0,0,0.75,1/4))
lines(density(CT_R6[,feature]), col=rgb(0,0,1,1/4))

lines(density(CT_L2[,feature]), col=rgb(0.25,0,0,1/4))
lines(density(CT_L3[,feature]), col=rgb(0.5,0,0,1/4))
lines(density(CT_L7[,feature]), col=rgb(0.75,0,0,1/4))
lines(density(CT_L8[,feature]), col=rgb(1,0,0,1/4))

dev.off()

#pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CT_bignodenum.pdf")

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CT_bignodenum_depth.pdf",width=10,height=5)
par(mfrow=c(1,2))

feature<-3
g<-density(CT_R1[,feature],bw=1)

plot(g, col=rgb(0,0,0.7,1/2), main="CT big node numbers", xlim=c(min(CT_R1[,feature],CT_R4[,feature],CT_R5[,feature],CT_R6[,feature],CT_L2[,feature],CT_L3[,feature],CT_L7[,feature],CT_L8[,feature]),max(CT_R1[,feature],CT_R4[,feature],CT_R5[,feature],CT_R6[,feature],CT_L2[,feature],CT_L3[,feature],CT_L7[,feature],CT_L8[,feature])),
		ylim=c(0,0.5),lwd=2)
lines(density(CT_R4[,feature],bw=1), col=rgb(0,0,0.8,1/2),lwd=2)
lines(density(CT_R5[,feature],bw=1), col=rgb(0,0,0.9,1/2),lwd=2)
lines(density(CT_R6[,feature],bw=1), col=rgb(0,0,1,1/2),lwd=2)

lines(density(CT_L2[,feature],bw=1), col=rgb(0.7,0,0,1/2),lwd=2)
lines(density(CT_L3[,feature],bw=1), col=rgb(0.8,0,0,1/2),lwd=2)
lines(density(CT_L7[,feature],bw=1), col=rgb(0.9,0,0,1/2),lwd=2)
lines(density(CT_L8[,feature],bw=1), col=rgb(1,0,0,1/2),lwd=2)

#dev.off()

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CT_depth.pdf")

feature<-5
band<-0.2
g<-density(CT_R1[,feature],bw=band)

plot(g, col=rgb(0,0.7,0.7,1/2), main="CT max depth", xlim=c(min(CT_R1[,feature],CT_R4[,feature],CT_R5[,feature],CT_R6[,feature],CT_L2[,feature],CT_L3[,feature],CT_L7[,feature],CT_L8[,feature]),max(CT_R1[,feature],CT_R4[,feature],CT_R5[,feature],CT_R6[,feature],CT_L2[,feature],CT_L3[,feature],CT_L7[,feature],CT_L8[,feature])),
		ylim=c(0,2),lwd=2, xaxt='n')
axis(1, at=c(2,3,4,5), labels=c(2,3,4,5))
lines(density(CT_R4[,feature],bw=band), col=rgb(0,0.8,0.8,1/2),lwd=2)
lines(density(CT_R5[,feature],bw=band), col=rgb(0,0.9,0.9,1/2),lwd=2)
lines(density(CT_R6[,feature],bw=band), col=rgb(0,1,1,1/2),lwd=2)

lines(density(CT_L2[,feature],bw=band), col=rgb(0.7,0,0.7,1/2),lwd=2)
lines(density(CT_L3[,feature],bw=band), col=rgb(0.8,0,0.8,1/2),lwd=2)
lines(density(CT_L7[,feature],bw=band), col=rgb(0.9,0,0.9,1/2),lwd=2)
lines(density(CT_L8[,feature],bw=band), col=rgb(1,0,1,1/2),lwd=2)


wilcox.test(c(var(CT_R1[,feature]), var(CT_R4[,feature]), var(CT_R5[,feature]), var(CT_R6[,feature])),c(var(CT_L2[,feature]),var(CT_L3[,feature]),var(CT_L7[,feature]),var(CT_L8[,feature])),
	alternative="greater")

dev.off()

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CT_base_std_rootbias.pdf",height=5,width=15)

par(mfrow=c(1,2))
feature<-6
band<-0.2
g<-density(CT_R1[,feature],bw=band)

plot(g, col=rgb(0,0,0.7,1/2), main="CT base value", xlim=c(min(CT_R1[,feature],CT_R4[,feature],CT_R5[,feature],CT_R6[,feature],CT_L2[,feature],CT_L3[,feature],CT_L7[,feature],CT_L8[,feature]),max(CT_R1[,feature],CT_R4[,feature],CT_R5[,feature],CT_R6[,feature],CT_L2[,feature],CT_L3[,feature],CT_L7[,feature],CT_L8[,feature])),
		ylim=c(0,1.6),lwd=2)
lines(density(CT_R4[,feature],bw=band), col=rgb(0,0,0.8,1/2),lwd=2)
lines(density(CT_R5[,feature],bw=band), col=rgb(0,0,0.9,1/2),lwd=2)
lines(density(CT_R6[,feature],bw=band), col=rgb(0,0,1,1/2),lwd=2)

lines(density(CT_L2[,feature],bw=band), col=rgb(0.7,0,0,1/2),lwd=2)
lines(density(CT_L3[,feature],bw=band), col=rgb(0.8,0,0,1/2),lwd=2)
lines(density(CT_L7[,feature],bw=band), col=rgb(0.9,0,0,1/2),lwd=2)
lines(density(CT_L8[,feature],bw=band), col=rgb(1,0,0,1/2),lwd=2)

feature<-7
band<-0.2
g<-density(CT_R1[,feature],bw=band)

plot(g, col=rgb(0,0,0.7,1/2), main="CT standard deviation", xlim=c(min(CT_R1[,feature],CT_R4[,feature],CT_R5[,feature],CT_R6[,feature],CT_L2[,feature],CT_L3[,feature],CT_L7[,feature],CT_L8[,feature]),max(CT_R1[,feature],CT_R4[,feature],CT_R5[,feature],CT_R6[,feature],CT_L2[,feature],CT_L3[,feature],CT_L7[,feature],CT_L8[,feature])),
		ylim=c(0,1.05),lwd=2)
lines(density(CT_R4[,feature],bw=band), col=rgb(0,0,0.8,1/2),lwd=2)
lines(density(CT_R5[,feature],bw=band), col=rgb(0,0,0.9,1/2),lwd=2)
lines(density(CT_R6[,feature],bw=band), col=rgb(0,0,1,1/2),lwd=2)

lines(density(CT_L2[,feature],bw=band), col=rgb(0.7,0,0,1/2),lwd=2)
lines(density(CT_L3[,feature],bw=band), col=rgb(0.8,0,0,1/2),lwd=2)
lines(density(CT_L7[,feature],bw=band), col=rgb(0.9,0,0,1/2),lwd=2)
lines(density(CT_L8[,feature],bw=band), col=rgb(1,0,0,1/2),lwd=2)

dev.off()


pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CT_mass_barplot.pdf",width=9,height=5)

feature<-23

massesR<-matrix(c(mean(CT_R1[,feature]),mean(CT_R1[,feature+1]),mean(CT_R1[,feature+2]),mean(CT_R1[,feature+3]),mean(CT_R1[,feature+4]),mean(CT_R4[,feature]),mean(CT_R4[,feature+1]),mean(CT_R4[,feature+2]),mean(CT_R4[,feature+3]),mean(CT_R4[,feature+4]),
		mean(CT_R5[,feature]),mean(CT_R5[,feature+1]),mean(CT_R5[,feature+2]),mean(CT_R5[,feature+3]),mean(CT_R5[,feature+4]),mean(CT_R6[,feature]),mean(CT_R6[,feature+1]),mean(CT_R6[,feature+2]),mean(CT_R6[,feature+3]),mean(CT_R6[,feature+4])),nrow=5)
massesL<-matrix(c(mean(CT_L2[,feature]),mean(CT_L2[,feature+1]),mean(CT_L2[,feature+2]),mean(CT_L2[,feature+3]),mean(CT_L2[,feature+4]),mean(CT_L3[,feature]),mean(CT_L3[,feature+1]),mean(CT_L3[,feature+2]),mean(CT_L3[,feature+3]),mean(CT_L3[,feature+4]),
				mean(CT_L7[,feature]),mean(CT_L7[,feature+1]),mean(CT_L7[,feature+2]),mean(CT_L7[,feature+3]),mean(CT_L7[,feature+4]),mean(CT_L8[,feature]),mean(CT_L8[,feature+1]),mean(CT_L8[,feature+2]),mean(CT_L8[,feature+3]),mean(CT_L8[,feature+4])),nrow=5)
masses<-cbind(massesL,massesR)
masses<-apply(masses,2,rev)
par(mar=c(5.1,4.1,4.1,6.1))
barplot(masses,col=rev(c("red","blue","green","pink","orange")),names.arg=c("CT_L2", "CT_L3", "CT_L7", "CT_L8", "CT_R1", "CT_R4", "CT_R5", "CT_R6"),main="mean posterior layer-wise tumor masses",
		args.legend = list(x=11.5, y=1), bty = "n")


x<-0
massesR<-massesR+0.00000001
massesL<-massesL+0.00000001
rightDiff<-rep(0,6)
leftDiff<-rep(0,6)
for (j in 1:3 ) {
	for (i in (j+1):4) {
		x<-x+1
			rightDiff[x]<-(sum(log(massesR[,j]/massesR[,i])*massesR[,j])+sum(log(massesR[,i]/massesR[,j])*massesR[,i]))/2
			leftDiff[x]<-(sum(log(massesL[,j]/massesL[,i])*massesL[,j])+sum(log(massesL[,i]/massesL[,j])*massesL[,i]))/2
		}
}
wilcox.test(rightDiff,leftDiff,alternative='less')

dev.off()



CT_R1branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/incorrespondent-6142/CT_IRX2P_R1.csv/traces-CT_IRX2P_R1.csv_branch_traces.csv",sep=",",header=F)
CT_R4branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/scolia-5868/CT_IRX2P_R4.csv/traces-CT_IRX2P_R4.csv_branch_traces.csv",sep=",",header=F)
CT_R5branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/neuromerous-9376/CT_IRX2P_R5.csv/traces-CT_IRX2P_R5.csv_branch_traces.csv",sep=",",header=F)
CT_R6branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/proaccelerin-3347/CT_IRX2P_R6.csv/traces-CT_IRX2P_R6.csv_branch_traces.csv",sep=",",header=F)

CT_L2branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/oceanity-4938/CT_IRX2P_L2.csv/traces-CT_IRX2P_L2.csv_branch_traces.csv",sep=",",header=F)
CT_L3branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/nonprolific-2900/CT_IRX2P_L3.csv/traces-CT_IRX2P_L3.csv_branch_traces.csv",sep=",",header=F)
CT_L7branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/buoy-3083/CT_IRX2P_L7.csv/traces-CT_IRX2P_L7.csv_branch_traces.csv",sep=",",header=F)
CT_L8branch<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/epistyle-3189/CT_IRX2P_L8.csv/traces-CT_IRX2P_L8.csv_branch_traces.csv",sep=",",header=F)

pdf("CT_total_branch_length.pdf")

band=0.5
tb1<-0
tb1list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb1list[i]<-sum(unique(unlist(CT_R1branch[i,])))
}
tb1<-sum(tb1list)
plot(density(tb1list,bw=band ), ylim=c(0,0.7), xlim=c(0,12), col=rgb(0,0.7,0.7,1/2),lwd=2)
axis(1,at=c(1,3,5,10),labels=c(1,3,5,10))

tb4<-0
tb4list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb4list[i]<-sum(unique(unlist(CT_R4branch[i,])))
}
tb4<-sum(tb4list)
lines(density(tb4list,bw=band ), col=rgb(0,0.8,0.8,1/2),lwd=2)

tb5<-0
tb5list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb5list[i]<-sum(unique(unlist(CT_R5branch[i,])))
}
tb5<-sum(tb5list)
lines(density(tb5list,bw=band ), col=rgb(0,0.9,0.9,1/2),lwd=2)

tb6<-0
tb6list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb6list[i]<-sum(unique(unlist(CT_R6branch[i,])))
}
tb6<-sum(tb6list)
lines(density(tb6list,bw=band ), col=rgb(0,1,1,1/2),lwd=2)

tb2<-0
tb2list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb2list[i]<-sum(unique(unlist(CT_L2branch[i,])))
}
tb2<-sum(tb2list)
lines(density(tb2list,bw=band ), col=rgb(0.7,0,0.7,1/2),lwd=2)

tb3<-0
tb3list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb3list[i]<-sum(unique(unlist(CT_L3branch[i,])))
}
tb3<-sum(tb3list)
lines(density(tb3list,bw=band ), col=rgb(0.8,0,0.8,1/2),lwd=2)

tb7<-0
tb7list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb7list[i]<-sum(unique(unlist(CT_L7branch[i,])))
}
tb7<-sum(tb7list)
lines(density(tb7list,bw=band ), col=rgb(0.9,0,0.9,1/2),lwd=2)

tb8<-0
tb8list<-rep(0,10000)
for ( i in 1:10000 ) {
	tb8list[i]<-sum(unique(unlist(CT_L8branch[i,])))
}
tb8<-sum(tb8list)
lines(density(tb8list,bw=band ), col=rgb(1,0,1,1/2),lwd=2)

t.test(c(mean(tb1list),mean(tb4list),mean(tb5list),mean(tb6list)),c(mean(tb2list),mean(tb3list),mean(tb7list),mean(tb8list)))

dev.off()

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CT_bignodes_nodenum.pdf")

feature1<-3
feature2<-5
plot(mean(CT_R1[,feature1]),mean(CT_R1[,feature2]),xlim=c(2,12),ylim=c(2,6))
points(mean(CT_R4[,feature1]),mean(CT_R4[,feature2]))
points(mean(CT_R5[,feature1]),mean(CT_R5[,feature2]))
points(mean(CT_R6[,feature1]),mean(CT_R6[,feature2]))

points(mean(CT_L2[,feature1]),mean(CT_L2[,feature2]))
points(mean(CT_L3[,feature1]),mean(CT_L3[,feature2]))
points(mean(CT_L8[,feature1]),mean(CT_L7[,feature2]))
points(mean(CT_L7[,feature1]),mean(CT_L8[,feature2]))

dev.off()
