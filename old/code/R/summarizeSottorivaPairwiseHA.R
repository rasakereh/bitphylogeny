HA_R5<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/backboards-6156/HA_IRX2P_R5.csv/traces-HA_IRX2P_R5.csv",sep=",",header=T)
HA_R7<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/uprose-8122/HA_IRX2P_R7.csv/traces-HA_IRX2P_R7.csv",sep=",",header=T)
HA_R8<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/prasinous-7644/HA_IRX2P_R8.csv/traces-HA_IRX2P_R8.csv",sep=",",header=T)
HA_R9<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/hiccough-9006/HA_IRX2P_R9.csv/traces-HA_IRX2P_R9.csv",sep=",",header=T)
HA_R10<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/slickers-4748/HA_IRX2P_R10.csv/traces-HA_IRX2P_R10.csv",sep=",",header=T)


HA_L1<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/toxoids-25/HA_IRX2P_L1.csv/traces-HA_IRX2P_L1.csv",sep=",",header=T)
HA_L2<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/guttier-7649/HA_IRX2P_L2.csv/traces-HA_IRX2P_L2.csv",sep=",",header=T)
HA_L3<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/nonirrigated-1537/HA_IRX2P_L3.csv/traces-HA_IRX2P_L3.csv",sep=",",header=T)
HA_L4<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/dismissingly-4485/HA_IRX2P_L4.csv/traces-HA_IRX2P_L4.csv",sep=",",header=T)
HA_L6<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/jacarandi-2246/HA_IRX2P_L6.csv/traces-HA_IRX2P_L6.csv",sep=",",header=T)


pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_HA_nodenum.pdf")

feature<-2
band<-0.7
g<-density(HA_R5[,feature],bw=band)

plot(g, col=rgb(0,0.6,0.6,1/2), lwd=2, main="HA node numbers", xlim=c(min(HA_R5[,feature],HA_R7[,feature],HA_R8[,feature],HA_R9[,feature],HA_R10[,feature],HA_L1[,feature],HA_L2[,feature],HA_L4[,feature],HA_L6[,feature]),max(HA_R5[,feature],HA_R7[,feature],HA_R9[,feature],HA_R10[,feature],HA_L1[,feature],HA_L2[,feature],HA_L3[,feature],HA_L4[,feature],HA_L6[,feature])),
	ylim=c(0,0.35))
lines(density(HA_R7[,feature],bw=band), col=rgb(0,0.7,0.7,1/2), lwd=2)
lines(density(HA_R8[,feature],bw=band), col=rgb(0,0.8,0.8,1/2), lwd=2)
lines(density(HA_R9[,feature],bw=band), col=rgb(0,0.9,0.9,1/2), lwd=2)
lines(density(HA_R10[,feature],bw=band), col=rgb(0,1,1,1/2), lwd=2)

lines(density(HA_L1[,feature],bw=band), col=rgb(0.6,0,0.6,1/2), lwd=2)
lines(density(HA_L2[,feature],bw=band), col=rgb(0.7,0,0.7,1/2), lwd=2)
lines(density(HA_L3[,feature],bw=band), col=rgb(0.8,0,0.8,1/2), lwd=2)
lines(density(HA_L4[,feature],bw=band), col=rgb(0.9,0,0.9,1/2), lwd=2)
lines(density(HA_L6[,feature],bw=band), col=rgb(1,0,0,1/2), lwd=2)

dev.off()

#pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_HA_bignodenum.pdf")

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_HA_bignodenum_depth.pdf",width=10,height=5)
par(mfrow=c(1,2))

feature<-3
band<-0.7
g<-density(HA_R5[,feature],bw=band)

plot(g, col=rgb(0,0.6,0.6,1/2), lwd=2, main="HA big node numbers", xlim=c(min(HA_R5[,feature],HA_R7[,feature],HA_R8[,feature],HA_R9[,feature],HA_R10[,feature],HA_L1[,feature],HA_L2[,feature],HA_L4[,feature],HA_L6[,feature]),max(HA_R5[,feature],HA_R7[,feature],HA_R9[,feature],HA_R10[,feature],HA_L1[,feature],HA_L2[,feature],HA_L3[,feature],HA_L4[,feature],HA_L6[,feature])))
lines(density(HA_R7[,feature],bw=band), col=rgb(0,0.7,0.7,1/2), lwd=2)
lines(density(HA_R8[,feature],bw=band), col=rgb(0,0.8,0.8,1/2), lwd=2)
lines(density(HA_R9[,feature],bw=band), col=rgb(0,0.9,0.9,1/2), lwd=2)
lines(density(HA_R10[,feature],bw=band), col=rgb(0,1,1,1/2), lwd=2)

lines(density(HA_L1[,feature],bw=band), col=rgb(0.6,0,0.6,1/2), lwd=2)
lines(density(HA_L2[,feature],bw=band), col=rgb(0.7,0,0.7,1/2), lwd=2)
lines(density(HA_L3[,feature],bw=band), col=rgb(0.8,0,0.8,1/2), lwd=2)
lines(density(HA_L4[,feature],bw=band), col=rgb(0.9,0,0.9,1/2), lwd=2)
lines(density(HA_L6[,feature],bw=band), col=rgb(1,0,0,1/2), lwd=2)

#dev.off()

#pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_HA_depth.pdf")

feature<-5
band<-0.2
g<-density(HA_R5[,feature],bw=band)

plot(g, col=rgb(0,0.6,0.6,1/2), lwd=2, main="HA max depth", xlim=c(min(HA_R5[,feature],HA_R7[,feature],HA_R8[,feature],HA_R9[,feature],HA_R10[,feature],HA_L1[,feature],HA_L2[,feature],HA_L4[,feature],HA_L6[,feature]),max(HA_R5[,feature],HA_R7[,feature],HA_R9[,feature],HA_R10[,feature],HA_L1[,feature],HA_L2[,feature],HA_L3[,feature],HA_L4[,feature],HA_L6[,feature])))
lines(density(HA_R7[,feature]), col=rgb(0,0.7,0.7,1/2), lwd=2)
lines(density(HA_R8[,feature]), col=rgb(0,0.8,0.8,1/2), lwd=2)
lines(density(HA_R9[,feature]), col=rgb(0,0.9,0.9,1/2), lwd=2)
lines(density(HA_R10[,feature]), col=rgb(0,1,1,1/2), lwd=2)

lines(density(HA_L1[,feature]), col=rgb(0.6,0,0.6,1/2), lwd=2)
lines(density(HA_L2[,feature]), col=rgb(0.7,0,0.7,1/2), lwd=2)
lines(density(HA_L3[,feature]), col=rgb(0.8,0,0.8,1/2), lwd=2)
lines(density(HA_L4[,feature]), col=rgb(0.9,0,0.9,1/2), lwd=2)
lines(density(HA_L6[,feature]), col=rgb(1,0,0,1/2), lwd=2)

dev.off()

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_HA_base_std_rootbias.pdf",height=5,width=15)

par(mfrow=c(1,2))
feature<-6
band<-0.2
g<-density(HA_R5[,feature],bw=band)

plot(g, col=rgb(0,0.6,0.6,1/2), lwd=2, main="HA base value", xlim=c(min(HA_R5[,feature],HA_R7[,feature],HA_R8[,feature],HA_R9[,feature],HA_R10[,feature],HA_L1[,feature],HA_L2[,feature],HA_L4[,feature],HA_L6[,feature]),max(HA_R5[,feature],HA_R7[,feature],HA_R9[,feature],HA_R10[,feature],HA_L1[,feature],HA_L2[,feature],HA_L3[,feature],HA_L4[,feature],HA_L6[,feature])))
lines(density(HA_R7[,feature]), col=rgb(0,0.7,0.7,1/2), lwd=2)
lines(density(HA_R8[,feature]), col=rgb(0,0.8,0.8,1/2), lwd=2)
lines(density(HA_R9[,feature]), col=rgb(0,0.9,0.9,1/2), lwd=2)
lines(density(HA_R10[,feature]), col=rgb(0,1,1,1/2), lwd=2)

lines(density(HA_L1[,feature]), col=rgb(0.6,0,0.6,1/2), lwd=2)
lines(density(HA_L2[,feature]), col=rgb(0.7,0,0.7,1/2), lwd=2)
lines(density(HA_L3[,feature]), col=rgb(0.8,0,0.8,1/2), lwd=2)
lines(density(HA_L4[,feature]), col=rgb(0.9,0,0.9,1/2), lwd=2)
lines(density(HA_L6[,feature]), col=rgb(1,0,0,1/2), lwd=2)

feature<-7
band<-0.2
g<-density(HA_R5[,feature],bw=band)

plot(g, col=rgb(0,0.6,0.6,1/2), lwd=2, main="HA standard deviation", xlim=c(min(HA_R5[,feature],HA_R7[,feature],HA_R8[,feature],HA_R9[,feature],HA_R10[,feature],HA_L1[,feature],HA_L2[,feature],HA_L4[,feature],HA_L6[,feature]),max(HA_R5[,feature],HA_R7[,feature],HA_R9[,feature],HA_R10[,feature],HA_L1[,feature],HA_L2[,feature],HA_L3[,feature],HA_L4[,feature],HA_L6[,feature])))
lines(density(HA_R7[,feature]), col=rgb(0,0.7,0.7,1/2), lwd=2)
lines(density(HA_R8[,feature]), col=rgb(0,0.8,0.8,1/2), lwd=2)
lines(density(HA_R9[,feature]), col=rgb(0,0.9,0.9,1/2), lwd=2)
lines(density(HA_R10[,feature]), col=rgb(0,1,1,1/2), lwd=2)

lines(density(HA_L1[,feature]), col=rgb(0.6,0,0.6,1/2), lwd=2)
lines(density(HA_L2[,feature]), col=rgb(0.7,0,0.7,1/2), lwd=2)
lines(density(HA_L3[,feature]), col=rgb(0.8,0,0.8,1/2), lwd=2)
lines(density(HA_L4[,feature]), col=rgb(0.9,0,0.9,1/2), lwd=2)
lines(density(HA_L6[,feature]), col=rgb(1,0,0,1/2), lwd=2)

dev.off()


pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_HA_mass_barplot.pdf",width=9,height=5)

feature<-23

massesR<-matrix(c(mean(HA_R5[,feature]),mean(HA_R5[,feature+1]),mean(HA_R5[,feature+2]),mean(HA_R5[,feature+3]),mean(HA_R5[,feature+4]),mean(HA_R7[,feature]),mean(HA_R7[,feature+1]),mean(HA_R7[,feature+2]),mean(HA_R7[,feature+3]),mean(HA_R7[,feature+4]),
				mean(HA_R8[,feature]),mean(HA_R8[,feature+1]),mean(HA_R8[,feature+2]),mean(HA_R8[,feature+3]),mean(HA_R8[,feature+4]),
		mean(HA_R9[,feature]),mean(HA_R9[,feature+1]),mean(HA_R9[,feature+2]),mean(HA_R9[,feature+3]),mean(HA_R9[,feature+4]),mean(HA_R10[,feature]),mean(HA_R10[,feature+1]),mean(HA_R10[,feature+2]),mean(HA_R10[,feature+3]),mean(HA_R10[,feature+4])),nrow=5)
massesL<-matrix(c(mean(HA_L1[,feature]),mean(HA_L1[,feature+1]),mean(HA_L1[,feature+2]),mean(HA_L1[,feature+3]),mean(HA_L1[,feature+4]),mean(HA_L2[,feature]),mean(HA_L2[,feature+1]),mean(HA_L2[,feature+2]),mean(HA_L2[,feature+3]),mean(HA_L2[,feature+4]),
				mean(HA_L3[,feature]),mean(HA_L3[,feature+1]),mean(HA_L3[,feature+2]),mean(HA_L3[,feature+3]),mean(HA_L3[,feature+4]),
				mean(HA_L4[,feature]),mean(HA_L4[,feature+1]),mean(HA_L4[,feature+2]),mean(HA_L4[,feature+3]),mean(HA_L4[,feature+4]),mean(HA_L6[,feature]),mean(HA_L6[,feature+1]),mean(HA_L6[,feature+2]),mean(HA_L6[,feature+3]),mean(HA_L6[,feature+4])),nrow=5)
masses<-cbind(massesL,massesR)
masses<-apply(masses,2,rev)
par(mar=c(5.1,4.1,4.1,6.1))
barplot(masses,col=rev(c("red","blue","green","pink","orange")),names.arg=c("HA_L1", "HA_L2", "HA_L3", "HA_L4", "HA_L6", "HA_R5", "HA_R7", "HA_R8", "HA_R9", "HA_R10"),main="mean posterior layer-wise tumor masses",
		args.legend = list(x=11.5, y=1), bty = "n")


dev.off()
