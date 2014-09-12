CU_R1<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/unfixedness-3111/CU_IRX2P_R1.csv/traces-CU_IRX2P_R1.csv",sep=",",header=T)
CU_R2<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/sigillistic-3254/CU_IRX2P_R2.csv/traces-CU_IRX2P_R2.csv",sep=",",header=T)
CU_R3<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/hypotrophy-9619/CU_IRX2P_R3.csv/traces-CU_IRX2P_R3.csv",sep=",",header=T)
CU_R5<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/elocular-9411/CU_IRX2P_R5.csv/traces-CU_IRX2P_R5.csv",sep=",",header=T)

CU_L4<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/steersmen-5353/CU_IRX2P_L4.csv/traces-CU_IRX2P_L4.csv",sep=",",header=T)
CU_L6<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/ungrouped-3338/CU_IRX2P_L6.csv/traces-CU_IRX2P_L6.csv",sep=",",header=T)
CU_L7<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/circumnavigations-9284/CU_IRX2P_L7.csv/traces-CU_IRX2P_L7.csv",sep=",",header=T)
CU_L8<-read.table("/mnt/home/final_mcmcrun_incomplete/Sottoriva/final/ransom-2122/CU_IRX2P_L8.csv/traces-CU_IRX2P_L8.csv",sep=",",header=T)

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CU_nodenum.pdf")

feature<-2
g<-density(CU_R1[,feature])

plot(g, col=rgb(0,0,0.25,1/4), main="CU node numbers", xlim=c(min(CU_R1[,feature],CU_R2[,feature],CU_R3[,feature],CU_R5[,feature],CU_L4[,feature],CU_L6[,feature],CU_L7[,feature],CU_L8[,feature]),max(CU_R1[,feature],CU_R2[,feature],CU_R3[,feature],CU_R5[,feature],CU_L4[,feature],CU_L6[,feature],CU_L7[,feature],CU_L8[,feature])))
lines(density(CU_R2[,feature]), col=rgb(0,0,0.5,1/4))
lines(density(CU_R3[,feature]), col=rgb(0,0,0.75,1/4))
lines(density(CU_R5[,feature]), col=rgb(0,0,1,1/4))

lines(density(CU_L4[,feature]), col=rgb(0.25,0,0,1/4))
lines(density(CU_L6[,feature]), col=rgb(0.5,0,0,1/4))
lines(density(CU_L7[,feature]), col=rgb(0.75,0,0,1/4))
lines(density(CU_L8[,feature]), col=rgb(1,0,0,1/4))

dev.off()

#pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CU_bignodenum.pdf")

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CU_bignodenum_depth.pdf",width=10,height=5)
par(mfrow=c(1,2))

feature<-3
g<-density(CU_R1[,feature],bw=1)

plot(g, col=rgb(0,0,0.7,1/2), main="CU big node numbers", xlim=c(min(CU_R1[,feature],CU_R2[,feature],CU_R3[,feature],CU_R5[,feature],CU_L4[,feature],CU_L6[,feature],CU_L7[,feature],CU_L8[,feature]),max(CU_R1[,feature],CU_R2[,feature],CU_R3[,feature],CU_R5[,feature],CU_L4[,feature],CU_L6[,feature],CU_L7[,feature],CU_L8[,feature])),
		ylim=c(0,0.5),lwd=2)
lines(density(CU_R2[,feature],bw=1), col=rgb(0,0,0.8,1/2),lwd=2)
lines(density(CU_R3[,feature],bw=1), col=rgb(0,0,0.9,1/2),lwd=2)
lines(density(CU_R5[,feature],bw=1), col=rgb(0,0,1,1/2),lwd=2)

lines(density(CU_L4[,feature],bw=1), col=rgb(0.7,0,0,1/2),lwd=2)
lines(density(CU_L6[,feature],bw=1), col=rgb(0.8,0,0,1/2),lwd=2)
lines(density(CU_L7[,feature],bw=1), col=rgb(0.9,0,0,1/2),lwd=2)
lines(density(CU_L8[,feature],bw=1), col=rgb(1,0,0,1/2),lwd=2)

#dev.off()

#pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CU_depth.pdf")

feature<-5
band<-0.2
g<-density(CU_R1[,feature],bw=band)

plot(g, col=rgb(0,0,0.7,1/2), main="CU max depth", xlim=c(min(CU_R1[,feature],CU_R2[,feature],CU_R3[,feature],CU_R5[,feature],CU_L4[,feature],CU_L6[,feature],CU_L7[,feature],CU_L8[,feature]),max(CU_R1[,feature],CU_R2[,feature],CU_R3[,feature],CU_R5[,feature],CU_L4[,feature],CU_L6[,feature],CU_L7[,feature],CU_L8[,feature])),
		ylim=c(0,2),lwd=2, xaxt='n')
axis(1, at=c(2,3,4,5), labels=c(2,3,4,5))
lines(density(CU_R2[,feature],bw=band), col=rgb(0,0,0.8,1/2),lwd=2)
lines(density(CU_R3[,feature],bw=band), col=rgb(0,0,0.9,1/2),lwd=2)
lines(density(CU_R5[,feature],bw=band), col=rgb(0,0,1,1/2),lwd=2)

lines(density(CU_L4[,feature],bw=band), col=rgb(0.7,0,0,1/2),lwd=2)
lines(density(CU_L6[,feature],bw=band), col=rgb(0.8,0,0,1/2),lwd=2)
lines(density(CU_L7[,feature],bw=band), col=rgb(0.9,0,0,1/2),lwd=2)
lines(density(CU_L8[,feature],bw=band), col=rgb(1,0,0,1/2),lwd=2)

dev.off()

pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CU_base_std_rootbias.pdf",height=5,width=15)

par(mfrow=c(1,2))
feature<-6
band<-0.2
g<-density(CU_R1[,feature],bw=band)

plot(g, col=rgb(0,0,0.7,1/2), main="CU base value", xlim=c(min(CU_R1[,feature],CU_R2[,feature],CU_R3[,feature],CU_R5[,feature],CU_L4[,feature],CU_L6[,feature],CU_L7[,feature],CU_L8[,feature]),max(CU_R1[,feature],CU_R2[,feature],CU_R3[,feature],CU_R5[,feature],CU_L4[,feature],CU_L6[,feature],CU_L7[,feature],CU_L8[,feature])),
		ylim=c(0,1.6),lwd=2)
lines(density(CU_R2[,feature],bw=band), col=rgb(0,0,0.8,1/2),lwd=2)
lines(density(CU_R3[,feature],bw=band), col=rgb(0,0,0.9,1/2),lwd=2)
lines(density(CU_R5[,feature],bw=band), col=rgb(0,0,1,1/2),lwd=2)

lines(density(CU_L4[,feature],bw=band), col=rgb(0.7,0,0,1/2),lwd=2)
lines(density(CU_L6[,feature],bw=band), col=rgb(0.8,0,0,1/2),lwd=2)
lines(density(CU_L7[,feature],bw=band), col=rgb(0.9,0,0,1/2),lwd=2)
lines(density(CU_L8[,feature],bw=band), col=rgb(1,0,0,1/2),lwd=2)

feature<-7
band<-0.2
g<-density(CU_R1[,feature],bw=band)

plot(g, col=rgb(0,0,0.7,1/2), main="CU standard deviation", xlim=c(min(CU_R1[,feature],CU_R2[,feature],CU_R3[,feature],CU_R5[,feature],CU_L4[,feature],CU_L6[,feature],CU_L7[,feature],CU_L8[,feature]),max(CU_R1[,feature],CU_R2[,feature],CU_R3[,feature],CU_R5[,feature],CU_L4[,feature],CU_L6[,feature],CU_L7[,feature],CU_L8[,feature])),
		ylim=c(0,1.05),lwd=2)
lines(density(CU_R2[,feature],bw=band), col=rgb(0,0,0.8,1/2),lwd=2)
lines(density(CU_R3[,feature],bw=band), col=rgb(0,0,0.9,1/2),lwd=2)
lines(density(CU_R5[,feature],bw=band), col=rgb(0,0,1,1/2),lwd=2)

lines(density(CU_L4[,feature],bw=band), col=rgb(0.7,0,0,1/2),lwd=2)
lines(density(CU_L6[,feature],bw=band), col=rgb(0.8,0,0,1/2),lwd=2)
lines(density(CU_L7[,feature],bw=band), col=rgb(0.9,0,0,1/2),lwd=2)
lines(density(CU_L8[,feature],bw=band), col=rgb(1,0,0,1/2),lwd=2)

dev.off()


pdf(file="git/phylo-tree/code/figures/sottoriva/Sottoriva_pairwise_CU_mass_barplot.pdf",width=9,height=5)

feature<-23

massesR<-matrix(c(mean(CU_R1[,feature]),mean(CU_R1[,feature+1]),mean(CU_R1[,feature+2]),mean(CU_R1[,feature+3]),mean(CU_R1[,feature+4]),mean(CU_R2[,feature]),mean(CU_R2[,feature+1]),mean(CU_R2[,feature+2]),mean(CU_R2[,feature+3]),mean(CU_R2[,feature+4]),
		mean(CU_R3[,feature]),mean(CU_R3[,feature+1]),mean(CU_R3[,feature+2]),mean(CU_R3[,feature+3]),mean(CU_R3[,feature+4]),mean(CU_R5[,feature]),mean(CU_R5[,feature+1]),mean(CU_R5[,feature+2]),mean(CU_R5[,feature+3]),mean(CU_R5[,feature+4])),nrow=5)
massesL<-matrix(c(mean(CU_L4[,feature]),mean(CU_L4[,feature+1]),mean(CU_L4[,feature+2]),mean(CU_L4[,feature+3]),mean(CU_L4[,feature+4]),mean(CU_L6[,feature]),mean(CU_L6[,feature+1]),mean(CU_L6[,feature+2]),mean(CU_L6[,feature+3]),mean(CU_L6[,feature+4]),
				mean(CU_L7[,feature]),mean(CU_L7[,feature+1]),mean(CU_L7[,feature+2]),mean(CU_L7[,feature+3]),mean(CU_L7[,feature+4]),mean(CU_L8[,feature]),mean(CU_L8[,feature+1]),mean(CU_L8[,feature+2]),mean(CU_L8[,feature+3]),mean(CU_L8[,feature+4])),nrow=5)
masses<-cbind(massesL,massesR)
masses<-apply(masses,2,rev)
par(mar=c(5.1,4.1,4.1,6.1))
barplot(masses,col=rev(c("red","blue","green","pink","orange")),names.arg=c("CU_L4", "CU_L6", "CU_L7", "CU_L8", "CU_R1", "CU_R2", "CU_R3", "CU_R5"),main="mean posterior layer-wise tumor masses",
		args.legend = list(x=11.5, y=1), bty = "n")


dev.off()
