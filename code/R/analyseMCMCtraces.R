R1<-read.table("git/phylo-tree/code/mcmc-traces/teatime_CT_IRX2P_R1.csv_traces.csv",sep=",",header=T)
R4<-read.table("git/phylo-tree/code/mcmc-traces/teatime_CT_IRX2P_R4.csv_traces.csv",sep=",",header=T)
R5<-read.table("git/phylo-tree/code/mcmc-traces/teatime_CT_IRX2P_R5.csv_traces.csv",sep=",",header=T)
R6<-read.table("git/phylo-tree/code/mcmc-traces/teatime_CT_IRX2P_R6.csv_traces.csv",sep=",",header=T)

RC1<-read.table("git/phylo-tree/code/mcmc-traces/saturnalia_CT_IRX2P_R1.csv_1234_traces.csv",sep=",",header=T)
RC4<-read.table("git/phylo-tree/code/mcmc-traces/saturnalia_CT_IRX2P_R4.csv_1234_traces.csv",sep=",",header=T)
RC5<-read.table("git/phylo-tree/code/mcmc-traces/saturnalia_CT_IRX2P_R5.csv_1234_traces.csv",sep=",",header=T)
RC6<-read.table("git/phylo-tree/code/mcmc-traces/saturnalia_CT_IRX2P_R6.csv_1234_traces.csv",sep=",",header=T)


pdf(file="git/phylo-tree/code/figures/sottoriva/CT_right_logistic_mixlaplace_methy_test.pdf")

range=c(5000:20000)

for( i in 1:dim(R1)[2]) {
	par(mfrow=c(2,2))
	plot(R1[range,i],ylab=colnames(R1)[i], type="l",ylim=c(min(R1[range,i],RC1[range,i]),max(R1[range,i],RC1[range,i])),xlab="iteration 5000 to 20000")
	lines(RC1[range,i],ylab=colnames(R1)[i], col="red")
	plot(R4[range,i],ylab=colnames(R1)[i], type="l",ylim=c(min(R4[range,i],RC4[range,i]),max(R4[range,i],RC4[range,i])),xlab="iteration 5000 to 20000")
	lines(RC4[range,i],ylab=colnames(R1)[i], col="red")
	plot(R5[range,i],ylab=colnames(R1)[i], type="l",ylim=c(min(R5[range,i],RC5[range,i]),max(R5[range,i],RC5[range,i])),xlab="iteration 5000 to 20000")
	lines(RC5[range,i],ylab=colnames(R1)[i], col="red")
	plot(R6[range,i],ylab=colnames(R1)[i], type="l",ylim=c(min(R6[range,i],RC6[range,i]),max(R6[range,i],RC6[range,i])),xlab="iteration 5000 to 20000")
	lines(RC6[range,i],ylab=colnames(R1)[i], col="red")
}

dev.off()