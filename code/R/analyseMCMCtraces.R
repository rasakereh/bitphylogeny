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

pdf(file="git/phylo-tree/code/figures/sottoriva/CT_R6_long_logistic_mixlaplace_methy_test.pdf")

Rx6_long<-read.table("git/phylo-tree/code/mcmc-traces/CT_R6_long_1234_traces.csv",sep="\t",header=T)

RC6_long<-read.table("git/phylo-tree/code/mcmc-traces/alphabetized_CT_IRX2P_R6.csv_1264_traces.csv",sep=",",header=F)

range=c(1000:49999)

for( i in 1:dim(Rx6_long)[2]) {
	par(mfrow=c(2,2))
	plot(Rx6_long[range,i],ylab=colnames(Rx6_long)[i], type="l",ylim=c(min(Rx6_long[range,i],RC6_long[range,i]),max(Rx6_long[range,i],RC6_long[range,i])),xlab="iteration 1000 to 49999")
	lines(RC6_long[range,i],ylab=colnames(Rx6_long)[i], col="red")
	hist(Rx6_long[range,i])
	hist(RC6_long[range,i])
}

dev.off()

pdf(file="git/phylo-tree/code/figures/sottoriva/CT_R6_long_logistic_mixlaplace_methy_test_noprior_sampling.pdf")

Rx6_long<-read.table("git/phylo-tree/code/mcmc-traces/quintessences_CT_IRX2P_R6.csv_1264_traces.csv",sep=",",header=T)

RC6_long<-read.table("git/phylo-tree/code/mcmc-traces/mangled_CT_IRX2P_R6.csv_1234_traces.csv",sep=",",header=T)

range=c(1000:49999)

for( i in 1:dim(Rx6_long)[2]) {
	par(mfrow=c(2,2))
	plot(Rx6_long[range,i],ylab=colnames(Rx6_long)[i], type="l",ylim=c(min(Rx6_long[range,i],RC6_long[range,i]),max(Rx6_long[range,i],RC6_long[range,i])),xlab="iteration 1000 to 49999")
	lines(RC6_long[range,i],ylab=colnames(Rx6_long)[i], col="red")
	hist(Rx6_long[range,i])
	hist(RC6_long[range,i])
}

dev.off()

