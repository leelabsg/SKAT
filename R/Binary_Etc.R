SKATPermu<-function(Z, y, nPermu=100000, epsilon=10^-6){


	nSNP = ncol(Z);
	nSample = length(y);

	#RSKATPermu(double *Z, int *Y, int *nSNP, int *nSample, int *nPermu,double * pval, double *pval_same,  double * epsilon){

	RE<-.C("RSKATPermu", as.double(as.vector(Z)), as.integer(y),  as.integer(nSNP),  as.integer(nSample)
		, as.integer(nPermu), double(1), double(1), as.double(epsilon));
		

	pval.re<-cbind(RE[[6]], RE[[7]])	
	pval1<-cbind(pval.re[,1], pval.re[,1] - pval.re[,2] /2, pval.re[,2])
	re=list(p.value=pval1[1,2], p.value.standard=pval1[1,1]);
	
	return(re);
		
}

#dyn.load("~/Project/Kernel_Machine/R_Package/SKAT/src/Binary_ComputeExact.so")

#nSample<-50000
#nSNP<-20
#Z<-matrix(rep(0,nSample*nSNP), ncol=nSNP)
#y<-rep(c(1,0), each=nSample/2)
#system.time((re=SKATPermu(Z, y, nPermu=100000, epsilon=10^-6)))





