
#
#	param should be n.r *3 matrix
#
#File, re.arr$resarray, re.arr$nres, re.arr$nres_k, Z0, Z1, test_Z0, pr$k, 
#			m, pr$n.total, pr$n.total.k, pr$prob_k, odds, p1_adj, pr$IsExact
			
Generate_TestData<-function(File, resarray, nres, nres_k, Z0, Z1, test_Z0, k, m, n.total, n.total.k, prob_k, odds, p1, IsExact, 
IsSKAT_0=FALSE, r.corr=NULL, param=NULL){


	cat("Run")
	write.table(cbind(nres, rbind(nres_k)), File, append=FALSE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	write.table(rbind(resarray), File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)

	write.table(k, File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	write.table(m, File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	write.table(n.total, File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	write.table(rbind(sprintf("%d", n.total.k)), File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	write.table(rbind(IsExact), File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	write.table(rbind(prob_k), File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	write.table(rbind(odds), File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	write.table(rbind(p1), File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	
	write.table(rbind(Z0), File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	write.table(rbind(Z1), File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)

	if(IsSKAT_0){
		# param example
		# a<-matrix(1:9, byrow=TRUE, ncol=3); as.vector(t(a))
		param1<-as.vector(t(param))
		write.table(rbind(r.corr), File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
		write.table(rbind(param1), File, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
	}
}



