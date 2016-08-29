Mk_SSD_ALL_Output<-function(out.list){
	
	N.Set<-length(out.list)
	OUT.Pvalue<-rep(NA,N.Set)
	OUT.Marker<-rep(NA,N.Set)
	OUT.Marker.Test<-rep(NA,N.Set)
	OUT.Error<-rep(-1,N.Set)
	OUT.Pvalue.Resampling<-NULL

	Is.Resampling=FALSE
	if(!is.null(out.list[[1]]$p.value.resampling)){
		Is.Resampling=TRUE
		n.Resampling<-length(out.list[[1]]$p.value.resampling)
		OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
	
	}
	
	for(i in 1:N.Set){
		re<-out.list[[i]]
		OUT.Pvalue[i]<-re$p.value
		OUT.Marker[i]<-re$param$n.marker
		OUT.Marker.Test[i]<-re$param$n.marker.test
		if(Is.Resampling){
			OUT.Pvalue.Resampling[i,]<-re$p.value.resampling
		}
	}
	
	out.tbl<-data.frame(SetID=1:N.Set, P.value=OUT.Pvalue, N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test)
	re<-list(results=out.tbl,P.value.Resampling=OUT.Pvalue.Resampling)
	class(re)<-"SKAT_SSD_ALL"

	return(re)

}
