
K1_Help= function(x,y){
  # Helper function for 2 way interaction kernel
  p = length(x)
  a = x*y
  b = cumsum(a)
  return(sum(a[-1]*b[-p]))
}

call_Kernel_IBS<-function(Z,n,p){

	#Kernel_IBS(double * Z, int * pn, int * pp, double * Kernel)
	K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
	temp<-.C("Kernel_IBS",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.double(as.vector(K)))[[4]]
	matrix(temp,nrow=n)
}

call_Kernel_IBS_Weight<-function(Z,n,p,weights){

	#Kernel_IBS_Weight(int * Z, int * pn, int * pp, int *UseGivenWeight ,  double * weight, double * Kernel)
	given_weight = 1;
	if( is.null(weights)){
		weights = rep(0,p);
		given_weight = 0;
	} else {
		# change!!
		weights<-weights^2;
	}
	K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
	temp<-.C("Kernel_IBS_Weight",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.integer(given_weight),
	as.double(weights),as.double(as.vector(K)))[[6]]
	matrix(temp,nrow=n)
}

call_Kernel_2wayIX<-function(Z,n,p){

	#Kernel_IBS(double * Z, int * pn, int * pp, double * Kernel)
	K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
	temp<-.C("Kernel_2wayIX",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.double(as.vector(K)))[[4]]
	matrix(temp,nrow=n)
}

lskmTest.GetKernel = function(Z, kernel, weights,n,m){

    	if (kernel == "quadratic") {
      		K = (Z%*%t(Z)+1)**2
    	}


	if (kernel == "IBS") {
      		K = call_Kernel_IBS(Z,n,m)
    	}
    	if (kernel == "IBS.weighted") {

      		K = call_Kernel_IBS_Weight(Z,n,m,weights)
    	}
  	if (kernel == "2wayIX") {
      		K = call_Kernel_2wayIX(Z,n,m)
    	}  
   	if (kernel == "IBS.weighted_OLD") {
      		#K = matrix(nrow = n, ncol = n)
      		if (is.null(weights)) {
        		qs = apply(Z, 2, mean)/(2)
        		weights = 1/sqrt(qs)
      		} else {
			weights<-weights^2
		}
      		K1 = matrix(nrow =n, ncol = n)
      		for (i in 1:n) {
        		K1[i,] = apply(abs(t(Z)-Z[i,])*weights,2, sum)
      		}
      		K= 1-(K1)/(2*sum(weights))
    	}

    	if (kernel == "IBS_OLD") {
      		K1=matrix(nrow=n,ncol=n)
      		for (i in 1:n) {
        		K1[i,] = apply(abs(t(Z)-Z[i,]),2, sum)
      		}
      		K = (2*m-K1)/(2*m)
    	}
   	if (kernel == "2wayIX_OLD") {
      		K = 1+Z%*%t(Z)
      		N1=  matrix(nrow = n, ncol = n)
      		for (i in 1:n){
        		for (j in i:n){
	    			N1[j,i] = N1[i,j] = K1_Help(Z[i,], Z[j,])
	  		}
      		}
      		K = K+N1
    	}
	return(K)

}

#.First.lib <- function(lib, pkg) { library.dynam('SKAT', pkg, lib) } 


