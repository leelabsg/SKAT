#include <R.h>
#include <Rmath.h>
#include <math.h>
/********************************************************
*
*	Kernel function ver 0.2
*
*
********************************************************/


/********************************************************
 *
 *  C style index (Row-wise)
 *
 *              K1=matrix(nrow=n,ncol=n)
      		for (i in 1:n) {
        		K1[i,] = apply(abs(t(Z)-Z[i,]),2, sum)
      		}
      		K = (2*m-K1)/(2*m)
 *
 *  Z : n by p matrix
 *  Kernel : output n by n matrix
 * ******************************************************/

void Kernel_IBS_1(int * Z, int * pn, int * pp, double * Kernel){

    int i,j,k,n,p,diff;
    double temp;

    n = *pn;
    p = *pp;

   /*Rprintf("Dimension: n [%d], p[%d] \n",n,p);*/

    for(i=0;i<n-1;i++){
        for(j=i+1;j<n;j++){
            temp = 0;
            for(k=0;k<p;k++){
                diff = 2 - abs(Z[i*p + k] - Z[j*p + k]);
		temp += diff ;
               /* Rprintf("ind1 [%d][%d],  ind2 [%d][%d]  \n",i,Z[i*p + k],j,Z[j*p + k]);
               */
		/*
                if(diff == 0)
                    temp = temp + 2;
                else if(diff ==1 || diff == -1)
                    temp = temp +1;
		*/
            }
            Kernel[i*n + j] = Kernel[j*n + i] = temp / 2/ p;
        }
    }

    for(i=0;i<n;i++){
        Kernel[i*n + i] = 1;
    }
    
}

/********************************************
 *
 *  IBS Weight

      		if (is.null(weights)) {
        		qs = apply(Z, 2, mean)/(2)
        		weights = 1/sqrt(qs)
      		}
      		K1 = matrix(nrow =n, ncol = n)
      		for (i in 1:n) {
        		K1[i,] = apply(abs(t(Z)-Z[i,])*weights,2, sum)
      		}
      		K= 1-(K1)/(2*sum(weights))
 *
 * *****************************************************/


void Kernel_IBS_Weight_1(int * Z, int * pn, int * pp, int *UseGivenWeight ,  double * weight, double * Kernel){

    int i,j,k,n,p,diff,temp1;
    double temp, w_total;


    n = *pn;
    p = *pp;

   /* Rprintf("Dimension: n [%d], p[%d] \n",n,p);*/

    if(*UseGivenWeight != 1){
        for(k=0;k<p;k++){
            temp1 = 0;
            for(i=0;i<n;i++){
                temp1 += Z[i*p + k] ;
            }
            weight[k] = sqrt(2.0 * p) / sqrt((double)temp1);
        }
    }
    
    w_total = 0;
    for(k=0;k<p;k++){
        w_total += weight[k];
    }


    for(i=0;i<n-1;i++){
        for(j=i+1;j<n;j++){
            temp = 0;
            for(k=0;k<p;k++){
                diff =  abs(Z[i*p + k] - Z[j*p + k]);
		temp += diff*weight[k] ;

               /* Rprintf("ind1 [%d][%d],  ind2 [%d][%d]  \n",i,Z[i*p + k],j,Z[j*p + k]);
               */
		/*
                if(diff == 2 || diff == -2)
                    temp += 2*weight[k];
                else if(diff ==1 || diff == -1)
                    temp += weight[k];
		*/
            }
            Kernel[i*n + j] = Kernel[j*n + i] = 1 - temp / 2 / w_total;
        }
    }

    for(i=0;i<n;i++){
        Kernel[i*n + i] = 1;
    }

}

/********************************************************
 *
 *  2wayIX Kernel
      		K = 1+Z%*%t(Z)
      		N1=  matrix(nrow = n, ncol = n)
      		for (i in 1:n){
        		for (j in i:n){
	    			N1[j,i] = N1[i,j] = K1(Z[i,], Z[j,])
	  		}
      		}
      		K = K+N1
 *
 * K1= function(x,y){
  p = length(x)
  a = x*y
  b = cumsum(a)
  return(sum(a[-1]*b[-p]))
}
 *
 *  Z : n by p matrix
 *  Kernel : output n by n matrix
 * ******************************************************/

void Kernel_2wayIX_1(int * Z, int * pn, int * pp, double * Kernel){

    int i,j,k,n,p;
    double temp, temp1, temp2;

    n = *pn;
    p = *pp;

   /*Rprintf("Dimension: n [%d], p[%d] \n",n,p);*/

    for(i=0;i<n;i++){
        for(j=i;j<n;j++){
            temp = 1;

            for(k=0;k<p;k++){
                temp2 = Z[i*p + k] * Z[j*p + k];
                temp += temp2;
                if(k==0){
                    temp1 = temp2;
                    continue ;
                }
                
                temp += temp1* Z[i*p + k] * Z[j*p + k];
                temp1 += temp2;

            }
            Kernel[i*n + j] = Kernel[j*n + i] = temp ;
        }
    }
}


