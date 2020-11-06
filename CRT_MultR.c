#include <stdlib.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define RAND_MAX_32 4294967295.0

int BinarySearch(double probrnd, double *prob_cumsum, int Ksize) {
    int k, kstart, kend;
    if (probrnd <=prob_cumsum[0])
        return(0);
    else {
        for (kstart=1, kend=Ksize-1; ; ) {
            if (kstart >= kend) {
                /*//k = kend;*/
                return(kend);
            }
            else {
                k = kstart+ (kend-kstart)/2;
                if (prob_cumsum[k-1] > probrnd && prob_cumsum[k]>probrnd)
                    kend = k-1;
                else if (prob_cumsum[k-1]<probrnd && prob_cumsum[k]<probrnd)
                    kstart = k+1;
                else
                    return(k);
            }
        }
    }
    return(k);
}
        
void CRT_MultR( double *y, double *Phi, double *Theta, int *V, int *K, int *J, int *lkj, int *lvk) {
    int k, j, v, token, total=0, table;
    double *prob_cumsum;
    double cum_sum, probrnd;
    
    prob_cumsum = (double *) calloc(K[0], sizeof(double));
    
    for(j=0;j<J[0];j++){
        for(v=0;v<V[0];v++){
            for(cum_sum=0,k=0;k<K[0];k++){
                cum_sum += Phi[v+ k*V[0]]*Theta[k + K[0]*j];
                prob_cumsum[k] = cum_sum;
            }
            if (y[total]<0.5)
                table = 0;
            else{
                for (token=1, table=1;token < (int) y[total];token++){
                    if  (((double) rand() / RAND_MAX) <= (cum_sum/(cum_sum+ token)))
                        table++;
                }
            }
            
            for (token=0;token< table;token++){
                probrnd = (double) rand() / RAND_MAX *cum_sum;
                k = BinarySearch(probrnd, prob_cumsum, K[0]);
                lkj[k+j*K[0]]++;
                lvk[v+k*V[0]]++;
            }
            total++;
        }
    }
    free(prob_cumsum);    
    
}

