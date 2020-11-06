#define RAND_MAX_32 4294967295.0
#include <stdlib.h>
/* L = CRT_sum_mex(x,r,rand(sum(x),1),max(x));   
 or L = CRT_sum_mex(x,r);   */
        
void CRT_sum( double *x, double *r, int *Lenx, int *Lsum) {
    int i, j;
    double *prob;
    double maxx; 
    
    for(i=0, maxx=0;i<Lenx[0];i++)
       if (maxx<x[i]) maxx = x[i];
    
    prob = (double *) calloc(maxx, sizeof(double));
    
    for(i=0;i<maxx;i++)
        prob[i] = r[0]/(r[0]+i);
    
    for(i=0, Lsum[0]=0;i<Lenx[0];i++)
        for(j=0;j<x[i];j++) {
            /*if ( ((double) randomMT() / (double) 4294967296.0) <prob[j])     Lsum[0]++;*/
            /*if  ((double) randomMT() <= prob[j]*RAND_MAX_32)     Lsum[0]++;*/
            if ( (double) rand() <= prob[j]*RAND_MAX)     Lsum[0]++;
        }
}
