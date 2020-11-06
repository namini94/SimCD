#include <stdlib.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define RAND_MAX_32 4294967295.0
        
void CRT_sum_matrix( double *x, double *r, int *K, int *J, int *Lsum) {
    int i, j, token;
    double *prob;
    double maxx; 
    
    for(j=0;j<J[0];j++){
        for(i=0, maxx=0;i<K[0];i++)
            maxx = MAX(maxx,x[i+j*K[0]]);
        prob = (double *) calloc(maxx, sizeof(double));
        for(i=0;i<maxx;i++)
            prob[i] = r[0]/(r[0]+i);
        for(Lsum[j]=0,i=0;i<K[0];i++)
            for(token=0;token<x[i+j*K[0]];token++) {
                if ( (double) rand() <= prob[i]*RAND_MAX)
                    Lsum[j]++; 
            }
        free(prob);    
    }
    
}
