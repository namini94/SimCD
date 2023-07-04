#define RAND_MAX_32 4294967295.0
#define MAX(a,b) ((a) > (b) ? a : b)
#include <stdlib.h>
        
void CRT_vector( double *x, double *r, int *Lenx, int *L) {
    int i, j;
    double *prob;

    for( i=0;i<Lenx[0];i++){
        prob = (double *) calloc(x[i], sizeof(double));   
        for(j=0, L[i]=0; j<x[i]; j++) {
            prob[j] = r[i]/(r[i]+j);
            if  ((double) rand() <= prob[j]*RAND_MAX)     L[i]++;
        }
    }
}
