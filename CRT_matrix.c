
#include <stdlib.h>
#define RAND_MAX_32 4294967295.0
#define MAX(a,b) ((a) > (b) ? a : b)
        

void CRT_matrix(double *x, double *r, int *Msize, int *Nsize, int *L) {
    
    int i, j;  
    
    for (j=0;j<Msize[0]*Nsize[0];j++) 
    {
        for (i=0; i< (int) x[j]; i++)
        {
            if (((double) rand() / RAND_MAX) <= (r[j]/(r[j]+ i)))
                L[j]++;
        }
    }

            
}
