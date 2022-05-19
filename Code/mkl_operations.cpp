/* mkl wrapper */

#include "mkl_operations.h"


void Matrix_mult(int N, int M, int L, double * A, double *B, double *C)
{
    CBLAS_TRANSPOSE TRANSA = CblasNoTrans;    
    CBLAS_TRANSPOSE TRANSB = CblasTrans;    
    double alpha = 1.;
    double beta = 0.;  
    cblas_dgemm (CblasColMajor, TRANSA, TRANSB, N, M, L, alpha, A, N, B, L, 
    beta, C, N); 
}