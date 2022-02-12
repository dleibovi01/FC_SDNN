#include <iostream>
#include "mkl.h"
#include "mkl_dfti.h"
#include <complex.h>
//#include "mat.h"
#include <fstream>
#include <string>
#include "time.h"
#include "FC.h"


void Print_Mat(double * A, int nrows, int ncols)
{
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            std::cout << A[j*nrows + i] << " ";
        }
        std::cout << std::endl;
    }   
}



double _Complex * Fcont_Gram_Blend(const double * fx, int N, int d, int C, double fourPts, const double * AQ, const double * FAQF, const DFTI_DESCRIPTOR_HANDLE &desc_handle)
{
    double * fr = new double[d];
    double * fl = new double[d];
    int incx = 1;
    int incy = 1;
    double alpha = 1.;
    double beta = 0.;    
    int lda = C;
    CBLAS_LAYOUT Layout = CblasColMajor;
    CBLAS_TRANSPOSE TRANS = CblasNoTrans;    
    double  *z = new double[C];
    double *fx0 = new double[N];
    double _Complex *f_ext = new double _Complex[N + C];
    MKL_LONG status;
    cblas_daxpy(N, (1.0/fourPts), fx, incx, fx0, incy);
    for(int i = 0; i < C; i++)
    {
        z[i] = 0.;
    }
    for(int i = 0; i < d; i++)
    {
        fr[i] = fx0[N - d + i];
        fl[i] = fx0[i];
    }
    cblas_dgemv (Layout, TRANS, C, d, alpha, AQ, lda, fr, incx, beta, z, incy);
    cblas_dgemv (Layout, TRANS, C, d, alpha, FAQF, lda, fl, incx, alpha, z, incy);
    std::copy(fx0, fx0 + N, f_ext);
    std::copy(z, z + C, f_ext + N);
    status = DftiComputeForward(desc_handle, f_ext); 

    delete [] fx0;
    delete [] fl;
    delete [] fr;
    delete [] z;
    return f_ext;

}



int main()
{ 

    int C = 27;
    int d = 5;
    int incx = 1;
    int incy = 1;
    double alpha = 1.;
    double beta = 0.;    
    int lda = C;
    CBLAS_LAYOUT Layout = CblasColMajor;
    CBLAS_TRANSPOSE TRANSA = CblasNoTrans;    
    CBLAS_TRANSPOSE TRANSQ = CblasTrans;  
    double  *z = new double[C];
    MKL_LONG status;
    DFTI_DESCRIPTOR_HANDLE desc_handle;  
    double *A = new double[C*d];
    double *Q = new double[d*d];
    double *F1 = new double [d*d];
    double *F2 = new double [C*C];

    std::string filename_A = "Ad5C27.txt";
    std::string filename_Q = "Qd5C27.txt";
    read_FC_Data(A, Q, d, C, filename_A, filename_Q);
    // int i = 0;
    // double data;
    // std::ifstream Adata ("Ad5C27.txt");
    // if(Adata.is_open())
    // {
    //     while(Adata >> data)
    //     {
    //         A[i] = data;
    //         i = i + 1;       
    //     }
    // }

    // i = 0;
    // std::ifstream Qdata ("Qd5C27.txt");
    // if(Qdata.is_open())
    // {
    //     while(Qdata >> data)
    //     {
    //         Q[i] = data;
    //         i = i + 1;       
    //     }
    // }

    for(int i = 0; i < d; i++)
    {
        for(int j = 0; j < d; j++)
        {
            if (j == d - 1 - i)
            {
                F1[i*d + j] = 1.0;
            }
            else
            {
                F1[i*d + j] = 0.0;
            }
            
        }
    }

    for(int i = 0; i < C; i++)
    {
        for(int j = 0; j < C; j++)
        {
            if (j == C - 1 - i)
            {
                F2[i*C + j] = 1.0;
            }
            else
            {
                F2[i*C + j] = 0.0;
            }
            
        }
    }

    double *AQ = new double[C*d];
    double *prod1 = new double[d*d];
    double *prod2 = new double[C*d];
    double *FAQF = new double[C*d];
    cblas_dgemm (CblasColMajor, TRANSA, TRANSQ, C, d, d, alpha, A, C, Q, d, beta, AQ, C);

    //Print_Mat(AQ, C, d);
    cblas_dgemm (CblasColMajor, TRANSQ, TRANSA, d, d, d, alpha, Q, d, F1, d, beta, prod1, d);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSA, C, d, d, alpha, A, C, prod1, d, beta, prod2, C);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSA, C, d, C, alpha, F2, C, prod2, C, beta, FAQF, C);
    // Print_Mat(FAQF, C, d);

    int N = 100;
    double *f = new double[N];
    double _Complex *f_ext = new double _Complex[N + C];
    for (int i = 0; i < N; i++)
    {
        f[i] = double(i);
    }

    
    status = DftiCreateDescriptor( &desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, N + C);
    status = DftiCommitDescriptor( desc_handle );

    int Ntrials = 1000000;
    double fourPts = double(N + C);
    clock_t t;
    t = clock();
    for(int i = 0; i < Ntrials; i++)
    {
        f_ext = Fcont_Gram_Blend(f, N, d, C, fourPts, AQ, FAQF, desc_handle);
        delete [] f_ext;
        double _Complex *f_ext = new double _Complex[N + C];

    }
    t = clock() - t;
    std::cout << double(t)/CLOCKS_PER_SEC << std::endl;
    status = DftiFreeDescriptor(&desc_handle); 


    delete [] A;
    delete [] Q;
    delete [] F1;
    delete [] F2;
    delete [] FAQF;
    delete [] AQ;  
    delete [] prod1;
    delete [] prod2;   


    return 0;   


}