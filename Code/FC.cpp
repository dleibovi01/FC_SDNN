/* FC related functions */

#include <iostream>
#include <complex>
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include "mkl.h"
#include "mkl_dfti.h"
#include <fstream>
#include <string>
#include "time.h"
#include "FC.h"
#include "printing.h"


void read_FC_Data(double *A, double *Q, int d, int C, std::string filename_A, 
    std::string filename_Q)
{
    int i = 0;
    double data;
    std::ifstream Adata (filename_A.c_str());
    if(Adata.is_open())
    {
        while(Adata >> data)
        {
            A[i] = data;
            i = i + 1;       
        }
    }
    i = 0;
    std::ifstream Qdata (filename_Q.c_str());
    if(Qdata.is_open())
    {
        while(Qdata >> data)
        {
            Q[i] = data;
            i = i + 1;       
        }
    }
}


 void build_Cont_Mat(const double *A, const double *Q, int d, int C, double *AQ, 
    double *FAQF)
 {
    double F1[d*d];
    double F2[C*C];
    double prod1[d*d];
    double prod2[C*d];
    CBLAS_TRANSPOSE TRANSA = CblasNoTrans;    
    CBLAS_TRANSPOSE TRANSQ = CblasTrans;    
    double alpha = 1.;
    double beta = 0.;      
    getF1_F2(d, C, F1, F2);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSQ, C, d, d, alpha, A, C, Q, d, 
        beta, AQ, C);
    cblas_dgemm (CblasColMajor, TRANSQ, TRANSA, d, d, d, alpha, Q, d, F1, d, 
        beta, prod1, d);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSA, C, d, d, alpha, A, C, prod1, d, 
        beta, prod2, C);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSA, C, d, C, alpha, F2, C, prod2, C, 
        beta, FAQF, C);        
 }


 void build_Cont_Mat_DN(const double *A, const double *Q, const double* Q_tilde,
    int d, int C, double *AQ, double *FAQF)
 {
    double F1[d*d];
    double F2[C*C];
    double prod1[d*d];
    double prod2[C*d];
    CBLAS_TRANSPOSE TRANSA = CblasNoTrans;    
    CBLAS_TRANSPOSE TRANSQ = CblasTrans;    
    double alpha = 1.;
    double beta = 0.;      
    getF1_F2(d, C, F1, F2);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSQ, C, d, d, alpha, A, C, Q_tilde, 
        d, beta, AQ, C);
    cblas_dgemm (CblasColMajor, TRANSQ, TRANSA, d, d, d, alpha, Q, d, F1, d, 
        beta, prod1, d);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSA, C, d, d, alpha, A, C, prod1, d, 
        beta, prod2, C);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSA, C, d, C, alpha, F2, C, prod2, C, 
        beta, FAQF, C);        
 }

 void build_Cont_Mat_ND(const double *A, const double *Q, const double* Q_tilde,
    int d, int C, double *AQ, double *FAQF)
 {
    double F1[d*d];
    double F2[C*C];
    double prod1[d*d];
    double prod2[C*d];
    CBLAS_TRANSPOSE TRANSA = CblasNoTrans;    
    CBLAS_TRANSPOSE TRANSQ = CblasTrans;    
    double alpha = 1.;
    double beta = 0.;      
    getF1_F2(d, C, F1, F2);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSQ, C, d, d, alpha, A, C, Q, d,
        beta, AQ, C);
    cblas_dgemm (CblasColMajor, TRANSQ, TRANSA, d, d, d, alpha, Q_tilde, d, F1,
        d, beta, prod1, d);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSA, C, d, d, alpha, A, C, prod1, d, 
        beta, prod2, C);
    cblas_dgemm (CblasColMajor, TRANSA, TRANSA, C, d, C, alpha, F2, C, prod2, C, 
        beta, FAQF, C);        
 }




void Fcont_Gram_Blend(const double * fx, std::complex<double> *f_ext, int N,
    int d, int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle)
{
    double fr[d];
    double fl[d];
    double scaling = 1.0/fourPts;
    static int incx = 1;
    static int incy = 1;
    static double alpha = 1.;
    static double beta = 0.;    
    static int lda = C;
    static CBLAS_LAYOUT Layout = CblasColMajor;
    static CBLAS_TRANSPOSE TRANS = CblasNoTrans;    
   
    MKL_LONG status;

    for(int i = 0; i < d; i++)
    {
        fr[i] = fx[N - d + i];
        fl[i] = fx[i];
    }
    double  z[C];
    for(int i = 0; i < C; i++)
    {
        z[i] = 0.;
    }
    cblas_dgemv (Layout, TRANS, C, d, alpha, AQ, lda, fr, incx, beta, z, incy);
    cblas_dgemv (Layout, TRANS, C, d, alpha, FAQF, lda, fl, incx, alpha, z, incy);
    std::copy(z, z + C, f_ext + N);
    std::copy(fx, fx + N , f_ext); 


    status = DftiComputeForward(desc_handle, f_ext); 
    #pragma omp simd
    for (int i = 0; i < N + C; i++)
    {
        f_ext[i] = f_ext[i]*scaling;
    }

}

void Fcont_Gram_Blend_DD(const double * fx, std::complex<double> *f_ext, int N,
    int d, int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle)
{
    double fr[d];
    double fl[d];
    double scaling = 1.0/fourPts;
    static int incx = 1;
    static int incy = 1;
    static double alpha = 1.;
    static double beta = 0.;    
    static int lda = C;
    static CBLAS_LAYOUT Layout = CblasColMajor;
    static CBLAS_TRANSPOSE TRANS = CblasNoTrans;    
   
    MKL_LONG status;

    for(int i = 0; i < d; i++)
    {
        fr[i] = fx[N - d + i];
        fl[i] = fx[i];
    }
    double  z[C];
    for(int i = 0; i < C; i++)
    {
        z[i] = 0.;
    }
    cblas_dgemv (Layout, TRANS, C, d, alpha, AQ, lda, fr, incx, beta, z, incy);
    cblas_dgemv (Layout, TRANS, C, d, alpha, FAQF, lda, fl, incx, alpha, z, incy);
    std::copy(z, z + C, f_ext + N);
    std::copy(fx, fx + N , f_ext); 


    status = DftiComputeForward(desc_handle, f_ext); 
    #pragma omp simd
    for (int i = 0; i < N + C; i++)
    {
        f_ext[i] = f_ext[i]*scaling;
    }

}

void Fcont_Gram_Blend_DN(const double * fx, std::complex<double> *f_ext, int N,
    int d, int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle, double h, double bc_r)
{
    double fr[d];
    double fl[d];
    static double h0 = 0.01;
    double scaling = 1.0/fourPts;
    static int incx = 1;
    static int incy = 1;
    static double alpha = 1.;
    static double beta = 0.;    
    static int lda = C;
    static CBLAS_LAYOUT Layout = CblasColMajor;
    static CBLAS_TRANSPOSE TRANS = CblasNoTrans;    
   
    MKL_LONG status;

    for(int i = 0; i < d - 1; i++)
    {
        fr[i] = fx[N - d + i];
        fl[i] = fx[i];
    }
    fl[d - 1] = fx[d - 1];
    fr[d - 1] = bc_r * h / h0;
    double  z[C];
    for(int i = 0; i < C; i++)
    {
        z[i] = 0.;
    }
    cblas_dgemv (Layout, TRANS, C, d, alpha, AQ, lda, fr, incx, beta, z, incy);
    cblas_dgemv (Layout, TRANS, C, d, alpha, FAQF, lda, fl, incx, alpha, z, incy);
    std::copy(z, z + C, f_ext + N);
    std::copy(fx, fx + N , f_ext); 


    status = DftiComputeForward(desc_handle, f_ext); 
    #pragma omp simd
    for (int i = 0; i < N + C; i++)
    {
        f_ext[i] = f_ext[i]*scaling;
    }

}

void Fcont_Gram_Blend_ND(const double * fx, std::complex<double> *f_ext, int N,
    int d, int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle, double h, double bc_l)
{
    double fr[d];
    double fl[d];
    static double h0 = 0.01;
    double scaling = 1.0/fourPts;
    static int incx = 1;
    static int incy = 1;
    static double alpha = 1.;
    static double beta = 0.;    
    static int lda = C;
    static CBLAS_LAYOUT Layout = CblasColMajor;
    static CBLAS_TRANSPOSE TRANS = CblasNoTrans;    
   
    MKL_LONG status;

    for(int i = 1; i < d; i++)
    {
        fr[i] = fx[N - d + i];
        fl[i] = fx[i];
    }
    fl[0] = - bc_l*h/h0;
    fr[0] = fx[0];
    double  z[C];
    for(int i = 0; i < C; i++)
    {
        z[i] = 0.;
    }
    cblas_dgemv (Layout, TRANS, C, d, alpha, AQ, lda, fr, incx, beta, z, incy);
    cblas_dgemv (Layout, TRANS, C, d, alpha, FAQF, lda, fl, incx, alpha, z, incy);
    std::copy(z, z + C, f_ext + N);
    std::copy(fx, fx + N , f_ext); 
    status = DftiComputeForward(desc_handle, f_ext); 
    #pragma omp simd
    for (int i = 0; i < N + C; i++)
    {
        f_ext[i] = f_ext[i]*scaling;
    }
}

void Fcont_Gram_Blend_NN(const double * fx, std::complex<double> *f_ext, int N,
    int d, int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle, double h, double bc_l, double bc_r)
{
    double fr[d];
    double fl[d];
    static double h0 = 0.01;
    double scaling = 1.0/fourPts;
    static int incx = 1;
    static int incy = 1;
    static double alpha = 1.;
    static double beta = 0.;    
    static int lda = C;
    static CBLAS_LAYOUT Layout = CblasColMajor;
    static CBLAS_TRANSPOSE TRANS = CblasNoTrans;    
   
    MKL_LONG status;

    for(int i = 1; i < d - 1; i++)
    {
        fr[i] = fx[N - d + i];
        fl[i] = fx[i];
    }
    fl[0] = - bc_l*h/h0;
    fr[0] = fx[N - d];
    fl[d - 1] = fx[d - 1];
    fr[d - 1] = bc_r * h / h0;


    double  z[C];
    for(int i = 0; i < C; i++)
    {
        z[i] = 0.;
    }
    cblas_dgemv (Layout, TRANS, C, d, alpha, AQ, lda, fr, incx, beta, z, incy);

    // std::cout << "z" << std::endl;
    // Print_Mat(z, C, 1);

    cblas_dgemv (Layout, TRANS, C, d, alpha, FAQF, lda, fl, incx, alpha, z, incy);
    std::copy(z, z + C, f_ext + N);
    std::copy(fx, fx + N , f_ext); 
    status = DftiComputeForward(desc_handle, f_ext); 
    #pragma omp simd
    for (int i = 0; i < N + C; i++)
    {
        f_ext[i] = f_ext[i]*scaling;
    }
}


double * FC_Der(const double * fx, std::complex<double> * der_coeffs, 
    std::complex<double> * filter_coeffs, int N, int d, int C, double fourPts, 
        const double * AQ, const double * FAQF, 
        const DFTI_DESCRIPTOR_HANDLE &desc_handle)
{
    MKL_LONG status;
    std::complex<double> f_ext[N + C];
    double * f_der = new double[N];
    Fcont_Gram_Blend(fx, f_ext, N, d, C, fourPts, AQ, FAQF, desc_handle);
    VectorMul(N + C, der_coeffs, f_ext, f_ext); 
    status = DftiComputeBackward(desc_handle, f_ext); 
    #pragma omp simd
    for (int j = 0; j < N; j++)
    {
        f_der[j] = f_ext[j].real();
    }
    return f_der;
}


// void FC_Der(const double * fx, double *f_der, std::complex<double> * der_coeffs, 
//     int N, int d, int C, double fourPts, 
//     const double * AQ, const double * FAQF, 
//     const DFTI_DESCRIPTOR_HANDLE &desc_handle)
void FC_Der(const double * fx, double *f_der, const double * der_coeffs, 
    int N, int d, int C, double fourPts, 
    const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle)
{
    MKL_LONG status;
    std::complex<double> f_ext[N + C];
    Fcont_Gram_Blend(fx, f_ext, N, d, C, fourPts, AQ, FAQF, desc_handle); 
    VectorMulImCmp(N + C, der_coeffs, f_ext, f_ext); 
    status = DftiComputeBackward(desc_handle, f_ext); 
    for (int j = 0; j < N; j++)
    {
        f_der[j] = f_ext[j].real();
    }
}

// void FC_Der(const double * fx, double *f_der, 
//     const std::complex<double> * der_coeffs, std::complex<double> * fft, int N,
//     int d, int C, double fourPts, const double * AQ, const double * FAQF, 
//     const DFTI_DESCRIPTOR_HANDLE &desc_handle, bool flag)
// {
//     MKL_LONG status;
//     Fcont_Gram_Blend(fx, fft, N, d, C, fourPts, AQ, FAQF, desc_handle); 
//     VectorMul(N + C, der_coeffs, fft, fft); 
//     status = DftiComputeBackward(desc_handle, fft); 
//     for (int j = 0; j < N; j++)
//     {
//         f_der[j] = fft[j].real();
//     }    
// }


// void FC_Der(double *f_der, const std::complex<double> * f_hat,
//     const std::complex<double> * der_coeffs, int N, int C, 
//     const DFTI_DESCRIPTOR_HANDLE &desc_handle)
void FC_Der(double *f_der, const std::complex<double> * f_hat,
    const double * der_coeffs, int N, int C, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle)
{

   std::complex<double> f_hat_temp[N+C];
   VectorMulImCmp(N + C, der_coeffs, f_hat, f_hat_temp);
   int status = DftiComputeBackward(desc_handle, f_hat_temp); 
   #pragma omp simd
   for (int j = 0; j < N; j++)
   {
       f_der[j] = f_hat_temp[j].real();
   } 
}

// void FC_Der(double *f_der, const std::complex<double> * f_hat,
//     const std::complex<double> * der_coeffs, int N, int C, 
//     const DFTI_DESCRIPTOR_HANDLE &desc_handle, bool flag)
void FC_Der(double *f_der, const std::complex<double> * f_hat,
    const double * der_coeffs, int N, int C, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle, bool flag)
{

   std::complex<double> f_hat_temp[N+C];
   VectorMulReCmp(N + C, der_coeffs, f_hat, f_hat_temp);
   int status = DftiComputeBackward(desc_handle, f_hat_temp); 
   #pragma omp simd
   for (int j = 0; j < N; j++)
   {
       f_der[j] = f_hat_temp[j].real();
   } 
}



double * Fcont(const double * fx, int N, int d, int C, double fourPts, 
    const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle)
{
    MKL_LONG status;
    std::complex<double> f_ext[N + C];
    double * f_der = new double[N];

    Fcont_Gram_Blend(fx, f_ext, N, d, C, fourPts, AQ, FAQF, desc_handle);
    status = DftiComputeBackward(desc_handle, f_ext); 
    for (int j = 0; j < N; j++)
    {
        f_der[j] = f_ext[j].real();
    }
    return f_der;
}

void Fcont(const double * fx, double * f_der, int N, int d, int C,
     double fourPts, const double * AQ, const double * FAQF, 
     const DFTI_DESCRIPTOR_HANDLE &desc_handle)
{
    MKL_LONG status;
    std::complex<double> f_ext[N + C];
    Fcont_Gram_Blend(fx, f_ext, N, d, C, fourPts, AQ, FAQF, desc_handle);
    status = DftiComputeBackward(desc_handle, f_ext); 
    #pragma omp simd
    for (int j = 0; j < N + C; j++)
    {
        f_der[j] = f_ext[j].real();
    }
}


void Fcont_shift(const double * fx, double * f_shift, 
    std::complex<double> * shift_coeffs, int N, int d, int C, double fourPts, 
    const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle)
{
    MKL_LONG status;
    std::complex<double> f_ext[N + C];
    Fcont_Gram_Blend(fx, f_ext , N, d, C, fourPts, AQ, FAQF, desc_handle);
    vzMul(N + C, shift_coeffs, f_ext, f_ext);
    status = DftiComputeBackward(desc_handle, f_ext); 
    for (int j = 0; j < N + C; j++)
    {
        f_shift[j] = f_ext[j].real();
    }
}


void getK(int * k, int fourPts)
{
   int l = 0;
   if (fourPts % 2 == 0)
      {
         for(int j = 0; j <= fourPts / 2; j++)
         {
               k[j] = j;
         }
         for(int j = fourPts / 2 + 1; j < fourPts; j++)
         {
               k[j] = -fourPts / 2 + 1 + l;
               l = l + 1;
         }   
      }     
   else
   {
      for(int j = 0; j <= (fourPts - 1)/2; j++)
      {
            k[j] = j;
      }
      for(int j = (fourPts + 1)/2 ; j < fourPts; j++)
      {
            k[j] = - (fourPts - 1)/2 + l;
            l = l + 1;
      }              
   }     
}


void getFiltCoeffs(std::complex<double> * filt_coeffs, int fourPts, 
    double alpha, double p)
{
    int k[fourPts];
    getK(k, fourPts);
    for(int i = 0; i < fourPts; i++)
    {
        filt_coeffs[i] = std::exp(- alpha*std::pow(2*k[i]/double(fourPts), p));
    }
}

void getFiltCoeffs(double * filt_coeffs, int fourPts, double alpha, double p)
{
    int k[fourPts];
    getK(k, fourPts);
    for(int i = 0; i < fourPts; i++)
    {
        filt_coeffs[i] = std::exp(- alpha*std::pow(2*k[i]/double(fourPts), p));
    }
}


void getF1_F2(const int d, const int C, double* F1, double* F2)
{
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
}