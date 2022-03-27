/* FC related functions */

#ifndef FC_H
#define FC_H

#include <iostream>
#include <complex>
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include "mkl.h"
#include "mkl_dfti.h"
#include <fstream>
#include <string>
#include <math.h>
#include "VectorOperations.h"


void read_FC_Data(double *A, double *Q, int d, int C, std::string filename_A, 
    std::string filename_Q);

void build_Cont_Mat(const double *A, const double *Q, int d, int C, double *AQ, 
    double *FAQF);

void build_Cont_Mat_DN(const double *A, const double *Q, const double* Q_tilde, 
    int d, int C, double *AQ, double *FAQF);

void build_Cont_Mat_ND(const double *A, const double *Q, const double* Q_tilde, 
    int d, int C, double *AQ, double *FAQF);

void Fcont_Gram_Blend(const double * fx, std::complex<double> *, int N, int d, 
    int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle);

void Fcont_Gram_Blend_DD(const double * fx, std::complex<double> *, int N, int d, 
    int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle);

void Fcont_Gram_Blend_DN(const double * fx, std::complex<double> *, int N, int d, 
    int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle, double h, double bc_r);

void Fcont_Gram_Blend_ND(const double * fx, std::complex<double> *, int N, int d, 
    int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle, double h, double bc_l);

void Fcont_Gram_Blend_NN(const double * fx, std::complex<double> *, int N, int d, 
    int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle, double h, double bc_l,
    double bc_r);

double * FC_Der(const double * fx, std::complex<double> * der_coeffs, 
    std::complex<double> * filter_coeffs, int N, int d, int C, double fourPts, 
    const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle);

void FC_Der(const double * fx, double *f_der, const double * der_coeffs, 
    int N, int d, int C, double fourPts, const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle);

void FC_Der(double *f_der, const std::complex<double> * f_hat, 
    const double * der_coeffs, int N, int C,
    const DFTI_DESCRIPTOR_HANDLE &desc_handle);

void FC_Der(double *f_der, const std::complex<double> * f_hat, 
    const double * der_coeffs, int N, int C,
    const DFTI_DESCRIPTOR_HANDLE &desc_handle, bool flag);

double * Fcont(const double * fx, int N, int d, int C, double fourPts, 
    const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle);

void Fcont(const double * fx, double * f_der, int N, int d, int C,
     double fourPts, const double * AQ, const double * FAQF, 
     const DFTI_DESCRIPTOR_HANDLE &desc_handle);

void Fcont_shift(const double * fx, double * f_shift, 
    std::complex<double> * shift_coeffs, int N, int d, int C, double fourPts, 
    const double * AQ, const double * FAQF, 
    const DFTI_DESCRIPTOR_HANDLE &desc_handle);

void getK(int *k, int fourPts);  

void getFiltCoeffs(std::complex<double>* filt_coeffs, int fourPts, double alpha,
    double p);

void getFiltCoeffs(double* filt_coeffs, int fourPts, double alpha,
    double p);

void getF1_F2(const int d, const int C, double* F1, double* F2);

#endif 