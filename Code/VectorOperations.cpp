#include <iostream>
#include <complex>
#include "VectorOperations.h"


void VectorMul(int N, const std::complex<double> * a, const std::complex<double> * b,
    std::complex<double> * c)
{
    // #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        c[i] = a[i] * b[i];
    }
}

void VectorMulReCmp(int N, const double * a, const std::complex<double> * b,
    std::complex<double> * c)
{
    // 
    double * data_b = const_cast<double* >(reinterpret_cast<const double *>(b));
    double * data_c = reinterpret_cast<double *>(c);
    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        data_c[2*i] = data_b[2*i]*a[i];
        data_c[2*i + 1] = data_b[2*i + 1]*a[i];
    }
}

void VectorMulImCmp(int N, const double * a, const std::complex<double> * b,
    std::complex<double> * c)
{
    // #pragma omp simd
    double * data_b = const_cast<double* >(reinterpret_cast<const double *>(b));
    double * data_c = reinterpret_cast<double *>(c);
    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        data_c[2*i] = - data_b[2*i + 1]*a[i];
        data_c[2*i + 1] = data_b[2*i]*a[i];
    }
}


void VectorMul(int N, const double * a, const double * b, double * c)
{
    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        c[i] = a[i] * b[i];
    }
}

void VectorAdd(int N, const double * a, const double * b, double * c)
{
    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        c[i] = a[i] + b[i];
    }
}

int VectorSum(int N, const bool* a)
{
    int sum = 0;
    for(int i = 0; i < N; i++)
    {
        sum += !a[i];
    }
    return sum;
}