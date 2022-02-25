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