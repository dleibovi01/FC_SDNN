#include <iostream>
#include <complex>
#include "VectorOperations.h"


void VectorMul(int N, const std::complex<double> * a, const std::complex<double> * b,
    std::complex<double> * c)
{
    for(int i = 0; i < N; i++)
    {
        c[i] = a[i] * b[i];
    }
}

void VectorMul(int N, const double * a, const double * b, double * c)
{
    for(int i = 0; i < N; i++)
    {
        c[i] = a[i] * b[i];
    }
}