/* Vector operations*/

#ifndef VECTOROPERATIONS_H
#define VECTOROPERATIONS_H


#include <iostream>
#include <complex>

#pragma omp declare simd
void VectorMul(int N, const std::complex<double> * a, 
    const std::complex<double> * b, std::complex<double> * c);

void VectorMulReCmp(int N, const double* a, const std::complex<double> * b,
    std::complex<double> * c);

void VectorMulImCmp(int N, const double* a, const std::complex<double> * b,
    std::complex<double> * c);

void VectorMul(int N, const double * a, const double * b, double * c);

void VectorAdd(int N, const double * a, const double * b, double *c);

int VectorSum(int N, const bool* a);

#endif

