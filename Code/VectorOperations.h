/* Vector operations*/

#ifndef VECTOROPERATIONS_H
#define VECTOROPERATIONS_H


#include <iostream>
#include <complex>


void VectorMul(int N, const std::complex<double> * a, 
    const std::complex<double> * b, std::complex<double> * c);

void VectorMul(int N, const double * a, const double * b, double * c);

#endif

