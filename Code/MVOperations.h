/* matrix Vector operations*/

#ifndef MVOPERATIONS_H
#define MVOPERATIONS_H

#include <iostream>

#pragma omp declare simd
void MVMult_coldom(const int N, const int M, const double * A, const double * B,
    double * C);

void MVMult_rowdom(const int N, const int M, const double * A, const double * B,
    double * C);

// First matrix is N*K and row dominant, second is K*N and col dominant.
// Output is row dominant
// void MMMult_rowdom(const int N, const int K, const int M, const double * A, 
//     const double * B, double * C);

#endif