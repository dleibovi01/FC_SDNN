

#include <iostream>
#include "MVOperations.h"
#include "printing.h"

void MVMult_coldom(const int N, const int M, const double * A, const double * B, double * C)
{
    
    double entries;
    // std::cout << "A" << std::endl;
    // Print_Mat(A, N, 1);
    // std::cout << std::endl;
    // std::cout << "B" << std::endl;
    // Print_Mat(B, M, 1);
    for(int i = 0; i < N; i++)
    {
        entries = 0.0;
        for(int j = 0; j < M; j++)
        {
            entries += A[j*N + i]*B[j];
            // std::cout << "entries = " << entries << std::endl;
            // std::cout << std::endl;
        }
        C[i] = entries;
    }
}


void MVMult_rowdom(const int N, const int M, const double * A, const double * B, double * C)
{
    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < M; j++)
        {
            C[i] += A[i*M + j]*B[j];
        }
    }
}

// void MMMult_rowdom(const int N, const int K, const int M, const double * A, 
//     const double * B, double * C)
// {
//     for(int j = 0; j < M; j++)
//     {
//         for(int i = 0; i < N; i++)
//         {

//         }
//     }
// }