/* Definition of ANN routines. */

#include "ANN.h"


ANN::ANN(std::string W1_filename, std::string W2_filename, 
    std::string W3_filename, std::string W4_filename, 
    std::string B1_filename, std::string B2_filename, 
    std::string B3_filename, std::string B4_filename, double _alpha) : 
    alpha{_alpha}
{
    set_NN_weights(W1_filename, B1_filename, W2_filename, B2_filename, 
        W3_filename, B3_filename, W4_filename, B4_filename);
}


void ANN::set_NN_weights(std::string W1_file, std::string B1_file, 
    std::string W2_file, std::string B2_file, std::string W3_file, 
    std::string B3_file, std::string W4_file, std::string B4_file)
{ 
    W1 = new double[7*16];
    B1 = new double[16];
    W2 = new double[16*16];
    B2 = new double[16];
    W3 = new double[16*16];
    B3 = new double[16];
    W4 = new double[4*16];
    B4 = new double[4*16];    
    readMatrix(W1, W1_file);  
    readMatrix(W2, W2_file); 
    readMatrix(W3, W3_file); 
    readMatrix(W4, W4_file); 
    readMatrix(B1, B1_file); 
    readMatrix(B2, B2_file); 
    readMatrix(B3, B3_file); 
    readMatrix(B4, B4_file);                 
} 


void ANN::set_NN_weights()
{
    W1 = new double[7*16];
    B1 = new double[16];
    W2 = new double[16*16];
    B2 = new double[16];
    W3 = new double[16*16];
    B3 = new double[16];
    W4 = new double[4*16];
    B4 = new double[4];   
    std::copy(W1_data.data(), W1_data.data() + 112, W1);
    std::copy(W2_data.data(), W2_data.data() + 256, W2);
    std::copy(W3_data.data(), W3_data.data() + 256, W3);
    std::copy(W4_data.data(), W4_data.data() + 64, W4);
    std::copy(B1_data.data(), B1_data.data() + 16, B1);
    std::copy(B2_data.data(), B2_data.data() + 16, B2);
    std::copy(B3_data.data(), B3_data.data() + 16, B3);
    std::copy(B4_data.data(), B4_data.data() + 4, B4);
} 


ANN & ANN::operator= (const ANN &ann)
{
    if(this == &ann)
    {
        return *this;
    }
    else
    {
        delete[] W1;
        delete[] W2;
        delete[] W3;
        delete[] W4;
        delete[] B1;
        delete[] B2;
        delete[] B3;
        delete[] B4;   
        set_NN_weights();  
        return *this;                    
    }   
}

void ANN::free_mem()
{
    delete[] W1;
    delete[] W2;
    delete[] W3;
    delete[] W4;
    delete[] B1;
    delete[] B2;
    delete[] B3;
    delete[] B4;       
}


void ANN::fwdClassif(int *tau, const double *regStencils, int N0)
{
    if(N0 > 0)
    {
        int M = 16;
        int output = 4;            
        double input_1[M*N0];
        double input_2[M*N0];
        double input_3[M*N0];
        double output_4[output*N0];
        double ones[N0];

        #pragma omp simd
        for(int i = 0; i < M*N0; i++)
        {
            input_1[i] = 0.0;
            input_2[i] = 0.0;
            input_3[i] = 0.0;
        }
        #pragma omp simd
        for(int i = 0; i < output*N0; i++)
        {
            output_4[i] = 0.0;
        }
        #pragma omp simd
        for(int i = 0; i < N0; i++)
        {
            ones[i] = 1.0;
        }

        cblas_dger (CblasColMajor, M, N0, 1.0, B1, 1, ones, 1, input_1, M);
        cblas_dger (CblasColMajor, M, N0, 1.0, B2, 1, ones, 1, input_2, M);
        cblas_dger (CblasColMajor, M, N0, 1.0, B3, 1, ones, 1, input_3, M);
        cblas_dger (CblasColMajor, output, N0, 1.0, B4, 1, ones, 1,
            output_4, output);

        cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, M, N0, s,
            1.0, W1, M, regStencils, s, 1.0, input_1, M);
        elu(M*N0, input_1, alpha);

        cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, M, N0, M, 1.0,
            W2, M, input_1, M, 1.0, input_2, M);
        elu(M*N0, input_2, alpha);

        cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, M, N0, M, 1.0,
            W3, M, input_2, M, 1.0, input_3, M);
        elu(M*N0, input_3, alpha);

        cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, output, N0, M,
            1.0, W4, output, input_3, M, 1.0, output_4, output);

        #pragma omp simd
        for(int i = 0; i < N0; i++)
        {
            tau[i] = std::distance(output_4 + output*i, 
                std::max_element(output_4 + output*i,
                    output_4 + output*(i + 1))) + 1;
        }
    }        
}

inline void ANN::elu(const int N, double* input, const double alpha) const
{
    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        if(input[i] <= 0)
        {
            input[i] =  alpha*(std::exp(input[i]) - 1.0);
        } 
    }
}

void ANN::softmax(int N, double* input) const
{
    double sum_exp = 0;
    for(int i = 0; i < N; i++)
    {
        input[i] = std::exp(input[i]);
        sum_exp +=  input[i]; 
    }
    for(int i = 0; i < N; i++)
    {
        input[i] = input[i] / sum_exp;
    }    
}

void ANN::readMatrix(double* A, std::string filename_A)
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
}